import os
import uuid
from tempfile import TemporaryDirectory
from pathlib import Path
from pymol import cmd
from psico import exporting

class PyMol:
    """PyMol environment
    """
    def __init__(self):
        self.group = uuid.uuid4().hex

    def __enter__(self):
        self._old_cwd = os.getcwd()
        # settings that cannot be applied to selections
        self._old_settings = {}
        # Private temp directory
        self._tempdir = TemporaryDirectory(prefix = "pymol_")
        self.tempdir = Path(self._tempdir.name)
        settings = {
            'pdb_conect_nodup': 0,
            'fetch_path': str(self.tempdir)
        }
        for setting in settings:
            self._old_settings[setting] = cmd.get(setting)
        cmd.group(self.group)
        for setting, value in settings.items():
            cmd.set(setting, value)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # delete temp directory
        try:
            self._tempdir.cleanup()
        except Exception:
            pass
        # delete our group and restore settings
        if cmd is not None:
            try:
                cmd.delete(self.group)
            except:
                pass
            try:
                for setting, value in self._old_settings.items():
                    cmd.set(setting, value)
            except:
                pass
        return False # never swallow exceptions

    def path(self, filename: str) -> Path:
        return self.tempdir / filename

    def save(self, output_file: Path | str | None):
        cmd.sort(self.group)
        if output_file:
            exporting.save_pdb(output_file, self.group, seqres = True)
        else:
            output_path = self.path('output.pdb')
            exporting.save_pdb(output_path, self.group, seqres = True)
            with open(output_path) as file:
                print(file.read(), end = '')

    @staticmethod
    def check_parens(sel_text: str) -> bool:
        balance = 0
        for char in sel_text:
            if char == '(':
                balance += 1
            elif char == ')':
                balance += -1
                if balance < 0:
                    return False
        return balance == 0

    def __getattr__(self, name):
        attr = getattr(cmd, name)
        if not callable(attr):
            return attr
        if attr.__name__ not in [ "set", "fetch", "load", "remove", "alter", "iterate", "sort", "delete" ]:
            raise ValueError(f"Unsupported command to pymol: {attr}")
        def wrapper(*args, **kwargs):
            if attr.__name__ in [ "fetch", "load" ]:
                if len(args) != 2:
                    raise ValueError(f"Implicit object names for {attr.__name__} are not allowed")
                args = (args[0], f"{self.group}.{args[1]}")
            elif attr.__name__ in [ "set" ]:
                if len(args) < 2:
                    raise ValueError(f"{attr.__name__} needs at least two arguments")
                elif len(args) == 3:
                    if not self.check_parens(args[2]):
                        raise ValueError(f"Invalid selection {args[2]}")
                    args = (args[0], args[1], f"{self.group} AND ({args[2]})")
                else:
                    args = (args[0], args[1], self.group)
            elif attr.__name__ in [ "remove", "alter", "iter", "sort", "delete" ]:
                if not self.check_parens(args[0]):
                    raise ValueError(f"Invalid selection {args[0]}")
                args = (f"{self.group} AND ({args[0]})",) + args[1:] if args else (self.group,)
            return attr(*args, **kwargs)
        return wrapper
