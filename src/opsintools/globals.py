import re, logging
from pathlib import Path
from typing import Final

logger = logging.getLogger(__name__)

PAD_N: Final[int]= 30
PAD_C: Final[int] = 30

THREADS: Final[int] = 1
N_REPS: Final[int] = 20
MAX_SEQ_ID: Final[float] = 0.9
METHODS: list[str] = [ 'sap_pair', 'mustang_pair', 't_coffee_msa', 'probcons_msa' ]

HMM_MIN_GAP: Final[int] = -20
HMM_MAX_GAP: Final[int] = 100
HMM_MIN_SCORE: Final[int] = 15

BUILD_REPO = 'BejaLab/opsintools-build'

class OpsinToolsLoggerFormatter(logging.Formatter):
    """Logger formatter class for CLIs
    """
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    fmt = "%(asctime)s - %(levelname)s - %(message)s"
    FORMATS = {
        logging.DEBUG: grey + fmt + reset,
        logging.INFO: grey + fmt + reset,
        logging.WARNING: yellow + fmt + reset,
        logging.ERROR: red + fmt + reset,
        logging.CRITICAL: bold_red + fmt + reset
    }
    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt = "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)

def read_database(data_dir: str, only_exptl: bool = False) -> tuple[dict, dict]:
    """Checks the opsintools database

    :param str data_dir: data directory
    :returns: a tuple of ref (parsed reference data), rep_dict (dict of { rep: its pdb path }, pdb_list (list of reps to be included always)
    """
    import json
    data_path = Path(data_dir)
    if not data_path.exists():
        raise FileExistsError(f"The data directory {data_dir} does not exist")
    if not data_path.is_dir():
        raise FileExistsError(f"{data_dir} is not a directory")
    ref_json = data_path / 'ref.json'
    if not ref_json.is_file():
        raise FileNotFoundError(f"File {ref_json} does not exist")
    with open(ref_json) as file:
        ref = json.load(file)

    rep_dir = data_path / 'reps'
    ref['filename'] = rep_dir / (ref['id'] + '.pdb')
    if not ref['filename'].is_file():
        raise FileNotFoundError(f"{ref['filename']} file not found in the data directory")

    exptl_txt = data_path / 'exptl.txt'
    exptl_dict = {}
    with open(exptl_txt) as file:
        for line in file:
            rep_name = line.rstrip()
            exptl_dict[rep_name] = True

    rep_dict = {}
    for rep_pdb in rep_dir.glob('*.pdb'):
        rep_name = rep_pdb.stem
        if rep_name in exptl_dict or not only_exptl:
            rep_dict[rep_name] = { 'filename': rep_pdb, 'pdb': exptl_dict.pop(rep_name, False), 'ref': False }
            rep_dict[rep_name]['ref'] = False

    rep_dict[ref['id']]['ref'] = True
    rep_dict[ref['id']]['pdb'] = True

    if sum(exptl_dict.values()) > 0:
        raise FileNotFoundError("Some of the PDBs listed in the data exptl.txt file are missing from the reps/ subdirectory")

    profile_path = data_path / 'profile.hmm'
    if not profile_path.is_file():
        raise FileNotFoundError("profile.hmm not found in the data directory")

    profile_fasta_path = data_path / 'profile.fasta'
    if not profile_fasta_path.is_file():
        raise FileNotFoundError("profile.fasta not found in the data directory")

    profile_ref_path = data_path / 'profile_ref.txt'
    if not profile_ref_path.is_file():
        raise FileNotFoundError("profile_ref.txt not found in the data directory")

    return { 'ref': ref, 'rep_dict': rep_dict, 'profile': profile_path, 'profile_fasta': profile_fasta_path, 'profile_ref': profile_ref_path }

def create_output_dir(dir_name, force = False):
    dir_path = Path(dir_name)
    if dir_path.exists():
        if not dir_path.is_dir():
            raise FileExistsError(f"The path {dir_name} exists and is not a directory")
        if not force:
            raise FileExistsError(f"Directory {dir_name} exists, not overriding")
    else:
        dir_path.mkdir(parents = True)

def run_with_logger(func, **args):

    from subprocess import CalledProcessError
    import sys
    # Define logger
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(OpsinToolsLoggerFormatter())
    logger.addHandler(ch)

    def _fatal(msg):
        logger.fatal(msg)
        sys.exit(1)

    try:
        func(**args)
    except CalledProcessError as e:
        if e.stderr:
            _fatal(f"Got return code {e.returncode}: {e.stderr.decode()}, check the log files for more details")
        else:
            _fatal(f"Got return code {e.returncode}, check the log files for more details")
    except (FileNotFoundError, FileExistsError) as e:
        _fatal(e)
    except ValueError as e:
        _fatal(f"Unexpected value: {e}")
    except AssertionError as e:
        _fatal(f"Assumption violated: {e}")
    except Exception as e:
        _fatal(f"Got exception: {e}")
