import subprocess
import shutil
from os import path

def us_align(ref_pdb, query_pdb, output):
    """Run US-align for each cluster representative on the query"""
    if shutil.which("USalign") is None:
        raise Exception("USalign not found in PATH")
    cmd = [ "USalign", ref_pdb, query_pdb, "-do" ]
    with open(output, 'w') as out_file, open(f"{output}.log", 'w') as err_file:
        process = subprocess.run(cmd, stdout = out_file, stderr = None)
        process.check_returncode()  # raise error if the command failed
    if path.getsize(output) == 0:
        raise ValueError(f"Output file {output} is empty.")
