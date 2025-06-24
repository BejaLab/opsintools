import subprocess
import shutil
from os import path

def us_align(ref_pdb, query_pdb, output):
    """Run US-align for each cluster representative on the query"""
    if shutil.which("USalign") is None:
        raise FileNotFoundError("USalign not found in PATH")
    cmd = [ "USalign", ref_pdb, query_pdb, "-do" ]
    with open(output, 'w') as out_file:
        result = subprocess.run(cmd, stdout = out_file, stderr = subprocess.PIPE)
    if path.getsize(output) < 100:
        raise ValueError(f"Alignment file {output} is too small, US-align returned: {result.stderr.decode()}")
