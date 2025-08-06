import re
from opsintools.classes.USalign import USalign
from pathlib import Path

def score_alignments(alignment_files, max_seq_id, preferred):

    # Create an empty dictionary for RMSD valuse per representative
    rmsd_dict = {}

    def get_value(line, prefix):
        if match := re.search(prefix + '=\\s*([0-9.]+)', line):
            return float(match.group(1))

    def sorting_fn(item):
        rep, rmsd = item
        return rep not in preferred, rmsd

    # Load all the input files
    for aln_file in alignment_files:
        # Read the aln file
        aln = USalign(aln_file)
        ref_pdb, query_pdb = aln.pdb_files
        ref_name = Path(ref_pdb).stem
        if aln.seq_id <= max_seq_id or ref_name in preferred:
            rmsd_dict[ref_name] = aln.rmsd

    # Sort the reps by the presence in preferred and rmsd score
    rmsd_sorted = sorted(rmsd_dict.items(), key = sorting_fn)
    output_list = [ rep for rep, rmsd in rmsd_sorted ]

    return output_list
