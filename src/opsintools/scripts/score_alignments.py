import re
from opsintools.classes.USalign import USalign
from pathlib import Path

def score_alignments(alignment_files, n_reps, max_seq_id, must_include):

    # Create an empty dictionary for RMSD valuse per representative
    rmsd_dict = {}

    def get_value(line, prefix):
        if match := re.search(prefix + '=\\s*([0-9.]+)', line):
            return float(match.group(1))

    # Load all the input files
    for aln_file in alignment_files:
        # Read the aln file
        aln = USalign(aln_file)
        if aln.seq_id <= max_seq_id:
            ref_pdb, query_pdb = aln.pdb_files
            ref_name = Path(ref_pdb).stem
            rmsd_dict[ref_name] = aln.rmsd

    # Sort the rmsd dictrionary by rmsd score
    sorted_rmsd = dict(sorted(rmsd_dict.items(), key = lambda x: x[1]))

    output_list = must_include.copy()
    n = len(must_include)
    for ref_name in sorted_rmsd.keys():
        if n >= n_reps: break
        if ref_name not in must_include:
            output_list.append(ref_name)
            n += 1

    return output_list
