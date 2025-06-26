from Bio import PDB
from Bio.Seq import Seq
from opsintools.classes.USalign import USalign
from opsintools.classes.TrimmedChain import TrimmedChain
from math import inf

def prot_trim_filter(
        aln_file, query_pdb_file, ref_pdb_file, trimmed_pdb_file,
        pad_n = 0, pad_c = 0, max_missing_n = inf, max_missing_c = inf, max_rmsd = inf, min_aln_len = 0, max_aln_len = inf, min_len = 0, ref_lys_pos = None
    ):
    aln = filter_aln_seq(aln_file, max_rmsd, min_len, min_aln_len, max_aln_len, ref_lys_pos)

    success = False

    if aln is not None:
        parser = PDB.PDBParser(QUIET = True)
        ref_structure = parser.get_structure('ref', ref_pdb_file)
        query_structure = parser.get_structure('query', query_pdb_file)

        ref_first, ref_last = get_first_and_last(ref_structure)

        aln_first, *_, aln_last = aln.alignment
        ref_first_aln_pos, query_first_aln_pos, *_ = aln_first
        ref_last_aln_pos,  query_last_aln_pos, *_  = aln_last

        missing_n = ref_first - ref_first_aln_pos
        missing_c = ref_last  - ref_last_aln_pos

        if missing_n <= max_missing_n and missing_c <= max_missing_c:
            # Trim the strcuture and save it to output pdb
            trim_struct(query_structure, trimmed_pdb_file, start = query_first_aln_pos - missing_n - pad_n, end = query_last_aln_pos + missing_c + pad_c)
            success = True

    if not success:
        with open(trimmed_pdb_file, 'w'):
            pass

# Function to check if query needs to be filtered out. Save aligned sequence for both referance and query
def filter_aln_seq(aln_file, max_rmsd, min_len, min_aln_len, max_aln_len, ref_lys_pos):
    aln = USalign(aln_file)
    
    assert aln.alignment, f"Alignment is empty in {aln_file}"

    # If criteria not met return empty
    if aln.rmsd > max_rmsd or len(aln.alignment) < min_aln_len or aln.seq_lens[1] < min_len:
        return None

    # Check if for the referance lysine position there is a lysine in the query 
    lys_reached = False
    for pos_r, pos_q, res_r, res_q, distance in aln.alignment:
        if pos_r == ref_lys_pos:
            assert res_r == 'K', f"Lysine not found at reference position {ref_lys_pos}"
            if res_q != "K":
                return None
            lys_reached = True

    if ref_lys_pos is not None and not lys_reached:
        return None

    return aln

# Function to trim the pdb file according to start and stop positions
def trim_struct(structure, trimmed_pdb, start, end):
    io = PDB.PDBIO()
    io.set_structure(structure)

    trimmed_chain = TrimmedChain(model = 0, chain = 'A', start = start, end = end)
    io.save(trimmed_pdb, trimmed_chain)

def get_first_and_last(structure):
    first, *_, last = structure[0]['A'].get_residues()
    return first.id[1], last.id[1]


