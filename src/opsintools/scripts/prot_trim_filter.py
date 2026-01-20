from Bio import PDB
from Bio.Seq import Seq
from opsintools.classes.USalign import USalign
from opsintools.classes.TrimmedChain import TrimmedChain
from opsintools.scripts import utils
from math import inf

def get_seqres(pdb_file, structure, model = 0, chain = 'A'):
    """Get SEQRES lines from a PDB file or generate
    new ones based on the sequence in the ATOMs
    """
    seqres_lines = []
    with open(pdb_file) as file:
        for line in file:
            if line.startswith('SEQRES'):
                seqres_lines.append(line)
    if not seqres_lines:
        seq = []
        first_pos, last_pos = utils.get_first_and_last(structure)
        for i in range(first_pos - 1):
            seq.append('XAA')
        for pos in range(first_pos, last_pos + 1):
            if pos in structure[model][chain]:
                res = structure[model][chain][pos]
                seq.append(res.resname)
            else:
                seq.append('XAA')
        serial = 1
        res_per_line = 13
        res_count = len(seq)
        for i in range(0, res_count, res_per_line):
            chunk = seq[i:i + res_per_line]
            residues_str = " ".join(chunk)
            line = f"SEQRES {serial:>3} {chain} {res_count:>4}  {residues_str} \n"
            seqres_lines.append(line)
            serial += 1

    return seqres_lines

def prot_trim_filter(
        aln_file, query_pdb_file, ref_pdb_file, trimmed_pdb_file,
        pad_n = 0, pad_c = 0, max_missing_n = inf, max_missing_c = inf, max_rmsd = inf, min_aln_len = 0, max_aln_len = inf, min_len = 0, ref_lys_pos = None
    ):
    """Trim protein based on an alignment
    """
    aln = filter_aln_seq(aln_file, max_rmsd, min_len, min_aln_len, max_aln_len, ref_lys_pos)

    success = False

    if aln is not None:
        parser = PDB.PDBParser(QUIET = True)
        ref_structure = parser.get_structure('ref', ref_pdb_file)
        query_structure = parser.get_structure('query', query_pdb_file)
        seqres = get_seqres(query_pdb_file, query_structure)

        ref_first, ref_last = utils.get_first_and_last(ref_structure)

        aln_first, *_, aln_last = aln.alignment
        ref_first_aln_pos, query_first_aln_pos, *_ = aln_first
        ref_last_aln_pos,  query_last_aln_pos, *_  = aln_last

        missing_n = ref_first - ref_first_aln_pos
        missing_c = ref_last  - ref_last_aln_pos

        if missing_n <= max_missing_n and missing_c <= max_missing_c:
            # Trim the structure and save it to output pdb
            start = query_first_aln_pos - missing_n - pad_n
            end = query_last_aln_pos + missing_c + pad_c
            trim_struct(query_structure, seqres, trimmed_pdb_file, start, end)
            success = True

    if not success:
        with open(trimmed_pdb_file, 'w'):
            pass

def filter_aln_seq(aln_file, max_rmsd, min_len, min_aln_len, max_aln_len, ref_lys_pos):
    """Function to check if query needs to be filtered out
    Save aligned sequence for both referance and query
    """
    aln = USalign(aln_file)
    
    if not aln.alignment:
        raise ValueError(f"Alignment is empty in {aln_file}")

    # If criteria not met return empty
    if aln.rmsd > max_rmsd or len(aln.alignment) < min_aln_len or aln.seq_lens[1] < min_len:
        return None

    # Check if for the referance lysine position there is a lysine in the query 
    lys_reached = False
    for pos_r, pos_q, res_r, res_q, distance in aln.alignment:
        if pos_r == ref_lys_pos:
            if res_r != 'K':
                raise ValueError(f"Lysine not found at reference position {ref_lys_pos}")
            if res_q != "K":
                return None
            lys_reached = True

    if ref_lys_pos is not None and not lys_reached:
        return None

    return aln

def trim_struct(structure, seqres, trimmed_pdb, start, end, model = 0, chain = 'A'):
    """Function to trim the pdb file according to start and stop positions
    """

    io = PDB.PDBIO()
    io.set_structure(structure)

    trimmed_chain = TrimmedChain(model = model, chain = chain, start = start, end = end)
    with open(trimmed_pdb, 'w') as file:
        for line in seqres:
            file.write(line)
        io.save(file, trimmed_chain)
