from Bio import SeqIO, PDB
from Bio.Seq import Seq
import warnings
import re

warnings.filterwarnings('ignore', message = ".*(can't determine PDB ID|Ignoring unrecognized record).*")

def get_pdb_record(pdb_file, seq_id = 'chain_A'):
    """Get record corresponding to (potentially trimmed) sequence
    defined in ATOM records but filling the gaps from SEQRES
    """
    parser = PDB.PDBParser(QUIET = True)
    structure = parser.get_structure("protein", pdb_file)
    if len(structure) != 1 or len(structure[0]) != 1:
        raise ValueError(f"{pdb_file} has multiple models or chains")
    first_pos, last_pos = get_first_and_last(structure)

    record = SeqIO.read(pdb_file, "pdb-atom")
    record.id = record.description = seq_id
    try:
        record_seqres = SeqIO.read(pdb_file, "pdb-seqres")
    except ValueError: # No SEQRES
        record_seqres = None
    if record_seqres is None:
        if first_pos > 1:
            raise ValueError(f"First residue in ATOM records > 1 ({first_pos}) but no SEQRES found in {pdb_file}")
        return record.seq, 0, len(record.seq)
    regex = str(record.seq).replace('X', '.')
    search = re.search(regex, str(record_seqres.seq))
    if not search:
        raise ValueError(f"SEQRES does not match ATOM records in {pdb_file}")
    start, stop = search.span()
    if start != first_pos - 1:
        # A rare case when SEQRES was trimmed in the original PDB file
        missing = first_pos - 1
        record_seqres.seq = 'X' * missing + record_seqres.seq
        start += missing
        stop += missing
    record_seqres.id = record_seqres.description = seq_id
    return record_seqres, start, stop

def get_first_and_last(structure, model = 0, chain = 'A'):
    """Get first and last aminoacid residues in a structure"""
    first = last = None
    for res in structure[model][chain]:
        het, pos, ins_code = res.id
        if not het.strip():
            if first is None:
                first = pos
            last = pos
    return first, last
