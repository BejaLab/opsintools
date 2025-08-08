from Bio import SeqIO
from Bio.Seq import Seq
import warnings
import re

warnings.filterwarnings('ignore', message = ".*(can't determine PDB ID|Ignoring unrecognized record).*")

def get_pdb_record(pdb_file, seq_id = 'chain_A'):
    """
    Get record corresponding to (potnetially trimmed) sequence
    defined in ATOM records but filling the gaps from SEQRES
    """
    record = SeqIO.read(pdb_file, "pdb-atom")
    record.id = record.description = seq_id
    if 'X' in record.seq:
        try:
            record_seqres = SeqIO.read(pdb_file, "pdb-seqres")
        except ValueError: # No SEQRES
            return record
        regex = str(record.seq).replace('X', '.')
        search = re.search(regex, str(record_seqres.seq))
        if not search:
            raise ValueError(f"SEQRES does not match ATOM records in {pdb_file}")
        record.seq = Seq(search.group())
    return record
