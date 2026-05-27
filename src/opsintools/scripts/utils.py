from Bio import SeqIO, PDB
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
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
    bfactors = get_ca_bfactors(structure)

    atom_record = SeqIO.read(pdb_file, "pdb-atom")
    try:
        seqres_record = SeqIO.read(pdb_file, "pdb-seqres")
    except ValueError: # No SEQRES
        if first_pos > 1:
            raise ValueError(f"First residue in ATOM records is >1 ({first_pos}) but no SEQRES found in {pdb_file}")
        record = atom_record
        start, stop = 0, len(record.seq)
    else:
        record = seqres_record
        regex = str(atom_record.seq).replace('X', '.')
        search = re.search(regex, str(record.seq))
        if not search:
            raise ValueError(f"SEQRES does not match ATOM records in {pdb_file}")
        start, stop = search.span()
        if start != first_pos - 1:
            # A rare case when SEQRES was trimmed in the original PDB file
            missing = first_pos - 1
            record.seq = 'X' * missing + record.seq
            start += missing
            stop += missing

    record.id = record.description = seq_id
    record.letter_annotations["b_factors"] = [ bfactors.get(i+1, "") for i in range(len(record.seq)) ]
    feature = SeqFeature(FeatureLocation(start, stop), type = "trimmed")
    record.features.append(feature)
    return record

def get_coords(record, target_type = "trimmed"):
    for feature in record.features:
        if feature.type == target_type:
            return feature.location.start, feature.location.end
    raise ValueError(f"Record {record.id} does not have any {type} feature")

def get_ca_bfactors(structure, model = 0, chain = 'A'):
    ca_bfactors = {}
    for res in structure[model][chain]:
        het, pos, ins_code = res.id
        if het == ' ' and 'CA' in res:
            ca_bfactors[pos] = res['CA'].get_bfactor()
    return ca_bfactors

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
