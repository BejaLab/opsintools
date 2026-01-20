from Bio import AlignIO, SeqIO, PDB
from opsintools.classes.Tcoffee import Tcoffee
from opsintools.scripts import utils

def get_pos_to_index(pdb_file):
    """Match residue positions (id as defined in ATOM) and 0-based indexes
    in the (potentially trimmed) sequence
    """
    record, start, stop = utils.get_pdb_record(pdb_file)
    pos_to_index = {}
    for index_in_trimmed, index_in_original in enumerate(range(start, stop)):
        pos_to_index[index_in_original + 1] = index_in_trimmed
    return pos_to_index, record.seq[:start], record.seq[start:stop], record.seq[stop:]

def pad_sequence(left, middle, right, pad_left, pad_right):
    return str(left).lower().rjust(pad_left, '.') + middle + str(right).ljust(pad_right, '.')

def aln_mapping(aln_file, query_pdb_file, ref_pdb_file, query, ref_data):
    ref = ref_data['id']
    tms = ref_data['tms']

    t_coffee = Tcoffee(aln_file)
    if t_coffee.aln_score < 0:
        raise ValueError("Something went wrong: check the output of t_coffee")

    query_pos_to_index, query_left, query_mid, query_right = get_pos_to_index(query_pdb_file)
    ref_pos_to_index,   ref_left,   ref_mid,   ref_right   = get_pos_to_index(ref_pdb_file)

    query_index_to_pos = { value: key for key, value in query_pos_to_index.items() }
    ref_index_to_pos   = { value: key for key, value in ref_pos_to_index.items() }

    lys_index = ref_pos_to_index[ref_data['lysine']]

    # Load the positions of the protein which are in the membrane
    ref_tms = {}
    for tm_label, (start, end) in tms.items():
        start_ix = ref_pos_to_index[start]
        end_ix = ref_pos_to_index[end]
        for i in range(start_ix, end_ix + 1):
            ref_tms[i] = tm_label

    ref_aln      = t_coffee.alignment[ref]
    query_aln    = t_coffee.alignment[query]
    ref_scores   = t_coffee.res_scores[ref]
    query_scores = t_coffee.res_scores[query]

    # Create a dictonary with a map of the aligned sequences with the correct residue numeration
    aln = {}
    ref_ix = query_ix = -1
    ref_trimmed = []
    query_trimmed = []
    query_score_trimmed = []
    ref_score_trimmed = []
    for aln_ix, (ref_res, query_res, ref_score, query_score) in enumerate(zip(ref_aln, query_aln, ref_scores, query_scores)):
        query_not_gap = query_res != '-'
        ref_not_gap = ref_res != '-'
        if query_not_gap:
            query_ix += 1
        if ref_not_gap:
            ref_ix += 1
        if query_not_gap and ref_not_gap:
            tm = ref_tms[ref_ix] if ref_ix in ref_tms else '-'
            aln[ref_ix] = {
                'ref_pos': ref_index_to_pos[ref_ix],
                'query_pos': query_index_to_pos[query_ix],
                'ref_res': ref_res,
                'query_res': query_res,
                'ref_score': int(ref_score),
                'query_score': int(query_score),
                'TM': tm
            }
        if query_not_gap or ref_not_gap:
            ref_trimmed.append(ref_res)
            query_trimmed.append(query_res)
            query_score_trimmed.append(query_score)
            ref_score_trimmed.append(ref_score)

    warnings = []

    # A couple of checks for the lysine position
    if lys_index in aln:
        ref_res, query_res = aln[lys_index]['ref_res'], aln[lys_index]['query_res']
        if ref_res != 'K':
            raise ValueError(f"Reference lysine poisition is {ref_res}")
        if query_res != 'K':
            warnings.append(f"Lysine position is occupied by {query_res}")
    else:
        warnings.append(f"Lysine position not identified in the query")

    pad_left = max(len(query_left), len(ref_left))
    pad_right = max(len(query_right), len(ref_right))

    query_seq = pad_sequence(query_left, ''.join(query_trimmed), query_right, pad_left, pad_right)
    ref_seq   = pad_sequence(ref_left,   ''.join(ref_trimmed),   ref_right,   pad_left, pad_right)
    query_score_seq = pad_sequence('', ''.join(query_score_trimmed), '', pad_left, pad_right)
    ref_score_seq   = pad_sequence('', ''.join(ref_score_trimmed),   '', pad_left, pad_right)

    # Add a description to the data that will be the json file
    out_data = {
        "map": list(aln.values()),
        "alignment": { "query": query_seq, "ref": ref_seq, "query_score": query_score_seq, "ref_score": ref_score_seq }
    }
    if warnings:
        out_data['warnings'] = warnings

    return out_data
