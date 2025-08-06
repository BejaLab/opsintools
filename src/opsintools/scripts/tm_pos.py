from Bio import AlignIO, SeqIO, PDB
from opsintools.classes.Tcoffee import Tcoffee

def tm_pos(aln_file, pdb_file, query, ref_data):
    ref = ref_data['id']
    tms = ref_data['tms']

    t_coffee = Tcoffee(aln_file)
    assert t_coffee.aln_score > -1, "Something went wrong: check the output of t_coffee"

    parser = PDB.PDBParser(QUIET = True)
    structure = parser.get_structure("protein", pdb_file)

    query_index_to_pos = {}
    for index, res in enumerate(structure[0]['A']):
        pos = res.id[1]
        query_index_to_pos[index] = pos

    ref_pos_to_index = {}
    ref_index_to_pos = {}
    for index, pos in enumerate(ref_data['res_pos']):
        ref_pos_to_index[pos] = index
        ref_index_to_pos[index] = pos

    lys_index = ref_pos_to_index[ref_data['lysine']]

    # Load the positions of the protein which are in the membrane
    membrane_indexes = []
    for tm_label, (start, end) in tms.items():
        start_ix = ref_pos_to_index[start]
        end_ix = ref_pos_to_index[end]
        for i in range(start_ix, end_ix + 1):
            membrane_indexes.append((i, tm_label))

    ref_aln      = t_coffee.alignment[ref]
    query_aln    = t_coffee.alignment[query]
    ref_scores   = t_coffee.res_scores[ref]
    query_scores = t_coffee.res_scores[query]

    # Create a dictonary with a map of the aligned sequences with the correct residue numeration
    aln = {}
    scores = {}
    ref_ix = query_ix = -1
    for aln_ix, (ref_res, query_res, ref_score, query_score) in enumerate(zip(ref_aln, query_aln, ref_scores, query_scores)):
        if query_res != '-':
            query_ix += 1
        if ref_res != '-':
            ref_ix += 1
            aln[ref_ix] = query_ix, ref_res, query_res, ref_score, query_score

    warnings = []

    query_ix_lys, ref_res_lys, query_res_lys, ref_score_lys, query_score_lys = aln[lys_index]

    assert ref_res_lys == 'K', f"Reference lysine poisition is {ref_res_lys}"
    if query_res_lys != 'K':
        warnings.append(f"Lysine position is occupied by {query_res_lys}")

    tm_map = {}
    for ref_ix, tm_label in membrane_indexes:
        ref_pos = ref_index_to_pos[ref_ix]
        query_ix, ref_res, query_res, ref_score, query_score = aln[ref_ix]
        query_pos = query_index_to_pos[query_ix] if query_ix >= 0 else -1
        score = int(query_score) if query_score != "-" else -1
        tm_map[ref_pos] = { "ref_residue": ref_res, "query_residue": query_res, "helix": tm_label, "query_pos": query_pos, "query_score": score }

    # Add a description to the data that will be the json file
    out_data = {
        "description": f"Mapping for membrane-embedded residues. Key is position number: (residue in {ref}, residue in {query} and the number of the TM helix",
        "alignment_map": tm_map
    }
    if warnings:
        out_data['warnings'] = warnings

    return out_data
