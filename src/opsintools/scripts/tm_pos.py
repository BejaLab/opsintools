from Bio import AlignIO, SeqIO, PDB

def tm_pos(aln_file, pdb_file, query, ref_data):
    ref = ref_data['id']
    tms = ref_data['tms']

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

    # Load the alignment file and filter out all irrelevant sequences
    alignment = SeqIO.to_dict(AlignIO.read(aln_file, "clustal"))

    query_aln = alignment[query]
    ref_aln = alignment[ref]

    # Create a dictonary with a map of the aligned sequences with the correct residue numeration
    aln = {}
    ref_ix = query_ix = -1
    for aln_ix, (ref_res, query_res) in enumerate(zip(ref_aln, query_aln)):
        if query_res != '-':
            query_ix += 1
        if ref_res != '-':
            ref_ix += 1
            aln[ref_ix] = query_ix, ref_res, query_res

    warnings = []

    query_ix_lys, ref_res_lys, query_res_lys = aln[lys_index]

    assert ref_res_lys == 'K', f"Reference lysine poisition is {ref_res_lys}"
    if query_res_lys != 'K':
        warnings.append(f"Lysine position is occupied {query_res_lys}")

    tm_map = {}
    for ref_ix, tm_label in membrane_indexes:
        ref_pos = ref_index_to_pos[ref_ix]
        query_ix, ref_res, query_res = aln[ref_ix]
        query_pos = query_index_to_pos[query_ix]
        tm_map[ref_pos] = { "ref_residue": ref_res, "query_residue": query_res, "helix": tm_label, "query_pos": query_pos }

    # Add a description to the data that will be the json file
    out_data = {
        "description": f"Mapping for rhodopsin residues inside the transmembrane. Key is position number: (residue in {ref}, residue in {query} and the number of the TM helix",
        "alignment_map": tm_map
    }
    if warnings:
        out_data['warnings'] = warnings

    return out_data
