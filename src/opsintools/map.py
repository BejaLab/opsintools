from .globals import run_with_logger, logger, create_output_dir, read_database
from .globals import N_REPS, METHODS, THREADS, PAD_N, PAD_C, MAX_SEQ_ID, HMM_MAX_GAP, HMM_MIN_SCORE
from opsintools.classes.Hmmer import Hmmer, HmmerContainer
from pathlib import Path
from timeit import default_timer as timer

def opsinmap3d(
        query_pdb: str,
        output_dir: str,
        data_dir: str,
        n_templates: int = N_REPS,
        methods: list = METHODS,
        threads: int = THREADS,
        force: bool = False,
        pad_n: int = PAD_N,
        pad_c: int = PAD_C,
        max_seq_id: float = MAX_SEQ_ID,
        only_exptl: bool = False,
        prefer_exptl: bool = False) -> dict | None:
    """Runs the opsinmap3d workflow

    :param query_pdb: input PDB file
    :param output_dir: output directory path
    :param data_dir: database directory
    :param n_templates: number of termplates to be included in the structural alignment
    :param methods: methods to be used for the structural alignment
    :param threads: number of threads
    :param force: whether to overwrite data in the output directory if exists
    :param pad_n: pad this number of residues to N-terminus when trimming the query
    :param pad_c: pad this number of residues to C-terminus when trimming the query
    :param max_seq_id: exclude templates with sequence identity higher than this
    :param only_exptl: only use experimental structures
    :param prefer_exptl: always prefer experimental structure over predictions
    """
    from opsintools.scripts.prot_trim_filter import prot_trim_filter
    from opsintools.scripts.t_coffee import t_coffee, check_t_coffee_methods
    from opsintools.scripts.us_align import us_align
    from opsintools.scripts.score_alignments import score_alignments
    from opsintools.scripts.aln_mapping import aln_mapping
    from pathlib import Path
    import json
    from multiprocessing import Pool

    checking_start = timer()
    logger.info("Checking the input")

    check_t_coffee_methods(methods)

    output_path = Path(output_dir)
    aln_dir = output_path / "alignments"

    if not Path(query_pdb).is_file():
        raise FileNotFoundError(f"Input file {query_pdb} not found")

    create_output_dir(output_dir, force)
    create_output_dir(aln_dir, force)

    ref: dict
    rep_dict: dict
    pdb_list: list
    database = read_database(data_dir, only_exptl)
    rep_dict = database['rep_dict']
    ref = database['ref']

    aln_to_ref = output_path / 'aln_to_ref.txt'
    trimmed_pdb = output_path / 'trimmed.pdb'
    t_coffee_aln = output_path / 't_coffee.aln'
    json_output = output_path / 'opsinmap.json'

    alignment_to_ref_start = timer()
    logger.info("Aligning the query to the reference")
    us_align(ref['filename'], query_pdb, aln_to_ref)

    trimming_start = timer()
    logger.info("Trimming the query")
    prot_trim_filter(aln_to_ref, query_pdb, ref['filename'], trimmed_pdb, pad_n = pad_n, pad_c = pad_c)

    templates_start = timer()
    chosen_templates: list = []
    if not n_templates or n_templates == 1:
        chosen_templates = [ ref['id'] ]
    elif n_templates >= len(rep_dict):
        chosen_templates = list(rep_dict.keys())
    else:
        logger.info("Picking representatives for structural alignment")
        rep_alns: dict = { rep_id: output_path / "alignments" / (rep_id + '.txt') for rep_id in rep_dict }
        preferred: list = [ rep_id for rep_id, rep in rep_dict.items() if rep['pdb'] ] if prefer_exptl else [ ref['id'] ]
        with Pool(threads) as pool:
            results: list = pool.starmap(us_align, [ (rep_dict[rep_id]['filename'], trimmed_pdb, rep_alns[rep_id]) for rep_id in rep_dict ])
        chosen_templates = score_alignments(rep_alns.values(), max_seq_id, preferred)
        chosen_templates = chosen_templates[:n_templates]

    chosen_pdbs: dict = { 'query': trimmed_pdb }
    for rep in chosen_templates:
        chosen_pdbs[rep] = rep_dict[rep]['filename']

    alignment_start = timer()
    logger.info("Doing the structural alignment")
    t_coffee(chosen_pdbs, t_coffee_aln, methods = methods, threads = threads)

    finishing_time = timer()
    output: dict = aln_mapping(t_coffee_aln, trimmed_pdb, ref['filename'], 'query', ref)
    output['params'] = {
        'methods': methods,
        'data_dir': str(Path(data_dir).resolve()),
        'n_templates': n_templates,
        'pad_n': pad_n,
        'pad_c': pad_c,
        'max_seq_id': max_seq_id,
        'only_exptl': only_exptl,
        'prefer_exptl': prefer_exptl
    }
    output['templates'] = chosen_templates
    output['timer'] = {
        'input checking': alignment_to_ref_start - checking_start,
        'alignment to reference': trimming_start - alignment_to_ref_start,
        'query trimming': templates_start - trimming_start,
        'picking templates': alignment_start - templates_start,
        'alignment to templates': finishing_time - alignment_start,
        'finishing': timer() - finishing_time
    }
    with open(json_output, 'w') as file:
        json.dump(output, file, indent = 2)
    logger.info("Finished")
    return output

def opsinmap3d_cli():
    from argparse import ArgumentParser
    from sys import argv

    parser = ArgumentParser(description = 'Opsin homology based on 3D structure.')

    main_group = parser.add_argument_group('Arguments')

    main_group.add_argument('-i', metavar = 'INPUT', required = True, help = 'query PDB structure (only chain A will be used)')
    main_group.add_argument('-d', metavar = 'DATADIR', required = True, help = 'opsin data directory')
    main_group.add_argument('-o', metavar = 'OUTPUT', required = True, help = 'output directory')
    main_group.add_argument('-f', action = 'store_true', help = 'whether to overwrite files in the output directory if it exists')
    main_group.add_argument('-n', metavar = 'N_REPS', type = int, default = N_REPS, help = f'maximum number of template structures to use (default: {N_REPS})')
    main_group.add_argument('-t', metavar = 'THREADS', type = int, default = THREADS, help = f'number of threads to use (default: {THREADS})')

    adv_group = parser.add_argument_group('Advanced')

    adv_group.add_argument('--only-exptl', action = 'store_true', help = f'only use experimental structures (default: use predictions as well)')
    adv_group.add_argument('--prefer-exptl', action = 'store_true', help = f'prefer experimental structure over predictions (default: choose by similarity only)')
    adv_group.add_argument('--pad-n', default = PAD_N, type = int, help = f'N-terminal padding for query trimming (default: {PAD_N})')
    adv_group.add_argument('--pad-c', default = PAD_C, type = int, help = f'C-terminal padding for query trimming (default: {PAD_C})')
    adv_group.add_argument('--max-seq-id', default = MAX_SEQ_ID, type = int, help = f'maximum sequence identity for predicted templates (default: {MAX_SEQ_ID})')
    adv_group.add_argument('--methods', default = ','.join(METHODS), help = f't-coffee methods to use (default: {",".join(METHODS)})')

    args = parser.parse_args()

    methods = args.methods.split(',')

    run_with_logger(opsinmap3d,
        query_pdb = args.i,
        output_dir = args.o,
        data_dir = args.d,
        n_templates = args.n,
        methods = methods,
        threads = args.t,
        force = args.f,
        pad_n = args.pad_n, pad_c = args.pad_c, max_seq_id = args.max_seq_id, only_exptl = args.only_exptl, prefer_exptl = args.prefer_exptl
    )

def opsinmaphmm(
        query_fasta: str,
        output_dir: str,
        data_dirs: list,
        pad_n: int = PAD_N,
        pad_c: int = PAD_C,
        max_gap: int = HMM_MAX_GAP,
        min_score: float = HMM_MIN_SCORE,
        threads: int = THREADS,
        force: bool = False) -> dict | None:
    """Runs the opsinhmm workflow

    :param query_fasta: input (multi)fasta file
    :param output_dir: output directory path
    :param data_dirs: list of database directories
    :param pad_n: pad this number of residues to N-terminus when trimming the queries
    :param pad_c: pad this number of residues to C-terminus when trimming the queries
    :param threads: number of threads
    :param force: whether to overwrite data in the output directory if exists
    """
    from collections import defaultdict
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import json
    from opsintools.scripts.prot_trim_filter import prot_trim_filter
    from opsintools.scripts.us_align import us_align
    from opsintools.scripts.score_alignments import score_alignments
    from opsintools.scripts.t_coffee import t_coffee, check_t_coffee_methods
    from opsintools.scripts.hmmer import hmmsearch, local_to_global, sync_pwas

    logger.info("Checking the input")

    if not Path(query_fasta).is_file():
        raise FileNotFoundError(f"Input file {query_fasta} not found or not a file")

    databases = {}

    for data_dir in data_dirs:
        database = read_database(data_dir)
        profile_path = database['profile']
        profile_name = profile_path.parent.stem
        if profile_name in databases:
            raise ValueError(f"Profile named '{profile_name}' specified multiple times")
        profile_ref_file = database['profile_ref']
        profile_fasta_file = database['profile_fasta']
        consens_seq = str(SeqIO.read(profile_fasta_file, 'fasta').seq)
        ref = database['ref']
        hmmer = Hmmer(profile_path, consens_seq, search = profile_ref_file)
        match = hmmer.matches[0]
        domain = match['domains'][0]
        database['ref_tms'] = {}
        for tm_label, (start, end) in ref['tms'].items():
            for pos in range(start, end + 1):
                database['ref_tms'][pos] = tm_label
        database['ref_domain'] = local_to_global(domain, consens_seq, ref['seq'])
        database['ref_id'] = ref['id']
        database['profile_cons'] = consens_seq
        database['profile_path'] = profile_path
        databases[profile_name] = database

    output_path = Path(output_dir)
    create_output_dir(output_dir, force)
    json_output = output_path / 'opsinmap.json'
    trimmed_fasta = output_path / 'trimmed.fasta'

    hmmsearch_start = timer()
    logger.info("Doing the searches")

    records = SeqIO.to_dict(SeqIO.parse(query_fasta, "fasta"))
    hmmers = HmmerContainer()
    for profile_name, database in databases.items():
        profile_path = database['profile_path']
        profile_cons = database['profile_cons']
        output_file, log_file = hmmsearch(query_fasta, profile_path, Path(output_dir) / profile_name, no_filters = True, threads = threads)
        hmmers += Hmmer(profile_path, profile_cons, search = output_file)

    hmmers.resolve()
    hmmers.chain_domains(records, max_gap = max_gap)

    logger.info("hmmsearch finished")

    output = []
    trimmed_records = []
    for match in hmmers.matches:
        profile_file = match["profile_file"]
        profile_name = Path(profile_file).parent.stem
        database = databases[profile_name]
        ref_dom = database['ref_domain']
        profile_cons = database['profile_cons']
        profile_len = len(profile_cons)
        seq_name = match["seq_name"]
        record = records[seq_name]
        has_domains = False
        for dom in match['domains']:
            dom_score = sum(dom['score']) if isinstance(dom['score'], list) else dom['score'][0]
            if dom_score >= min_score:
                has_domains = True
                dom = local_to_global(dom, profile_cons, record.seq)
                hmm_seq, ref_seq, ref_pp, query_seq, query_pp = sync_pwas(ref_dom['hmm']['seq'], ref_dom['ali']['seq'], ref_dom['PP'], dom['hmm']['seq'], dom['ali']['seq'], dom['PP'])
                hmm_pos = ref_pos = query_pos = 0
                aln_map = []
                first_hmm_pos = first_query_pos = None
                last_hmm_pos  = last_query_pos = None
                for hmm_res, ref_res, ref_pp_res, query_res, query_pp_res in zip(hmm_seq, ref_seq, ref_pp, query_seq, query_pp):
                    hmm_not_gap = hmm_res != '-'
                    ref_not_gap = ref_res != '-'
                    query_not_gap = query_res != '-'
                    hmm_pos += hmm_not_gap
                    ref_pos += ref_not_gap
                    query_pos += query_not_gap
                    if hmm_not_gap and query_not_gap:
                        if first_hmm_pos is None:
                            first_hmm_pos, first_query_pos = hmm_pos, query_pos
                        else:
                            last_hmm_pos, last_query_pos = hmm_pos, query_pos
                        if ref_not_gap:
                            aln_map.append({
                                "ref_pos": ref_pos, "query_pos": query_pos, "hmm_pos": hmm_pos,
                                "ref_res": ref_res, "query_res": query_res, "hmm_res": hmm_res,
                                "ref_score": Hmmer.decode_prob(ref_pp_res), "query_score": Hmmer.decode_prob(query_pp_res),
                                "TM": database['ref_tms'][ref_pos] if ref_pos in database['ref_tms'] else '-'
                            })
                trim_start = max(0, first_query_pos - first_hmm_pos - pad_n)
                trim_end = min(len(record.seq), last_query_pos + profile_len - last_hmm_pos + pad_c)
                trimmed_name = f"{seq_name}/{trim_start+1}-{trim_end}"
                trimmed_seq = query_seq.replace('-', '')[trim_start:trim_end]
                trimmed_rec = SeqRecord(Seq(trimmed_seq), id = trimmed_name, description = record.description)
                trimmed_records.append(trimmed_rec)
                output.append({
                    "profile": profile_name,
                    "query": seq_name,
                    "trimmed": trimmed_name,
                    "ref": database['ref_id'],
                    "domain": dom['num'],
                    "alignment": {
                        "query": query_seq,
                        "ref": ref_seq,
                        "query_score": query_pp,
                        "ref_score": ref_pp,
                        "hmm_consensus": hmm_seq
                    },
                    "map": aln_map
                })
        if not has_domains:
            logger.warning(f"No domains found in {seq_name}")
    with open(json_output, 'w') as file:
        json.dump(output, file, indent = 2)
    with open(trimmed_fasta, 'w') as file:
        SeqIO.write(trimmed_records, file, "fasta")
    logger.info("Finished") 
    return output

def opsinmaphmm_cli():
    from argparse import ArgumentParser
    from sys import argv

    parser = ArgumentParser(description = 'Opsin alignments to reference HMM profiles.')

    main_group = parser.add_argument_group('Arguments')

    main_group.add_argument('-i', metavar = 'INPUT', required = True, help = 'query sequences in fasta format')
    main_group.add_argument('-d', metavar = 'DATADIR', nargs='+', help = 'opsin data directories')
    main_group.add_argument('-o', metavar = 'OUTPUT', required = True, help = 'output directory')
    main_group.add_argument('--pad-n', default = PAD_N, type = int, help = f'N-terminal padding for query trimming (default: {PAD_N})')
    main_group.add_argument('--pad-c', default = PAD_C, type = int, help = f'C-terminal padding for query trimming (default: {PAD_C})')
    main_group.add_argument('--max-gap', default = HMM_MAX_GAP, type = int, help = f'Maximum gap for chaining (default: {HMM_MAX_GAP})')
    main_group.add_argument('--min-score', default = HMM_MIN_SCORE, type = float, help = f'Minimum score (default: {HMM_MIN_SCORE})')
    main_group.add_argument('-f', action = 'store_true', help = 'whether to overwrite files in the output directory if it exists')
    main_group.add_argument('-t', metavar = 'THREADS', type = int, default = THREADS, help = f'number of threads to use (default: {THREADS})')

    args = parser.parse_args()

    run_with_logger(opsinmaphmm,
        query_fasta = args.i,
        output_dir = args.o,
        data_dirs = args.d,
        pad_n = args.pad_n,
        pad_c = args.pad_c,
        max_gap = args.max_gap,
        min_score = args.min_score,
        threads = args.t,
        force = args.f
    )
