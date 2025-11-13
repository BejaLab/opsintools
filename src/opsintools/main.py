from glob import glob
from os import path, makedirs
import json, re, logging
from pathlib import Path
from subprocess import CalledProcessError
from typing import Final
from opsintools.classes.Tcoffee import Tcoffee
from timeit import default_timer as timer

logger = logging.getLogger(__name__)

# Constants
PAD_N: Final[int]= 30
PAD_C: Final[int] = 30
THREADS: Final[int] = 1
N_REPS: Final[int] = 20
MAX_SEQ_ID: Final[float] = 0.9
METHODS: list[str] = [ 'sap_pair', 'mustang_pair', 't_coffee_msa', 'probcons_msa' ]

class OpsinToolsLoggerFormatter(logging.Formatter):
    """Logger formatter class for CLIs
    """
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    fmt = "%(asctime)s - %(levelname)s - %(message)s"
    FORMATS = {
        logging.DEBUG: grey + fmt + reset,
        logging.INFO: grey + fmt + reset,
        logging.WARNING: yellow + fmt + reset,
        logging.ERROR: red + fmt + reset,
        logging.CRITICAL: bold_red + fmt + reset
    }
    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt = "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)

def read_database(data_dir: str, only_exptl: bool) -> tuple[dict, dict]:
    """Checks the opsinmap3d database

    :param str data_dir: data directory
    :returns: a tuple of ref (parsed reference data), rep_dict (dict of { rep: its pdb path }, pdb_list (list of reps to be included always)
    """
    if not path.exists(data_dir):
        raise FileExistsError(f"The data directory {data_dir} does not exist")
    ref_json: str = path.join(data_dir, 'ref.json')
    if not path.isfile(ref_json):
        raise FileNotFoundError(f"File {ref_json} does not exist")
    with open(ref_json) as file:
        ref: dict = json.load(file)

    rep_dir: str = path.join(data_dir, 'reps')
    ref['filename'] = path.join(rep_dir, ref['id'] + '.pdb')
    if not path.isfile(ref['filename']):
        raise FileNotFoundError(f"{ref['filename']} file not found in the data directory")

    exptl_txt: str = path.join(data_dir, 'exptl.txt')
    exptl_dict: dict = {}
    with open(exptl_txt) as file:
        for line in file:
            rep_name = line.rstrip()
            exptl_dict[rep_name] = True

    rep_dict: dict = {}
    for rep_pdb in glob(path.join(rep_dir, '*.pdb')):
        rep_name = Path(rep_pdb).stem
        if rep_name in exptl_dict or not only_exptl:
            rep_dict[rep_name] = { 'filename': rep_pdb, 'pdb': exptl_dict.pop(rep_name, False), 'ref': False }
            rep_dict[rep_name]['ref'] = False

    rep_dict[ref['id']]['ref'] = True
    rep_dict[ref['id']]['pdb'] = True

    if sum(exptl_dict.values()) > 0:
        raise FileNotFoundError("Some of the PDBs listed in the data exptl.txt file are missing from the reps/ subdirectory")

    return ref, rep_dict

def create_output_dir(dir_name, force = False):
    if path.exists(dir_name):
        if not path.isdir(dir_name):
            raise FileExistsError(f"The path {dir_name} exists and is not a directory")
        if not force:
            raise FileExistsError(f"Directory {dir_name} exists, not overriding")
    else:
        makedirs(dir_name)

def opsinalign3d(
        query_pdbs: list[str],
        output_dir: str,
        methods: list = METHODS,
        threads: int = THREADS,
        force: bool = False) -> Tcoffee | None:
    """Runs the opsinalign3d workflow

    :param query_pdbs: input PDB files
    :param output_dir: output directory path
    :param methods: methods to be used for the structural alignment
    :param threads: number of threads
    :param force: whether to overwrite data in the output directory if exists
    """
    from opsintools.scripts.t_coffee import t_coffee, check_t_coffee_methods
    from opsintools.scripts.tm_pos import tm_pos
    from multiprocessing import Pool

    check_t_coffee_methods(methods)

    create_output_dir(output_dir, force)
    t_coffee_aln: str = path.join(output_dir, 't_coffee.aln')
    pdb_dict: dict = {}
    for query_pdb in query_pdbs:
        pdb_stem: str = Path(query_pdb).stem
        if pdb_stem in pdb_dict:
            raise ValueError(f"Duplicate query name: {pdb_stem}")
        pdb_dict[pdb_stem] = query_pdb

    logger.info("Doing the structural alignment")
    t_coffee(pdb_dict, t_coffee_aln, methods, threads)
    output: Tcoffee | None = Tcoffee(t_coffee_aln)
    logger.info("Finished")
    return output

def opsinpdb(
    file_or_accession: str,
    output_file: str,
    is_file: bool,
    chains: list[str] = [],
    non_cov_lig: list[str] = [],
    dont_remove_w: bool = False,
    dont_remove_h: bool = False,
    dont_remove_alt: bool = False,
    no_remap: bool = False) -> None:
    """Read a structure file or fetch from PDB, cleaunp and save as PDB file

    :param accession: PDB accession
    :param output_file: output file path
    :no_remap: whether to remap LYR atoms
    """
    from pathlib import Path
    import importlib.resources

    if not file_or_accession:
        raise ValueError("No file or accession provided")
    if is_file and not Path(file_or_accession).is_file():
        raise ValueError(f"File {file_or_accession} does not exist")
    
    from opsintools.scripts.pdb import process_pdb

    atom_map_file = None if no_remap else importlib.resources.files('opsintools.resources').joinpath('atom_map.json')

    process_pdb(
        file_or_accession = file_or_accession,
        is_file = is_file,
        output_file = output_file,
        chains = chains,
        atom_map_file = atom_map_file,
        non_cov_lig = non_cov_lig,
        dont_remove_w = dont_remove_w,
        dont_remove_h = dont_remove_h,
        dont_remove_alt = dont_remove_alt,
        no_remap = no_remap)

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
    from opsintools.scripts.us_align import us_align
    from opsintools.scripts.score_alignments import score_alignments
    from opsintools.scripts.t_coffee import t_coffee, check_t_coffee_methods
    from opsintools.scripts.tm_pos import tm_pos
    from multiprocessing import Pool

    checking_start = timer()
    logger.info("Checking the input")

    check_t_coffee_methods(methods)

    aln_dir: str = path.join(output_dir, "alignments")

    if not path.isfile(query_pdb):
        raise FileNotFoundError(f"Input file {query_pdb} not found")

    create_output_dir(output_dir, force)
    create_output_dir(aln_dir, force)

    ref: dict
    rep_dict: dict
    pdb_list: list
    ref, rep_dict = read_database(data_dir, only_exptl)

    aln_to_ref: str = path.join(output_dir, 'aln_to_ref.txt')
    trimmed_pdb: str = path.join(output_dir, 'trimmed.pdb')
    t_coffee_aln: str = path.join(output_dir, 't_coffee.aln')
    json_output: str = path.join(output_dir, 'opsinmap.json')

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
        rep_alns: dict = { rep_id: path.join(output_dir, "alignments", rep_id + '.txt') for rep_id in rep_dict }
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
    output: dict = tm_pos(t_coffee_aln, trimmed_pdb, ref['filename'], 'query', ref)
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

def opsinalign3d_cli() -> None:
    from argparse import ArgumentParser
    from sys import argv

    parser = ArgumentParser(description = 'Opsin alignment with t-coffee.')

    main_group = parser.add_argument_group('Arguments')

    main_group.add_argument('-i', nargs='+', metavar = 'INPUT', required = True, help = 'input PDB structures (only chain A will be used)')
    main_group.add_argument('-o', metavar = 'OUTPUT', required = True, help = 'output directory')
    main_group.add_argument('-f', action = 'store_true', help = 'whether to overwrite files in the output directory if it exists')
    main_group.add_argument('-t', metavar = 'THREADS', type = int, default = THREADS, help = f'number of threads to use (default: {THREADS})')

    adv_group = parser.add_argument_group('Advanced')

    adv_group.add_argument('--methods', default = ','.join(METHODS), help = f't-coffee methods to use (default: {",".join(METHODS)})')

    args = parser.parse_args()

    methods = args.methods.split(',')
    run_with_logger(opsinalign3d,
        query_pdbs = args.i,
        output_dir = args.o,
        methods = methods,
        threads = args.t,
        force = args.f,
    )

def mustang_msa_cli():
    from argparse import ArgumentParser
    from sys import argv
    from opsintools.scripts.t_coffee import run_mustang_msa

    parser = ArgumentParser(description = 'A wrapper for mustang MSA to be used with t_coffee.')

    parser.add_argument('-i', metavar = 'INPUT', required = True, help = 'concatenated PDB structures')
    parser.add_argument('-o', metavar = 'OUTPUT', required = True, help = 'output file')

    args = parser.parse_args()

    run_mustang_msa(args.i, args.o)

def mtm_align_msa_cli():
    from argparse import ArgumentParser
    from sys import argv
    from opsintools.scripts.t_coffee import run_mtm_align_msa

    parser = ArgumentParser(description = 'A wrapper for mTM-align MSA to be used with t_coffee.')

    parser.add_argument('-i', metavar = 'INPUT', required = True, help = 'concatenated PDB structures')
    parser.add_argument('-o', metavar = 'OUTPUT', required = True, help = 'output file')

    args = parser.parse_args()

    run_mtm_align_msa(args.i, args.o)

def opsinpdb_cli():
    from argparse import ArgumentParser
    from sys import argv

    parser = ArgumentParser(description = 'Read a structure or fetch from PDB, cleanup and save as a PDB file.')

    main_group = parser.add_argument_group('Arguments')

    main_group.add_argument('-i', metavar = 'FILE', help = 'Input file')
    main_group.add_argument('-a', metavar = 'FILE', help = 'PDB accession (optionally including chain)')

    main_group.add_argument('-o', metavar = 'OUTPUT', required = False, help = 'output path (default: stdout)')
    main_group.add_argument('-c', metavar = 'CHAINS', required = False, help = 'comma-separated list of chains (default: all chains)')
    main_group.add_argument('-L', action = 'store_true', help = 'do not remap LYR atoms (default: do remap)')
    main_group.add_argument('-H', action = 'store_true', help = 'do not remove hydrogens (default: remove)')
    main_group.add_argument('-W', action = 'store_true', help = 'do not remove water molecules (default: remove)')
    main_group.add_argument('-A', action = 'store_true', help = 'do not remove alternative conformations (default: remove)')
    main_group.add_argument('--ligands', metavar = 'LIGANDS', help = 'comma-separated list of non-covalent ligands to retain (default: remove all)')

    args = parser.parse_args()

    chains = [] if not args.c else args.c.split(',')
    ligands = [] if not args.ligands else args.ligands.split(',')
    if bool(args.i) == bool(args.a):
        raise ValueError("Must specify either -a or -i")

    if args.i:
        is_file = True
        file_or_accession = args.i
    else:
        is_file = False
        file_or_accession = args.a
    output_file = args.o if args.o else '/dev/stdout'


    run_with_logger(opsinpdb,
        file_or_accession = file_or_accession,
        output_file = args.o,
        is_file = is_file,
        chains = chains,
        no_remap = args.L,
        non_cov_lig = ligands,
        dont_remove_w = args.W,
        dont_remove_h = args.H,
        dont_remove_alt = args.A
    )

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

def run_with_logger(func, **args):

    # Define logger
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(OpsinToolsLoggerFormatter())
    logger.addHandler(ch)

    try:
        func(**args)
    except CalledProcessError as e:
        if e.stderr:
            logger.fatal(f"Got return code {e.returncode}: {e.stderr.decode()}, check the log files for more details")
        else:
            logger.fatal(f"Got return code {e.returncode}, check the log files for more details")
    except FileNotFoundError as e:
        logger.fatal(e)
    except FileExistsError as e:
        logger.fatal(e)
    except ValueError as e:
        logger.fatal(f"Unexpected value: {e}")
    except AssertionError as e:
        logger.fatal(f"Assumption violated: {e}")
    #except Exception as e:
    #    logger.fatal(f"An error occurred: {e}")
