from glob import glob
from os import path, makedirs
import json, re, logging
from pathlib import Path
from subprocess import CalledProcessError
from typing import Final
from opsintools.classes.Tcoffee import Tcoffee

# Constants
PAD_N: Final[int]= 30
PAD_C: Final[int] = 30
THREADS: Final[int] = 1
N_REPS: Final[int] = 20
MAX_SEQ_ID: Final[float] = 0.9
METHODS: list[str] = [ 'sap_pair', 'mustang_pair', 't_coffee_msa', 'probcons_msa' ]

class OpsinToolsLoggerFormatter(logging.Formatter):
    """
    logger formatter class for CLIs
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

def read_database(data_dir: str) -> tuple[dict, dict, list]:
    """
    checks the opsinmap3d database

    :param data_dir: data directory
    :return: tuple of ref (parsed reference data), rep_dict (dict of { rep: its pdb path }, must_list (list of reps to be included always)
    """
    if not path.exists(data_dir):
        raise FileExistsError(f"The data directory {data_dir} does not exist")
    ref_json: str = path.join(data_dir, 'ref.json')
    if not path.isfile(ref_json):
        raise FileNotFoundError(f"File {ref_json} does not exist")
    with open(ref_json) as file:
        ref: dict = json.load(file)
    rep_dir: str = path.join(data_dir, 'reps')
    ref['pdb'] = path.join(rep_dir, ref['id'] + '.pdb')
    rep_pdbs: list[str] = glob(path.join(rep_dir, '*.pdb'))
    rep_dict: dict = { Path(rep_pdb).stem: rep_pdb for rep_pdb in rep_pdbs }
    must_txt: str = path.join(data_dir, 'must.txt')
    must_list: list = []
    with open(must_txt) as file:
        for line in file:
            must_list.append(line.rstrip())

    if not path.isfile(ref['pdb']):
        raise FileNotFoundError(f"{ref['pdb']} file not found in the data directory")
    if len(rep_pdbs) == 0:
        raise FileNotFoundError("Expected several files in reps/ in the data directory, found zero instead")

    return ref, rep_dict, must_list

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
    """
    runs the opsinalign3d workflow

    :param query_pdbs: input PDB files
    :param output_dir: output directory path
    :param methods: methods to be used for the structural alignment
    :param threads: number of threads
    :param force: whether to overwrite data in the output directory if exists
    """
    from opsintools.scripts.t_coffee import t_coffee, check_t_coffee_methods
    from opsintools.scripts.tm_pos import tm_pos
    from multiprocessing import Pool

    logger = logging.getLogger(__name__)
    logger.propagate = False

    check_t_coffee_methods(methods)

    create_output_dir(output_dir, force)
    t_coffee_aln = path.join(output_dir, 't_coffee.aln')
    pdb_dict = {}
    for query_pdb in query_pdbs:
        pdb_stem = Path(query_pdb).stem
        if pdb_stem in pdb_dict:
            raise ValueError(f"Duplicate query name: {pdb_stem}")
        pdb_dict[pdb_stem] = query_pdb

    logger.info("Doing the structural alignment")
    t_coffee(pdb_dict, t_coffee_aln, methods, threads)
    output = Tcoffee(t_coffee_aln)
    logger.info("Finished")
    return output

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
        max_seq_id: float = MAX_SEQ_ID) -> dict | None:
    """
    runs the opsinmap3d workflow

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
    """
    from opsintools.scripts.prot_trim_filter import prot_trim_filter
    from opsintools.scripts.us_align import us_align
    from opsintools.scripts.score_alignments import score_alignments
    from opsintools.scripts.t_coffee import t_coffee, check_t_coffee_methods
    from opsintools.scripts.tm_pos import tm_pos
    from multiprocessing import Pool

    logger = logging.getLogger(__name__)
    logger.propagate = False

    logger.info("Checking the input")

    check_t_coffee_methods(methods)

    aln_dir: str = path.join(output_dir, "alignments")

    if not path.isfile(query_pdb):
        raise FileNotFoundError(f"Input file {query_pdb} not found")

    create_output_dir(output_dir, force)
    create_output_dir(aln_dir, force)

    ref: dict
    rep_dict: dict
    must_list: list
    ref, rep_dict, must_list = read_database(data_dir)

    aln_to_ref: str = path.join(output_dir, 'aln_to_ref.txt')
    trimmed_pdb: str = path.join(output_dir, 'trimmed.pdb')
    t_coffee_aln: str = path.join(output_dir, 't_coffee.aln')
    json_output: str = path.join(output_dir, 'opsinmap.json')

    logging.info("Aligning the query to the reference")
    us_align(ref['pdb'], query_pdb, aln_to_ref)

    rep_alns: dict = { rep_id: path.join(output_dir, "alignments", rep_id + '.txt') for rep_id in rep_dict }

    logging.info("Trimming the query")
    prot_trim_filter(aln_to_ref, query_pdb, ref['pdb'], trimmed_pdb, pad_n = pad_n, pad_c = pad_c)

    logger.info("Picking representatives for structural alignment")
    with Pool(threads) as pool:
        results: list = pool.starmap(us_align, [ (rep_dict[rep_id], trimmed_pdb, rep_alns[rep_id]) for rep_id in rep_dict ])

    chosen_templates: list = score_alignments(rep_alns.values(), n_templates, max_seq_id, must_list)
    chosen_pdbs: dict = { 'query': trimmed_pdb }
    for rep in chosen_templates:
        chosen_pdbs[rep] = rep_dict[rep]

    logger.info("Doing the structural alignment")
    t_coffee(chosen_pdbs, t_coffee_aln, methods = methods, threads = threads)

    output: dict = tm_pos(t_coffee_aln, trimmed_pdb, 'query', ref)

    with open(json_output, 'w') as file:
        json.dump(output, file, indent = 2)
    logger.info("Finished")
    return output

def opsinalign3d_cli() -> None:
    from argparse import ArgumentParser
    from sys import argv

    parser = ArgumentParser(description = 'Opsin alignment with t-coffee.')

    main_group = parser.add_argument_group('Arguments')

    main_group.add_argument('-i', nargs='+', metavar = 'FILENAME', required = True, help = 'input PDB structures (only chain A will be used)')
    main_group.add_argument('-o', metavar = 'DIRNAME', required = True, help = 'output directory')
    main_group.add_argument('-f', action = 'store_true', help = 'whether to overwrite files in the output directory if it exists')
    main_group.add_argument('-t', metavar = 'INT', type = int, default = THREADS, help = f'number of threads to use (default: {THREADS})')

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

def opsinmap3d_cli():
    from argparse import ArgumentParser
    from sys import argv

    parser = ArgumentParser(description = 'Opsin homology based for a query 3D structure.')

    main_group = parser.add_argument_group('Arguments')

    main_group.add_argument('-i', metavar = 'FILENAME', required = True, help = 'query PDB structure (only chain A will be used)')
    main_group.add_argument('-d', metavar = 'FILENAME', required = True, help = 'opsin data directory')
    main_group.add_argument('-o', metavar = 'DIRNAME', required = True, help = 'output directory')
    main_group.add_argument('-f', action = 'store_true', help = 'whether to overwrite files in the output directory if it exists')
    main_group.add_argument('-n', metavar = 'INT', type = int, default = THREADS, help = f'number of template structures to use (default: {N_REPS})')
    main_group.add_argument('-t', metavar = 'INT', type = int, default = THREADS, help = f'number of threads to use (default: {THREADS})')

    adv_group = parser.add_argument_group('Advanced')

    adv_group.add_argument('--pad-n', default = PAD_N, type = int, help = f'N-terminal padding for query trimming (default: {PAD_N})')
    adv_group.add_argument('--pad-c', default = PAD_C, type = int, help = f'C-terminal padding for query trimming (default: {PAD_C})')
    adv_group.add_argument('--max-seq-id', default = MAX_SEQ_ID, type = int, help = f'Maximum sequence identity for templates (default: {MAX_SEQ_ID})')
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
        pad_n = args.pad_n, pad_c = args.pad_c, max_seq_id = args.max_seq_id
    )

def run_with_logger(func, **args):

    # Define logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(OpsinToolsLoggerFormatter())
    logger.addHandler(ch)

    try:
        func(**args)
    except CalledProcessError as e:
        logger.fatal(f"Got return code {e.returncode}: {e.stderr.decode()}")
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
