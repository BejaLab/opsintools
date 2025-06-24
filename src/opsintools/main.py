from pathlib import Path
from glob import glob
from os import path, makedirs
import json, re, logging
from subprocess import CalledProcessError
from typing import Final

# Constants
PAD_N: Final[int]= 30
PAD_C: Final[int] = 30
THREADS: Final[int] = 1
N_REPS: Final[int] = 20
MIN_SEQ_ID: Final[int] = 90
METHODS: list[str] = [ 'sap_pair', 'mustang_pair', 't_coffee_msa', 'probcons_msa' ]
ALL_METHODS: list[str] = [ 'test_pair', 'fast_pair', 'exon3_pair', 'exon2_pair', 'exon_pair', 'blastr_pair', 'promo_pair', 'clean_slow_pair', 'slow_pair', 'hash_pair', 'biphasic_pair', 'proba_prfpair', 'fs_pair', 'proba_pair', 'best_pair4prot', 'best_pair4rna', 'lalign_id_pair', 'fs_lalign_id_pair', 'seq_pair', 'externprofile_pair', 'hh_pair', 'co_pair', 'cwprofile_pair', 'cdna_fast_pair', 'cdna_cfast_pair', 'old_clustalo_pair', 'mafftsparsecore_pair', 'dynamic_pair', '3dcoffee_pair', 'expresso_pair', 'accurate_pair', 'psicoffee_pair', 'clustaloNF_pair', 'clustalw2_pair', 'clustalw_pair', 'uppNF_pair', 'upp_pair', 'msa_pair', 'dca_pair', 'dialigntx_pair', 'dialignt_pair', 'poa_pair', 'msaprobs_pair', 'probcons_pair', 'probconsRNA_pair', 'muscle_pair', 'mus4_pair', 't_coffee_pair', 'pcma_pair', 'kalign_pair', 'amap_pair', 'proda_pair', 'prank_pair', 'fsa_pair', 'consan_pair', 'famsa_pair', 'align_pdbpair', 'lalign_pdbpair', 'extern_pdbpair', 'thread_pair', 'fugue_pair', 'pdb_pair', 'sap_pair', 'sara_pair', 'daliweb_pair', 'dali_pair', 'mustang_pair', 'TMalign_pair', 'ktup_msa', 'blastp_msa', 'old_clustalo_msa', 'dynamic_msa', '3dcoffee_msa', 'expresso_msa', 'accurate_msa', 'psicoffee_msa', 'famsa_msa', 'clustalo_msa', 'mafft_msa', 'mafftginsi_msa', 'mafftfftns1_msa', 'mafftfftnsi_msa', 'mafftnwnsi_msa', 'mafftsparsecore_msa', 'mafftsparsecore_msa', 'mafftlinsi_msa', 'maffteinsi_msa', 'maffteinsi_pair', 'clustaloNF_msa', 'clustalw2_msa', 'clustalw_msa', 'uppNF_msa', 'upp_msa', 'msa_msa', 'dca_msa', 'dialigntx_msa', 'dialignt_msa', 'poa_msa', 'msaprobs_msa', 'probcons_msa', 'probconsRNA_msa', 'muscle_msa', 'mus4_msa', 't_coffee_msa', 'pcma_msa', 'kalign_msa', 'amap_msa', 'proda_msa', 'fsa_msa', 'tblastx_msa', 'tblastpx_msa', 'plib_msa', 'famsa_msa' ] 

# Define logger
class CustomFormatter(logging.Formatter):
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
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(CustomFormatter())
logger.addHandler(ch)

def read_database(data_dir: str) -> tuple[dict, dict, str, list]:
    """
    checks the opsinmap3d database

    :param data_dir: data directory
    :return: tuple of ref (parsed reference data), rep_dict (dict of { rep: its pdb path }, rep_fasta (fasta of representatives), must_list (list of reps to be included always)
    """
    ref_json: str = path.join(data_dir, 'ref.json')
    if not path.isfile(ref_json):
        raise FileNotFoundError(f"File {ref_json} does not exist")
    with open(ref_json) as file:
        ref: dict = json.load(file)
    rep_dir: str = path.join(data_dir, 'reps')
    ref['pdb'] = path.join(rep_dir, ref['id'] + '.pdb')
    rep_pdbs: list[str] = glob(path.join(rep_dir, '*.pdb'))
    rep_dict: dict = { Path(rep_pdb).stem: rep_pdb for rep_pdb in rep_pdbs }
    rep_fasta: str = path.join(data_dir, 'reps.fasta')
    must_txt: str = path.join(data_dir, 'must.txt')
    must_list: list = []
    with open(must_txt) as file:
        for line in file:
            must_list.append(line.rstrip())

    if not path.isfile(rep_fasta):
        raise FileNotFoundError("reps.fasta is not found in the data directory")
    if not path.isfile(ref['pdb']):
        raise FileNotFoundError(f"{ref['pdb']} file not found in the data directory")
    if len(rep_pdbs) == 0:
        raise FileNotFoundError("Expected several files in reps/ in the data directory, found zero instead")

    return ref, rep_dict, rep_fasta, must_list

def opsinmap3d(
        query_pdb: str,
        output_dir: str,
        data_dir: str,
        n_reps: int = N_REPS,
        methods: list = METHODS,
        threads: int = THREADS,
        force: bool = False,
        pad_n: int = PAD_N,
        pad_c: int = PAD_C,
        min_seq_id: int = MIN_SEQ_ID) -> None:
    """
    runs the opsinmap3d workflow

    :param query_pdb: input PDB file
    :param output_dir: output directory path
    :param data_dir: database directory
    :param n_reps: number of termplates to be included in the structural alignment
    :param methods: methods to be used for the structural alignment
    :param threads: number of threads
    :param force: whether to overwrite data in the output directory if exists
    :param pad_n: pad this number of residues to N-terminus when trimming the query
    :param pad_c: pad this number of residues to C-terminus when trimming the query
    :param min_seq_id: exclude templates with sequence identity higher than this
    """
    from opsintools.scripts.prot_trim_filter import prot_trim_filter
    from opsintools.scripts.us_align import us_align
    from opsintools.scripts.score_alignments import score_alignments
    from opsintools.scripts.t_coffee import t_coffee
    from opsintools.scripts.tm_pos import tm_pos
    from importlib import resources
    from multiprocessing import Pool

    try:
        for method in methods:
            if method not in ALL_METHODS:
                raise ValueError(f"{method} is not supported")

        logging.info("Checking the input")

        query_id: str = Path(query_pdb).stem
        aln_dir: str = path.join(output_dir, "alignments")

        if not path.isfile(query_pdb):
            raise FileNotFoundError(f"Input file {query_pdb} not found")

        if not re.search('^[\\w+-]+$', query_id):
            raise ValueError(f"Query name should only contain alphanumeric symbols or [_+-]: '{query_id}'")

        if path.exists(output_dir):
            if not path.isdir(output_dir):
                raise FileExistsError(f"The path {output_dir} exists and is not a directory")
            if not force:
                raise FileExistsError(f"Directory {output_dir} exists, not overriding")
        else:
            makedirs(output_dir)

        if path.exists(aln_dir):
            if not path.isdir(aln_dir):
                raise FileExistsError(f"The path {aln_dir} exists and is not a directory")
            if not force:
                raise FileExistsError(f"Directory {aln_dir} exists, not overriding")
        else:
            makedirs(aln_dir)

        ref: dict
        rep_dict: dict
        rep_fasta: str
        must_list: list
        ref, rep_dict, rep_fasta, must_list = read_database(data_dir)

        aln_to_ref: str = path.join(output_dir, 'aln_to_ref.txt')
        trimmed_pdb: str = path.join(output_dir, 'trimmed.pdb')
        trimmed_fasta: str = path.join(output_dir, 'trimmed.fasta')
        t_coffee_aln: str = path.join(output_dir, 't_coffee.fasta')
        t_coffee_log: str = path.join(output_dir, 't_coffee.log')
        json_output: str = path.join(output_dir, 'opsinmap.json')

        logging.info("Aligning the query to the reference")
        us_align(ref['pdb'], query_pdb, aln_to_ref)

        rep_alns: dict = { rep_id: path.join(output_dir, "alignments", rep_id + '.txt') for rep_id in rep_dict }

        logging.info("Trimming the query")
        prot_trim_filter(aln_to_ref, query_pdb, ref['pdb'], trimmed_fasta, trimmed_pdb, query_id, pad_n = pad_n, pad_c = pad_c)

        logging.info("Picking representatives for multiple structural alignment")
        with Pool(threads) as pool:
            results: list = pool.starmap(us_align, [ (rep_dict[rep_id], trimmed_pdb, rep_alns[rep_id]) for rep_id in rep_dict ])

        chosen_reps: list = score_alignments(rep_alns.values(), n_reps, min_seq_id, must_list)
        chosen_pdbs: dict = { rep: rep_dict[rep] for rep in chosen_reps }

        logging.info("Doing the structural alignment")
        t_coffee(query_id, trimmed_pdb, trimmed_fasta, chosen_pdbs, rep_fasta, t_coffee_aln, t_coffee_log, methods, threads)

        logging.info("Writing the output")
        tm_pos(t_coffee_aln, trimmed_pdb, json_output, query_id, ref)

        logging.info("Finished")

        logging.info(f"""Check out the output:
        * {aln_to_ref} - alignment of the query to the reference
        * {trimmed_pdb} - trimmed query structure
        * {trimmed_fasta} - trimmed query sequence
        * {t_coffee_aln} - structural alignment
        * {t_coffee_aln} - structural alignment log
        * {json_output} - mapping between the query and the reference""")
    except CalledProcessError as e:
        logging.fatal(f"Got return code {e.returncode}: {e.stderr.decode()}")
    except FileNotFoundError as e:
        logging.fatal(e)
    except FileExistsError as e:
        logging.fatal(e)
    except ValueError as e:
        logging.fatal(f"Unexpected value: {e}")
    except AssertionError as e:
        logging.fatal(f"Assumption violated: {e}")
    except Exception as e:
        logging.fatal(f"An error occurred: {e}")

def opsinmap3d_cli() -> None:
    from argparse import ArgumentParser
    from sys import argv

    parser = ArgumentParser(description = 'Opsin homology based for a query 3D structure.')

    main_group = parser.add_argument_group('Arguments')

    main_group.add_argument('-i', metavar = 'FILENAME', required = True, help = 'query PDB structure (only chain A will be used)')
    main_group.add_argument('-d', metavar = 'FILENAME', required = False, help = 'opsin data directory (the default is microbial rhodopsins)')
    main_group.add_argument('-o', metavar = 'DIRNAME', required = True, help = 'output directory')
    main_group.add_argument('-f', action = 'store_true', help = 'whether to overwrite files in the output directory if it exists')
    main_group.add_argument('-n', metavar = 'INT', type = int, default = THREADS, help = f'number of representative structures to use (default: {N_REPS})')
    main_group.add_argument('-t', metavar = 'INT', type = int, default = THREADS, help = f'number of threads to use (default: {THREADS})')

    adv_group = parser.add_argument_group('Advanced')

    adv_group.add_argument('--pad-n', default = PAD_N, type = int, help = f'N-terminal padding for query trimming (default: {PAD_N})')
    adv_group.add_argument('--pad-c', default = PAD_C, type = int, help = f'C-terminal padding for query trimming (default: {PAD_C})')
    adv_group.add_argument('--min-seq-id', default = PAD_C, type = int, help = f'Minimum sequence identity for templates (default: {MIN_SEQ_ID})')
    adv_group.add_argument('--methods', default = ','.join(METHODS), help = f't-coffee methods to use (default: {",".join(METHODS)})')

    args = parser.parse_args()

    methods = args.methods.split(',')
    opsinmap3d(args.i, args.o, args.d, n_reps = args.n, methods = methods, threads = args.t, force = args.f, pad_n = args.pad_n, pad_c = args.pad_c, min_seq_id = args.min_seq_id)
