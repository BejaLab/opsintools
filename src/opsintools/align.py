from .globals import run_with_logger, logger, create_output_dir
from .globals import METHODS, THREADS
from opsintools.classes.Tcoffee import Tcoffee

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
    from pathlib import Path
    from multiprocessing import Pool

    check_t_coffee_methods(methods)

    create_output_dir(output_dir, force)
    output_path = Path(output_dir)

    t_coffee_aln = output_path / 't_coffee.aln'
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
