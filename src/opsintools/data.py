from .globals import run_with_logger

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

    if not file_or_accession:
        raise ValueError("No file or accession provided")
    if is_file and not Path(file_or_accession).is_file():
        raise ValueError(f"File {file_or_accession} does not exist")
    
    import importlib.resources
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
