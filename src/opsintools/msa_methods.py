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
