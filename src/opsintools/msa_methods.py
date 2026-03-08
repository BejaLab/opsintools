def parse_args(description):
    from argparse import ArgumentParser

    parser = ArgumentParser(description = 'A wrapper for mustang MSA to be used with t_coffee.')

    parser.add_argument('input', nargs = '?', metavar = 'INPUT', help = 'concatenated PDB structures')
    parser.add_argument('-i', metavar = 'INPUT', help = 'concatenated PDB structures')
    parser.add_argument('-o', metavar = 'OUTPUT', required = True, help = 'output file')

    args = parser.parse_args()
    if not args.input and not args.i:
        parser.error("An input must be provided either via -i or as a positional argument.")
    if args.input and args.i:
        parser.error("Ambiguous input: provided both -i and a positional argument.")

    from pathlib import Path

    full_input = args.input if args.input else args.i
    output_file = args.o

    check_files = [ full_input ]
    if len(full_input) > 1 and full_input.startswith('P'):
        check_files.append(full_input[1:])

    for input_file in check_files:
        if Path(input_file).is_file():
            return input_file, output_file

    parser.error("Input file does not exist")

def mustang_msa_cli():
    input_file, output_file = parse_args('A wrapper for mustang MSA for TC_LINE_SEPARATOR-concatenated input.')
    from opsintools.scripts.t_coffee import run_mustang_msa
    run_mustang_msa(input_file, output_file)

def mtm_align_msa_cli():
    input_file, output_file = parse_args('A wrapper for mTM-align MSA for TC_LINE_SEPARATOR-concatenated input.')
    from opsintools.scripts.t_coffee import run_mtm_align_msa
    run_mtm_align_msa(input_file, output_file)
