from tempfile import TemporaryDirectory
import subprocess
import shutil
from os import path
import os
from Bio import SeqIO
from Bio.Seq import Seq
import warnings
import re
from pathlib import Path
from opsintools.scripts import utils

warnings.filterwarnings('ignore', message = ".*(can't determine PDB ID|Ignoring unrecognized record).*")

NAME_MAX_LEN = 20

def check_pdb_file(pdb_name, pdb_file):
    if not path.exists(pdb_file):
        raise FileNotFoundError(f"File '{pdb_file}' does not exist")
    if re.search('[^\w.+-]', pdb_name):
        raise ValueError(f"Name '{pdb_name}' contains characters [^\w.+-]")
    if len(pdb_name) >= NAME_MAX_LEN:
        raise ValueError(f"Name '{pdb_name}' too long, shorten it to {NAME_MAX_LEN} characters")
    if pdb_name == 'cons':
        raise ValueError(f"Name 'cons' has a special meaning and cannot be used")
    return pdb_name

def write_inputs(pdb_dict, fasta_file, template_file):
    with open(fasta_file, 'w') as fasta_fh, open(template_file, 'w') as templ_fh:
        for seq_id, pdb_file in pdb_dict.items():
            record = utils.get_pdb_record(pdb_file, seq_id)
            SeqIO.write(record, fasta_fh, "fasta")
            pdb_file_name = Path(pdb_file).name
            templ_fh.write(f">{seq_id} _P_ {pdb_file_name}\n")

def split_t_coffee_pdb_input(cat_pdb_file: os.PathLike, dir_path: os.PathLike) -> [os.PathLike]:
    pdb_data = [[]]
    with open(cat_pdb_file) as file:
        for line in file:
            if line.startswith('#TC_LINE_SEPARATOR'):
                pdb_data.append([])
            elif line.startswith('TARGET_SEQ_NAME:'):
                names = line.split(':')[1].split()
            else:
                pdb_data[-1].append(line)
    if not pdb_data or pdb_data[-1] or pdb_data[-2]:
        raise ValueError(f"Expected t-coffee concatenated PDB file, got something else")
    split_paths = []
    for pdb, name in zip(pdb_data, names):
        split_path = dir_path / name
        split_paths.append(split_path.resolve())
        with open(split_path, 'w') as out_file:
            out_file.write(''.join(pdb))
    return split_paths

def run_mtm_align_msa(cat_pdb_file: str, output_alignment: str) -> None:
    """Run mTM-align MSA with the input structures in the form of a concatenated PDB file

    :param str cat_pdb_file: file name of the concatenated PDB
    :param str output_alignment: file name of the output FASTA file
    """
    output_alignment_path = Path(output_alignment).resolve()
    output_alignment_path.parent.mkdir(exist_ok = True)
    cat_pdb_path = Path(cat_pdb_file).resolve()

    with TemporaryDirectory() as work_dir:
        work_path = Path(work_dir)
        input_list_file = work_path / "_mtm_align_list_of_inputs.txt"
        with open(input_list_file, "w") as file:
            for pdb_path in split_t_coffee_pdb_input(cat_pdb_path, work_path):
                file.write(f"{pdb_path}\n")
        outdir = work_path / "_mtm_align_results"
        cmd = [ 'mtm-align', '-i', input_list_file, '-outdir', outdir ]
        process = subprocess.run(cmd, cwd = work_dir)
        process.check_returncode()
        fasta_path = outdir / "result.fasta"
        fasta_path.rename(output_alignment_path)

# TODO: Not ready yet, outputs a2m
def run_learnmsa_msa(input_file: str, output_alignment: str) -> None:
    """Run learnMSA

    :param str input_file: file name of the input fasta file
    :param str output_alignment: file name of the output FASTA file
    """
    cmd = [ 'learnMSA', '-i', input_file, '-o', output_alignment ]
    process = subprocess.run(cmd)
    process.check_returncode()

def run_mustang_msa(cat_pdb_file: str, output_alignment: str) -> None:
    """Run mustang MSA with the input structures in the form of a concatenated PDB file

    :param str cat_pdb_file: file name of the concatenated PDB
    :param str output_alignment: file name of the output FASTA file
    """
    output_alignment_path = Path(output_alignment).resolve()
    output_alignment_path.parent.mkdir(exist_ok = True)
    cat_pdb_path = Path(cat_pdb_file).resolve()

    with TemporaryDirectory() as work_dir:
        work_path = Path(work_dir)
        prefix = work_path / "_mustang_results"
        cmd = [ 'mustang', '-F', 'fasta', '-o', prefix, '-i' ]
        cmd += split_t_coffee_pdb_input(cat_pdb_path, work_path)
        process = subprocess.run(cmd, cwd = work_dir)
        process.check_returncode()
        fasta_path = prefix.with_suffix(".afasta")
        fasta_path.rename(output_alignment_path)

def run_t_coffee(work_dir, fasta_file, template_file, output_alignment, log_file, methods, threads):
    output_alignment_real = path.realpath(output_alignment)
    log_file_real = path.realpath(log_file)
    args = {
        'outfile': output_alignment_real,
        'output': 'aln,score_ascii',
        'method': ','.join(methods),
        'template_file': template_file,
        'pdb_min_sim': 90,
        'pdb_min_cov': 0,
        'thread': threads
    }
    cmd = [ 't_coffee', fasta_file ]
    for key, value in args.items():
        cmd += [ '-' + key, str(value) ]
    with open(log_file_real, 'w') as log_fh:
        process = subprocess.run(cmd, stdout = log_fh, stderr = log_fh, cwd = work_dir)
        process.check_returncode()

def write_tc_method(
        work_dir: str,
        executable: str,
        aln_mode: str,
        in_flag: str = "-i",
        out_flag: str = "-o",
        out_mode: str = "aln",
        seq_type: str = "P"
    ) -> str:
    """Writes a custom .tc_method file
    
    :param str executable: executable (= method)
    :param str aln_mode:   alignment mode (multiple or pairwise)
    :param str in_flag:    flag for the input file(s)
    :param str out_flag:   flag for the output file
    :param str out_mode:   output mode
    :param str seq_type:   sequence mode (P for PDB)
    :returns: output file name
    :rtype: str
    """
    file_path = Path(work_dir) / (executable + ".tc_method")
    cfg = {
        'EXECUTABLE': executable,
        'ALN_MODE': aln_mode,
        'IN_FLAG': in_flag,
        'OUT_FLAG': out_flag,
        'OUT_MODE': out_mode,
        'SEQ_TYPE': seq_type
    }
    with open(file_path, 'w') as file:
        file.write("*TC_METHOD_FORMAT_01\n")
        for key, val in cfg.items():
            file.write(f"{key} {val}\n")
    return str(file_path)

def check_t_coffee_methods(methods):
    ALL_METHODS = [ 'sap_pair', 'mustang_pair', 't_coffee_msa', 'probcons_msa', 'mustang_msa', 'mtm_align_msa', 'maffteinsi_msa', 'mafftfftns1_msa', 'mafftfftnsi_msa', 'mafftginsi_msa', 'mafftlinsi_msa', 'mafft_msa', 'mafftnwnsi_msa', 'mafftsparsecore_msa' ]
    if not methods or not isinstance(methods, list):
        raise ValueError(f"'methods' should be a non-empty list, got {methods} instead")
    for method in methods:
        if method not in ALL_METHODS:
            raise ValueError(f"Method '{method}' is not supported")

def write_custom_methods(methods: [str], work_dir: str) -> [str]:
    """Check list of methods for custom methods, write the corresponding
    .tc_method files and replace the method names with file names

    :param [str] methods: list of the methods
    :param str work_dir:  directory where to write the .tc_method files
    :returns: potentially updated list of methods
    :rtype: [str]
    """
    CUSTOM_METHODS = { 'mustang_msa': 'P', 'mtm_align_msa': 'P', 'learnmsa_msa': 'S' }
    methods_out = []
    for method in methods:
        if method in CUSTOM_METHODS:
            aln_mode = 'multiple' if method.endswith('msa') else 'pairwise'
            seq_type = CUSTOM_METHODS[method]
            method = write_tc_method(work_dir, method, aln_mode, seq_type)
        methods_out.append(method)
    return methods_out

def t_coffee(pdb_dict, output_prefix, methods, threads):
    if len(pdb_dict) < 2:
        raise ValueError("Cannot align less than two structures")
    if shutil.which("t_coffee") is None:
        raise FileNotFoundError("t_coffee not found in PATH")
    check_t_coffee_methods(methods)
    custom_methods = [ "mustang_msa" ]
    with TemporaryDirectory() as work_dir:
        methods = write_custom_methods(methods, work_dir)
        fasta_file = path.join(work_dir, 'sequences.fasta')
        template_file = path.join(work_dir, 'template.txt')
        for pdb_name, pdb_file in pdb_dict.items():
            check_pdb_file(pdb_name, pdb_file)
            source = path.realpath(pdb_file)
            dest = path.join(work_dir, Path(pdb_file).name)
            os.symlink(source, dest)

        log_file = output_prefix + '.log'
        write_inputs(pdb_dict, fasta_file, template_file)
        run_t_coffee(work_dir, fasta_file, template_file, output_prefix, log_file, methods, threads)
