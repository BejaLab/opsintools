from tempfile import TemporaryDirectory
import subprocess
import shutil
from os import path
import os
from Bio import SeqIO
import warnings
import re
from pathlib import Path

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
            record = next(SeqIO.parse(pdb_file, "pdb-atom"))
            record.id = record.description = seq_id
            SeqIO.write(record, fasta_fh, "fasta")
            pdb_file_name = Path(pdb_file).name
            templ_fh.write(f">{seq_id} _P_ {pdb_file_name}\n")

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

def check_t_coffee_methods(methods):
    ALL_METHODS = [ 'test_pair', 'fast_pair', 'exon3_pair', 'exon2_pair', 'exon_pair', 'blastr_pair', 'promo_pair', 'clean_slow_pair', 'slow_pair', 'hash_pair', 'biphasic_pair', 'proba_prfpair', 'fs_pair', 'proba_pair', 'best_pair4prot', 'best_pair4rna', 'lalign_id_pair', 'fs_lalign_id_pair', 'seq_pair', 'externprofile_pair', 'hh_pair', 'co_pair', 'cwprofile_pair', 'cdna_fast_pair', 'cdna_cfast_pair', 'old_clustalo_pair', 'mafftsparsecore_pair', 'dynamic_pair', '3dcoffee_pair', 'expresso_pair', 'accurate_pair', 'psicoffee_pair', 'clustaloNF_pair', 'clustalw2_pair', 'clustalw_pair', 'uppNF_pair', 'upp_pair', 'msa_pair', 'dca_pair', 'dialigntx_pair', 'dialignt_pair', 'poa_pair', 'msaprobs_pair', 'probcons_pair', 'probconsRNA_pair', 'muscle_pair', 'mus4_pair', 't_coffee_pair', 'pcma_pair', 'kalign_pair', 'amap_pair', 'proda_pair', 'prank_pair', 'fsa_pair', 'consan_pair', 'famsa_pair', 'align_pdbpair', 'lalign_pdbpair', 'extern_pdbpair', 'thread_pair', 'fugue_pair', 'pdb_pair', 'sap_pair', 'sara_pair', 'daliweb_pair', 'dali_pair', 'mustang_pair', 'TMalign_pair', 'ktup_msa', 'blastp_msa', 'old_clustalo_msa', 'dynamic_msa', '3dcoffee_msa', 'expresso_msa', 'accurate_msa', 'psicoffee_msa', 'famsa_msa', 'clustalo_msa', 'mafft_msa', 'mafftginsi_msa', 'mafftfftns1_msa', 'mafftfftnsi_msa', 'mafftnwnsi_msa', 'mafftsparsecore_msa', 'mafftsparsecore_msa', 'mafftlinsi_msa', 'maffteinsi_msa', 'maffteinsi_pair', 'clustaloNF_msa', 'clustalw2_msa', 'clustalw_msa', 'uppNF_msa', 'upp_msa', 'msa_msa', 'dca_msa', 'dialigntx_msa', 'dialignt_msa', 'poa_msa', 'msaprobs_msa', 'probcons_msa', 'probconsRNA_msa', 'muscle_msa', 'mus4_msa', 't_coffee_msa', 'pcma_msa', 'kalign_msa', 'amap_msa', 'proda_msa', 'fsa_msa', 'tblastx_msa', 'tblastpx_msa', 'plib_msa', 'famsa_msa' ]
    if not methods or not isinstance(methods, list):
        raise ValueError(f"'methods' should be a non-empty list, got {methods} instead")
    for method in methods:
        if method not in ALL_METHODS:
            raise ValueError(f"Method '{method}' is not supported")

def t_coffee(pdb_dict, output_prefix, methods, threads):
    if len(pdb_dict) < 2:
        raise ValueError("Cannot align less than two structures")
    if shutil.which("t_coffee") is None:
        raise FileNotFoundError("t_coffee not found in PATH")
    check_t_coffee_methods(methods)
    with TemporaryDirectory() as work_dir:
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
