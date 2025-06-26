from tempfile import TemporaryDirectory
from pyfaidx import Fasta
import subprocess
import shutil
from os import path
from Bio import SeqIO

def query_reps_fasta(chosen_pdbs, query_fasta, rep_fasta, output_fasta):
    with open(output_fasta, 'w') as out_fh:
        query = SeqIO.read(query_fasta, "fasta")
        out_fh.write(f">{query.id}\n{query.seq}\n")
        fasta = Fasta(rep_fasta)
        for rep in chosen_pdbs.keys():
            rep_seq = fasta[rep]
            out_fh.write(f">{rep}\n{rep_seq}\n")

# Create a template for query alignment
def query_reps_template(chosen_pdbs, query_pdb, query_id, template_file):
    with open(template_file, 'w') as out_fh:
        query_pdb_real = path.realpath(query_pdb)
        out_fh.write(f">{query_id}  _P_ {query_pdb_real}\n")
        for rep, rep_pdb in chosen_pdbs.items():
            rep_pdb_real = path.realpath(rep_pdb)
            out_fh.write(f">{rep}  _P_ {rep_pdb_real}\n")

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
        cmd +=[ '-' + key, str(value) ]
    with open(log_file_real, 'w') as log_fh:
        process = subprocess.run(cmd, stdout = log_fh, stderr = log_fh, cwd = work_dir)
        process.check_returncode()

def t_coffee(query_id, query_pdb, query_fasta, chosen_pdbs, rep_fasta, output_alignment, log_file, methods, threads):
    if shutil.which("t_coffee") is None:
        raise FileNotFoundError("t_coffee not found in PATH")
    with TemporaryDirectory() as work_dir:
        combined_fasta = path.join(work_dir, 'combined.fasta')
        template_file = path.join(work_dir, 'template.txt')
        query_reps_fasta(chosen_pdbs, query_fasta, rep_fasta, combined_fasta)
        query_reps_template(chosen_pdbs, query_pdb, query_id, template_file)
        run_t_coffee(work_dir, combined_fasta, template_file, output_alignment, log_file, methods, threads)
