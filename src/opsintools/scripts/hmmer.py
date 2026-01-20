import subprocess
import shutil
from pathlib import Path
from Bio import AlignIO
from Bio.Align import Alignment
from io import StringIO
from os import PathLike
from typing import Any

def hmmsearch(fasta_file: PathLike, hmm_file: PathLike, output_dir: PathLike, no_filters: bool = False, no_bias: bool = False, threads: int = 1):
    """hmmsearch for a fasta file and an hmm file"""
    if shutil.which("hmmsearch") is None:
        raise FileNotFoundError("hmmsearch not found in PATH")
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / "hmmsearch.txt"
    log_file = output_path / "hmmsearch.log"
    cmd = [ "hmmsearch", "-o", output_file, "--cpu", str(threads) ]
    if no_filters:
        cmd.append("--max")
    if no_bias:
        cmd.append("--nobias")
    cmd += [ hmm_file, fasta_file ]
    with open(log_file, 'w') as fh:
        process = subprocess.run(cmd, stdout = fh, stderr = fh)
        process.check_returncode()
    return output_file, log_file

def hmmalign(fasta_file, hmm_file, output_dir: PathLike):
    """hmmalign for a fasta file and an hmm file"""
    if shutil.which("hmmalign") is None:
        raise FileNotFoundError("hmmalign not found in PATH")
    output_path = Path(output_dir)
    output_file = output_path / "hmmalign.sto"
    log_file = output_path / "hmmalign.log"
    cmd = [ "hmmalign", "-o", output_file, hmm_file, fasta_file ]
    with open(log_file, 'w') as fh:
        process = subprocess.run(cmd, stdout = fh, stderr = fh)
        process.check_returncode()
    return output_file, log_file

def local_to_global(dom, hmm_seq, query_seq):
    hmm_missing_l = dom['hmm']['from'] - 1
    ali_missing_l = dom['ali']['from'] - 1
    hmm_missing_r = len(hmm_seq) - dom['hmm']['to']
    ali_missing_r = len(query_seq) - dom['ali']['to']
    dom['hmm'] = {
        'from': 1,
        'to': len(hmm_seq),
        'seq': hmm_seq[:hmm_missing_l] + '-' * ali_missing_l + dom['hmm']['seq'].replace('.','-').upper() + hmm_seq[dom['hmm']['to']:] + '-' * ali_missing_r
    }
    dom['ali'] = {
        'from': 1,
        'to': len(query_seq),
        'seq': '-' * hmm_missing_l + query_seq[:ali_missing_l].lower() + dom['ali']['seq'] + '-' * hmm_missing_r + query_seq[dom['ali']['to']:].lower()
    }
    dom['PP'] = '.' * (hmm_missing_l + ali_missing_l) + dom['PP'] + '.' * (hmm_missing_r + ali_missing_r)
    if 'matches' in dom:
        dom['matches'] = ' ' * (hmm_missing_l + ali_missing_l) + dom['matches'] + ' ' * (hmm_missing_r + ali_missing_r)
    return dom

def sync_pwas(ref1, seq1, pp1, ref2, seq2, pp2):
    fasta1 = '\n'.join([ ">ref", ref1, ">query1", str(seq1), ">pp1", pp1 ])
    fasta2 = '\n'.join([ ">ref", ref2, ">query2", str(seq2), ">pp2", pp2 ])
    aln1 = AlignIO.read(StringIO(fasta1), "fasta").alignment
    aln2 = AlignIO.read(StringIO(fasta2), "fasta").alignment
    return Alignment.from_alignments_with_same_reference([aln1, aln2])

def chain_local_alignments(
    query_seq: str,
    ref_seq: str,
    alignments: list[dict[str, Any]],
    min_gap: int = -20,
    max_gap: int = 100
) -> list[dict[str, Any]]:
    class Align:
        def __init__(self, data: dict[str, Any]):
            q = data['ali']
            r = data['hmm']
            self.q_from = q['from']
            self.q_to = q['to']
            self.q_seq = q['seq']
            self.q_aligned_len = len(q['seq'])
            self.r_from = r['from']
            self.r_to = r['to']
            self.r_seq = r['seq']
            self.pp_scores = data['PP']  # Posterior probabilities string
            self.length = len(self.pp_scores)
            
            # Metadata for the domain
            self.score = data.get('score', 0.0)
            self.c_evalue = data.get('c_evalue', 1.0)
            self.i_evalue = data.get('i_evalue', 1.0)
            
            self.r_res_count = sum(1 for c in self.r_seq if c not in '-.')
            self.r_end = self.r_from + self.r_res_count - 1

    blocks = [Align(a) for a in alignments]
    if not blocks:
        return []

    # 1. Filter: Remove alignments contained within others on the reference
    n = len(blocks)
    sorted_by_score = sorted(range(n), key=lambda i: blocks[i].length, reverse=True)
    keep = [True] * n
    for i in sorted_by_score:
        if not keep[i]: continue
        bi = blocks[i]
        for j in range(n):
            if i == j or not keep[j]: continue
            bj = blocks[j]
            if (bi.r_from <= bj.r_from and bj.r_end <= bi.r_end):
                keep[j] = False

    blocks = [blocks[i] for i in range(n) if keep[i]]
    blocks.sort(key=lambda b: b.q_from)
    
    # 2. Dynamic Programming with Collinearity Check
    n = len(blocks)
    dp = [0.0] * n
    prev = [-1] * n

    for i in range(n):
        dp[i] = blocks[i].q_aligned_len
        for j in range(i):
            gap_q = blocks[i].q_from - blocks[j].q_to - 1
            gap_r = blocks[i].r_from - blocks[j].r_end - 1
            
            if gap_r >= min_gap and gap_q <= max_gap and gap_r <= max_gap:
                candidate = dp[j] + blocks[i].q_aligned_len
                if candidate > dp[i]:
                    dp[i] = candidate
                    prev[i] = j

    # 3. Greedy Extraction of Chains
    used = [False] * n
    results = []
    while True:
        best_i = -1
        best_score = -float('inf')
        for i in range(n):
            if not used[i] and dp[i] > best_score:
                best_score = dp[i]
                best_i = i
        if best_i == -1:
            break
            
        chain = []
        curr = best_i
        while curr != -1:
            chain.append(curr)
            used[curr] = True
            curr = prev[curr]
        chain.reverse()
        chain_blocks = [blocks[k] for k in chain]

        # 4. Merging Logic and Metadata Collection
        merged_q_parts = [chain_blocks[0].q_seq]
        merged_r_parts = [chain_blocks[0].r_seq]
        merged_pp_parts = [chain_blocks[0].pp_scores]
        
        # Lists for the requested metrics
        chain_scores = [b.score for b in chain_blocks]
        chain_c_evals = [b.c_evalue for b in chain_blocks]
        chain_i_evals = [b.i_evalue for b in chain_blocks]

        for i in range(1, len(chain_blocks)):
            prev_b = chain_blocks[i - 1]
            b = chain_blocks[i]

            gap_q = b.q_from - prev_b.q_to - 1
            gap_r = b.r_from - prev_b.r_end - 1

            if gap_r < 0:
                overlap_residues = -gap_r
                chars_to_strip = 0
                res_found = 0
                for char in b.r_seq:
                    chars_to_strip += 1
                    if char not in '-.':
                        res_found += 1
                    if res_found == overlap_residues:
                        break
                
                b.r_seq = b.r_seq[chars_to_strip:]
                b.q_seq = b.q_seq[chars_to_strip:]
                b.pp_scores = b.pp_scores[chars_to_strip:]
                gap_r = 0 

            gap_len = max(gap_q, gap_r, 0)
            if gap_len > 0:
                q_gap_seq = query_seq[prev_b.q_to : prev_b.q_to + max(0, gap_q)].lower()
                r_gap_seq = ref_seq[prev_b.r_end : prev_b.r_end + max(0, gap_r)].lower()
                merged_q_parts.append(q_gap_seq + '-' * (gap_len - len(q_gap_seq)))
                merged_r_parts.append(r_gap_seq + '-' * (gap_len - len(r_gap_seq)))
                merged_pp_parts.append('.' * gap_len)

            merged_q_parts.append(b.q_seq)
            merged_r_parts.append(b.r_seq)
            merged_pp_parts.append(b.pp_scores)

        results.append({
            'ali': {'seq': ''.join(merged_q_parts), 'from': chain_blocks[0].q_from, 'to': chain_blocks[-1].q_to},
            'hmm': {'seq': ''.join(merged_r_parts), 'from': chain_blocks[0].r_from, 'to': chain_blocks[-1].r_end},
            'PP': ''.join(merged_pp_parts),
            'components': len(chain_blocks),
            'score': chain_scores,
            'c_evalue': chain_c_evals,
            'i_evalue': chain_i_evals
        })

    results.sort(key=lambda x: x['ali']['from'])
    for i, res in enumerate(results):
        res['num'] = i + 1
    return results
