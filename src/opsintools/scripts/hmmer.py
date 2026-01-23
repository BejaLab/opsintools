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

from typing import Any

def chain_local_alignments(
    query_seq: str,
    ref_seq: str,
    alignments: list[dict[str, Any]],
    max_gap: int = 150
) -> list[dict[str, Any]]:
    class Align:
        def __init__(self, data: dict[str, Any]):
            q, r = data['ali'], data['hmm']
            self.q_from, self.q_to, self.q_seq = q['from'], q['to'], q['seq']
            self.r_from, self.r_to, self.r_seq = r['from'], r['to'], r['seq']
            self.pp_scores = data['PP']
            self.score = data.get('score', 0.0)
            self.c_evalue = data.get('c_evalue', 1.0)
            self.i_evalue = data.get('i_evalue', 1.0)
            self.r_end = self.r_from + sum(1 for c in self.r_seq if c not in '-.') - 1

        def truncate_end_to_query(self, target_q_to: int):
            if target_q_to >= self.q_to: return
            new_q, new_r, new_pp = [], [], []
            curr_q = self.q_from
            for q_c, r_c, pp_c in zip(self.q_seq, self.r_seq, self.pp_scores):
                new_q.append(q_c); new_r.append(r_c); new_pp.append(pp_c)
                if q_c not in '-.':
                    if curr_q == target_q_to: break
                    curr_q += 1
            self.q_seq, self.r_seq, self.pp_scores = "".join(new_q), "".join(new_r), "".join(new_pp)
            self.q_to = target_q_to
            self.r_end = self.r_from + sum(1 for c in self.r_seq if c not in '-.') - 1

        def strip_start_by_hmm_overlap(self, overlap_count: int):
            if overlap_count <= 0: return
            strip_idx, r_found = 0, 0
            for r_c in self.r_seq:
                strip_idx += 1
                if r_c not in '-.': r_found += 1
                if r_found == overlap_count: break
            q_removed = sum(1 for c in self.q_seq[:strip_idx] if c not in '-.')
            self.q_from += q_removed
            self.r_from += overlap_count
            self.q_seq, self.r_seq, self.pp_scores = self.q_seq[strip_idx:], self.r_seq[strip_idx:], self.pp_scores[strip_idx:]

    blocks = [Align(a) for a in alignments]
    if not blocks: return []

    # 1. DP Pathfinding: Must move forward on Query AND HMM
    n = len(blocks)
    blocks.sort(key=lambda b: b.q_from)
    dp, prev = [float(len(b.pp_scores)) for b in blocks], [-1] * n
    for i in range(n):
        for j in range(i):
            # Collinearity Constraint: block i must be after block j on both sequences
            if blocks[i].q_from > blocks[j].q_to and blocks[i].r_from > blocks[j].r_end:
                gap_q = blocks[i].q_from - blocks[j].q_to - 1
                gap_r = blocks[i].r_from - blocks[j].r_end - 1
                if gap_q <= max_gap and gap_r <= max_gap:
                    score = dp[j] + len(blocks[i].pp_scores)
                    if score > dp[i]:
                        dp[i], prev[i] = score, j

    # 2. Chain Extraction
    used, results = [False] * n, []
    while True:
        best_i, best_val = -1, -1.0
        for i in range(n):
            if not used[i] and dp[i] > best_val:
                best_val, best_i = dp[i], i
        if best_i == -1: break
            
        chain_idx = []
        curr = best_i
        while curr != -1:
            chain_idx.append(curr); used[curr] = True; curr = prev[curr]
        chain_idx.reverse()
        chain = [blocks[k] for k in chain_idx]

        # 3. Stitching with literal filling
        res_q, res_r, res_pp = [chain[0].q_seq], [chain[0].r_seq], [chain[0].pp_scores]
        
        for i in range(1, len(chain)):
            prev_b, curr_b = chain[i-1], chain[i]
            
            # Resolve overlaps (should be minimal with the DP constraint)
            if curr_b.q_from <= prev_b.q_to:
                prev_b.truncate_end_to_query(curr_b.q_from - 1)
                res_q[-1], res_r[-1], res_pp[-1] = prev_b.q_seq, prev_b.r_seq, prev_b.pp_scores

            # Fill the physical gap using original sequences
            q_bridge = query_seq[prev_b.q_to : curr_b.q_from - 1].lower()
            r_bridge = ref_seq[prev_b.r_end : curr_b.r_from - 1].lower()
            
            bridge_len = max(len(q_bridge), len(r_bridge))
            if bridge_len > 0:
                res_q.append(q_bridge + '-' * (bridge_len - len(q_bridge)))
                res_r.append(r_bridge + '-' * (bridge_len - len(r_bridge)))
                res_pp.append('.' * bridge_len)

            res_q.append(curr_b.q_seq); res_r.append(curr_b.r_seq); res_pp.append(curr_b.pp_scores)

        # Final string assembly
        final_q, final_r, final_pp = "".join(res_q), "".join(res_r), "".join(res_pp)

        # Calculate coordinates by inspecting final sequence (most reliable method)
        def get_coords(seq, start_val):
            first_idx = next((i for i, c in enumerate(seq) if c not in '-.'), 0)
            last_idx = next((i for i, c in enumerate(reversed(seq)) if c not in '-.'), 0)
            residues_before = sum(1 for c in seq[:first_idx] if c not in '-.')
            total_residues = sum(1 for c in seq if c not in '-.')
            return start_val + residues_before, start_val + residues_before + total_residues - 1

        # We use the chain's original first block as the coordinate anchor
        q_start, q_end = get_coords(final_q, chain[0].q_from)
        r_start, r_end = get_coords(final_r, chain[0].r_from)

        results.append({
            'ali': {'seq': final_q, 'from': q_start, 'to': q_end},
            'hmm': {'seq': final_r, 'from': r_start, 'to': r_end},
            'PP': final_pp, 'components': len(chain),
            'score': [b.score for b in chain], 'c_evalue': [b.c_evalue for b in chain]
        })

    results.sort(key=lambda x: x['ali']['from'])
    for i, res in enumerate(results): res['num'] = i + 1
    return results
