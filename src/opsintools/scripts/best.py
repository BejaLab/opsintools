from typing import List, Dict, Any

def chain_local_alignments(
    query_seq: str,
    ref_seq: str,
    alignments: List[Dict[str, Any]],
    max_gap: int = 100
) -> List[Dict[str, Any]]:
    class Align:
        def __init__(self, data: Dict[str, Any]):
            q = data['query']
            r = data['ref']
            self.q_from = q['from']
            self.q_to = q['to']
            self.q_seq = q['seq']
            self.q_aligned_len = len(q['seq'])
            self.r_from = r['from']
            self.r_to = r.get('to', r['from'] + len(r['seq']) - 1)
            self.r_seq = r['seq']
            self.scores = data['scores']
            self.length = len(self.scores)
            self.r_res_count = sum(1 for c in self.r_seq if c not in '-.')
            self.r_end = self.r_from + self.r_res_count - 1

    blocks = [Align(a) for a in alignments]
    if not blocks:
        return []

    # Remove contained lower-score alignments (contained on reference)
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
    n = len(blocks)
    dp = [0.0] * n
    prev = [-1] * n
    for i in range(n):
        dp[i] = blocks[i].q_aligned_len
        for j in range(i):
            gap_q = blocks[i].q_from - blocks[j].q_to - 1
            gap_r = blocks[i].r_from - blocks[j].r_end - 1
            gap = max(gap_q, gap_r)
            if gap <= max_gap:
                candidate = dp[j] + blocks[i].q_aligned_len
                if candidate > dp[i]:
                    dp[i] = candidate
                    prev[i] = j

    # Extract all maximal chains greedily
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
        i = best_i
        while i != -1:
            chain.append(i)
            used[i] = True
            i = prev[i]
        chain.reverse()
        chain_blocks = [blocks[k] for k in chain]

        merged_q_parts = [chain_blocks[0].q_seq]
        merged_r_parts = [chain_blocks[0].r_seq]
        merged_s_parts = [chain_blocks[0].scores]

        for i in range(1, len(chain_blocks)):
            prev_b = chain_blocks[i - 1]
            b = chain_blocks[i]

            gap_q = b.q_from - prev_b.q_to - 1
            gap_r = b.r_from - prev_b.r_end - 1

            if gap_r < 0:
                overlap = -gap_r
                b.r_seq = ('-' * overlap) + b.r_seq[overlap:]
                b.scores = ('.' * overlap) + b.scores[overlap:]

            gap_len = max(gap_q, gap_r, 0)
            if gap_len > 0:
                q_gap = query_seq[prev_b.q_to : prev_b.q_to + gap_q].lower()
                r_gap = ref_seq[prev_b.r_end : prev_b.r_end + gap_r].lower()
                q_pad = '-' * (gap_len - len(q_gap))
                r_pad = '-' * (gap_len - len(r_gap))
                merged_q_parts.append(q_gap + q_pad)
                merged_r_parts.append(r_gap + r_pad)
                merged_s_parts.append('.' * gap_len)

            merged_q_parts.append(b.q_seq)
            merged_r_parts.append(b.r_seq)
            merged_s_parts.append(b.scores)

        merged_q = ''.join(merged_q_parts)
        merged_r = ''.join(merged_r_parts)
        merged_s = ''.join(merged_s_parts)

        first = chain_blocks[0]
        last = chain_blocks[-1]
        results.append({
            'query': {'seq': merged_q, 'from': first.q_from, 'to': last.q_to},
            'ref': {'seq': merged_r, 'from': first.r_from, 'to': last.r_end},
            'scores': merged_s,
            'components': len(chain_blocks)
        })

    results.sort(key=lambda x: x['query']['from'])
    return results
