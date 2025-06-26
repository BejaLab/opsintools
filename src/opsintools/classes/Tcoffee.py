from Bio import SeqIO, AlignIO
from os import path

class Tcoffee:
    @staticmethod
    def parse_scores(scores_file: str) -> tuple[int, dict, dict]:
        aln_score: int = -1
        seq_scores: dict[int] = {}
        res_scores: dict[str] = {}
        with open(scores_file) as file:
            # find the line with the alignment score
            for line in file:
               if line.startswith('SCORE='):
                    aln_score = int(line.split('=')[1])
                    break
            assert aln_score > -1, f"Alignment score not found in file {scores_file}"
            # parse the block of the sequence scores
            for line in file:
                if not line.strip():
                    break
                if ':' in line:
                    parts: list[str] = line.split(':')
                    seq_name: str = parts[0].strip()
                    seq_score: int = int(parts[1])
                    seq_scores[seq_name] = seq_score
                    res_scores[seq_name] = ''
            # parse the block of the residue scores
            for line in file:
                if line.strip():
                    parts: list[str] = line.split()
                    seq_name: str = parts[0]
                    assert seq_name in res_scores, f"Unexpected format: no sequence score was found to {seq_name}"
                    res_scores[seq_name] += parts[1]
        return aln_score, seq_scores, res_scores

    def __init__(self, aln_file: str):
        self.aln_score: int = -1
        self.seq_scores: dict[int] = {}
        self.res_scores: dict[str] = {}
        scores_file: str = aln_file + '.score_ascii'
        if path.isfile(scores_file):
            self.aln_score, self.seq_scores, self.res_scores = Tcoffee.parse_scores(scores_file)
        self.alignment = SeqIO.to_dict(AlignIO.read(aln_file, "clustal"))
