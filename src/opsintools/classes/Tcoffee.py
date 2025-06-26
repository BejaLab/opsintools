from Bio import SeqIO, AlignIO
from os import path

class Tcoffee():
    @staticmethod
    def parse_scores(scores_file: str) -> tuple[int, dict, dict, str]:
        aln_score: int = -1
        seq_scores: dict[str, int] = {}
        res_scores: dict[str, str] = {}
        cons_scores: str = ''
        with open(scores_file) as file:
            # find the line with the alignment score
            for line in file:
               if line.startswith('SCORE='):
                    aln_score = int(line.split('=')[1])
                    break
            assert aln_score > -1, f"Alignment score not found in file {scores_file}"
            # parse the block of the sequence scores
            parts: list[str]
            seq_name: str
            for line in file:
                if not line.strip():
                    break
                if ':' in line:
                    parts = line.split(':')
                    seq_name = parts[0].strip()
                    seq_score: int = int(parts[1])
                    seq_scores[seq_name] = seq_score
                    res_scores[seq_name] = ''
            # parse the block of the residue scores
            for line in file:
                if line.strip():
                    parts = line.split()
                    seq_name = parts[0]
                    assert seq_name in res_scores, f"Unexpected format: no sequence score was found to {seq_name}"
                    if seq_name == 'cons':
                        cons_scores += parts[1]
                    else:
                        res_scores[seq_name] += parts[1]
        return aln_score, seq_scores, res_scores, cons_scores

    def __init__(self, aln_file: str):
        self.aln_file = aln_file
        self.scores_file: str = aln_file + '.score_ascii'
        self.aln_score: int = -1
        self.seq_scores: dict[str, int] = {}
        self.res_scores: dict[str, str] = {}
        self.cons_scores: str
        if path.isfile(self.scores_file):
            self.aln_score, self.seq_scores, self.res_scores, self.cons_scores = Tcoffee.parse_scores(self.scores_file)
        self.alignment = SeqIO.to_dict(AlignIO.read(aln_file, "clustal"))
