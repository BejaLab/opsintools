import re
from Bio.SeqUtils import IUPACData

class USalign:
    """
    a class for parsing the output of US-align
    """
    @staticmethod
    def get_value(line: str, prefix: str, sep: str = '=') -> str | None: 
        """
        parses a key=value pair
        """
        if match := re.search(f"{prefix}{sep}" + '\\s*([0-9.]+)', line):
            return match.group(1)
        return None

    @staticmethod
    def get_file_name_chain(line: str) -> tuple[str, str]:
        """
        parses a line containing file name and chain
        """
        prefix, value = line.split(':', maxsplit = 1)
        file_name, chain_comments = value.rsplit(':', maxsplit = 1)
        chain, *comments = chain_comments.split()
        return file_name.lstrip(), chain

    @staticmethod
    def three_to_one(three_letter_code: str) -> str:
        """
        converts IUPAC three-letter amino acids to one-letter
        """
        return IUPACData.protein_letters_3to1.get(three_letter_code.title(), 'X')

    def __init__(self, aln_file: str):
        seq_id: float | None = None
        rmsd: float   | None = None
        tm_score_1: float | None = None
        tm_score_2: float | None = None
        seq_len_1: int = -1
        seq_len_2: int = -1
        alignment: list = []
        pdb_file_1: str = ""
        pdb_file_2: str = ""
        with open(aln_file) as file:
            for line in file:
                if line.startswith('#Aligned'):
                    break
                if line.startswith('Name of Structure_1'):
                    pdb_file_1, self.chain_1 = USalign.get_file_name_chain(line)
                elif line.startswith('Name of Structure_2'):
                    pdb_file_2, self.chain_2 = USalign.get_file_name_chain(line)
                if val := USalign.get_value(line, 'Length of Structure_1', sep = ':'):
                    seq_len_1 = int(val)
                if val := USalign.get_value(line, 'Length of Structure_2', sep = ':'):
                    seq_len_2 = int(val)
                if val := USalign.get_value(line, 'RMSD'):
                    rmsd = float(val)
                if val := USalign.get_value(line, 'Seq_ID=n_identical/n_aligned'):
                    seq_id = float(val)
                if line.startswith('TM-score') and 'Structure_1' in line:
                    if val := USalign.get_value(line, 'TM-score'):
                        tm_score_1 = float(val)
                if line.startswith('TM-score') and 'Structure_2' in line:
                    if val := USalign.get_value(line, 'TM-score'):
                        tm_score_2 = float(val)
            for line in file:    
                if not line.startswith('#'):
                    res1: str = USalign.three_to_one(line[5:8])
                    res2: str = USalign.three_to_one(line[21:24])
                    pos1, pos2 = int(line[10:14]), int(line[26:30])
                    distance = float(line[31:])
                    alignment.append((pos1, pos2, res1, res2, distance))
        assert pdb_file_1 and pdb_file_2 and seq_len_1 >= 0 and seq_len_2 >= 0, f"input file names or lenghts not found in file {aln_file}"
        assert seq_id is not None, f"Sequence identity not found in file {aln_file}"
        assert rmsd is not None, f"RMSD not found in file {aln_file}"
        assert tm_score_1 is not None or tm_score_2 is not None, f"TM-scores not found in file {aln_file}"

        self.alignment: list = alignment
        self.seq_id: float = seq_id
        self.rmsd: float = rmsd
        self.pdb_files: tuple[str, str] = (pdb_file_1, pdb_file_2)
        self.seq_lens: tuple[int, int] = (seq_len_1, seq_len_2)
        self.tm_scores: tuple[float | None, float | None] = (tm_score_1, tm_score_2)
