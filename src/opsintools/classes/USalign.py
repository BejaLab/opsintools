import re
from Bio.SeqUtils import IUPACData

class USalign:

    @staticmethod
    def get_value(line, prefix, var = None, sep = '='):
        if var is None:
            if match := re.search(f"{prefix}{sep}" + '\\s*([0-9.]+)', line):
                return match.group(1)

    def get_file_name_chain(line):
        prefix, value = line.split(':', maxsplit = 1)
        file_name, chain_comments = value.rsplit(':', maxsplit = 1)
        chain, *comments = chain_comments.split()
        return file_name.lstrip(), chain

    @staticmethod
    def three_to_one(three_letter_code):
        return IUPACData.protein_letters_3to1.get(three_letter_code.title(), 'X')

    @staticmethod
    def get_structure_and_chain(line):
        prefix, file_name_and_chain = line.split(': ', maxsplit = 1)
        file_name, chain_and_comments = value.rsplit(':', maxsplit = 1)
        chain, comments = chain_and_comments.split(maxsplit = 1)
        return file_name, chain

    def __init__(self, aln_file):
        seq_id = rmsd = None
        tm_score_1 = tm_score_2 = None
        seq_len_1 = seq_len_2 = None
        alignment = []
        pdb_file_1 = pdb_file_2 = None
        with open(aln_file) as file:
            for line in file:
                if line.startswith('#Aligned'):
                    break
                if line.startswith('Name of Structure_1'):
                    pdb_file_1, self.chain_1 = USalign.get_file_name_chain(line)
                elif line.startswith('Name of Structure_2'):
                    pdb_file_2, self.chain_2 = USalign.get_file_name_chain(line)
                if val := USalign.get_value(line, 'Length of Structure_1', seq_len_1, sep = ':'):
                    seq_len_1 = int(val)
                if val := USalign.get_value(line, 'Length of Structure_2', seq_len_2, sep = ':'):
                    seq_len_2 = int(val)
                if val := USalign.get_value(line, 'RMSD', rmsd):
                    rmsd = float(val)
                if val := USalign.get_value(line, 'Seq_ID=n_identical/n_aligned', seq_id):
                    seq_id = float(val)
                if line.startswith('TM-score') and 'Structure_1' in line:
                    tm_score_1 = float(USalign.get_value(line, 'TM-score'))
                if line.startswith('TM-score') and 'Structure_2' in line:
                    tm_score_2 = float(USalign.get_value(line, 'TM-score'))
            for line in file:    
                if not line.startswith('#'):
                    res1, res2 = USalign.three_to_one(line[5:8]), USalign.three_to_one(line[21:24])
                    pos1, pos2 = int(line[10:14]), int(line[26:30])
                    distance = float(line[31:])
                    alignment.append((pos1, pos2, res1, res2, distance))
        assert pdb_file_1 is not None and pdb_file_2 is not None and seq_len_1 is not None and seq_len_2 is not None, f"input file names or lenghts not found in file {aln_file}"
        assert seq_id is not None, f"Sequence identity not found in file {aln_file}"
        assert rmsd is not None, f"RMSD not found in file {aln_file}"
        assert tm_score_1 is not None or tm_score_2 is not None, f"TM-scores not found in file {aln_file}"

        self.alignment = alignment
        self.seq_id = seq_id
        self.rmsd = rmsd
        self.pdb_files = (pdb_file_1, pdb_file_2)
        self.seq_lens = (seq_len_1, seq_len_2)
