from Bio import PDB
from Bio.SeqUtils import IUPACData

class TrimmedChain(PDB.Select):
    """
    allows trimming a specific chain in PDB
    """
    def __init__(self, model: int, chain: str, start: int, end: int):
        self.residues: list = []
        self.model: int = model
        self.chain: str = chain
        self.start: int = start
        self.end:   int = end
        super(TrimmedChain, self).__init__()

    def accept_model(self, model) -> bool:
        return model.id == self.model

    def accept_chain(self, chain) -> bool:
        return chain.id == self.chain

    def accept_residue(self, residue) -> bool:
        allow = self.start <= residue.id[1] <= self.end
        if allow: self.residues.append(residue)
        return allow

    def get_seq(self):
        """
        get the trimmed sequence after processing
        """
        seq: list = [ IUPACData.protein_letters_3to1.get(res.resname.title(), 'X') for res in self.residues ]
        return ''.join(seq)
