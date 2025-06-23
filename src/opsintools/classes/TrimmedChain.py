from Bio import PDB
from Bio.SeqUtils import IUPACData

class TrimmedChain(PDB.Select):    
    def __init__(self, model, chain, start, end):
        self.residues = []
        self.model = model
        self.chain = chain
        self.start = start
        self.end   = end
        super(TrimmedChain, self).__init__()
    def accept_model(self, model):
        return model.id == self.model
    def accept_chain(self, chain):
        return chain.id == self.chain
    def accept_residue(self, residue):
        if self.start <= residue.id[1] <= self.end: 
            self.residues.append(residue)
            return True
    def get_seq(self):
        seq = [ IUPACData.protein_letters_3to1.get(res.resname.title(), 'X') for res in self.residues ]
        return ''.join(seq)
