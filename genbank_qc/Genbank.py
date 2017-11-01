import os
from genbank_qc import Species


class Genbank:
    def __init__(self, genbank):
        """Genbank"""
        self.genbank = genbank

    @property
    def species(self):
        for i in os.listdir(self.genbank):
            p = os.path.join(self.genbank, i)
            if os.path.isdir(p):
                yield Species(p)

    def qc(self):
        for i in self.species:
            i.qc()
