import os.path as path

from genbank_qc import Species


class Genbank:
    def __init__(self, genbank):
        """Genbank"""
        self.genbank = genbank
        self.species = self.get_species()

    @property
    def species(self):
        for i in path.listdir(self.genbank):
            p = path.join(self.genbank, i)
            if path.isdir(p):
                yield Species(p)
