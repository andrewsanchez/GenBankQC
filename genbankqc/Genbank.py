import os
import traceback

from genbankqc import Species


class Genbank:
    def __init__(self, genbank):
        """Genbank"""
        self.genbank = genbank

    @property
    def species(self):
        dirs = (os.path.join(self.genbank, d)
                for d in os.listdir(self.genbank))
        dirs = (d for d in dirs if os.path.isdir(d))
        for d in dirs:
            fastas = len([f for f in os.listdir(d)
                          if f.endswith('fasta')])
            if fastas > 5:
                try:
                    yield Species(d)
                except:
                    print('Skipping ', d)
                    traceback.print_exc()

    def qc(self):
        for i in self.species:
            try:
                i.qc()
            except:
                print('Failed ', i.species)
                traceback.print_exc()
