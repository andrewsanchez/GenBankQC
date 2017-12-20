import os, traceback
from genbankqc import Species


class Genbank:
    def __init__(self, genbank):
        """Genbank"""
        self.genbank = genbank

    @property
    def species(self):
        for root, dirs, files in os.walk(self.genbank):
            for d in dirs:
                d = os.path.join(root, d)
                try:
                    yield Species(p)
                except:
                    print('Skipping ', d)
                    traceback.print_exc()

    def qc(self):
            for i in self.species:
                try:
                    i.qc()
                    print('Completed ', i.species)
                except:
                    print('Failed ', i.species)
                    traceback.print_exc()
