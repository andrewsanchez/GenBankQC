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
            try:
                yield Species(d)
                print("Instantiated ", d)
            except:
                print('Skipping ', d)
                traceback.print_exc()
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
