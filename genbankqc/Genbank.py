import os
import pandas as pd
import traceback

from genbankqc import Species


class Genbank:
    def __init__(self, genbank):

        """Genbank"""

        self.genbank = genbank
        self.assembly_summary = os.path.join(
            self.genbank,
            ".info/assembly_summary.txt",
        )
        try:
            self.assembly_summary = pd.read_csv(
                self.assembly_summary,
                sep="\t",
                index_col=0,
            )
        except FileNotFoundError:
            pass


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
                except Exception:
                    print('Skipping ', d)
                    traceback.print_exc()

    def qc(self):
        for i in self.species:
            try:
                i.qc()
            except Exception:
                print('Failed ', i.species)
                traceback.print_exc()
