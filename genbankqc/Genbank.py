import os
import traceback
import pandas as pd
from logbook import Logger

from genbankqc import Species


class Genbank:
    def __init__(self, path):
        """
        GenBank
        """
        self.genbank = path
        self.path = path
        self.assembly_summary = os.path.join(self.genbank,
                                             ".info/assembly_summary.txt")
        try:
            self.assembly_summary = pd.read_csv(self.assembly_summary,
                                                sep="\t", index_col=0)
        except FileNotFoundError:
            self.log.error("No assembly summary")
        self.log = Logger("GenBank")
        self.log.info("Instantiated")

    @property
    def species(self):
        for d in os.path.listdir(self.genbank):
            if d.startswith('.'):
                continue
            path = os.path.join(self.genbank, d)
            if os.path.isdir(path):
                yield Species(path, self.assembly_summary)

    def qc(self):
        for i in self.species:
            i.qc()
