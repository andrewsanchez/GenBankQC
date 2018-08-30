import os
import pandas as pd
from logbook import Logger

from genbankqc import Species

taxdump_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"


class Genbank:
    def __init__(self, path):
        """
        GenBank
        """
        self.path = path
        self.assembly_summary_path = os.path.join(self.path, ".info/assembly_summary.txt")
        try:
            self.assembly_summary = pd.read_csv(self.assembly_summary, sep="\t", index_col=0)
        except FileNotFoundError:
            self.assembly_summary = pd.read_csv(assembly_summary_url, sep="\t", index_col=0, skiprows=1)
            self.assembly_summary.to_csv(self.assembly_summary_path, sep="\t")
        self.log = Logger("GenBank")
        self.log.info("Instantiated")

    @property
    def species(self):
        for d in os.listdir(self.path):
            if d.startswith('.'):
                continue
            path = os.path.join(self.path, d)
            if os.path.isdir(path):
                yield Species.Species(path, assembly_summary=self.assembly_summary)

    def init(self):
        """
        Create GenBank skeleton and get resources
        """
        if not os.path.isdir(self.path):
           os.mkdir(self.path)
            os.mkdir(os.path.join(self.path), '.info')

    def qc(self):
        for i in self.species:
            i.qc()
