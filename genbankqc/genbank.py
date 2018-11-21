import os
import attr
import pandas as pd
from logbook import Logger

from genbankqc import Paths
from genbankqc import Species

taxdump_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
assembly_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"


@attr.s
class Genbank(object):
    log = Logger("GenBank")
    root = attr.ib(default=os.getcwd())
    def __attrs_post_init__(self):
        self.paths = Paths(root=self.root, subdirs=['metadata'])

    def assembly_summary(self):
        self.paths.assembly_summary = os.path.join(self.paths.metadata, "assembly_summary.txt")
        try:
            self.assembly_summary = pd.read_csv(self.paths.assembly_summary, sep="\t", index_col=0)
        except FileNotFoundError:
            self.assembly_summary = pd.read_csv(assembly_summary_url, sep="\t",
                                                index_col=0, skiprows=1)
            self.assembly_summary.to_csv(self.paths.assembly_summary, sep="\t")
            self.log.info("Downloaded assembly_summary.txt")

    @property
    def species(self):
        """Iterate through all directories under self.path, yielding those
        that contain > 10 fastas.
        """
        for d in os.listdir(self.path):
            species_path = os.path.join(self.path, d)
            if not os.path.isdir(species_path):
                continue
            fastas = len([f for f in os.listdir(species_path) if f.endswith('fasta')])
            if fastas < 10:
                self.log.info("Not enough genomes for {}".format(d))
                continue
            yield Species(species_path, assembly_summary=self.assembly_summary)

    def qc(self):
        for species in self.species:
            species.qc()
