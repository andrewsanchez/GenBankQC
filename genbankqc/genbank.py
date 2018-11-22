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

    @property
    def species_directories(self):
        for dir_ in os.listdir(self.root):
            species_dir = os.path.join(self.root, dir_)
            if not os.path.isdir(species_dir):
                continue
            if species_dir.startswith('.'):
                continue
            yield species_dir

    def species(self):
        """Iterate through all directories under self.root, yielding those
        that contain > 10 fastas.
        """
        for dir_ in os.listdir(self.species_directories):
            fastas = len([f for f in os.listdir(dir_) if f.endswith('fasta')])
            if fastas < 10:
                self.log.info("Not enough genomes for {}".format(dir_))
                continue
            yield Species(dir_, assembly_summary=self.assembly_summary)

    def qc(self):
        for species in self.species():
            species.qc()
