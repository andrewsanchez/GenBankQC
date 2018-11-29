import os
from pathlib import Path

import attr
from logbook import Logger

from genbankqc import config, Species, metadata

taxdump_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


@attr.s
class Genbank(object):
    log = Logger("GenBank")
    root = attr.ib(default=Path())

    def __attrs_post_init__(self):
        self.paths = config.Paths(root=self.root, subdirs=["metadata", ".logs"])

    @property
    def species_directories(self):
        for dir_ in os.listdir(self.root):
            species_dir = os.path.join(self.root, dir_)
            if not os.path.isdir(species_dir):
                continue
            if species_dir.startswith("."):
                continue
            yield species_dir

    def species(self, assembly_summary=None):
        """Iterate through all directories under self.root, yielding those
        that contain > 10 fastas."""
        self.assembly_summary = metadata.AssemblySummary(self.paths.metadata)
        for dir_ in self.species_directories:
            fastas = len([f for f in os.listdir(dir_) if f.endswith("fasta")])
            if fastas < 10:
                self.log.info("Not enough genomes for {}".format(dir_))
                continue
            yield Species(dir_, assembly_summary=self.assembly_summary.df)

    def qc(self):
        for species in self.species():
            species.qc()

    def metadata(self):
        assembly_summary = metadata.AssemblySummary(self.paths.metadata)
        biosample = metadata.BioSample(self.paths.metadata)
        biosample.generate()
