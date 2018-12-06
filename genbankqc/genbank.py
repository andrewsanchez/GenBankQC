import os
from pathlib import Path

import attr
import pandas as pd
from logbook import Logger

from genbankqc import config, Species, metadata

taxdump_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


@attr.s
class Genbank(object):
    log = Logger("GenBank")
    root = attr.ib(default=Path())

    def __attrs_post_init__(self):
        if not isinstance(self.root, Path):
            self.root = Path(self.root)
        self.paths = config.Paths(root=self.root, subdirs=["metadata", ".logs"])

    @property
    def species_directories(self):
        """Generator of `Path` objects for directories under `self.root`.
        Only species with more than ten FASTAs are included."""
        for item in self.root.iterdir():
            if not item.is_dir():
                continue
            dir_ = item.absolute()
            if len(list(dir_.glob("*fasta"))) < 10:
                self.log.info("Not enough FASTAs in {}".format(dir_))
                continue
            yield dir_

    def species(self, assembly_summary=None):
        """Generator of Species objects for directories returned by `species_directories`."""
        self.assembly_summary = metadata.AssemblySummary(self.paths.metadata)
        for dir_ in self.species_directories:
            yield Species(dir_, assembly_summary=self.assembly_summary.df)

    def qc(self):
        for species in self.species():
            species.qc()

    def biosample_metadata(self):
        biosample = metadata.BioSample(self.path.metadata)
        biosample.generate()

    def metadata(self):
        """Assumes existence of metadata files"""
        biosample = metadata.BioSample(self.path.metadata, read_existing=True)
        runs = metadata.SRA(self.paths.metadata / "sra_runs.tsv")
        for species in self.species():
            species.metadata(biosample=biosample.df, runs=sra_runs)
