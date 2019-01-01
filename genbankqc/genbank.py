import re
from pathlib import Path

import attr
from logbook import Logger

from genbankqc import config, Species, metadata

taxdump_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


@attr.s
class Genbank(object):
    log = Logger("GenBank")
    root = attr.ib(default=Path(), converter=Path)

    def __attrs_post_init__(self):
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
        for dir_ in self.species_directories:
            yield Species(dir_, assembly_summary=assembly_summary)

    def qc(self):
        for species in self.species():
            species.qc()

    def prune(self):
        """Prune all files that aren't latest assembly versions."""
        # patterns for matching accession IDs
        p_id = re.compile("GCA_[0-9]*.[0-9]")
        p_glob = "GCA_[0-9]*.[0-9]_*[fasta|msh|csv]"

        # IDs and associated files
        d_local = {}

        def glob_local():
            """Yield all files that match `p_glob`"""
            for path in self.root.rglob(p_glob):
                id_ = p_id.match(path.name).group()
                d_local[id_] = []
                yield path

        # Update `d_local` with a list containing paths for all matches
        for path in glob_local():
            id_ = p_id.match(path.name).group()
            d_local[id_].append(path)

        # Remove local files that aren't latest assembly versions
        assumbly_summary = metadata.AssemblySummary(self.paths.metadata)
        for i in set(d_local.keys()) - set(assumbly_summary.ids):
            for f in d_local[i]:
                f.unlink()

    def biosample_metadata(self, email):
        biosample = metadata.BioSample(self.paths.metadata, email=email)
        biosample.with_runs()

    def species_metadata(self, email):
        """Assumes existence of metadata files"""
        biosample = metadata.BioSample(self.paths.metadata, email, read_existing=True)
        runs = metadata.SRA(self.paths.metadata / "sra_runs.tsv")
        all_metadata = biosample.df.join(runs.df)
        assembly_summary = metadata.AssemblySummary(self.paths.metadata)
        for species in self.species(assembly_summary):
            species.metadata(all_metadata)
