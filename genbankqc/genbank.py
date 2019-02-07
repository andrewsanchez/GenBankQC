import re
from pathlib import Path

import attr
from logbook import Logger

from genbankqc import config, Species, Metadata, AssemblySummary

taxdump_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


@attr.s
class Genbank(object):
    log = Logger("GenBank")
    root = attr.ib(default=Path(), converter=Path)

    def __attrs_post_init__(self):
        self.paths = config.Paths(root=self.root, subdirs=["metadata", ".logs"])

    def info(self):
        patterns = ["*/*/dmx.csv", "*/*/*/tree.svg", "*/*/stats.csv"]
        info = []
        for pattern in patterns:
            import pdb; pdb.set_trace()
            files = self.root.glob(pattern)
            count = 0
            empty_files = []
            for file_ in files:
                count += 1
                if not file_.stat().st_size:
                    empty_files.append(file_.name)
            info.append(f"{pattern.split('/')[-1]}:")
            info.append(f"{count:>8} existing files")
            info.append(f"{len(empty_files):>8} empty files")
            for empty in empty_files:
                info.append(f"Empty:  {empty:>8}")
        return "\n".join(info)

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
        self.prune()
        for species in self.species():
            try:
                species.qc()
            except Exception:
                self.log.error(f"qc command failed for {species.name}")
                self.log.exception()

    def prune(self):
        """Prune all files that aren't latest assembly versions."""
        p_id = re.compile("GCA_[0-9]*.[0-9]")  # patterns for matching accession IDs
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
        assumbly_summary = AssemblySummary(self.paths.metadata)
        previous_versions = set(d_local.keys()) - set(assumbly_summary.ids)
        for i in previous_versions:
            for f in d_local[i]:
                f.unlink()
                self.log.info(f"Removed {f}")

    def metadata(self, email, sample=False, update=True):
        """Download and join all metadata and write out .csv for each species"""
        metadata_ = Metadata(
            self.paths.metadata, email=email, sample=sample, update=update
        )
        return metadata_

    def species_metadata(self, metadata):
        for species in self.species(metadata.assembly_summary):
            species.select_metadata(metadata)
