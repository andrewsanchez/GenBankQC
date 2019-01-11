import subprocess
import xml.etree.cElementTree as ET
from pathlib import Path
from tempfile import mkdtemp

import attr
import pandas as pd
from logbook import Logger

from Bio import Entrez
from genbankqc import config
from tenacity import retry, stop_after_attempt, wait_fixed


@attr.s
class AssemblySummary(object):
    """Read in existing file or download latest assembly summary."""

    path = attr.ib(default=mkdtemp(), converter=Path)
    read = attr.ib(default=False)
    url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"

    def __attrs_post_init__(self):
        self.file_ = self.path / "assembly_summary.txt"
        if self.read:
            self.df = self._read()
        else:
            self.df = self._download()
        self.ids = self.df.index.tolist()

    @retry(stop=stop_after_attempt(3), wait=wait_fixed(2))
    def _download(self):
        print(f"Downloading {self.url}")
        df = pd.read_csv(self.url, sep="\t", index_col=0, skiprows=1)
        df.to_csv(self.file_, sep="\t")
        return df

    def _read(self):
        try:
            return pd.read_csv(self.file_, sep="\t", index_col=0)
        except FileNotFoundError:
            return self._download()


@attr.s
class BioSample(object):
    """Download and parse BioSample metadata for GenBank bacteria genomes."""

    outdir = attr.ib(converter=Path)
    email = attr.ib()
    sample = attr.ib(default=False)
    read_existing = attr.ib(default=False)

    attributes = [
        "BioSample",
        "SRA",
        "geo_loc_name",
        "collection_date",
        "strain",
        "isolation_source",
        "host",
        "collected_by",
        "sample_type",
        "sample_name",
        "host_disease",
        "isolate",
        "host_health_state",
        "serovar",
        "env_biome",
        "env_feature",
        "ref_biomaterial",
        "env_material",
        "isol_growth_condt",
        "num_replicons",
        "sub_species",
        "host_age",
        "genotype",
        "host_sex",
        "serotype",
        "host_disease_outcome",
    ]

    def __attrs_post_init__(self):
        self.paths = config.Paths(root=self.outdir)
        self.df = pd.DataFrame(index=["BioSample"], columns=self.attributes)
        self.data = [self.df]
        if self.read_existing:
            self.df = self.read()

    log = Logger("BioSample")

    # @retry(stop_max_attempt_number=3, stop_max_delay=10000, wait_fixed=100)
    def _esearch(
        self, db="biosample", term="bacteria[orgn] AND biosample_assembly[filter]"
    ):
        """Use NCBI's esearch to make a query"""
        Entrez.email = self.email
        esearch_handle = Entrez.esearch(db=db, term=term, usehistory="y")
        self.esearch_results = Entrez.read(esearch_handle)

    def _efetch(self):
        """Use NCBI's efetch to download esearch results"""
        db_xp = 'Ids/Id/[@db="{}"]'
        # Tuples for building XPath patterns
        xp_tups = [("SRA", "db", db_xp), ("BioSample", "db", db_xp)]
        for attrib in self.attributes:
            xp_tups.append(
                (
                    attrib,
                    "harmonized_name",
                    'Attributes/Attribute/[@harmonized_name="{}"]',
                )
            )

        def parse_record(xml):
            data = {}
            tree = ET.fromstring(xml)
            for attrib, key, xp in xp_tups:
                e = tree.find(xp.format(attrib))
                if e is not None:
                    name = e.get(key)
                    attribute = e.text
                    data[name] = attribute
            self.data.append(pd.DataFrame(data, index=[data["BioSample"]]))

        web_env = self.esearch_results["WebEnv"]
        query_key = self.esearch_results["QueryKey"]
        if self.sample:
            count = self.sample
            batch_size = self.sample
            group = range(0, count, batch_size)
        else:
            count = int(self.esearch_results["Count"])
            batch_size = 10000
            group = range(0, count, batch_size)
        # Return a list of handles for multiprocessing
        for start in group:
            end = min(count, start + batch_size)
            with Entrez.efetch(
                db="biosample",
                rettype="docsum",
                webenv=web_env,
                query_key=query_key,
                retstart=start,
                retmax=batch_size,
            ) as handle:
                try:
                    print("Downloading records {} to {}".format(start + 1, end))
                    efetch_record = Entrez.read(handle, validate=False)
                    for xml in efetch_record["DocumentSummarySet"]["DocumentSummary"]:
                        xml = xml["SampleData"]
                        parse_record(xml)
                except Entrez.Parser.CorruptedXMLError:
                    continue
                    # log here

    @property
    def sra_ids(self):
        ids = self.df[self.df.SRA.notnull()].SRA.tolist()
        return ids

    def _DataFrame(self):
        self.df = pd.concat(self.data)
        self.df.set_index("BioSample", inplace=True)
        self.paths.raw = self.outdir / "_biosample_raw.csv"
        self.df.to_csv(self.paths.raw)

    def generate(self):
        self._esearch()
        self._efetch()
        self._DataFrame()
        self.paths.sra_ids = self.outdir / "sra_ids.txt"
        with open(self.paths.sra_ids, "w") as f:
            f.write("\n".join(self.sra_ids))

    def with_runs(self):
        self.generate()
        summary = AssemblySummary(self.outdir)
        # Add accession summary accession IDs to BioSample DataFrame
        summary.df.reset_index(inplace=True)
        ids = summary.df[["biosample", "# assembly_accession"]]
        ids.set_index("biosample", inplace=True)
        self.df = self.df.join(ids)
        self.paths.csv = self.outdir / "biosample.csv"
        self.df.to_csv(self.paths.csv)

    def read(self):
        return pd.read_csv(self.paths.root / "biosample.csv", index_col=0)


@attr.s
class SRA:
    """Runs from the SRA database"""

    path = attr.ib()

    def __attrs_post_init__(self):
        self.paths = config.Paths(root=self.path)
        self.paths.runs = self.path / "sra_runs.tsv"
        self.id_files = list(self.path.glob("_sra_ids_*"))
        self.runs = pd.read_csv(
            self.paths.runs,
            index_col=0,
            names=["biosample", "runs"],
            sep="\t",
            error_bad_lines=False,
            warn_bad_lines=False,
        )


@attr.s
class Metadata:
    """Methods for accessing and joining all the available Metadata we are interested in."""

    path = attr.ib(converter=Path)
    email = attr.ib()
    sample = attr.ib(default=False)

    def __attrs_post_init__(self):
        self.csv = self.path / "metadata.csv"

    def update(self):
        self.summary = AssemblySummary()
        self.biosample = BioSample(
            email=self.email, outdir=self.path, sample=self.sample
        )
        self.biosample.generate()
        subprocess.run(["bash", "./scripts/efetch_sra_runs.sh", f"{self.path}/"])
        self.sra = SRA(self.path)

    def join_all(self):
        self.summary.df.reset_index(inplace=True)
        accession_ids = self.summary.df[["biosample", "# assembly_accession"]]
        accession_ids.set_index("biosample", inplace=True)
        self.df = self.biosample.df.join(accession_ids)
        self.df = self.biosample.df.join([self.sra.runs])
        self.df.to_csv(self.csv)
