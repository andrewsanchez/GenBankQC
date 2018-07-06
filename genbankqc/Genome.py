import os.path
import re
import subprocess

from retrying import retry
import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError
from collections import defaultdict
from genbankqc import Metadata

import pandas as pd
from Bio import SeqIO


class Genome:
    def __init__(self, genome, assembly_summary=None):
        """
        :param genome: Path to genome
        :returns: Path to genome and name of the genome
        :rtype:
        """
        self.path = genome
        self.basename = os.path.splitext(self.path)[0]
        self.name = self.basename.split('/')[-1]
        try:
            self.accession_id = re.match('GCA_.*\.\d', self.name).group()
        except AttributeError:
            self.accession_id = self.name
        self.metadata = defaultdict(
            lambda: 'missing',
            accession=self.accession_id,
        )
        if '/' not in self.path:
            self.species_dir = '.'
        else:
            self.species_dir = os.path.split(self.path)[0]
        self.qc_dir = os.path.join(self.species_dir, "qc")
        if assembly_summary is not None:
            self.assembly_summary = assembly_summary
            self.metadata["biosample_id"] = assembly_summary.loc[
                self.accession_id].biosample
        self.xml = {}
        self.msh = os.path.join(self.qc_dir, self.name + ".msh")
        self.stats_path = os.path.join(self.qc_dir, self.name + '.csv')
        if os.path.isfile(self.stats_path):
            self.stats_df = pd.read_csv(self.stats_path, index_col=0)
        else:
            self.stats_df = None
        # TODO: Maybe include the species_mean_distance here

    @retry(stop_max_attempt_number=7, stop_max_delay=10000, wait_fixed=2000)
    def efetch(self, db):
        # TODO: Figure out how to retry if the process doesn't complete
        # after 20 seconds
        if db == "biosample":
            db_id = db + "_id"
        elif db == "sra":
            db_id = db + "_id"
        cmd = ("esearch -db {} -query {} | "
               "efetch -format docsum".format(db, self.metadata[db_id]))
        p = subprocess.Popen(cmd, shell="True", stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL)
        xml, err = p.communicate()
        # Maybe this should return something else of xml is empty
        # or if err is not 0
        self.xml[db] = xml

    def parse_biosample(self):
        """
        Get what we need to get out of the xml returned by efetch("biosample")
        Including the SRA ID and fields of interest as defined in
        Metadata.biosample_fields
        """
        try:
            tree = ET.fromstring(self.xml["biosample"])
            sra = tree.find('DocumentSummary/SampleData/'
                            'BioSample/Ids/Id/[@db="SRA"]')
            try:
                self.metadata["sra_id"] = sra.text
            except AttributeError:
                self.metadata["sra_id"] = "missing"
            for name in Metadata.Metadata.biosample_fields:
                xp = ('DocumentSummary/SampleData/BioSample/Attributes/'
                      'Attribute/[@harmonized_name="{}"]'.format(name))
                attrib = tree.find(xp)
                try:
                    self.metadata[name] = attrib.text
                except AttributeError:
                    self.metadata[name] = "missing"
        except ParseError:
            self.metadata["sra_id"] = "missing"

    def parse_sra(self):
        # before running, make sure xml exists
        try:
            tree = ET.fromstring(self.xml["sra"])
            elements = tree.iterfind("DocumentSummary/Runs/Run/[@acc]")
            srs_accessions = []
            for el in elements:
                    items = el.items()
                    acc = [i[1] for i in items if i[0] == 'acc']
                    acc = acc[0]
                    srs_accessions.append(acc)
            self.metadata["srs_accessions"] = ','.join(srs_accessions)
        except ParseError:
            self.metadata["srs_accessions"] = "missing"

    def get_contigs(self):
        """Return a list of of Bio.Seq.Seq objects for fasta and calculate
        the total the number of contigs.
        """
        try:
            self.contigs = [seq.seq for seq in SeqIO.parse(self.path, "fasta")]
            self.count_contigs = len(self.contigs)
        except UnicodeDecodeError:
            self.contigs = UnicodeDecodeError

    def get_assembly_size(self):
        """Calculate the sum of all contig lengths"""
        # TODO: map or reduce might be more elegant here
        self.assembly_size = sum((len(str(seq)) for seq in self.contigs))

    def get_unknowns(self):
        """Count the number of unknown bases, i.e. not [ATCG]"""
        # TODO: Would it be useful to allow the user to define p?
        p = re.compile("[^ATCG]")
        self.unknowns = sum((len(re.findall(p, str(seq)))
                             for seq in self.contigs))

    def get_distance(self, dmx_mean):
        self.distance = dmx_mean.loc[self.name]

    def sketch(self):
        cmd = "mash sketch '{}' -o '{}'".format(self.path, self.msh)
        if not os.path.isfile(self.msh):
            subprocess.Popen(
                cmd, shell="True", stderr=subprocess.DEVNULL).wait()

    def get_stats(self, dmx_mean):
        from pandas import DataFrame
        if not os.path.isfile(self.stats_path):
            self.get_contigs()
            self.get_assembly_size()
            self.get_unknowns()
            self.get_distance(dmx_mean)
            data = {"contigs": self.count_contigs,
                    "assembly_size": self.assembly_size,
                    "unknowns": self.unknowns,
                    "distance": self.distance}
            self.stats = DataFrame(data, index=[self.name])
            self.stats.to_csv(self.stats_path)
