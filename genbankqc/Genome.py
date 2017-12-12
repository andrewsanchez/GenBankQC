import os.path
import re
import subprocess

import pandas as pd
from Bio import SeqIO


class Genome:
    def __init__(self, genome):
        """
        :param genome: Path to genome
        :returns: Path to genome and name of the genome
        :rtype:
        """
        self.path = genome
        p = re.compile('.*(GCA_\d+\.\d.*)(.fasta)')
        self.name = re.match(p, genome).group(1)
        self.basename = os.path.splitext(self.path)[0]
        if '/' not in self.path:
            self.species_dir = '.'
        else:
            self.species_dir = os.path.split(self.path)[0]
        self.qc_dir = os.path.join(self.species_dir, "qc")
        if not os.path.isdir(self.qc_dir):
            os.mkdir(self.qc_dir)
        self.msh = os.path.join(self.qc_dir, self.name + ".msh")
        self.stats_path = os.path.join(self.qc_dir, self.name+'.csv')
        if os.path.isfile(self.stats_path):
            self.stats_df = pd.read_csv(self.stats_path, index_col=0)
        else:
            self.stats_df = None
        # TODO: Maybe include the species_mean_distance here

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
