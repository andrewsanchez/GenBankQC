import os.path
import re
import subprocess

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
            self.species_dir = self.path
        else:
            self.species_dir = '/'.join(self.path.split('/')[:-1])
        self.qc_dir = os.path.join(self.species_dir, "qc")
        if not os.path.isdir(self.qc_dir):
            os.mkdir(self.qc_dir)
        self.msh = os.path.join(self.qc_dir, self.name + ".msh")
        if not os.path.isfile(self.msh):
            self.msh = None

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
        dst = os.path.join(self.qc_dir, self.name + '.msh')
        cmd = "mash sketch '{}' -o '{}'".format(self.path, dst)
        subprocess.Popen(cmd, shell="True", stdout=subprocess.DEVNULL).wait()
        self.msh = dst
        if not os.path.isfile(dst):
            self.msh = None

    def get_stats(self):
        from pandas import DataFrame
        self.get_contigs()
        self.get_assembly_size()
        self.get_unknowns()
        data = {"contigs": self.count_contigs,
                "assembly_size": self.assembly_size,
                "unknowns": self.unknowns}
        self.stats = DataFrame(data, index=[self.name])
        dst = os.path.join(self.qc_dir, self.name+'.csv')
        self.stats.to_csv(dst, index=0)
