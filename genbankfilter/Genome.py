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
        # TODO: Check if it has MASH file
        self.path = genome
        p = re.compile('.*(GCA_\d+\.\d.*)(.fasta)')
        self.name = re.match(p, genome).group(1)
        self.basename = os.path.splitext(self.path)[0]
        if os.path.isfile(self.basename + ".msh"):
            self.msh = self.basename + ".msh"
        else:
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

    def sketch(self):
        cmd = "mash sketch '{}' -o '{}.msh'".format(self.path, self.basename)
        subprocess.Popen(cmd, shell="True", stdout=subprocess.DEVNULL).wait()
        if os.path.isfile(self.basename + ".msh"):
            self.msh = self.basename + ".msh"

    def get_stats(self):
        from pandas import DataFrame
        self.get_contigs()
        self.get_assembly_size()
        self.get_unknowns()
        data = {"contigs": self.count_contigs,
                "assembly_size": self.assembly_size,
                "unknowns": self.unknowns}
        self.stats = DataFrame(data, index=[self.name])
