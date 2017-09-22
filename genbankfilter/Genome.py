import re

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

    def get_contigs(self):
        """ Return a list of of Bio.Seq.Seq objects for fasta and calculate
        the total the number of contigs.
        """
        try:
            self.contigs = [seq.seq for seq in SeqIO.parse(self.path, "fasta")]
        except UnicodeDecodeError:
            self.contigs = UnicodeDecodeError

    def get_assembly_size(self):
        """Calculate the sum of all contig lengths"""
        # TODO: map or reduce might be more elegant here
        self.assembly_size = sum((len(str(seq)) for seq in self.contigs))


    # def get_N_Count(contigs, n_counts):
    #     """
    #     Count the number of unknown bases, i.e. all bases that are not in [ATCG]
    #     """
    #     N_Count = sum([len(re.findall("[^ATCG]", str(seq))) for seq in contigs])
    #     n_counts.append(N_Count)
    #     return N_Count
