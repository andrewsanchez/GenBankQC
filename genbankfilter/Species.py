import os
import re

import pandas as pd
from Bio import Phylo, SeqIO

from ete3 import Tree


class Species:
    """Represents a collection of genomes in `species_dir`
    :Parameters:
        species_dir : str
            The path to the directory of related genomes you wish to analyze.
    """

    def __init__(self, species_dir):
        self.species_dir = species_dir
        self.species = species_dir
        if '/' in self.species:
            self.species = species_dir.strip('/').split('/')[-1]
        stats = os.path.join(self.species_dir, 'stats.csv')
        nw_file = os.path.join(species_dir, 'tree.nw')
        dmx = os.path.join(species_dir, 'dmx.txt')
        # TODO: What to do when these files don't exist?
        if os.path.isfile(stats):
            self.stats = pd.read_csv(stats, index_col=0)
        if os.path.isfile(nw_file):
            self.tree = Tree(nw_file, 1)
        if os.path.isfile(dmx):
            self.dmx = pd.read_csv(dmx, index_col=0, sep="\t")

    def genomes(self, ext="fasta"):
        # TODO: Maybe this should return a tuple (genome-path, genome-id)
        """Returns a generator for every file ending with `ext`

        :param ext: File extension of genomes in species directory
        :returns: Path for all genomes in species directory
        :rtype: generator
        """
        for f in os.listdir(self.species_dir):
            if f.endswith(ext):
                yield os.path.join(self.species_dir, f)


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
        """
        Return a list of of Bio.Seq.Seq objects for fasta and calculate
        the total the number of contigs.
        """
        try:
            self.contigs = [seq.seq for seq in SeqIO.parse(self.path, "fasta")]
        except UnicodeDecodeError:
            self.contigs = UnicodeDecodeError


