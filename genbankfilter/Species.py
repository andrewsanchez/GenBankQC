import os

import pandas as pd

from ete3 import Tree
from genbankfilter.Genome import Genome


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
        else:
            self.stats = None
        if os.path.isfile(nw_file):
            self.tree = Tree(nw_file, 1)
        else:
            self.tree = None
        # TODO: Throw error here if dmx.index and stats.index
        if os.path.isfile(dmx):
            self.dmx = pd.read_csv(dmx, index_col=0, sep="\t")
        else:
            self.tree = None
        # self.genomes = (Genome(os.path.join(self.species_dir, f)) for
        #         f in os.listdir(self.species_dir) if f.endswith('fasta'))

    def genomes(self, ext="fasta"):
        # TODO: Maybe this should return a tuple (genome-path, genome-id)
        """Returns a generator for every file ending with `ext`

        :param ext: File extension of genomes in species directory
        :returns: Path for all genomes in species directory
        :rtype: generator
        """
        # for f in os.listdir(self.species_dir):
        #     if f.endswith(ext):
        #         yield Genome(os.path.join(self.species_dir, f))

        genomes = (Genome(os.path.join(self.species_dir, f)) for
                   f in os.listdir(self.species_dir) if f.endswith(ext))
        return genomes

    def get_stats(self):
        pass
