import os
import pandas as pd
from ete3 import Tree


class Species:
    """Represents a collection of genomes in `species_dir`
    :Parameters:
        species_dir : str
            The path to the directory of related genomes you wish to analyze.
    """

    def __init__(self, species_dir):
        self.species_dir = species_dir
        stats = os.path.join(self.species_dir, 'stats.csv')
        if os.path.isfile(stats):
            self.stats = pd.read_csv(stats, index_col=0)
        nw_file = os.path.join(species_dir, 'tree.nw')
        if os.path.isfile(nw_file):
            self.tree = Tree(nw_file, 1)
