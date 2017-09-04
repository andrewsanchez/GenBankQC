import os
import re
import shutil
import pandas as pd
import numpy as np
from ete3 import Tree
from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

from collections import namedtuple


class Species:
    """Represents a collection of genomes in `species_dir`
    :Parameters:
        species_dir : str
            The path to the directory of related genomes you wish to analyze.
    """

    color_map = {
        "N_Count": "red",
        "Contigs": "green",
        "MASH": "blue",
        "Assembly_Size": "purple"
    }

    def __init__(self, species_dir, max_n_count=200, c_range=3.0, s_range=3.0,
                 m_range=3.0):
        self.species_dir = species_dir
        stats = os.path.join(self.species_dir, 'stats.csv')
        if os.path.isfile(stats):
            self.stats = pd.read_csv(stats, index_col=0)
        nw_file = os.path.join(species_dir, 'tree.nw')
        if os.path.isfile(nw_file):
            self.tree = Tree(nw_file, 1)
        self.max_n_count = max_n_count
        self.c_range = c_range
        self.s_range = s_range
        self.m_range = m_range
        # probably won't need this
        self.filter_ranges = [max_n_count, c_range, s_range, m_range]
