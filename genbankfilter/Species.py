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
        self.qc_dir = os.path.join(self.species_dir, "qc")
        if not os.path.isdir(self.qc_dir):
            os.mkdir(self.qc_dir)
        self.stats = os.path.join(self.qc_dir, 'stats.csv')
        self.nw_file = os.path.join(self.qc_dir, 'tree.nw')
        self.dmx = os.path.join(self.qc_dir, 'dmx.csv')
        if os.path.isfile(self.stats):
            self.stats = pd.read_csv(self.stats, index_col=0)
        else:
            self.stats = None
        if os.path.isfile(self.nw_file):
            self.tree = Tree(self.nw_file, 1)
        else:
            self.tree = None
        # TODO: Throw error here if dmx.index and stats.index
        if os.path.isfile(self.dmx):
            self.dmx = pd.read_csv(self.dmx, index_col=0, sep="\t")
        else:
            self.dmx = None

    def genomes(self, ext="fasta"):
        # TODO: Maybe this should return a tuple (genome-path, genome-id)
        """Returns a generator for every file ending with `ext`

        :param ext: File extension of genomes in species directory
        :returns: Generator of Genome objects for all genomes in species dir
        :rtype: generator
        """
        genomes = (Genome(os.path.join(self.species_dir, f)) for
                   f in os.listdir(self.species_dir) if f.endswith(ext))
        return genomes

    def sketches(self):
        return (i.msh for i in self.genomes())

    def genome_ids(self):
        ids = [i.name for i in self.genomes()]
        return pd.Index(ids)

    def get_stats(self):
        pass

    def sketch(self):
        for genome in self.genomes():
            if genome.msh is None:
                genome.sketch()

    def mash_paste(self):
        self.paste_file = os.path.join(self.qc_dir, 'all.msh')
        sketches = os.path.join(self.qc_dir, "GCA*msh")
        cmd = "mash paste {} {}".format(self.paste_file, sketches)
        Popen(cmd, shell="True", stdout=DEVNULL).wait()
        if not os.path.isfile(self.paste_file):
            self.paste_file = None

    def mash_dist(self):
        import re
        dst = os.path.join(self.qc_dir, 'dmx.csv')
        cmd = "mash dist -t '{}' '{}' > '{}'".format(
            self.paste_file, self.paste_file, dst)
        Popen(cmd, shell="True", stdout=DEVNULL).wait()
        self.dmx = pd.read_csv(dst, index_col=0, sep="\t")
        # Make distance matrix more readable
        p = re.compile('.*(GCA_\d+\.\d.*)(.fasta)')
        names = [re.match(p, i).group(1) for i in self.dmx.index]
        self.dmx.index = names
        self.dmx.columns = names
        self.dmx.to_csv(dst, sep="\t")

    def get_stats(self):
        from pandas import DataFrame
        dmx_mean = self.dmx.mean()
        stats = []
        for i in self.genomes():
            i.get_contigs()
            i.get_assembly_size()
            i.get_unknowns()
            i.get_distance(dmx_mean)
            data = {"contigs": i.count_contigs,
                    "assembly_size": i.assembly_size,
                    "unknowns": i.unknowns,
                    "distance": i.distance}
            i.stats = DataFrame(data, index=[i.name])
            dst = os.path.join(i.qc_dir, i.name+'.csv')
            stats.append(i.stats)
            i.stats.to_csv(dst)
        self.stats = pd.concat(stats)
        self.stats.to_csv(os.path.join(self.qc_dir, 'stats.csv'))

