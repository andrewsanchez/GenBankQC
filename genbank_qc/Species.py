import os
import re
from subprocess import DEVNULL, Popen

import pandas as pd

from ete3 import Tree


class Species:
    """Represents a collection of genomes in `path`
    :Parameters:
        path : str
            The path to the directory of related genomes you wish to analyze.
    """

    def __init__(self, path,
                 max_unknowns=200,
                 contigs=3.0,
                 assembly_size=3.0,
                 mash=3.0):
        self.path = path
        self.species = path
        self.max_unknowns = max_unknowns
        self.contigs = contigs
        self.assembly_size = assembly_size
        self.mash = mash
        self.qc_dir = os.path.join(self.path, "qc")
        self.label = '{}-{}-{}-{}'.format(
            max_unknowns, contigs, assembly_size, mash)
        self.qc_results_dir = os.path.join(self.qc_dir, self.label)
        self.stats_path = os.path.join(self.qc_dir, 'stats.csv')
        self.nw_path = os.path.join(self.qc_dir, 'tree.nw')
        self.dmx_path = os.path.join(self.qc_dir, 'dmx.csv')
        self.failed_path = os.path.join(self.qc_results_dir, "failed.csv")
        self.tree_img = os.path.join(self.qc_results_dir, "tree.svg")
        self.summary_path = os.path.join(self.qc_results_dir, "summary.txt")
        self.paste_file = os.path.join(self.qc_dir, 'all.msh')
        self.tree = None
        self.stats = None
        self.dmx = None
        if '/' in self.species:
            self.species = path.strip('/').split('/')[-1]
        if not os.path.isdir(self.qc_dir):
            os.mkdir(self.qc_dir)
        if not os.path.isdir(self.qc_results_dir):
            os.mkdir(self.qc_results_dir)
        if os.path.isfile(self.stats_path):
            self.stats = pd.read_csv(self.stats_path, index_col=0)
        if os.path.isfile(self.nw_path):
            self.tree = Tree(self.nw_path, 1)
        if os.path.isfile(self.dmx_path):
            self.dmx = pd.read_csv(self.dmx_path, index_col=0, sep="\t")
        self.criteria = ["unknowns", "contigs", "assembly_size", "distance"]
        self.tolerance = {"unknowns": max_unknowns,
                          "contigs": contigs,
                          "assembly_size": assembly_size,
                          "distance": mash}
        self.passed = self.stats
        self.failed = {}
        self.med_abs_devs = {}
        self.dev_refs = {}
        self.allowed = {"unknowns": max_unknowns}
        # Enable user defined colors
        self.colors = {"unknowns": "red",
                       "contigs": "green",
                       "distance": "purple",
                       "assembly_size": "orange"}
        self.complete = self.assess()

    def __str__(self):
        self.message = [
            "Species: {}".format(self.species), "Tolerance Levels:",
            "Unknown bases:  {}".format(self.max_unknowns),
            "Contigs: {}".format(self.contigs),
            "Assembly Size: {}".format(self.assembly_size),
            "MASH: {}".format(self.mash)]
        return '\n'.join(self.message)

    def assess(self):
        from pandas.util.testing import assert_index_equal
        try:
            assert self.stats is not None
            assert_index_equal(self.genome_ids().sort_values(),
                               self.stats.index.sort_values())
            return True
        except AssertionError:
            return False

    def genomes(self, ext="fasta"):
        # TODO: Maybe this should return a tuple (genome-path, genome-id)
        """Returns a generator for every file ending with `ext`

        :param ext: File extension of genomes in species directory
        :returns: Generator of Genome objects for all genomes in species dir
        :rtype: generator
        """
        from genbank_qc import Genome
        genomes = (Genome(os.path.join(self.path, f)) for
                   f in os.listdir(self.path) if f.endswith(ext))
        return genomes

    def sketches(self):
        return (i.msh for i in self.genomes())

    def genome_ids(self):
        ids = [i.name for i in self.genomes()]
        return pd.Index(ids)

    def sketch(self):
        for genome in self.genomes():
            genome.sketch()

    def mash_paste(self):
        if os.path.isfile(self.paste_file):
            os.remove(self.paste_file)
        sketches = os.path.join(self.qc_dir, "GCA*msh")
        cmd = "mash paste {} {}".format(self.paste_file, sketches)
        Popen(cmd, shell="True", stdout=DEVNULL).wait()
        if not os.path.isfile(self.paste_file):
            self.paste_file = None

    def mash_dist(self):
        cmd = "mash dist -t '{}' '{}' > '{}'".format(
            self.paste_file, self.paste_file, self.dmx_path)
        Popen(cmd, shell="True", stdout=DEVNULL).wait()
        self.dmx = pd.read_csv(self.dmx_path, index_col=0, sep="\t")
        # Make distance matrix more readable
        p = re.compile('.*(GCA_\d+\.\d.*)(.fasta)')
        names = [re.match(p, i).group(1) for i in self.dmx.index]
        self.dmx.index = names
        self.dmx.columns = names
        self.dmx.to_csv(self.dmx_path, sep="\t")

    def run_mash(self):
        """Run all mash related functions."""
        self.sketch()
        self.mash_paste()
        self.mash_dist()

    def assess_tree(self):
        try:
            assert (sorted(self.tree.get_leaf_names()) ==
                    sorted(self.stats.index.tolist()) ==
                    sorted(self.genome_ids().tolist()))
            return True
        except:
            return False

    def get_tree(self):
        if not self.assess_tree():
            import numpy as np
            from skbio.tree import TreeNode
            from scipy.cluster.hierarchy import weighted
            ids = self.dmx.index.tolist()
            triu = np.triu(self.dmx.as_matrix())
            hclust = weighted(triu)
            t = TreeNode.from_linkage_matrix(hclust, ids)
            nw = t.__str__().replace("'", "")
            self.tree = Tree(nw)
            # midpoint root tree
            self.tree.set_outgroup(self.tree.get_midpoint_outgroup())
            self.tree.write(outfile=self.nw_path)

    def get_stats(self):
        """Get stats for all genomes. Concat the results into a DataFrame
        """
        dmx_mean = self.dmx.mean()
        for genome in self.genomes():
            genome.get_stats(dmx_mean)
        species_stats = [genome.stats_df for genome in self.genomes()]
        self.stats = pd.concat(species_stats)
        self.stats.to_csv(self.stats_path)

    def filter_unknown_bases(self):
        """Filter out genomes with too many unknown bases."""
        self.failed["unknowns"] = self.stats.index[
            self.stats["unknowns"] > self.tolerance["unknowns"]]
        self.passed = self.stats.drop(self.failed["unknowns"])

    def filter_contigs(self):
        # Only look at genomes with > 10 contigs to avoid throwing off the
        # median absolute deviation
        # Extract genomes with < 10 contigs to add them back in later.
        eligible_contigs = self.passed.contigs[self.passed.contigs > 10]
        not_enough_contigs = self.passed.contigs[self.passed.contigs <= 10]
        # Median absolute deviation - Average absolute difference between
        # number of contigs and the median for all genomes
        # TODO Define separate function for this
        med_abs_dev = abs(eligible_contigs - eligible_contigs.median()).mean()
        self.med_abs_devs["contigs"] = med_abs_dev
        # Define separate function for this
        # The "deviation reference"
        # Multiply
        dev_ref = med_abs_dev * self.contigs
        self.dev_refs["contigs"] = dev_ref
        self.allowed["contigs"] = eligible_contigs.median() + dev_ref
        # self.passed["contigs"] = eligible_contigs[
        #     abs(eligible_contigs - eligible_contigs.median()) <= dev_ref]
        self.failed["contigs"] = eligible_contigs[
            abs(eligible_contigs - eligible_contigs.median()) > dev_ref].index
        eligible_contigs = eligible_contigs[
            abs(eligible_contigs - eligible_contigs.median()) <= dev_ref]
        # Add genomes with < 10 contigs back in
        eligible_contigs = pd.concat([eligible_contigs, not_enough_contigs])
        # We only need the index of passed genomes at this point
        eligible_contigs = eligible_contigs.index
        self.passed = self.passed.loc[eligible_contigs]

    def filter_med_abs_dev(self, criteria):
        """Filter based on median absolute deviation."""
        # Get the median absolute deviation
        med_abs_dev = abs(self.passed[criteria] -
                          self.passed[criteria].median()).mean()
        dev_ref = med_abs_dev * self.tolerance[criteria]
        self.allowed[criteria] = "{:.4f}".format(self.passed[criteria].median() + dev_ref)
        self.failed[criteria] = self.passed[
            abs(self.passed[criteria] -
                self.passed[criteria].median()) > dev_ref].index
        self.passed = self.passed[
            abs(self.passed[criteria] -
                self.passed[criteria].median()) <= dev_ref]
        # lower = self.passed[criteria].median() - dev_ref
        # upper = self.passed[criteria].median() + dev_ref

    def base_node_style(self):
        from ete3 import NodeStyle, AttrFace
        nstyle = NodeStyle()
        nstyle["shape"] = "sphere"
        nstyle["size"] = 2
        nstyle["fgcolor"] = "black"
        for n in self.tree.traverse():
            n.set_style(nstyle)
            if re.match('^GCA', n.name):
                nf = AttrFace('name', fsize=8)
                nf.margin_right = 100
                nf.margin_left = 3
                n.add_face(nf, column=0)

    def color_clade(self, criteria):
        """Color nodes using ete3 """
        from ete3 import NodeStyle
        for genome in self.failed[criteria]:
            n = self.tree.get_leaves_by_name(genome).pop()
            nstyle = NodeStyle()
            nstyle["fgcolor"] = self.colors[criteria]
            nstyle["size"] = 6
            n.set_style(nstyle)

    # Might be better in a layout function
    def style_and_render_tree(self, file_types=["svg", "pdf"]):
        from ete3 import TreeStyle, TextFace, CircleFace
        ts = TreeStyle()
        title_face = TextFace(self.species, fsize=20)
        ts.title.add_face(title_face, column=0)
        ts.branch_vertical_margin = 10
        # ts.show_leaf_name = False
        # Legend
        # TODO: use ordered dictionary
        for k, v in self.colors.items():
            failures = "Filtered: {}".format(len(self.failed[k]))
            failures = TextFace(failures, fgcolor=v)
            failures.margin_bottom = 5
            tolerance = "Tolerance: {}".format(self.tolerance[k])
            tolerance = TextFace(tolerance, fgcolor=v)
            tolerance.margin_bottom = 5
            allowed = "Allowed: {}".format(self.allowed[k])
            allowed = TextFace(allowed, fgcolor=v)
            allowed.margin_bottom = 5
            f = TextFace(k, fgcolor=v)
            f.margin_bottom = 5
            f.margin_right = 40
            cf = CircleFace(3, v, style="sphere")
            cf.margin_bottom = 5
            cf.margin_right = 5
            ts.legend.add_face(f, column=1)
            ts.legend.add_face(cf, column=2)
            ts.legend.add_face(failures, 1)
            ts.legend.add_face(TextFace(""), 2)
            ts.legend.add_face(tolerance, 1)
            ts.legend.add_face(TextFace(""), 2)
            ts.legend.add_face(allowed, 1)
            ts.legend.add_face(TextFace(""), 2)
        for f in file_types:
            out_tree = os.path.join(self.qc_results_dir, 'tree.{}'.format(f))
            self.tree.render(out_tree, tree_style=ts)

    def color_tree(self):
        self.base_node_style()
        for i in self.criteria:
            self.color_clade(i)
        self.style_and_render_tree()

    def filter(self):
        self.filter_unknown_bases()
        if check_df_len(self.passed, "unknowns"):
            self.filter_contigs()
        if check_df_len(self.passed, "assembly_size"):
            self.filter_med_abs_dev("assembly_size")
        if check_df_len(self.passed, "distance"):
            self.filter_med_abs_dev("distance")
        self.summary()
        self.failed_report()

    def failed_report(self):
        if os.path.isfile(self.failed_path):
            os.remove(self.failed_path)
        with open(self.failed_path, "a") as f:
            for criteria in self.failed.keys():
                ixs = self.failed[criteria]
                for genome in ixs:
                    value = str(self.stats.loc[genome, criteria])
                    f.write('\t'.join([genome, criteria, value]))
                    f.write('\n')

    def summary(self):
        summary = [
            "Filtered genomes",
            "Unknown Bases: {}".format(len(self.failed["unknowns"])),
            "Contigs: {}".format(len(self.failed["contigs"])),
            "Assembly Size: {}".format(len(self.failed["assembly_size"])),
            "MASH: {}".format(len(self.failed["distance"]))]
        summary = '\n'.join(summary)
        with open(os.path.join(self.summary_path), "w") as f:
            f.write(summary)
        return summary


def check_df_len(df, criteria, num=5):
    """
    Verify that df has > than num genomes
    """
    if len(df) > num:
        return True
    else:
        # TODO: Just pass and return false here.
        # info in this print statement will be apparent in summary
        print("Filtering based on {} resulted in less than 5 genomes.")
        return False
