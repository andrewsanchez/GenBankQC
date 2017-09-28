import os
import shutil
import tempfile
import unittest

from pandas import DataFrame, Index

import genbankfilter.filter as gbf


def test_init(aphidicola_multi):
    params, aphidicola = aphidicola_multi
    a, b, c, d = params
    assert aphidicola.max_unknowns == a
    assert aphidicola.contigs == b
    assert aphidicola.assembly_size == c
    assert aphidicola.mash == d
    assert aphidicola.tolerance["unknowns"] == a
    assert aphidicola.tolerance["contigs"] == b
    assert aphidicola.tolerance["Assembly_Size"] == c
    assert aphidicola.tolerance["MASH"] == d
    assert aphidicola.label == "-".join(map(str, params))
    assert id(aphidicola.stats) == id(aphidicola.passed)


def test_filter_unknowns(unknowns):
    aphidicola, expected_failures = unknowns
    aphidicola.filter_unknown_bases()
    passed_and_failed = sum(map(len, [aphidicola.failed["unknowns"],
                                      aphidicola.passed]))
    assert len(aphidicola.stats) == passed_and_failed
    assert isinstance(aphidicola.passed, DataFrame)
    assert isinstance(aphidicola.failed["unknowns"], Index)
    assert(id(aphidicola.stats) != id(aphidicola.passed))
    assert expected_failures == aphidicola.failed["unknowns"].tolist()


def test_filter_contigs(species):
    species.filter_contigs()
    total_genomes = len(species.passed) + len(species.failed["contigs"])
    assert total_genomes == len(species.stats)
    assert isinstance(species.med_abs_devs["contigs"], float)
    assert isinstance(species.dev_refs["contigs"], float)
    assert isinstance(species.failed["contigs"], Index)
    assert isinstance(species.allowed["contigs"], float)
    assert isinstance(species.passed, DataFrame)


def test_filter_med_abs_dev(species):
    for criteria in ["MASH", "Assembly_Size"]:
        genomes_before_filtering = len(species.passed)
        species.filter_med_abs_dev(criteria)
        assert type(species.passed) == DataFrame
        assert type(species.failed[criteria]) == Index
        passed_and_failed = sum(map(len, [species.failed[criteria],
                                          species.passed]))
        assert passed_and_failed == genomes_before_filtering


class TestFilteredSpecies(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('test/resources/', self.genbank)
        self.species = 'Buchnera_aphidicola'
        self.path = os.path.join(self.genbank, self.species)
        self.B_aphidicola = gbf.FilteredSpecies(self.path)

    def test_str(self):
        print(self.B_aphidicola)

    def test_summary(self):
        pass

    def test_filter_all(self):
        import subprocess
        species_dir = os.path.join(self.genbank, "Acinetobacter_baumannii")
        baumannii = gbf.FilteredSpecies(species_dir)
        gbf.filter_all(baumannii)
        tree_svg = os.path.join(baumannii.species_dir, "tree_200-3.0-3.0-3.0.svg")
        home = subprocess.Popen(
            "echo $HOME",
            shell=True,
            stdout=subprocess.PIPE).communicate()[0].decode().strip()
        shutil.move(tree_svg, "{}/scratch/test_tree.svg".format(home))
        # tree_svg = "/Users/andrew/scratch/test_tree.svg"
        # subprocess.Popen("open {}".format(tree_svg), shell=True)

    def tearDown(self):
        shutil.rmtree(self.genbank)
