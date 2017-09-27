import os
import shutil
import tempfile
import unittest

import genbankfilter.filter as gbf


def test_SpeciesQC_init(aphidicola_multi):
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


def test_filter_contigs(baumannii):
    baumannii.filter_contigs()
    total_genomes = len(baumannii.passed) + len(baumannii.failed["contigs"])
    assert total_genomes == len(baumannii.stats)
    assert type(baumannii.med_abs_devs["contigs"]) == gbf.np.float64
    assert type(baumannii.dev_refs["contigs"]) == gbf.np.float64
    assert type(baumannii.allowed["contigs"]) == gbf.np.float64
    assert type(baumannii.passed) == gbf.pd.DataFrame
    assert type(baumannii.failed["contigs"]) == gbf.pd.Index
    assert type(baumannii.allowed["contigs"]) == gbf.np.float64


def test_filter_med_abs_dev(baumannii):
    for criteria in ["MASH", "Assembly_Size"]:
        genomes_before_filtering = len(baumannii.passed)
        baumannii.filter_med_abs_dev(criteria)
        assert type(baumannii.passed) == gbf.pd.DataFrame
        assert type(baumannii.failed[criteria]) == gbf.pd.Index
        total_genomes = len(baumannii.passed) + len(baumannii.failed[criteria])
        assert (total_genomes == genomes_before_filtering)


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

    def test_filter_unknown_bases(self):
        self.B_aphidicola.filter_unknown_bases()
        self.assertIsInstance(self.B_aphidicola.passed, gbf.pd.DataFrame)
        self.assertIsInstance(self.B_aphidicola.failed["unknowns"],
                              gbf.pd.Index)
        # Set all rows in column N_Count to 0
        self.B_aphidicola.stats.iloc[:, 0] = 0
        # Make sure the 0th column is in fact N_Count
        self.assertEqual(self.B_aphidicola.stats.iloc[:, 0].name, "N_Count")
        self.B_aphidicola.stats.iloc[:10, 0] = 300
        expected_failures = self.B_aphidicola.stats.iloc[:10, 0].index.tolist()
        self.B_aphidicola.filter_unknown_bases()
        self.assertNotEqual(id(self.B_aphidicola.stats),
                            id(self.B_aphidicola.passed))
        self.assertIsInstance(self.B_aphidicola.failed["unknowns"],
                              gbf.pd.Index)
        self.assertEqual(len(self.B_aphidicola.stats),
                         len(self.B_aphidicola.passed) +
                         len(self.B_aphidicola.failed["unknowns"]))
        self.assertEqual(expected_failures,
                         self.B_aphidicola.failed["unknowns"].tolist())

    def test_filter_mash(self):
        baumannii = gbf.FilteredSpecies(
            "test/resources/Acinetobacter_baumannii")
        baumannii.passed = baumannii.stats
        baumannii.filter_med_abs_dev("MASH")
        self.assertEqual(len(baumannii.passed) +
                         len(baumannii.failed["MASH"]),
                         len(baumannii.stats))
        self.assertIsInstance(baumannii.passed, gbf.pd.DataFrame)
        self.assertIsInstance(baumannii.failed["MASH"],
                              gbf.pd.Index)

    def test_filter_assembly_size(self):
        baumannii = gbf.FilteredSpecies(
            "test/resources/Acinetobacter_baumannii")
        baumannii.passed = baumannii.stats
        baumannii.filter_med_abs_dev("Assembly_Size")
        self.assertEqual(len(baumannii.passed) +
                         len(baumannii.failed["Assembly_Size"]),
                         len(baumannii.stats))
        self.assertIsInstance(baumannii.passed, gbf.pd.DataFrame)
        self.assertIsInstance(baumannii.failed["Assembly_Size"],
                              gbf.pd.Index)

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
