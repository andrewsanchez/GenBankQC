import os
import shutil
import tempfile
import unittest

import pytest

import genbankfilter.filter as gbf


class TestFilter(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('test/resources/', self.genbank)
        self.species = 'Buchnera_aphidicola'
        self.species_dir = os.path.join(self.genbank, self.species)
        self.baumannii = os.path.join(self.genbank, "Acinetobacter_baumannii")
        self.B_aphidicola = gbf.FilteredSpecies(self.species_dir)
        self.fastas = gbf.get_all_fastas(self.species_dir)
        self.stats = gbf.pd.read_csv(
            os.path.join(self.species_dir, 'stats.csv'), index_col=0)
        self.filter_ranges = [200, 3, 3, 3]
        self.criteria_and_franges = gbf.criteria_dict(self.filter_ranges)
        self.nw_file = os.path.join(self.species_dir, 'tree.nw')

    def test_stats_functions(self):
        from Bio.Seq import Seq
        fasta = next(self.fastas)
        contigs, contig_count = gbf.get_contigs(fasta, [])
        self.assertEqual(type(contigs[0]), Seq)
        self.assertEqual(type(contig_count), int)
        assembly_size = gbf.get_assembly_size(contigs, [])
        self.assertEqual(type(assembly_size), int)
        self.assertNotEqual(type(assembly_size), 0)
        N_count = gbf.get_N_Count(contigs, [])
        self.assertEqual(type(N_count), int)

    def test_generate_stats(self):
        dmx = gbf.read_dmx(self.species_dir)
        stats = gbf.generate_stats(self.species_dir, dmx)
        self.assertTrue(type(stats) == gbf.pd.DataFrame)

    def test_filter_Ns(self):
        passed_N_count, failed_N_count = gbf.filter_Ns(self.stats, 200)
        self.assertTrue(type(passed_N_count) == gbf.pd.DataFrame)
        self.assertTrue(type(failed_N_count) == gbf.pd.DataFrame)

    def test_filter_contigs(self):
        passed = self.stats
        summary = {}
        filter_results = gbf.filter_contigs(self.stats, passed, 3.0, summary)
        passed, failed = filter_results
        self.assertTrue(type(passed) == gbf.pd.DataFrame)
        self.assertTrue(type(failed) == list)

    def test_filter_med_abs_dev(self):
        results = gbf.filter_med_abs_dev(self.stats, {}, "MASH",
                                         self.criteria_and_franges)
        self.assertTrue(type(results.passed) == gbf.pd.DataFrame)

    def test_read_tree(self):
        tree = gbf.read_nw_tree(self.nw_file)
        self.assertTrue(type(tree) == gbf.Tree)

    def test_read_dmx(self):
        dmx = gbf.read_dmx(self.species_dir)
        self.assertTrue(type(dmx) == gbf.pd.DataFrame)

    def test_nested_dmx(self):
        dmx = gbf.read_dmx(self.species_dir)
        nested_mx = gbf.nested_matrix(dmx)
        self.assertTrue((type(nested_mx) == list))
        self.assertTrue(len(nested_mx) != 0)

    def tearDown(self):
        shutil.rmtree(self.genbank)


class TestFilteredSpecies(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('test/resources/', self.genbank)
        self.species = 'Buchnera_aphidicola'
        self.species_dir = os.path.join(self.genbank, self.species)
        self.B_aphidicola = gbf.FilteredSpecies(self.species_dir)

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


def test_FilteredSpecies_init(provide_aphidicola_multi):
    params, aphidicola = provide_aphidicola_multi
    a, b, c, d = params
    assert type(aphidicola.stats) == gbf.pd.DataFrame
    assert type(aphidicola.tree) == gbf.Tree
    # Check default values
    assert aphidicola.max_unknowns == a
    assert aphidicola.contigs == b
    assert aphidicola.assembly_size == c
    assert aphidicola.mash == d
    assert aphidicola.tolerance["unknowns"] == aphidicola.max_unknowns
    assert aphidicola.tolerance["contigs"] == aphidicola.contigs
    assert aphidicola.tolerance["Assembly_Size"] == aphidicola.assembly_size
    assert aphidicola.tolerance["MASH"] == aphidicola.mash
    assert aphidicola.label == "-".join(map(str, params))


def test_filter_contigs(provide_baumannii):
    baumannii = provide_baumannii
    baumannii.filter_contigs()
    total_genomes = len(baumannii.passed) + len(baumannii.failed["contigs"])
    assert total_genomes == len(baumannii.stats)
    assert type(baumannii.med_abs_devs["contigs"]) == gbf.np.float64
    assert type(baumannii.dev_refs["contigs"]) == gbf.np.float64
    assert type(baumannii.allowed["contigs"]) == gbf.np.float64
    assert type(baumannii.passed) == gbf.pd.DataFrame
    assert type(baumannii.failed["contigs"]) == gbf.pd.Index
    assert type(baumannii.allowed["contigs"]) == gbf.np.float64


def test_filter_med_abs_dev(provide_baumannii):
    baumannii = provide_baumannii
    for criteria in ["MASH", "Assembly_Size"]:
        genomes_before_filtering = len(baumannii.passed)
        baumannii.filter_med_abs_dev(criteria)
        assert type(baumannii.passed) == gbf.pd.DataFrame
        assert type(baumannii.failed[criteria]) == gbf.pd.Index
        total_genomes = len(baumannii.passed) + len(baumannii.failed[criteria])
        assert (total_genomes == genomes_before_filtering)


if __name__ == '__main__':
    unittest.main()
