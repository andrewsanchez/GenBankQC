import os
import shutil
import unittest
import tempfile
import subprocess
import genbankfilter.filter as gbf


class TestFilter(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('test/resources/', self.genbank)
        self.species = 'Buchnera_aphidicola'
        self.species_dir = os.path.join(self.genbank, self.species)
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
        self.assertTrue(type(contigs[0]) == Seq)
        self.assertTrue(type(contig_count) == int)
        assembly_size = gbf.get_assembly_size(contigs, [])
        self.assertTrue(type(assembly_size) == int)
        self.assertTrue(type(assembly_size) != 0)
        N_count = gbf.get_N_count(contigs, [])
        self.assertTrue(type(N_count == int))

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

    def test_style_and_render_tree(self):
        from glob import glob
        file_types = ["*png", "*svg"]
        tree = gbf.read_nw_tree(self.nw_file)
        gbf.style_and_render_tree(self.species_dir, tree, self.filter_ranges)
        print(os.listdir(self.species_dir))
        for f in file_types:
            img = os.path.join(self.species_dir, f)
            p = glob(img)
            print(p)
            self.assertTrue(os.path.isfile(p[0]))
            self.assertTrue(len(p) == 1)

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

    def test_init(self):
        self.assertEqual(type(self.B_aphidicola.stats), gbf.pd.DataFrame)
        self.assertEqual(type(self.B_aphidicola.tree), gbf.Tree)
        # Check default values
        self.assertEqual(self.B_aphidicola.max_unknowns["tolerance"], 200)
        self.assertEqual(self.B_aphidicola.contigs["tolerance"], 3.0)
        self.assertEqual(self.B_aphidicola.assembly_size["tolerance"], 3.0)
        self.assertEqual(self.B_aphidicola.mash["tolerance"], 3.0)
        self.assertEqual(self.B_aphidicola.label, "200-3.0-3.0-3.0")
        # Check different values are set properly
        self.B_aphidicola = gbf.FilteredSpecies(self.species_dir, 300,
                                                2.0, 2.0, 2.0)
        self.assertEqual(self.B_aphidicola.max_unknowns["tolerance"], 300)
        self.assertEqual(self.B_aphidicola.contigs["tolerance"], 2.0)
        self.assertEqual(self.B_aphidicola.assembly_size["tolerance"], 2.0)
        self.assertEqual(self.B_aphidicola.mash["tolerance"], 2.0)
        self.assertEqual(self.B_aphidicola.label, "300-2.0-2.0-2.0")

    def test_str(self):
        print(self.B_aphidicola)

    def test_filter_unknown_bases(self):
        self.B_aphidicola.filter_unknown_bases()
        self.assertEqual(self.B_aphidicola.max_unknowns["passed"].tolist(),
                         self.B_aphidicola.passed.index.tolist())
        self.assertIsInstance(self.B_aphidicola.passed, gbf.pd.DataFrame)
        self.assertIsInstance(self.B_aphidicola.max_unknowns["failed"],
                              gbf.pd.Index)
        # Set all rows in column N_Count to 0
        self.B_aphidicola.stats.iloc[:, 0] = 0
        # Make sure the 0th column is in fact N_Count
        self.assertEqual(self.B_aphidicola.stats.iloc[:, 0].name, "N_Count")
        self.B_aphidicola.stats.iloc[:10, 0] = 300
        expected_failures = self.B_aphidicola.stats.iloc[:10, 0].index.tolist()
        self.B_aphidicola.filter_unknown_bases()
        self.assertEqual(self.B_aphidicola.max_unknowns["passed"].tolist(),
                         self.B_aphidicola.passed.index.tolist())
        self.assertEqual(expected_failures,
                         self.B_aphidicola.max_unknowns["failed"].tolist())

    def test_filter_contigs(self):
        baumanii = gbf.FilteredSpecies("test/resources/Acinetobacter_baumanii")
        baumanii.passed = baumanii.stats
        baumanii.filter_contigs()
        self.assertEqual(len(baumanii.contigs["passed"]) +
                         len(baumanii.contigs["failed"]), len(baumanii.stats))
        self.assertIsInstance(baumanii.passed, gbf.pd.DataFrame)
        self.assertIsInstance(baumanii.contigs["passed"], gbf.pd.Index)
        self.assertIsInstance(baumanii.contigs["failed"], gbf.pd.Index)
        self.assertEqual(baumanii.contigs["passed"].tolist(),
                         baumanii.passed.index.tolist())
        self.assertIsInstance(baumanii.contigs["allowed"], float)

    def test_filter_med_ad(self):
        for criteria in ["MASH", "Assembly_Size"]:
            self.B_aphidicola.filter_med_ad(criteria)
            self.assertIsInstance(self.B_aphidicola.passed, gbf.pd.DataFrame)
            self.assertIsInstance(self.B_aphidicola.failed, list)

    def test_filter_all(self):
        gbf._filter_all(self.B_aphidicola)
        tree_svg = os.path.join(self.species_dir, "tree_200-3.0-3.0-3.0.svg")
        shutil.move(tree_svg, "/Users/andrew/scratch/test_tree.svg")
        tree_svg = "/Users/andrew/scratch/test_tree.svg"
        subprocess.Popen("open {}".format(tree_svg), shell=True)

    def tearDown(self):
        shutil.rmtree(self.genbank)


if __name__ == '__main__':
    unittest.main()
