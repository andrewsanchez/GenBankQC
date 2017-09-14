import os
import shutil
import unittest
import tempfile
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
        self.assertEqual(self.B_aphidicola.max_unknowns, 200)
        self.assertEqual(self.B_aphidicola.contigs, 3.0)
        self.assertEqual(self.B_aphidicola.assembly_size, 3.0)
        self.assertEqual(self.B_aphidicola.mash, 3.0)
        self.assertEqual(self.B_aphidicola.tolerance["unknowns"],
                         self.B_aphidicola.max_unknowns)
        self.assertEqual(self.B_aphidicola.tolerance["contigs"],
                         self.B_aphidicola.contigs)
        self.assertEqual(self.B_aphidicola.tolerance["Assembly_Size"],
                         self.B_aphidicola.assembly_size)
        self.assertEqual(self.B_aphidicola.tolerance["MASH"],
                         self.B_aphidicola.mash)
        self.assertEqual(self.B_aphidicola.label, "200-3.0-3.0-3.0")
        # Check different values are set properly
        self.B_aphidicola = gbf.FilteredSpecies(self.species_dir, 300,
                                                2.0, 2.0, 2.0)
        self.assertEqual(self.B_aphidicola.max_unknowns, 300)
        self.assertEqual(self.B_aphidicola.contigs, 2.0)
        self.assertEqual(self.B_aphidicola.assembly_size, 2.0)
        self.assertEqual(self.B_aphidicola.mash, 2.0)
        self.assertEqual(self.B_aphidicola.tolerance["unknowns"],
                         self.B_aphidicola.max_unknowns)
        self.assertEqual(self.B_aphidicola.tolerance["contigs"],
                         self.B_aphidicola.contigs)
        self.assertEqual(self.B_aphidicola.tolerance["Assembly_Size"],
                         self.B_aphidicola.assembly_size)
        self.assertEqual(self.B_aphidicola.tolerance["MASH"],
                         self.B_aphidicola.mash)
        self.assertEqual(self.B_aphidicola.label, "300-2.0-2.0-2.0")

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

    def test_filter_contigs(self):
        baumannii = gbf.FilteredSpecies(
            "test/resources/Acinetobacter_baumannii")
        baumannii.passed = baumannii.stats
        baumannii.filter_contigs()
        self.assertEqual(len(baumannii.passed) +
                         len(baumannii.failed["contigs"]),
                         len(baumannii.stats))
        self.assertIsInstance(baumannii.med_abs_devs["contigs"], float)
        self.assertIsInstance(baumannii.dev_refs["contigs"], float)
        self.assertIsInstance(baumannii.allowed["contigs"], float)
        self.assertIsInstance(baumannii.passed, gbf.pd.DataFrame)
        self.assertIsInstance(baumannii.failed["contigs"], gbf.pd.Index)
        self.assertIsInstance(baumannii.allowed["contigs"], float)

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

    def test_filter_med_abs_dev(self):
        baumannii = gbf.FilteredSpecies(
            "test/resources/Acinetobacter_baumannii")
        baumannii.passed = baumannii.stats
        for criteria in ["MASH", "Assembly_Size"]:
            genomes_before_filtering = len(baumannii.passed)
            baumannii.filter_med_abs_dev(criteria)
            self.assertIsInstance(baumannii.passed, gbf.pd.DataFrame)
            self.assertIsInstance(baumannii.failed[criteria],
                                  gbf.pd.Index)
            self.assertEqual(len(baumannii.passed) +
                             len(baumannii.failed[criteria]),
                             genomes_before_filtering)

    def test_filter_all(self):
        import subprocess
        baumannii = gbf.FilteredSpecies(
            "test/resources/Acinetobacter_baumannii")
        gbf.filter_all(baumannii)
        print(baumannii)
        print(baumannii.summary())
        print(os.listdir(baumannii.species_dir))
        tree_svg = os.path.join(baumannii.species_dir, "tree_200-3.0-3.0-3.0.svg")
        shutil.move(tree_svg, "/Users/andrew/scratch/test_tree.svg")
        tree_svg = "/Users/andrew/scratch/test_tree.svg"
        subprocess.Popen("open {}".format(tree_svg), shell=True)

    def tearDown(self):
        shutil.rmtree(self.genbank)


if __name__ == '__main__':
    unittest.main()
