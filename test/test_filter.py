import os
import shutil
import tempfile
import unittest

import genbankfilter.filter as gbf


class TestFilter(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('test/resources/', self.genbank)
        self.species = 'Buchnera_aphidicola'
        self.path = os.path.join(self.genbank, self.species)
        self.baumannii = os.path.join(self.genbank, "Acinetobacter_baumannii")
        self.B_aphidicola = gbf.FilteredSpecies(self.path)
        self.fastas = gbf.get_all_fastas(self.path)
        self.stats = gbf.pd.read_csv(
            os.path.join(self.path, 'stats.csv'), index_col=0)
        self.filter_ranges = [200, 3, 3, 3]
        self.criteria_and_franges = gbf.criteria_dict(self.filter_ranges)
        self.nw_file = os.path.join(self.path, 'tree.nw')

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
        dmx = gbf.read_dmx(self.path)
        stats = gbf.generate_stats(self.path, dmx)
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
        dmx = gbf.read_dmx(self.path)
        self.assertTrue(type(dmx) == gbf.pd.DataFrame)

    def test_nested_dmx(self):
        dmx = gbf.read_dmx(self.path)
        nested_mx = gbf.nested_matrix(dmx)
        self.assertTrue((type(nested_mx) == list))
        self.assertTrue(len(nested_mx) != 0)

    def tearDown(self):
        shutil.rmtree(self.genbank)
if __name__ == '__main__':
    unittest.main()
