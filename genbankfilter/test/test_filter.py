import os
import unittest
import genbankfilter.filter as gbf


class TestFilter(unittest.TestCase):
    def setUp(self):
        self.species_dir = "genbankfilter/test/resources/Buchnera_aphidicola"
        self.fastas = gbf.get_all_fastas(self.species_dir)
        stats = os.path.join(self.species_dir, 'stats.csv')
        self.stats = gbf.pd.read_csv(
            os.path.join(self.species_dir, 'stats.csv'), index_col=0)

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


if __name__ == '__main__':
    unittest.main()
