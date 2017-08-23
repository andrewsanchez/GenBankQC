import unittest
import os
import genbankfilter.filter as gbf


class TestFilter(unittest.TestCase):
    def setUp(self):
        self.species_dir = "genbankfilter/test/resources/Buchnera_aphidicola"
        self.fastas = (f for f in os.listdir(self.species_dir))

    def test_stats_function(self):
        from Bio.Seq import Seq
        fasta = os.path.join(self.species_dir, next(self.fastas))
        contigs, contig_count = gbf.get_contigs(fasta)
        self.assertTrue(type(contigs[0]) == Seq)
        self.assertTrue(type(contig_count) == int)
        assembly_size = gbf.get_assembly_size(contigs)
        self.assertTrue(type(assembly_size) == int)
        self.assertTrue(type(assembly_size) != 0)


if __name__ == '__main__':
    unittest.main()
