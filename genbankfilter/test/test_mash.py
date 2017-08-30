import tempfile
import unittest
import os
import shutil
import pandas as pd
import genbankfilter.mash as mash


class TestMash(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('genbankfilter/test/resources/', self.genbank)
        # self.assembly_summary = pd.read_csv(
        #     'genbankfilter/test/resources/updated_assembly_summary.txt',
        #     sep="\t",
        #     index_col=0)
        self.species = 'Buchnera_aphidicola'
        self.species_dir = os.path.join(self.genbank, self.species)
        self.genomes = mash.find_all_genome_paths(self.genbank)

    def test_find_all_genome_paths(self):
        for genome in self.genomes:
            self.assertTrue(os.path.isfile(genome))

    def test_sketch_genome(self):
        for genome in self.genomes:
            mash.sketch_genome(genome)
        sketches = mash.find_all_sketches(self.species_dir)
        self.assertTrue(len(sketches) == len(self.genomes))
        for sketch in sketches:
            self.assertTrue(os.path.isfile(sketch))

    def test_sketch_species_dir(self):
        mash.sketch_dir(self.species_dir)
        sketches = mash.find_all_sketches(self.species_dir)
        self.assertTrue(len(sketches) == len(self.genomes))
        for sketch in sketches:
            self.assertTrue(os.path.isfile(sketch))

    def test_sketch_genbank(self):
        mash.sketch_dir(self.genbank)
        sketches = mash.find_all_sketches(self.genbank)
        self.assertTrue(len(sketches) == len(self.genomes))
        for sketch in sketches:
            self.assertTrue(os.path.isfile(sketch))

    def test_paste_and_dist(self):
        mash.sketch_dir(self.genbank)
        all_msh = mash.paste(self.species_dir)
        dst_mx = mash.dist(self.species_dir)
        self.assertTrue(os.path.isfile(all_msh))
        self.assertIsInstance(dst_mx, pd.DataFrame)

    def tearDown(self):
        shutil.rmtree(self.genbank)


if __name__ == '__main__':
    unittest.main()
