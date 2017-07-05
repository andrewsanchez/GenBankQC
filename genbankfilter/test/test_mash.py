import tempfile
import unittest
import os
import re
import shutil
import pandas as pd

from genbankfilter import get_resources
from genbankfilter import config
from genbankfilter import curate
from genbankfilter import mash
from genbankfilter import filter


class TestMash(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('genbankfilter/test/resources/',
                        self.genbank)
        # self.assembly_summary = pd.read_csv(
        #     'genbankfilter/test/resources/updated_assembly_summary.txt',
        #     sep="\t",
        #     index_col=0)
        self.species = 'Buchnera_aphidicola'
        self.species_dir = os.path.join(self.genbank, self.species)
        self.genomes = mash.find_all_genome_paths(self.genbank)

        for genome in self.genomes:
            mash.sketch(self.genbank, self.assembly_summary, genome)
        self.all_msh = mash.paste(self.genbank, self.assembly_summary, self.species)
        self.dst_mx = mash.dist(self.genbank,
                                self.assembly_summary,
                                self.species)
        self.stats = filter.generate_stats(self.species_dir, self.dst_mx)

    def test_mash_and_stats(self):
        for genome in self.genomes:
            sketch = os.path.join(self.species_dir, '{}.msh'.format(genome))
            self.assertTrue(os.path.isfile(sketch))
        self.assertTrue(os.path.isfile(self.all_msh))
        self.assertIsInstance(self.dst_mx, pd.DataFrame)
        self.assertIsInstance(self.stats, pd.DataFrame)

    def test_filter_and_link(self):
        filter_summary, failed, passed_final = filter.filter_med_ad(self.species_dir, self.stats)
        passed_dir = filter.check_passed_dir(self.species_dir)
        filter.link_passed_genomes(self.species_dir, passed_final, passed_dir)
        self.assertTrue(passed_dir)
        self.assertIsInstance(filter_summary, pd.DataFrame)
    
    def tearDown(self):
        shutil.rmtree(self.genbank)


if __name__ == '__main__':
    unittest.main()
