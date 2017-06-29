import tempfile
import unittest
import os
import shutil
import pandas as pd

from genbankfilter import get_resources
from genbankfilter import config
from genbankfilter import curate
from genbankfilter import mash


class TestMash(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        self.assembly_summary = pd.read_csv(
            'genbankfilter/test/resources/updated_assembly_summary.txt',
            sep="\t",
            index_col=0)
        self.species = 'Buchnera_aphidicola'
        self.species_dir = os.path.join(self.tmp, 'genbank',
                                        'Buchnera_aphidicola')
        self.genomes = ['GCA_000007365.1', 'GCA_000007725.1']
        shutil.copytree('genbankfilter/test/resources/', self.genbank)

    def test_mash(self):
        for genome in self.genomes:
            sketch_file = mash.sketch(self.genbank, self.assembly_summary, genome)
            self.assertTrue(os.path.isfile(sketch_file))
        master_sketch = mash.paste(self.genbank, self.assembly_summary, self.species)
        dist_matrix = mash.dist(self.genbank,
                                self.assembly_summary,
                                self.species)
        self.assertTrue(os.path.isfile(master_sketch))
        self.assertTrue(os.path.isfile(dist_matrix))

    def tearDown(self):
        shutil.rmtree(self.genbank)


if __name__ == '__main__':
    unittest.main()
