import tempfile
import unittest
import os
import shutil
import pandas as pd

from genbankfilter import get_resources
from genbankfilter import config
from genbankfilter import curate
from genbankfilter import mash


class TestCurate(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        self.assembly_summary = pd.read_csv(
            'genbankfilter/test/resources/updated_assembly_summary.txt',
            sep="\t",
            index_col=0)
        self.species_dir = os.path.join(self.tmp, 'genbank',
                                        'Francisella_tularensis')
        self.test_fasta = 'genbankfilter/test/resources/GCA_000009245.1_Francisella_tularensis_holarctica_LVS_Complete_Genome.fasta'
        self.mash_file = os.path.join(self.species_dir, 'GCA_000009245.1.msh')
        self.genome = 'GCA_000009245.1'
        shutil.copytree('genbankfilter/test/resources/', self.genbank)

    def test_sketch(self):
        mash.sketch(self.genbank, self.assembly_summary, self.genome)
        print(os.listdir(self.genbank))
        print(os.listdir(self.species_dir))
        self.assertTrue(os.path.isfile(self.mash_file))

    def tearDown(self):
        shutil.rmtree(self.genbank)


if __name__ == '__main__':
    unittest.main()
