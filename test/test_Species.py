import unittest
import shutil
import os
import tempfile
from genbankfilter.Species import Species


class TestSpecies(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('genbankfilter/test/resources/', self.genbank)
        self.species = 'Buchnera_aphidicola'
        self.species_dir = os.path.join(self.genbank, self.species)
        self.B_aphidicola = Species(self.species_dir)
        self.assertEqual(type(self.B_aphidicola), Species)

    def test_species_init(self):
        from pandas import DataFrame
        from ete3 import Tree
        self.assertEqual(type(self.B_aphidicola.stats), DataFrame)
        self.assertEqual(type(self.B_aphidicola.tree), Tree)

    def tearDown(self):
        shutil.rmtree(self.genbank)


if __name__ == '__main__':
    unittest.main()
