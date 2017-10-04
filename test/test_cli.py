import os
import shutil
import tempfile
import unittest

from click.testing import CliRunner

from genbank_qc.__main__ import cli


class TestCli(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('test/resources/',
                        self.genbank)
        self.species = 'Buchnera_aphidicola'
        self.species_dir = os.path.join(self.genbank, self.species)

    def test_levels(self):
        result = self.runner.invoke(cli, ['-n', 1000, '-m', 3, '-s', 1, '-c',
                                          2.0, '-d', self.species_dir])
        print(result.output)
        print(result.exit_code)
        self.assertEqual(result.exit_code, 0)
        self.assertIn(str(1000), result.output)
        self.assertIn(str(3.0), result.output)
        self.assertIn(str(1.0), result.output)
        self.assertIn(str(2.0), result.output)

    def tearDown(self):
        shutil.rmtree(self.tmp)
