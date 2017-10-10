<<<<<<< HEAD
import os
import shutil
import tempfile
import unittest

=======
>>>>>>> develop
from click.testing import CliRunner

from genbank_qc.__main__ import cli


<<<<<<< HEAD
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

    def test_stats_and_filter(self):
        result = self.runner.invoke(cli, [self.species_dir])
        print(os.listdir(self.species_dir))
        self.assertEqual(result.exit_code, 0)

    def tearDown(self):
        shutil.rmtree(self.tmp)
=======
def test_cli(path):
    runner = CliRunner()
    result = runner.invoke(cli, [path])
    assert result.exit_code == 0
>>>>>>> develop
