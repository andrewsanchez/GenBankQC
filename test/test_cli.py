import unittest
import tempfile
import shutil
import os
import click
from click.testing import CliRunner

from genbankfilter.__main__ import cli

class TestCli(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self.tmp = tempfile.mkdtemp()
        self.genbank = os.path.join(self.tmp, 'genbank')
        shutil.copytree('genbankfilter/test/resources/',
                        self.genbank)
        self.species = 'Buchnera_aphidicola'
        self.species_dir = os.path.join(self.genbank, self.species)

    def test_levels(self):
        level = 2.0
        result = self.runner.invoke(cli, ['--filter-level', level, self.species_dir])
        print(result.exit_code)
        self.assertEqual(result.exit_code, 0)
        self.assertTrue('Contigs:  {}'.format(level) in result.output)
        self.assertTrue('Assembly size:  {}'.format(level) in result.output)
        self.assertTrue('MASH distances:  {}'.format(level) in result.output)

    def tearDown(self):
        shutil.rmtree(self.tmp)
