import unittest
import tempfile
import shutil
import click
from click.testing import CliRunner

from genbankfilter.__main__ import cli

class TestCli(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()
        self.tempdir = tempfile.mkdtemp(prefix='gbfilter-cli-test-temp-') 

    def test_levels(self):
        level = 2.0
        result = self.runner.invoke(cli, ['--filter-level', level, self.tempdir])
        self.assertEqual(result.exit_code, 0)
        self.assertTrue('Contigs:  {}'.format(level) in result.output)
        self.assertTrue('Assembly size:  {}'.format(level) in result.output)
        self.assertTrue('MASH distances:  {}'.format(level) in result.output)

    def tearDown(self):
        shutil.rmtree(self.tempdir)
