from click.testing import CliRunner

from genbankqc.__main__ import cli


def test_cli():
    runner = CliRunner()
    result = runner.invoke(cli)
    assert result.exit_code == 0


def test_help():
    runner = CliRunner()
    result = runner.invoke(cli, ['--help'])
    assert result.exit_code == 0


def test_species(genbank, aphidicola):
    runner = CliRunner()
    result = runner.invoke(cli, [genbank.root, 'species', aphidicola.path])
    assert result.exit_code == 0


def test_genome(genbank, genome):
    genome, handler = genome
    runner = CliRunner()
    result = runner.invoke(cli, [genbank.root, 'genome', genome.path])
    assert result.exit_code == 0
