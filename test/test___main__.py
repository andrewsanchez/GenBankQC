from click.testing import CliRunner

from genbankqc.__main__ import cli


def test_cli(aphidicola_bare):
    runner = CliRunner()
    result = runner.invoke(cli)
    assert result.exit_code == 0


def test_help():
    runner = CliRunner()
    result = runner.invoke(cli, ['--help'])
    assert result.exit_code == 0


def test_species(genbank, aphidicola):
    runner = CliRunner()
    result = runner.invoke(
        cli, [
            genbank.path,
            'species',
            aphidicola.path,
            ]
    )
    assert result.exit_code == 0


def test_genome(genbank, genome):
    genome, handler = genome
    runner = CliRunner()
    result = runner.invoke(
        cli, [
            genbank.path,
            'genome',
            genome.path,
            ]
    )
    assert result.exit_code == 0
