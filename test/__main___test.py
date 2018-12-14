from click.testing import CliRunner

from genbankqc.__main__ import cli


def test_cli_noargs():
    runner = CliRunner()
    # Should only display help
    result = runner.invoke(cli)
    assert result.exit_code == 0


def test_cli_nosubcommand(genbank):
    runner = CliRunner()
    result = runner.invoke(cli, [genbank.root.as_posix()])
    assert result.exit_code == 0


def test_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0


def test_genbank():
    pass


def test_metadata():
    pass


def test_species(genbank, aphidicola):
    runner = CliRunner()
    result = runner.invoke(cli, [genbank.root.as_posix(), "species", aphidicola.path])
    assert result.exit_code == 0


def test_genome(genbank, genome):
    genome, handler = genome
    runner = CliRunner()
    result = runner.invoke(cli, [genbank.root.as_posix(), "genome", genome.path])
    assert result.exit_code == 0
