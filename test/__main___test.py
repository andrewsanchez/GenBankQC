import os

import pytest
from click.testing import CliRunner

from genbankqc.__main__ import cli


def test_cli_noargs():
    runner = CliRunner()
    result = runner.invoke(cli)  # Should only display help
    assert result.exit_code == 0


@pytest.mark.skipif(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason="Fails because the prune function downloads the assembly summary",
)
def test_cli_no_subcommand(genbank):
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
    result = runner.invoke(
        cli, [genbank.root.as_posix(), "species", aphidicola.path.as_posix()]
    )
    assert result.exit_code == 0


def test_genome(genbank, genome):
    genome, handler = genome
    runner = CliRunner()
    result = runner.invoke(cli, [genbank.root.as_posix(), "genome", genome.path])
    assert result.exit_code == 0
