from click.testing import CliRunner

from genbankqc.__main__ import cli


def test_cli(aphidicola_bare):
    runner = CliRunner()
    result = runner.invoke(cli, ['--species', aphidicola_bare.path])
    assert result.exit_code == 0


def test_help():
    runner = CliRunner()
    result = runner.invoke(cli, ['--help'])
    print(result.output)
    assert result.exit_code == 0
