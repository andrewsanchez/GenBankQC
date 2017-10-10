from click.testing import CliRunner

from genbank_qc.__main__ import cli


def test_cli(path):
    runner = CliRunner()
    result = runner.invoke(cli, [path])
    assert result.exit_code == 0
