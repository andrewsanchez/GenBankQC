import os
import shutil
import tempfile
import unittest

from click.testing import CliRunner

from genbank_qc.__main__ import cli


def test_cli(aphidicola_bare):
    aphidicola = aphidicola_bare
    runner = CliRunner()
    result = runner.invoke(cli, [aphidicola.path])
    print(result.output)
    print(aphidicola.path)
    assert result.exit_code == 0
