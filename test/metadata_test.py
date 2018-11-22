import tempfile
from pathlib import Path

import pandas as pd
from genbankqc import metadata


def test_existing_assembly_summary():
    path = Path().cwd() / "resources" / ".info" / "assembly_summary.txt"
    summary = metadata.AssemblySummary(path)
    assert isinstance(summary.df, pd.DataFrame)


def test_download_assembly_summary():
    summary = metadata.AssemblySummary(tempfile.mktemp())
    assert isinstance(summary.df, pd.DataFrame)
