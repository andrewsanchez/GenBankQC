import os
import tempfile
from pathlib import Path

import pandas as pd
from genbankqc import metadata


assembly_summary = Path("test/resources/metadata/assembly_summary.txt")


def test_existing_assembly_summary():
    summary = metadata.AssemblySummary(assembly_summary, latest=False)
    assert isinstance(summary.df, pd.DataFrame)


def test_download_assembly_summary():
    summary = metadata.AssemblySummary(tempfile.mktemp())
    assert os.path.isfile(summary.path)
    assert isinstance(summary.df, pd.DataFrame)
