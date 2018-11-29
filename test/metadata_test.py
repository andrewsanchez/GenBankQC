import os
import tempfile
from pathlib import Path

import pandas as pd
from genbankqc import metadata


assembly_summary = Path("test/resources/metadata/assembly_summary.txt")


def test_existing_assembly_summary():
    summary = metadata.AssemblySummary(assembly_summary, read=True)
    assert isinstance(summary.df, pd.DataFrame)


def test_download_assembly_summary():
    summary = metadata.AssemblySummary(tempfile.mkdtemp())
    assert os.path.isfile(summary.outfile.as_posix())
    assert isinstance(summary.df, pd.DataFrame)
