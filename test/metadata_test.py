import os
import shutil
import tempfile
from pathlib import Path

import pytest
import pandas as pd
from genbankqc import metadata


def test_existing_assembly_summary():
    summary = metadata.AssemblySummary("test/resources/metadata", read=True)
    assert isinstance(summary.df, pd.DataFrame)


def test_download_assembly_summary():
    summary = metadata.AssemblySummary(tempfile.mkdtemp())
    assert os.path.isfile(summary.file_.as_posix())
    assert isinstance(summary.df, pd.DataFrame)


@pytest.fixture()
def biosample():
    temp = Path(tempfile.mkdtemp())
    biosample = metadata.BioSample(temp, "inbox.asanchez@gmail.com", sample=100)
    yield biosample
    shutil.rmtree(temp)


def test_biosample(biosample):
    biosample.generate()
    assert biosample.paths.csv.is_file()
    assert (biosample.outdir / "sra_ids.txt").is_file()
