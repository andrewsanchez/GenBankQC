import os
import shutil
import tempfile
from pathlib import Path

import pytest
import pandas as pd
from genbankqc import AssemblySummary, BioSample


def test_existing_assembly_summary():
    summary = AssemblySummary("test/resources/metadata", update=False)
    assert isinstance(summary.df, pd.DataFrame)


@pytest.mark.xfail(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason="Fails on Travis for unknown reason.",
)
def test_download_assembly_summary():
    summary = AssemblySummary(tempfile.mkdtemp())
    assert os.path.isfile(summary.file_.as_posix())
    assert isinstance(summary.df, pd.DataFrame)


@pytest.fixture()
def biosample():
    temp = Path(tempfile.mkdtemp())
    biosample = BioSample(temp, "inbox.asanchez@gmail.com", sample=100)
    yield biosample
    shutil.rmtree(temp)


def test_biosample(biosample):
    biosample.generate()
    assert biosample.paths.raw.is_file()
    assert biosample.paths.sra_ids.is_file()


@pytest.mark.xfail(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason="Fails on Travis for unknown reason.",
)
def test_Metadata(metadata):
    assert metadata.biosample.paths.raw.is_file()
    assert metadata.biosample.paths.sra_ids.is_file()
    assert metadata.sra.paths.runs.is_file()
    assert metadata.sra.id_files
    assert isinstance(metadata.sra.runs, pd.DataFrame)
    assert metadata.csv.is_file()
    assert isinstance(metadata.joined, pd.DataFrame)
