import os
import shutil
import tempfile
from pathlib import Path

import pytest
import pandas as pd
from genbankqc import AssemblySummary, BioSample, Metadata


def test_existing_assembly_summary():
    summary = AssemblySummary("test/resources/metadata", read=True)
    assert isinstance(summary.df, pd.DataFrame)


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
    biosample.with_runs()
    assert biosample.paths.raw.is_file()
    assert biosample.paths.sra_ids.is_file()
    # Needs to be tested in with_runs()
    assert "# assembly_accession" in biosample.df.columns


@pytest.fixture()
def metadata():
    temp = tempfile.mkdtemp()
    metadata = Metadata(temp, "inbox.asanchez@gmail.com", sample=100)
    yield metadata
    shutil.rmtree(temp)


def test_Metadata(metadata):
    metadata.update()
    assert metadata.biosample.paths.raw.is_file()
    assert metadata.biosample.paths.sra_ids.is_file()
    assert metadata.sra.paths.runs.is_file()
    assert isinstance(metadata.sra.runs, pd.DataFrame)
    assert metadata.sra.id_files
    metadata.join_all()
    assert metadata.csv.is_file()
    assert isinstance(metadata.df, pd.DataFrame)
