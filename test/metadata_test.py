import os
import shutil
import tempfile
from pathlib import Path

import pytest
import pandas as pd
from genbankqc import metadata


assembly_summary = Path("test/resources/metadata/assembly_summary.txt")


def test_existing_assembly_summary():
    summary = metadata.AssemblySummary(assembly_summary, read=True)
    assert isinstance(summary.df, pd.DataFrame)


def test_download_assembly_summary():
    summary = metadata.AssemblySummary(tempfile.mkdtemp())
    assert os.path.isfile(summary.file_.as_posix())
    assert isinstance(summary.df, pd.DataFrame)


@pytest.fixture()
def biosample():
    temp = Path(tempfile.mkdtemp())
    biosample = metadata.BioSample("inbox.asanchez@gmail.com", temp)
    yield biosample
    shutil.rmtree(temp)


def test_biosample(biosample):
    biosample.generate()
    assert biosample.paths.csv.is_file()
    assert len(list(biosample.paths.sra_ids.iterdir())) > 5
    assert (biosample.paths.sra_ids / "sra_ids_0.txt").is_file()
