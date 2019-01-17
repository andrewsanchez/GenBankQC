import os
import shutil
import pytest
import tempfile
from pathlib import Path

import pandas as pd

from genbankqc import Genbank
from genbankqc import Species


@pytest.fixture()
def genbank_bare():
    temp_dir = Path(tempfile.mkdtemp())
    yield Genbank(temp_dir)
    shutil.rmtree(temp_dir)


def test_genbank_init(genbank):
    assert isinstance(genbank, Genbank)
    for i in genbank.species():
        assert isinstance(i, Species)


@pytest.mark.skip(reason="Fails because our test genomes aren't in first 1000k samples")
def test_species_metadata(genbank):
    metadata = genbank.metadata(email="inbox.asanchez@gmail.com", sample=1000)
    genbank.species_metadata(metadata)
    for species in genbank.species():
        assert isinstance(species.metadata, pd.DataFrame)
        assert os.path.isfile(species.metadata_path)
