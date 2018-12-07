import shutil
import pytest
import tempfile
from pathlib import Path

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
