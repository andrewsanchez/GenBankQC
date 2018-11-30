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


def test_metadata(genbank):
    genbank.metadata()
    for species in genbank.species():
        assert os.path.isfile(species.metadata_path)

# def test_genbank_bare(genbank_bare):
#     genbank = genbank_bare
#     assert os.path.isdir(genbank.path)
#     assert os.path.isdir(genbank.info_dir)
#     assert os.path.isfile(genbank.assembly_summary_path)
#     assert isinstance(genbank.assembly_summary, DataFrame)
