import os
import shutil
import pytest
import tempfile

from genbankqc import Genbank
from genbankqc import Species


def test_genbank_init(genbank):
    assert isinstance(genbank, Genbank)
    species = genbank.species
    for i in species:
        assert isinstance(i, Species)
