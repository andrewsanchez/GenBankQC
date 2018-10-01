import os.path
import shutil
import tempfile

import pytest
import pandas as pd
from logbook import TestHandler

from genbankqc import Genbank
from genbankqc import Species
from genbankqc import Genome
from genbankqc import Metadata


assembly_summary = pd.read_csv('test/resources/.info/assembly_summary.txt', sep="\t", index_col=0)


@pytest.fixture(scope="module")
def genbank():
    resources = 'test/resources'
    genbank = tempfile.mkdtemp()
    for resource in os.listdir(resources):
        source = os.path.join(resources, resource)
        target = os.path.join(genbank, resource)
        if os.path.isdir(source):
            shutil.copytree(source, target)
    yield Genbank(genbank)
    shutil.rmtree(genbank)


@pytest.fixture(scope="module")
def aphidicola():
    tmp = tempfile.mkdtemp()
    aphidicola = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', aphidicola)
    yield Species(aphidicola)
    shutil.rmtree(tmp)


@pytest.fixture(scope="module")
def metadata(genbank):
    yield Metadata(genbank.path)
