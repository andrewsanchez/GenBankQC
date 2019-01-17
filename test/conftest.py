import os.path
import shutil
import tempfile
from pathlib import Path

import pytest
import pandas as pd
from logbook import TestHandler

from genbankqc import Genome
from genbankqc import Species
from genbankqc import Genbank
from genbankqc import Metadata


assembly_summary = pd.read_csv(
    "test/resources/metadata/assembly_summary.txt", sep="\t", index_col=0
)


@pytest.fixture(scope="module")
def genome(aphidicola):
    handler = TestHandler()
    genome = (
        "GCA_000521565.1_Buchnera_aphidicola_G002_Myzus_persicae_Complete_Genome.fasta"
    )
    genome = os.path.join(aphidicola.path, genome)
    with handler:
        genome = Genome(genome, assembly_summary)
        genome.sketch()
        genome.get_contigs()
        genome.get_assembly_size()
        genome.get_unknowns()
        yield genome, handler


@pytest.fixture(scope="module")
def aphidicola():
    """A species object for Buchnera aphidicola"""
    tmp = tempfile.mkdtemp()
    aphidicola = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree("test/resources/Buchnera_aphidicola", aphidicola)
    yield Species(aphidicola)
    shutil.rmtree(tmp)


@pytest.fixture(scope="module")
def genbank():
    resources = Path("test/resources").absolute()
    tmp = Path(tempfile.mkdtemp())
    for resource in resources.iterdir():
        target = tmp / resource.name
        shutil.copytree(resource, target)
    yield Genbank(tmp)
    shutil.rmtree(tmp)


@pytest.fixture()
def metadata():
    temp = tempfile.mkdtemp()
    metadata = Metadata(temp, "inbox.asanchez@gmail.com", sample=1000)
    yield metadata
    shutil.rmtree(temp)
