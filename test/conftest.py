import os.path
import shutil
import tempfile

import pytest
import pandas as pd
from logbook import TestHandler

from genbankqc import Genome
from genbankqc import Species
from genbankqc import Genbank
from genbankqc import Metadata


assembly_summary = pd.read_csv('test/resources/.info/assembly_summary.txt', sep="\t", index_col=0)


@pytest.fixture(scope="module")
def genome(aphidicola):
    handler = TestHandler()
    genome = ("GCA_000521565.1_Buchnera_aphidicola_G002_Myzus_persicae_Complete_Genome.fasta")
    genome = os.path.join(aphidicola.path, genome)
    with handler:
        genome = Genome(genome, assembly_summary)
        genome.sketch()
        genome.get_contigs()
        genome.get_assembly_size()
        genome.get_unknowns()
        yield genome, handler


@pytest.fixture(scope="module")
def genbank():
    resources = 'test/resources'
    genbank = tempfile.mkdtemp()
    for resource in os.listdir(resources):
        source = os.path.join(resources, resource)
        target = os.path.join(genbank, resource)
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
