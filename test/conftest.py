import os.path
import shutil
import tempfile

import pytest

import genbankfilter.filter as gbf
from genbankfilter.Genome import Genome


@pytest.fixture(scope="module",
                params=["Buchnera_aphidicola", "Acinetobacter_baumannii"])
def provide_Species(request):
    species = "test/resources/" + request.param
    species = gbf.Species(species)
    yield species


@pytest.fixture()
def provide_baumannii(request):
    baumannii = "test/resources/Acinetobacter_baumannii"
    baumannii = gbf.FilteredSpecies(baumannii)
    # Initialize the otherwise empty `passed` DataFrame
    baumannii.passed = baumannii.stats
    yield baumannii


@pytest.fixture(scope="module")
def provide_aphidicola(request):
    tmp = tempfile.mkdtemp()
    aphidicola = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', aphidicola)
    aphidicola = gbf.Species(aphidicola)
    yield aphidicola
    shutil.rmtree(aphidicola.species_dir)


@pytest.fixture(params=[[200, 3.0, 3.0, 3.0], [300, 2.0, 2.0, 2.0]])
def provide_aphidicola_multi(request):
    a, b, c, d = request.param
    aphidicola = "test/resources/Buchnera_aphidicola"
    aphidicola = gbf.FilteredSpecies(aphidicola, a, b, c, d)
    yield request.param, aphidicola


@pytest.fixture(scope="module")
def genome(request, provide_aphidicola):
    aphidicola = provide_aphidicola
    genome = next(aphidicola.genomes())
    genome.get_contigs()
    genome.get_assembly_size()
    genome.get_unknowns()
    yield genome
