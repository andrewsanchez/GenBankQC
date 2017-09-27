import os.path
import shutil
import tempfile

import pytest

import genbankfilter.filter as gbf


@pytest.fixture(scope="module",
                params=["Buchnera_aphidicola", "Acinetobacter_baumannii"])
@pytest.fixture()
def provide_baumannii(request):
    baumannii = "test/resources/Acinetobacter_baumannii"
    baumannii = gbf.FilteredSpecies(baumannii)
    # Initialize the otherwise empty `passed` DataFrame
    baumannii.passed = baumannii.stats
    yield baumannii


@pytest.fixture(scope="module")
def aphidicola(request):
    tmp = tempfile.mkdtemp()
    aphidicola = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', aphidicola)
    aphidicola = gbf.Species(aphidicola)
    yield aphidicola
    shutil.rmtree(tmp)


@pytest.fixture(scope="module")
def aphidicola_bare(request, aphidicola):
    os.remove(aphidicola.stats_path)
    os.remove(aphidicola.dmx_path)
    os.remove(aphidicola.nw_path)
    del aphidicola.dmx
    del aphidicola.tree
    del aphidicola.stats
    yield aphidicola


@pytest.fixture(params=[[200, 3.0, 3.0, 3.0], [300, 2.0, 2.0, 2.0]])
def aphidicola_multi(request):
    a, b, c, d = request.param
    aphidicola = "test/resources/Buchnera_aphidicola"
    aphidicola = gbf.FilteredSpecies(aphidicola, a, b, c, d)
    yield request.param, aphidicola


@pytest.fixture(scope="module")
def genome(request, aphidicola):
    genome = next(aphidicola.genomes())
    genome.get_contigs()
    genome.get_assembly_size()
    genome.get_unknowns()
    yield genome
