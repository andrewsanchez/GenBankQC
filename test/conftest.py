import pytest

import genbankfilter.filter as gbf


@pytest.fixture(params=["Buchnera_aphidicola", "Acinetobacter_baumannii"])
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


@pytest.fixture()
def provide_aphidicola(request):
    aphidicola = "test/resources/Buchnera_aphidicola"
    aphidicola = gbf.Species(aphidicola)
    yield aphidicola


@pytest.fixture(params=[[200, 3.0, 3.0, 3.0], [300, 2.0, 2.0, 2.0]])
def provide_aphidicola_multi(request):
    a, b, c, d = request.param
    aphidicola = "test/resources/Buchnera_aphidicola"
    aphidicola = gbf.FilteredSpecies(aphidicola, a, b, c, d)
    yield request.param, aphidicola
