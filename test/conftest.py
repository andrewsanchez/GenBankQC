import os.path
import shutil
import tempfile

import pytest

from genbank_qc import QC, Species


@pytest.fixture(scope="module",
                params=["Buchnera_aphidicola", "Acinetobacter_baumannii"])
@pytest.fixture()
def species(request):
    species = "test/resources/{}".format(request.param)
    species = QC(species)
    yield species


@pytest.fixture(scope="module")
def aphidicola(request):
    tmp = tempfile.mkdtemp()
    aphidicola = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', aphidicola)
    yield Species(aphidicola)
    shutil.rmtree(tmp)


@pytest.fixture(params=[[200, 3.0, 3.0, 3.0], [300, 2.0, 2.0, 2.0]])
def aphidicola_multi(request):
    a, b, c, d = request.param
    aphidicola = "test/resources/Buchnera_aphidicola"
    aphidicola = QC(aphidicola, a, b, c, d)
    yield request.param, aphidicola


def altered_unknowns():
    aphidicola = QC("test/resources/Buchnera_aphidicola")
    expected_failures = []
    yield aphidicola, expected_failures
    aphidicola = QC("test/resources/Buchnera_aphidicola")
    aphidicola.stats.iloc[:, 0] = 0
    aphidicola.stats.iloc[:10, 0] = 300
    expected_failures = aphidicola.stats.iloc[:10, 0].index.tolist()
    yield aphidicola, expected_failures


@pytest.fixture(scope="module",
                params=altered_unknowns())
def unknowns(request):
    aphidicola, failures = request.param
    yield aphidicola, failures


@pytest.fixture(scope="module")
def aphidicolaQC():
    tmp = tempfile.mkdtemp()
    aphidicola = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', aphidicola)
    yield QC(aphidicola)
    shutil.rmtree(tmp)


@pytest.fixture(scope="module")
def aphidicola_bare(aphidicolaQC):
    aphidicola = aphidicolaQC
    shutil.rmtree(aphidicola.qc_dir)
    aphidicola.dmx = None
    aphidicola.tree = None
    aphidicola.stats = None
    yield aphidicola


@pytest.fixture(scope="module")
def filtered(aphidicolaQC):
    aphidicola = aphidicolaQC
    aphidicola.filter()
    print(aphidicola.path)
    yield aphidicola


@pytest.fixture(scope="module")
def genome(aphidicola):
    genome = next(aphidicola.genomes())
    genome.get_contigs()
    genome.get_assembly_size()
    genome.get_unknowns()
    yield genome


@pytest.fixture()
def path():
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', path)
    yield path
    shutil.rmtree(tmp)
