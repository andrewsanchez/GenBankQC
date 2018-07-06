import os.path
import shutil
import tempfile

import pytest

from genbankqc import Genbank
from genbankqc import Species
from genbankqc import Genome
from genbankqc import Metadata


@pytest.fixture(scope="module")
def genbank():
    tmp = tempfile.mkdtemp()
    genbank = os.path.join(tmp, 'genbank')
    shutil.copytree('test/resources/.info', os.path.join(genbank, ".info"))
    shutil.copytree('test/resources/Acinetobacter_baumannii',
                    os.path.join(genbank, "Acinetobacter_baumannii"))
    shutil.copytree('test/resources/Buchnera_aphidicola',
                    os.path.join(genbank, "Buchnera_aphidicola"))
    shutil.copytree('test/resources/Escherichia_coli',
                    os.path.join(genbank, "Escherichia_coli"))
    yield Genbank(genbank)
    shutil.rmtree(tmp)


@pytest.fixture(scope="module",
                params=["Buchnera_aphidicola", "Acinetobacter_baumannii"])
@pytest.fixture()
def species(request):
    species = "test/resources/{}".format(request.param)
    species = Species(species)
    yield species


@pytest.fixture(scope="module")
def aphidicola():
    tmp = tempfile.mkdtemp()
    aphidicola = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', aphidicola)
    yield Species(aphidicola)
    shutil.rmtree(tmp)


@pytest.fixture(params=[[200, 3.0, 3.0, 3.0], [300, 2.0, 2.0, 2.0]])
def aphidicola_multi(request):
    a, b, c, d = request.param
    aphidicola = "test/resources/Buchnera_aphidicola"
    aphidicola = Species(aphidicola, a, b, c, d)
    yield request.param, aphidicola


@pytest.fixture(scope="module")
def altered_unknowns():
    aphidicola = Species("test/resources/Buchnera_aphidicola")
    aphidicola.stats.iloc[:, 3] = 0
    aphidicola.stats.iloc[:3, 3] = 300
    expected_failures = aphidicola.stats.iloc[:3, 3].index.tolist()
    yield aphidicola, expected_failures


@pytest.fixture(scope="module")
def aphidicola_bare():
    tmp = tempfile.mkdtemp()
    aphidicola = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', aphidicola)
    shutil.rmtree(os.path.join(aphidicola, 'qc'))
    aphidicola = Species(aphidicola)
    yield aphidicola
    shutil.rmtree(tmp)


@pytest.fixture(scope="module")
def genome(genbank, aphidicola):
    genome = ("GCA_000521565.1_Buchnera_aphidicola_G002_"
              "Myzus_persicae_Complete_Genome.fasta")
    genome = os.path.join(aphidicola.path, genome)
    genome = Genome(genome, genbank.assembly_summary)
    genome.sketch()
    genome.get_contigs()
    genome.get_assembly_size()
    genome.get_unknowns()
    yield genome


@pytest.fixture(scope="module")
def ecoli(genbank):
    genome = ("GCA_002012025.1_Escherichia_coli_"
              "Ecol_542_Complete_Genome.fasta")
    genome = os.path.join(genbank.path, "Escherichia_coli", genome)
    genome = Genome(genome, genbank.assembly_summary)
    yield genome


@pytest.fixture(scope="module")
def five_genomes(aphidicola):
    shutil.rmtree(aphidicola.qc_dir)
    for genome in list(aphidicola.genomes())[:5]:
        os.remove(genome.path)
    yield aphidicola


@pytest.fixture(scope="module")
def metadata(genbank):
    yield Metadata(genbank.path)
