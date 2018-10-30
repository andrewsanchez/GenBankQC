import os
import pytest
import shutil
import tempfile

import pandas as pd
from pandas.util.testing import assert_index_equal

from genbankqc import Species
from genbankqc import Genome


assembly_summary = pd.read_csv('test/resources/.info/assembly_summary.txt', sep="\t", index_col=0)


@pytest.fixture()
def species():
    """
    Provides a Species object for B. aphidicola, which contains one misnamed genome, a qc
    directory, one sketch file, dmx.csv, stats.csv, tree.nw and two filtered directories.
    """
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', path)
    species = Species(path, assembly_summary=assembly_summary)
    yield species
    shutil.rmtree(tmp)


def test_filter_contigs(species):
    species.filter_contigs('contigs')
    total_genomes = len(species.passed) + len(species.failed["contigs"])
    assert total_genomes == len(species.stats)
    assert isinstance(species.med_abs_devs["contigs"], float)
    assert isinstance(species.dev_refs["contigs"], float)
    assert isinstance(species.failed["contigs"], pd.Index)
    assert isinstance(species.allowed["contigs"], float)
    assert isinstance(species.passed, pd.DataFrame)


def test_filter_MAD(species):
    genomes_before_filtering = len(species.passed)
    species.filter_MAD_range('assembly_size')
    assert type(species.passed) == pd.DataFrame
    assert type(species.failed['assembly_size']) == pd.Index
    passed_and_failed = sum(map(len, [species.failed['assembly_size'], species.passed]))
    assert passed_and_failed == genomes_before_filtering
    genomes_before_filtering = len(species.passed)
    species.filter_MAD_upper('distance')
    assert type(species.passed) == pd.DataFrame
    assert type(species.failed['distance']) == pd.Index
    passed_and_failed = sum(
        map(len, [species.failed['distance'], species.passed]))
    assert passed_and_failed == genomes_before_filtering


@pytest.fixture(params=[[200, 3.0, 3.0, 3.0], [300, 2.0, 2.0, 2.0]])
def aphidicola_multi(request):
    a, b, c, d = request.param
    aphidicola = "test/resources/Buchnera_aphidicola"
    aphidicola = Species(aphidicola, a, b, c, d)
    yield request.param, aphidicola


def test_init(aphidicola_multi):
    from ete3 import Tree
    params, aphidicola = aphidicola_multi
    a, b, c, d = params
    assert type(aphidicola) == Species
    assert type(aphidicola.stats) == pd.DataFrame
    assert type(aphidicola.tree) == Tree
    assert type(aphidicola.dmx) == pd.DataFrame
    assert aphidicola.total_genomes == 10
    assert (sorted(aphidicola.dmx.index.tolist()) == sorted(aphidicola.stats.index.tolist()))
    assert (sorted(aphidicola.dmx.mean().index.tolist()) ==
            sorted(aphidicola.stats.index.tolist()))
    assert aphidicola.max_unknowns == a
    assert aphidicola.contigs == b
    assert aphidicola.assembly_size == c
    assert aphidicola.mash == d
    assert aphidicola.tolerance["unknowns"] == a
    assert aphidicola.tolerance["contigs"] == b
    assert aphidicola.tolerance["assembly_size"] == c
    assert aphidicola.tolerance["distance"] == d
    assert aphidicola.label == "{}-{}-{}-{}".format(a, b, c, d)
    assert id(aphidicola.stats) == id(aphidicola.passed)
    assert aphidicola.tree_complete is True


@pytest.fixture(scope="module")
def altered_unknowns():
    aphidicola = Species("test/resources/Buchnera_aphidicola")
    aphidicola.stats.iloc[:, 3] = 0
    aphidicola.stats.iloc[:3, 3] = 300
    expected_failures = aphidicola.stats.iloc[:3, 3].index.tolist()
    yield aphidicola, expected_failures


def test_genomes(aphidicola):
    assert len(list(aphidicola.genomes)) == 10
    assert isinstance(aphidicola.genomes[0], Genome)


def test_genome_ids(aphidicola):
    assert_index_equal(aphidicola.genome_ids().sort_values(), aphidicola.stats.index.sort_values())


def test_sketches(aphidicola):
    from os.path import basename
    aphidicola_sketches = aphidicola.sketches()
    for i in aphidicola_sketches:
        assert isinstance(i, str)
        assert basename(i).replace('.msh', '') in aphidicola.genome_ids()


def test_filter(aphidicola):
    aphidicola.filter()
    # this won't work because Species.complete is set by the assess
    # decorator which only wraps the qc() method
    # maybe assessment should occur in instantiation
    # assert aphidicola.complete is False
    total_failed = sum(map(len, aphidicola.failed.values()))
    assert os.path.isfile(aphidicola.summary_path)
    assert sum([total_failed, len(aphidicola.passed)]) \
        == len(aphidicola.stats)
    aphidicola.filter()
    # assert aphidicola.complete is True
    assert isinstance(aphidicola.allowed, dict)


def test_link_genomes(aphidicola):
    aphidicola.link_genomes()
    passed = ["{}.fasta".format(i)
              for i in aphidicola.passed.index.tolist()]
    assert os.listdir(aphidicola.passed_dir)
    assert sorted(os.listdir(aphidicola.passed_dir)) == \
        sorted(passed)


def test_failed_report(aphidicola):
    assert os.path.isfile(aphidicola.failed_path)


def test_color_tree(aphidicola):
    aphidicola = aphidicola
    aphidicola.color_tree()
    import subprocess
    subprocess.call("open {}".format(aphidicola.tree_img), shell=True)
    assert os.path.isfile(aphidicola.tree_img)


@pytest.fixture(scope="module")
def species_bare():
    tmp = tempfile.mkdtemp()
    aphidicola = os.path.join(tmp, "Buchnera_aphidicola")
    shutil.copytree('test/resources/Buchnera_aphidicola', aphidicola)
    shutil.rmtree(os.path.join(aphidicola, 'qc'))
    aphidicola = Species(aphidicola)
    yield aphidicola
    shutil.rmtree(tmp)


@pytest.mark.usefixtures("species_bare")
class TestBare:
    def test_mash(self, species_bare):
        aphidicola = species_bare
        aphidicola.run_mash()
        assert os.path.isfile(aphidicola.paste_file)
        assert os.path.isfile(aphidicola.dmx_path)
        assert type(aphidicola.dmx) == pd.DataFrame
        aphidicola_sketches = aphidicola.sketches()
        for i in aphidicola_sketches:
            assert i is not None
            assert os.path.isfile(i)

    def test_get_stats(self, species_bare):
        aphidicola = species_bare
        aphidicola.get_stats()
        assert os.path.isfile(aphidicola.stats_path)
        assert type(aphidicola.stats) == pd.DataFrame

    def test_get_tree(self, species_bare):
        from ete3 import Tree
        aphidicola = species_bare
        aphidicola.get_tree()
        assert type(aphidicola.tree) == Tree
        assert os.path.isfile(aphidicola.nw_path)
        os.remove(aphidicola.nw_path)
        aphidicola.assess_tree()
        aphidicola.get_tree()
        assert not os.path.isfile(aphidicola.nw_path)
        assert type(aphidicola.tree) == Tree


def test_filter_unknowns(altered_unknowns):
    aphidicola, expected_failures = altered_unknowns
    aphidicola.filter_unknown_bases()
    passed_and_failed = sum(map(len, [aphidicola.failed["unknowns"],
                                      aphidicola.passed]))
    assert len(aphidicola.stats) == passed_and_failed
    assert isinstance(aphidicola.passed, pd.DataFrame)
    assert isinstance(aphidicola.failed["unknowns"], pd.Index)
    assert(id(aphidicola.stats) != id(aphidicola.passed))
    assert expected_failures == aphidicola.failed["unknowns"].tolist()


@pytest.fixture()
def five_genomes(aphidicola):
    shutil.rmtree(aphidicola.qc_dir)
    for genome in list(aphidicola.genomes)[:5]:
        os.remove(genome.path)
    yield aphidicola


def test_min_genomes(five_genomes):
    # with pytest.raises(FileNotFoundError):
    five_genomes.qc()
    assert not os.path.isdir(five_genomes.qc_dir)


# def test_metadata(species):
#     species.metadata()
#     assert isinstance(species.metadata_df, pd.DataFrame)
#     assert os.path.isfile(species.metadata_path)
