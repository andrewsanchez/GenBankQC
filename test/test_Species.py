import os
import pytest

import pandas as pd
from pandas.util.testing import assert_index_equal

from genbankqc import Genome, Species


def test_init(aphidicola_multi):
    from ete3 import Tree
    params, aphidicola = aphidicola_multi
    a, b, c, d = params
    assert type(aphidicola) == Species
    assert type(aphidicola.stats) == pd.DataFrame
    assert type(aphidicola.tree) == Tree
    assert type(aphidicola.dmx) == pd.DataFrame
    assert aphidicola.total_genomes == 10
    assert aphidicola.dmx.index.tolist() == aphidicola.stats.index.tolist()
    assert (aphidicola.dmx.mean().index.tolist() ==
            aphidicola.stats.index.tolist())
    assert os.path.isdir(aphidicola.qc_dir)
    assert os.path.isdir(aphidicola.qc_results_dir)
    assert aphidicola.max_unknowns == a
    assert aphidicola.contigs == b
    assert aphidicola.assembly_size == c
    assert aphidicola.mash == d
    assert aphidicola.tolerance["unknowns"] == a
    assert aphidicola.tolerance["contigs"] == b
    assert aphidicola.tolerance["assembly_size"] == c
    assert aphidicola.tolerance["distance"] == d
    assert aphidicola.label == "-".join(map(str, params))
    assert id(aphidicola.stats) == id(aphidicola.passed)
    assert aphidicola.tree_complete is True


def test_genomes(aphidicola):
    assert len(list(aphidicola.genomes())) == 10
    assert isinstance(next(aphidicola.genomes()), Genome)


def test_genome_ids(aphidicola):
    assert_index_equal(aphidicola.genome_ids().sort_values(),
                       aphidicola.stats.index.sort_values())


def test_sketches(aphidicola):
    from re import match
    from os.path import basename
    aphidicola_sketches = aphidicola.sketches()
    for i in aphidicola_sketches:
        assert isinstance(i, str)
        assert basename(i).replace('.msh', '') in aphidicola.genome_ids()


def test_filter(aphidicola):
    aphidicola.filter()
    assert aphidicola.complete is False
    total_failed = sum(map(len, aphidicola.failed.values()))
    assert os.path.isfile(aphidicola.summary_path)
    assert sum([total_failed, len(aphidicola.passed)]) \
        == len(aphidicola.stats)
    aphidicola.filter()
    assert aphidicola.complete is True
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


@pytest.mark.usefixtures("aphidicola_bare")
class TestBare:

    def test_sketch(self, aphidicola_bare):
        aphidicola = aphidicola_bare
        aphidicola.sketch()
        aphidicola_sketches = aphidicola.sketches()
        for i in aphidicola_sketches:
            assert i is not None
            assert os.path.isfile(i)

    def test_mash_paste(self, aphidicola_bare):
        aphidicola = aphidicola_bare
        aphidicola.mash_paste()
        assert os.path.isfile(aphidicola.paste_file)

    def test_mash_dist(self, aphidicola_bare):
        aphidicola = aphidicola_bare
        aphidicola.mash_dist()
        assert os.path.isfile(aphidicola.dmx_path)
        assert type(aphidicola.dmx) == pd.DataFrame

    def test_mash(self, aphidicola_bare):
        aphidicola = aphidicola_bare
        aphidicola.run_mash()
        assert os.path.isfile(aphidicola.paste_file)
        assert os.path.isfile(aphidicola.dmx_path)
        assert type(aphidicola.dmx) == pd.DataFrame
        aphidicola_sketches = aphidicola.sketches()
        for i in aphidicola_sketches:
            assert i is not None
            assert os.path.isfile(i)

    def test_get_stats(self, aphidicola_bare):
        aphidicola = aphidicola_bare
        aphidicola.get_stats()
        assert os.path.isfile(aphidicola.stats_path)
        assert type(aphidicola.stats) == pd.DataFrame

    def test_get_tree(self, aphidicola_bare):
        from ete3 import Tree
        aphidicola = aphidicola_bare
        aphidicola.get_tree()
        assert type(aphidicola.tree) == Tree
        assert os.path.isfile(aphidicola.nw_path)
        os.remove(aphidicola.nw_path)
        aphidicola.assess_tree()
        aphidicola.get_tree()
        assert not os.path.isfile(aphidicola.nw_path)
        assert type(aphidicola.tree) == Tree


def test_filter_unknowns(unknowns):
    aphidicola, expected_failures = unknowns
    aphidicola.filter_unknown_bases()
    passed_and_failed = sum(map(len, [aphidicola.failed["unknowns"],
                                      aphidicola.passed]))
    assert len(aphidicola.stats) == passed_and_failed
    assert isinstance(aphidicola.passed, pd.DataFrame)
    assert isinstance(aphidicola.failed["unknowns"], pd.Index)
    assert(id(aphidicola.stats) != id(aphidicola.passed))
    assert expected_failures == aphidicola.failed["unknowns"].tolist()


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
    passed_and_failed = sum(
        map(len, [species.failed['assembly_size'], species.passed]))
    assert passed_and_failed == genomes_before_filtering
    genomes_before_filtering = len(species.passed)
    species.filter_MAD_upper('distance')
    assert type(species.passed) == pd.DataFrame
    assert type(species.failed['distance']) == pd.Index
    passed_and_failed = sum(
        map(len, [species.failed['distance'], species.passed]))
    assert passed_and_failed == genomes_before_filtering


def test_min_genomes(five_genomes):
    five_genomes.qc()
    assert not os.listdir(five_genomes.qc_dir)
