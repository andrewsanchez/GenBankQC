import os.path

import genbankfilter.filter as gbf
from genbankfilter.Genome import Genome


def test_init(provide_aphidicola):
    from ete3 import Tree
    species = provide_aphidicola
    assert type(species) == gbf.Species
    assert type(species.stats) == gbf.pd.DataFrame
    assert type(species.tree) == Tree
    assert type(species.dmx) == gbf.pd.DataFrame
    assert species.dmx.index.tolist() == species.stats.index.tolist()
    assert species.dmx.mean().index.tolist() == species.stats.index.tolist()
    assert os.path.isdir(species.qc_dir)


def test_genomes(provide_aphidicola):
    aphidicola = provide_aphidicola
    assert len(list(aphidicola.genomes())) == 10
    assert isinstance(next(aphidicola.genomes()), Genome)


def test_genome_ids(provide_aphidicola):
    aphidicola = provide_aphidicola
    genome_ids = aphidicola.genome_ids()
    assert genome_ids.tolist() == aphidicola.stats.index.tolist()


def test_sketches(provide_aphidicola):
    aphidicola = provide_aphidicola
    aphidicola_sketches = aphidicola.sketches()
    for i in aphidicola_sketches:
        assert i is None or "GCA_000007365.1" in i


def test_sketch(provide_aphidicola):
    aphidicola = provide_aphidicola
    aphidicola.sketch()
    aphidicola_sketches = aphidicola.sketches()
    for i in aphidicola_sketches:
        assert i is not None


def test_mash_paste(provide_aphidicola):
    aphidicola = provide_aphidicola
    aphidicola.mash_paste()
    assert os.path.isfile(aphidicola.paste_file)


# TODO: Update conftest.py to provide this fixture scenario
def test_mash_dist(provide_aphidicola):
    aphidicola = provide_aphidicola
    os.remove(os.path.join(aphidicola.qc_dir, 'dmx.csv'))
    del aphidicola.dmx
    aphidicola.mash_dist()
    assert os.path.isfile(os.path.join(aphidicola.qc_dir, 'dmx.csv'))
    assert type(aphidicola.dmx) == gbf.pd.DataFrame


def test_get_stats(provide_aphidicola):
    aphidicola = provide_aphidicola
    del aphidicola.stats
    os.remove(os.path.join(aphidicola.qc_dir, 'stats.csv'))
    aphidicola.get_stats()
    assert os.path.isfile(os.path.join(aphidicola.qc_dir, 'stats.csv'))
    assert type(aphidicola.stats) == gbf.pd.DataFrame
    print(aphidicola.stats)
