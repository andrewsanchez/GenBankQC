import os.path

import genbankfilter.filter as gbf
from genbankfilter.Genome import Genome


def test_init(aphidicola):
    from ete3 import Tree
    assert type(aphidicola) == gbf.Species
    assert type(aphidicola.stats) == gbf.pd.DataFrame
    assert type(aphidicola.tree) == Tree
    assert type(aphidicola.dmx) == gbf.pd.DataFrame
    assert aphidicola.dmx.index.tolist() == aphidicola.stats.index.tolist()
    assert (aphidicola.dmx.mean().index.tolist() ==
            aphidicola.stats.index.tolist())
    assert os.path.isdir(aphidicola.qc_dir)


def test_genomes(aphidicola):
    assert len(list(aphidicola.genomes())) == 10
    assert isinstance(next(aphidicola.genomes()), Genome)


def test_genome_ids(aphidicola):
    genome_ids = aphidicola.genome_ids()
    assert genome_ids.tolist() == aphidicola.stats.index.tolist()


def test_sketches(aphidicola):
    aphidicola_sketches = aphidicola.sketches()
    for i in aphidicola_sketches:
        assert i is None or "GCA_000007365.1" in i


def test_sketch(aphidicola):
    aphidicola.sketch()
    aphidicola_sketches = aphidicola.sketches()
    for i in aphidicola_sketches:
        assert i is not None
        assert os.path.isfile(i)


def test_mash_paste(aphidicola):
    aphidicola.mash_paste()
    assert os.path.isfile(aphidicola.paste_file)


def test_mash_dist(aphidicola_bare):
    aphidicola = aphidicola_bare
    aphidicola.mash_dist()
    assert os.path.isfile(aphidicola.dmx_path)
    assert type(aphidicola.dmx) == gbf.pd.DataFrame


def test_mash(aphidicola_bare):
    aphidicola = aphidicola_bare
    assert os.path.isfile(aphidicola.paste_file)
    assert os.path.isfile(aphidicola.dmx_path)
    assert type(aphidicola.dmx) == gbf.pd.DataFrame
    aphidicola_sketches = aphidicola.sketches()
    aphidicola_sketches = aphidicola.sketches()
    aphidicola_sketches = aphidicola.sketches()
    for i in aphidicola_sketches:
        assert i is not None
        assert os.path.isfile(i)


def test_get_stats(aphidicola_bare):
    aphidicola = aphidicola_bare
    aphidicola.get_stats()
    assert os.path.isfile(aphidicola.stats_path)
    assert type(aphidicola.stats) == gbf.pd.DataFrame
