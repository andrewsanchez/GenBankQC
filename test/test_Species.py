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
    print(os.listdir(aphidicola.species_dir))
    for i in aphidicola_sketches:
        assert i is not None


