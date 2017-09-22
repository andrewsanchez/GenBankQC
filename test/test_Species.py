import genbankfilter.filter as gbf


def test_init(provide_Species):
    from ete3 import Tree
    species = provide_Species
    assert type(species) == gbf.Species
    assert type(species.stats) == gbf.pd.DataFrame
    assert type(species.tree) == Tree
    assert type(species.dmx) == gbf.pd.DataFrame
    assert species.dmx.index.tolist() == species.stats.index.tolist()
    assert species.dmx.mean().index.tolist() == species.stats.index.tolist()


def test_genomes(provide_aphidicola):
    from types import GeneratorType
    aphidicola = provide_aphidicola
    genomes = aphidicola.genomes()
    assert type(genomes) == GeneratorType
