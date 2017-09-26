import genbankfilter.filter as gbf
from genbankfilter.Genome import Genome


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
    aphidicola = provide_aphidicola
    assert len(list(aphidicola.genomes())) == 10
    assert isinstance(next(aphidicola.genomes()), Genome)


def test_genome_ids(provide_aphidicola):
    aphidicola = provide_aphidicola
    genome_ids = aphidicola.genome_ids()
    assert genome_ids.tolist() == aphidicola.stats.index.tolist()
