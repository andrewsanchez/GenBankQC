from genbankfilter.Genome import Genome


def test_Genome_init(provide_aphidicola):
    aphidicola = provide_aphidicola
    genomes = aphidicola.genomes()
    genome = Genome(next(genomes))
    assert genome.name == "GCA_000007365.1_Buchnera_aphidicola_Sg_" + \
        "Schizaphis_graminum_Complete_Genome"


def test_Genome_get_contigs(provide_aphidicola):
    from Bio.Seq import Seq
    aphidicola = provide_aphidicola
    genomes = aphidicola.genomes()
    genome = Genome(next(genomes))
    genome.get_contigs()
    assert type(genome.contigs) == list
    assert type(genome.contigs[0]) == Seq
