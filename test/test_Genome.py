def test_init(genome):
    from genbankfilter.Genome import Genome
    expected_name = ("GCA_000007365.1_Buchnera_aphidicola_Sg_"
                     "Schizaphis_graminum_Complete_Genome")
    assert isinstance(genome, Genome)
    assert genome.name == expected_name


def test_get_contigs(genome):
    from Bio.Seq import Seq
    assert type(genome.contigs) is list
    assert type(genome.contigs[0]) is Seq


def test_assembly_size(genome):
    assert type(genome.assembly_size) is int
