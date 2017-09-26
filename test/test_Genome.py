import os.path


def test_init(genome):
    from genbankfilter.Genome import Genome
    expected_name = ("GCA_000007365.1_Buchnera_aphidicola_Sg_"
                     "Schizaphis_graminum_Complete_Genome")
    assert isinstance(genome, Genome)
    assert genome.name == expected_name
    assert not genome.msh


def test_get_contigs(genome):
    from Bio.Seq import Seq
    assert type(genome.contigs) is list
    assert type(genome.contigs[0]) is Seq
    assert isinstance(genome.count_contigs, int)


def test_assembly_size(genome):
    assert type(genome.assembly_size) is int


def test_unknowns(genome):
    assert type(genome.unknowns) is int


def test_sketch(genome):
    genome.sketch()
    assert os.path.isfile(genome.msh)
