import re
import pytest
import os.path

import pandas as pd
from genbankqc import Genome

assembly_summary = pd.read_csv('test/resources/.info/assembly_summary.txt', sep="\t", index_col=0)


@pytest.fixture(scope="module")
def ecoli_genome(genbank):
    genome = ("GCA_002012025.1_Escherichia_coli_Ecol_542_Complete_Genome.fasta")
    genome = os.path.join(genbank.path, "Escherichia_coli", genome)
    genome = Genome(genome, assembly_summary)
    yield genome


def test_init(genome):
    genome, handler = genome
    expected_name = ("GCA_000521565.1_Buchnera_aphidicola_G002_Myzus_persicae_Complete_Genome")
    expected_path = os.path.join(genome.species_dir, expected_name+".fasta")
    assert genome.path == expected_path
    assert isinstance(genome, Genome)
    assert genome.name == expected_name
    assert genome.name in genome.msh
    assert handler.has_info(re.compile('Instantiated'))


def test_get_contigs(genome):
    genome, handler = genome
    from Bio.Seq import Seq
    assert type(genome.contigs) is list
    assert type(genome.contigs[0]) is Seq
    assert isinstance(genome.count_contigs, int)


def test_assembly_size(genome):
    genome, handler = genome
    assert type(genome.assembly_size) is int


def test_unknowns(genome):
    genome, handler = genome
    assert type(genome.unknowns) is int


def test_get_distance(aphidicola, genome):
    genome, handler = genome
    genome.get_distance(aphidicola.dmx.mean())
    assert isinstance(genome.distance, float)


def test_sketch(genome):
    genome, handler = genome
    genome.sketch()
    assert os.path.isfile(genome.msh)


def test_get_stats(genome, aphidicola):
    genome, handler = genome
    from pandas import DataFrame
    dmx_mean = aphidicola.dmx.mean()
    genome.get_stats(dmx_mean)
    assert isinstance(genome.stats, DataFrame)
    assert os.path.isfile(genome.stats_path)


# def test_parse_empty_biosample(ecoli_genome):
#     ecoli_genome.parse_biosample()
#     assert ecoli_genome.metadata["sra_id"] == "missing"


# def test_efetch_biosample(ecoli_genome, genome):
#     genome, handler = genome
#     ecoli_genome.efetch("biosample")
#     assert ecoli_genome.xml["biosample"] is not None
#     genome.efetch("biosample")
#     assert genome.xml["biosample"] is not None


# def test_parse_biosample(ecoli_genome):
#     ecoli_genome.parse_biosample()
#     assert ecoli_genome.metadata.items() is not None


# def test_efetch_sra(ecoli_genome, genome):
#     genome, handler = genome
#     ecoli_genome.efetch("sra")
#     assert ecoli_genome.xml["sra"] is not None
#     genome.efetch("sra")
#     assert genome.xml["sra"] is 'missing'


# # Test for genome that doesn't have SRA info
# def test_efetch_fail(genome):
#     genome, handler = genome
#     pass


# def test_parse_sra(ecoli_genome):
#     ecoli_genome.parse_sra()
#     assert ecoli_genome.metadata["srs_accessions"] is not None


# def test_get_metadata(ecoli_genome, genome):
#     genome, handler = genome
#     ecoli_genome.get_metadata()
#     genome.get_metadata()
