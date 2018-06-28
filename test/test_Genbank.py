from genbankqc import Genbank, Species


def test_genbank_init(genbank):
    assert isinstance(genbank, Genbank)
    species = genbank.species
    for i in species:
        assert isinstance(i, Species)
