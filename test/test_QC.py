import os

from pandas import DataFrame, Index

import genbankfilter.filter as gbf


def test_init(aphidicola_multi):
    params, aphidicola = aphidicola_multi
    a, b, c, d = params
    assert aphidicola.max_unknowns == a
    assert aphidicola.contigs == b
    assert aphidicola.assembly_size == c
    assert aphidicola.mash == d
    assert aphidicola.tolerance["unknowns"] == a
    assert aphidicola.tolerance["contigs"] == b
    assert aphidicola.tolerance["assembly_size"] == c
    assert aphidicola.tolerance["distance"] == d
    assert aphidicola.label == "-".join(map(str, params))
    assert id(aphidicola.stats) == id(aphidicola.passed)


def test_filter_unknowns(unknowns):
    aphidicola, expected_failures = unknowns
    aphidicola.filter_unknown_bases()
    passed_and_failed = sum(map(len, [aphidicola.failed["unknowns"],
                                      aphidicola.passed]))
    assert len(aphidicola.stats) == passed_and_failed
    assert isinstance(aphidicola.passed, DataFrame)
    assert isinstance(aphidicola.failed["unknowns"], Index)
    assert(id(aphidicola.stats) != id(aphidicola.passed))
    assert expected_failures == aphidicola.failed["unknowns"].tolist()


def test_filter_contigs(species):
    species.filter_contigs()
    total_genomes = len(species.passed) + len(species.failed["contigs"])
    assert total_genomes == len(species.stats)
    assert isinstance(species.med_abs_devs["contigs"], float)
    assert isinstance(species.dev_refs["contigs"], float)
    assert isinstance(species.failed["contigs"], Index)
    assert isinstance(species.allowed["contigs"], float)
    assert isinstance(species.passed, DataFrame)


def test_filter_med_abs_dev(species):
    for criteria in ["distance", "assembly_size"]:
        genomes_before_filtering = len(species.passed)
        species.filter_med_abs_dev(criteria)
        assert type(species.passed) == DataFrame
        assert type(species.failed[criteria]) == Index
        passed_and_failed = sum(map(len, [species.failed[criteria],
                                          species.passed]))
        assert passed_and_failed == genomes_before_filtering


