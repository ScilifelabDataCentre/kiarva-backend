# Tests for the pure query-building helpers in repositories/filters.py and the display
# order lookup in services/frequencies.py. No database or app fixture: these take a string
# and return a column or a list.

import pytest

import repositories.filters
from repositories.filters import allele_column, locus_of
from services.frequencies import population_display_order


def test_unknown_plot_type_raises_a_useful_error():
    # Four call sites chose their column with an if/elif and no else, so an unrecognised
    # plot_type left the variable unassigned and failed one line later with
    # UnboundLocalError - a 500 complaining about a local rather than about the argument.
    for bad in ("genomicc", "", None):
        with pytest.raises(ValueError, match="unknown plot_type"):
            allele_column(bad)
        with pytest.raises(ValueError, match="unknown population_type"):
            population_display_order(bad)

    assert allele_column("genomic") is not None
    assert allele_column("aminoacid") is not None
    assert population_display_order("superpopulation")[0] == "AFR"
    assert population_display_order("population")[0] == "ACB"


def test_locus_of_matches_known_prefixes():
    # Matched against KNOWN_LOCI rather than sliced at a fixed width, so a gene name shorter
    # than a prefix returns None instead of a truncated string that is in no list.
    assert locus_of("IGHV1-2") == "IGHV"
    assert locus_of("IGHD5-18/5-5") == "IGHD"
    assert locus_of("TRGV9") == "TRGV"
    assert locus_of("IGHV3-30,IGHV3-30-5") == "IGHV"
    assert locus_of("TRGJ1") is None
    assert locus_of("IGH") is None
    assert locus_of("") is None


def test_locus_of_does_not_assume_every_locus_is_four_characters(monkeypatch):
    # Slicing gene[:4] and matching KNOWN_LOCI by prefix are the same thing while every
    # locus happens to be four characters, so only a locus of another length tells them
    # apart. That is exactly the case this check exists to survive: a new locus arriving in
    # source data should be named in the crash, not turned into a truncated string that
    # matches nothing.
    monkeypatch.setattr(repositories.filters, "KNOWN_LOCI",
                        ("IGHV", "IGHD", "IGHJ", "TRGV", "IGKV1", "TRB"))

    assert locus_of("IGKV1-5") == "IGKV1"      # longer than four
    assert locus_of("TRB2-1") == "TRB"         # shorter than four
    assert locus_of("IGHV1-2") == "IGHV"       # unchanged
    assert locus_of("TRGJ1") is None
