# Tests for the source data validation in loaders/validation.py.
#
# validate_loaded_data() runs at the end of load_tsv_to_db(), so the mock dataset is
# already validated by the time any fixture yields - every other test in the suite would
# fail if these assumptions stopped holding for it. What is tested here is that each check
# actually catches its violation, and says which rows caused it.

import pytest

from db import db
from loaders.validation import SourceDataError, validate_loaded_data
from models.immunediscoverdata import ImmuneDiscoverDataModel

TRANSLATED_ALLELE = "IGHV1-8*01"
DELETION_ALLELE = "IGHV1-8*DEL"


def add_row(**overrides):
    """Insert a row, defaulting every non-nullable column to something valid."""
    fields = {
        "cohort": "1KGP",
        "case": "case_GBR_EUR",
        "db_name": "IGHD1-99*99",
        "gene": "IGHD1-99",
        "allele": "99",
        "flank_index": 1,
        "db_name_AA": None,
        "superpopulation": "EUR",
        "population": "GBR",
        "loaded_from_tsv": "test.tsv",
        "loaded_at": "2026-08-14",
    }
    fields.update(overrides)
    db.session.add(ImmuneDiscoverDataModel(**fields))
    db.session.commit()


def set_amino_acid_name(db_name, value):
    ImmuneDiscoverDataModel.query.filter(
        ImmuneDiscoverDataModel.db_name == db_name
    ).update({"db_name_AA": value})
    db.session.commit()


def test_mock_data_passes_validation(app):
    # The fixture already loaded and validated it; assert it explicitly so a failure here
    # points at the data rather than surfacing as every other test breaking at once.
    validate_loaded_data()


def test_rejects_translated_allele_without_amino_acid_name(app):
    # An IGHV allele that is not a flanking variant or a deletion must have db_name_AA.
    # Without it, it silently disappears from the amino acid plots and the translated FASTA.
    set_amino_acid_name(TRANSLATED_ALLELE, None)

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    assert TRANSLATED_ALLELE in str(raised.value)
    assert "should be translated" in str(raised.value)


def test_rejects_amino_acid_name_on_untranslated_row(app):
    # The reverse: a *DEL row with db_name_AA set would show up in the amino acid plots,
    # and in the translated FASTA download, which relies on the null to exclude it.
    set_amino_acid_name(DELETION_ALLELE, "IGHV1-8*01")

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    assert DELETION_ALLELE in str(raised.value)
    assert "never translated" in str(raised.value)


def test_rejects_unknown_population(app):
    # A population missing from the display order in services/frequencies.py is left out
    # of every plot rather than reported, so it has to fail loudly here instead.
    add_row(population="ZZZ")

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    assert "ZZZ" in str(raised.value)
    assert "population" in str(raised.value)


def test_rejects_unknown_superpopulation(app):
    add_row(superpopulation="QQQ")

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    assert "QQQ" in str(raised.value)


def test_rejects_ambiguous_gene_and_allele(app):
    # Two non-flanking allele names under one gene/allele pair makes the db_name that
    # get_db_name_from_options resolves a selection to arbitrary.
    add_row(db_name="IGHD1-1*01", gene="IGHD1-1", allele="01", case="case_IBS_EUR")
    add_row(db_name="IGHD1-1*01_duplicate", gene="IGHD1-1", allele="01", case="case_IBS_EUR")

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    assert "IGHD1-1,01" in str(raised.value)
    assert "more than one allele name" in str(raised.value)


def test_reports_every_problem_at_once(app):
    # One crash should list everything that needs fixing, so the research group gets a
    # complete report rather than one problem per restart.
    set_amino_acid_name(TRANSLATED_ALLELE, None)
    add_row(population="ZZZ")

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    message = str(raised.value)
    assert "breaks 2 assumption(s)" in message
    assert TRANSLATED_ALLELE in message
    assert "ZZZ" in message


def test_flanking_variants_are_allowed_to_lack_amino_acid_data(app):
    # The mock data includes IGHV1-8*01_F1 with no db_name_AA. That is correct, not a
    # violation, and this pins it so the check cannot be tightened into rejecting it.
    flanking = ImmuneDiscoverDataModel.query.filter(
        ImmuneDiscoverDataModel.db_name.contains("_F", autoescape=True)
    ).first()

    assert flanking is not None, "mock data no longer has a flanking row to cover this"
    assert flanking.db_name_AA is None
    validate_loaded_data()


def test_rejects_an_unrecognised_locus(app):
    # Deriving "never translated" from PLOT_LOCI means an unknown locus is treated as
    # untranslated and quietly left out of the plots. That is consistent but silent, so the
    # locus itself is checked - and reported as an unknown locus rather than mis-reported as
    # a translation problem, which is what the hardcoded IGHD/IGHJ list used to do.
    add_row(db_name="TRGJ1*01", gene="TRGJ1", allele="01", db_name_AA=None)

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    message = str(raised.value)
    assert "Unknown locus/loci TRGJ" in message
    assert "should be translated" not in message, "wrong diagnosis for a J locus"


def test_a_translated_locus_is_still_checked(app):
    # The derivation has not simply stopped checking: a V gene of a known locus without
    # amino acid data is still a problem.
    set_amino_acid_name(TRANSLATED_ALLELE, None)

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    assert "should be translated" in str(raised.value)
