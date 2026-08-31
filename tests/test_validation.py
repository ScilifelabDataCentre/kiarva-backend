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
    # get_db_name_from_options resolves a selection to arbitrary. A plotted locus, because
    # the check covers exactly the rows the resolver can reach - see the test below.
    add_row(db_name="IGHV9-99*01", gene="IGHV9-99", allele="01", case="case_IBS_EUR",
            db_name_AA="IGHV9-99*01")
    add_row(db_name="IGHV9-99*01_duplicate", gene="IGHV9-99", allele="01",
            case="case_TSI_EUR", db_name_AA="IGHV9-99*01")

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    assert "IGHV9-99,01" in str(raised.value)
    assert "more than one allele name" in str(raised.value)


def test_ambiguity_outside_the_plotted_loci_is_not_a_boot_failure(app):
    # get_db_name_from_options filters with the full plot_selection_criteria, locus
    # restriction included, so it can never resolve to an IGHD or IGHJ row. Checking those
    # made the boot check stricter than the function it protects: an ambiguous pair the app
    # never reaches and no plot ever shows would have stopped the service starting.
    # A gene not already in the mock data, so the only thing under test is the ambiguity -
    # reusing IGHD1-1*01 would add a row with no IgSNPer values beside one that has them and
    # trip that check instead.
    add_row(db_name="IGHD9-99*01", gene="IGHD9-99", allele="01", case="case_IBS_EUR")
    add_row(db_name="IGHD9-99*01_duplicate", gene="IGHD9-99", allele="01", case="case_TSI_EUR")

    validate_loaded_data()


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


def test_rejects_a_gene_from_an_unrecognised_locus(app):
    # Deriving "never translated" from PLOT_LOCI means an unknown locus is treated as
    # untranslated and quietly left out of the plots. That is consistent but silent, so the
    # locus is checked itself - and reported as an unknown locus rather than mis-reported as
    # a translation problem, which is what the hardcoded IGHD/IGHJ list used to do.
    add_row(db_name="TRGJ1*01", gene="TRGJ1", allele="01", db_name_AA=None)

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    message = str(raised.value)
    assert "no known locus: TRGJ1" in message
    assert "should be translated" not in message, "wrong diagnosis for a J locus"


def test_an_unrecognised_locus_is_reported_once(app):
    # A row from an unknown locus with db_name_AA set satisfies the amino acid check too, so
    # it used to be reported twice - "breaks 2 assumption(s)" about one row, under two
    # headings, one of which is the wrong diagnosis. The locus is the true problem.
    add_row(db_name="TRGJ1*01", gene="TRGJ1", allele="01", db_name_AA="TRGJ1*01")

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    message = str(raised.value)
    assert "breaks 1 assumption(s)" in message
    assert "no known locus: TRGJ1" in message
    assert "never translated" not in message


def test_a_gene_name_shorter_than_a_locus_prefix_is_not_misdiagnosed(app):
    # The locus is matched against KNOWN_LOCI rather than sliced at a fixed width. Sliced, a
    # short gene name yielded a truncated string that is in no list, so the crash cited a
    # locus that does not exist.
    add_row(db_name="IGH*01", gene="IGH", allele="01", db_name_AA=None)

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    assert "no known locus: IGH" in str(raised.value)


def test_rejects_divergent_igsnper_values_for_one_allele(app):
    # get_igSNPer_data answers from the first row a .distinct() query returns, with no
    # ORDER BY. That is correct because the researchers look IgSNPer values up per allele
    # name, so cohort cannot change them - but nothing enforced it, and one divergent row
    # would make a real score disappear depending on the query plan.
    add_row(db_name=TRANSLATED_ALLELE, gene="IGHV1-8", allele="01", case="case_XTRA_EUR",
            db_name_AA=TRANSLATED_ALLELE)

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    message = str(raised.value)
    assert TRANSLATED_ALLELE in message
    assert "more than one set of IgSNPer values" in message


def test_an_unknown_locus_does_not_hide_other_problems(app):
    # The de-duplication is scoped to the unknown locus's own rows. Suppressing the amino
    # acid checks outright whenever any unknown locus existed hid a genuine missing
    # db_name_AA on a known allele until the unrelated row was dealt with - the extra round
    # trip that reporting everything at once exists to avoid.
    add_row(db_name="TRGJ1*01", gene="TRGJ1", allele="01", db_name_AA=None)
    set_amino_acid_name(TRANSLATED_ALLELE, None)

    with pytest.raises(SourceDataError) as raised:
        validate_loaded_data()

    message = str(raised.value)
    assert "no known locus: TRGJ1" in message
    assert TRANSLATED_ALLELE in message, "the real amino acid problem was suppressed"
    assert "should be translated" in message


def test_a_deletion_is_recognised_by_either_spelling(app):
    # The homozygous deletions were identified two ways: allele == 'DEL' here and
    # db_name NOT LIKE '%*DEL' in the FASTA and MSA queries. A row spelled the other way was
    # reported as "should be translated" and stopped the service booting, while the FASTA
    # queries excluded it correctly.
    add_row(db_name="IGHV9-99*DEL", gene="IGHV9-99", allele="deletion", db_name_AA=None)

    validate_loaded_data()
