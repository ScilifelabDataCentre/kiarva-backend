# Tests for loaders/load_tsv_to_db.py.
#
# The mock dataset is 36 rows, too small to reach the 1000-row batch boundary, so these
# generate their own TSVs. That boundary is the point of the test: rows used to be committed
# as each batch filled, and a duplicate arriving after one of those commits left the earlier
# rows in the database with the file recorded in loaded_from_tsv - so every later startup
# skipped it and the data stayed silently truncated.

import csv

import pytest

from app import create_app, db
from config import TestConfig
from loaders import load_tsv_to_db, SourceDataError
from models.immunediscoverdata import ImmuneDiscoverDataModel

# The loader reads these by name; the rest of the source columns are ignored.
COLUMNS = [
    "cohort", "case", "db_name", "gene", "allele", "sequence", "flank_index",
    "IgSNPer_uncommon", "IgSNPer_SNPs", "db_name_AA", "db_name_AA_list", "sequence_AA",
]

BATCH_SIZE = 1000


def row(index):
    """A valid row, unique by case, which is what the table's unique constraint keys on."""
    return {
        "cohort": "1KGP",
        # The loader splits case on "_" to get population and superpopulation.
        "case": f"case{index:05d}_GBR_EUR",
        "db_name": "IGHV1-2*01",
        "gene": "IGHV1-2",
        "allele": "01",
        "sequence": "ACGT",
        "flank_index": "1",
        "IgSNPer_uncommon": "",
        "IgSNPer_SNPs": "",
        "db_name_AA": "IGHV1-2*01",
        "db_name_AA_list": "IGHV1-2*01",
        "sequence_AA": "AC",
    }


def write_tsv(path, indices):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for index in indices:
            writer.writerow(row(index))


@pytest.fixture
def loader_app(tmp_path):
    data_dir = tmp_path / "data"
    (data_dir / "in").mkdir(parents=True)
    (data_dir / "compressed").mkdir()

    app = create_app(TestConfig)
    app.config["DATA_DIR"] = str(data_dir) + "/"
    with app.app_context():
        db.create_all()
        yield app, data_dir / "in"
        db.session.remove()
        db.drop_all()


def rows_from(file_name):
    return ImmuneDiscoverDataModel.query.filter_by(loaded_from_tsv=file_name).count()


def test_loads_a_file_larger_than_one_batch(loader_app):
    _, in_dir = loader_app
    write_tsv(in_dir / "first.tsv", range(BATCH_SIZE + 100))

    load_tsv_to_db()

    assert rows_from("first.tsv") == BATCH_SIZE + 100


def test_a_file_failing_after_a_full_batch_leaves_nothing_behind(loader_app):
    _, in_dir = loader_app
    write_tsv(in_dir / "first.tsv", range(BATCH_SIZE))
    load_tsv_to_db()
    assert rows_from("first.tsv") == BATCH_SIZE

    # Enough new rows to fill and flush a batch, then rows duplicating the first file.
    write_tsv(in_dir / "second.tsv",
              list(range(BATCH_SIZE, BATCH_SIZE * 2 + 100)) + list(range(50)))

    with pytest.raises(SourceDataError) as raised:
        load_tsv_to_db()
    assert "second.tsv" in str(raised.value)

    # Nothing of the failed file survives, so loaded_from_tsv does not record it and the
    # next startup reads it again instead of skipping it as already loaded.
    assert rows_from("second.tsv") == 0
    assert rows_from("first.tsv") == BATCH_SIZE

    loaded = {name for (name,) in ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.loaded_from_tsv).distinct().all()}
    assert loaded == {"first.tsv"}


def test_missing_tsv_files_raises(loader_app):
    # Was a bare quit(), which said nothing to whoever read the pod log.
    with pytest.raises(SourceDataError) as raised:
        load_tsv_to_db()

    assert "No .tsv files found" in str(raised.value)
