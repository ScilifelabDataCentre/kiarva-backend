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


def test_a_missing_data_directory_still_names_the_archives(loader_app):
    # data/in/ is gitignored and dockerignored, so it exists only because extractall creates
    # it. With no archive to extract there was nothing to create it, and listing it raised
    # FileNotFoundError before the message that says to check data/compressed/ - unreachable
    # in exactly the case it describes.
    _, in_dir = loader_app
    in_dir.rmdir()

    with pytest.raises(SourceDataError) as raised:
        load_tsv_to_db()

    message = str(raised.value)
    assert "No .tsv files found" in message
    assert "data/compressed/" in message


def test_a_file_duplicating_its_own_rows_says_so(loader_app):
    # One commit per file means a duplicate inside a single file is not caught until that
    # commit, so it arrives as the same IntegrityError as a duplicate of an already-loaded
    # file. The message named only the second cause, which sent whoever read the pod log
    # looking for a file that had been loaded before rather than at the file in front of them.
    _, in_dir = loader_app
    write_tsv(in_dir / "self.tsv", list(range(10)) + list(range(5)))

    with pytest.raises(SourceDataError) as raised:
        load_tsv_to_db()

    message = str(raised.value)
    assert "self.tsv" in message
    assert "within itself" in message
    assert rows_from("self.tsv") == 0


def test_missing_tsv_files_raises(loader_app):
    # Was a bare quit(), which said nothing to whoever read the pod log.
    with pytest.raises(SourceDataError) as raised:
        load_tsv_to_db()

    assert "No .tsv files found" in str(raised.value)


def test_a_row_the_loader_cannot_parse_rolls_back(loader_app):
    # Only IntegrityError used to be rolled back, so a malformed row escaped with rows
    # already written into the open transaction. Nothing was committed, so no bad data
    # survived - but the next user of the session got a failed transaction instead of a
    # database.
    _, in_dir = loader_app
    # Loaded and committed on its own first, so the assertion below distinguishes a rollback
    # of the failed file from a rollback of everything.
    write_tsv(in_dir / "first.tsv", range(10))
    load_tsv_to_db()
    assert rows_from("first.tsv") == 10

    # flank_index is read with int(float(...)), so a non-numeric value raises ValueError
    # from inside load_one_tsv, after earlier rows have been written.
    # The bad row has to land after a full batch, or bulk_save_objects never runs and there
    # is nothing in the transaction to roll back.
    with open(in_dir / "broken.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for index in range(BATCH_SIZE, BATCH_SIZE * 2 + 100):
            writer.writerow(row(index))
        bad = row(BATCH_SIZE * 3)
        bad["flank_index"] = "not-a-number"
        writer.writerow(bad)

    with pytest.raises(ValueError):
        load_tsv_to_db()

    # Without the rollback the flushed rows are still visible inside the open transaction,
    # so this reads back the partial load rather than nothing.
    assert rows_from("broken.tsv") == 0
    assert rows_from("first.tsv") == 10


def test_create_app_does_not_return_an_app_when_there_is_no_data(tmp_path):
    """create_app() must propagate a load failure instead of returning a routed app.

    This is the guarantee the branch exists for: catching the error here logged it and
    returned a fully routed app anyway, so the pod passed its readiness check and answered
    requests against an empty database. Nothing else in the suite covers it - under
    TestConfig the in-memory database has no tables when create_app runs, so it takes the
    'table not found' branch and never reaches the loader, and the fixtures call
    load_tsv_to_db() directly. Hence the file-backed database: the table has to survive from
    the first create_app to the second.
    """
    data_dir = tmp_path / "data"
    (data_dir / "compressed").mkdir(parents=True)
    (data_dir / "in").mkdir()

    class FileConfig(TestConfig):
        SQLALCHEMY_DATABASE_URI = "sqlite:///" + str(tmp_path / "t.db")
        DATA_DIR = str(data_dir) + "/"

    # First call stands in for a fresh container before 'flask db upgrade': no table, so the
    # loader is skipped and an app comes back.
    first = create_app(FileConfig)
    with first.app_context():
        db.create_all()
        db.session.remove()

    # Second call is the real startup. The table is there now, the source data is not.
    with pytest.raises(SourceDataError) as raised:
        create_app(FileConfig)

    assert "No .tsv files found" in str(raised.value)
