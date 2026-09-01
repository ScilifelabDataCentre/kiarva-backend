# Setup for pytest fixtures, sets up Flask for pytest

import os
import pytest
from app import create_app, db
from config import TestConfig
from tests import generate_mock_data
from loaders import load_tsv_to_db
from constants import ROOT_DIR

@pytest.fixture
def app(monkeypatch):
    mock_data_file_path = ROOT_DIR + "/tests/mock_data/in/mock_allele_data.tsv"
    if not os.path.exists(mock_data_file_path):
        generate_mock_data.generate_mock_data()
    # create_app() refuses to return an app when the table is missing, which is the state
    # every test starts from - the schema is created below, immediately after. SKIP_DATA_LOAD
    # is the same flag docker/entrypoint.sh sets for 'flask db upgrade', for the same reason:
    # this call is not the server. It stays set for the rest of the test, so a create_app()
    # inside one does not try to load either.
    monkeypatch.setenv("SKIP_DATA_LOAD", "1")
    app = create_app(TestConfig)
    with app.app_context():
        db.create_all()
        load_tsv_to_db()
        yield app
        db.session.remove()
        db.drop_all()

@pytest.fixture
def client(app):
    return app.test_client()

@pytest.fixture
def db_session(app):
    yield db.session
    db.session.rollback()