# Setup for pytest fixtures, sets up Flask for pytest

import os
import pytest
from app import create_app, db
from config import TestConfig
from tests import generate_mock_data
from loaders import load_tsv_to_db
from constants import ROOT_DIR

@pytest.fixture
def app():
    mock_data_file_path = ROOT_DIR + "/tests/mock_data/in/mock_allele_data.tsv"
    if not os.path.exists(mock_data_file_path):
        generate_mock_data.generate_mock_data()
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