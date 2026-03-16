# Different configurations for prod/local dev (Config) or Pytest (TestConfig)

import os

from constants import ROOT_DIR

class Config:
    TESTING = False
    PROPAGATE_EXCEPTIONS = True
    SQLALCHEMY_DATABASE_URI = os.getenv("DATABASE_URL", "sqlite:///data.db")
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    DATA_DIR = ROOT_DIR + '/data/'
    API_TITLE = "Precision Medicine Portal REST API"
    API_VERSION = "v1"
    OPENAPI_VERSION = "3.0.3"
    OPENAPI_URL_PREFIX = "/"
    OPENAPI_SWAGGER_UI_PATH = "/swagger-ui"
    OPENAPI_SWAGGER_UI_URL = "https://cdn.jsdelivr.net/npm/swagger-ui-dist/"

class TestConfig:
    TESTING = True
    PROPAGATE_EXCEPTIONS = True
    SQLALCHEMY_DATABASE_URI = 'sqlite:///:memory:'
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    DATA_DIR = ROOT_DIR + '/tests/mock_data/'
    API_TITLE = "Precision Medicine Portal REST API"
    API_VERSION = "v1"
    OPENAPI_VERSION = "3.0.3"
    OPENAPI_URL_PREFIX = "/"
    OPENAPI_SWAGGER_UI_PATH = "/swagger-ui"
    OPENAPI_SWAGGER_UI_URL = "https://cdn.jsdelivr.net/npm/swagger-ui-dist/"