import os
from dotenv import load_dotenv

from flask import Flask
from flask_smorest import Api
from flask_migrate import Migrate
from flask_cors import CORS
from sqlalchemy import column
from sqlalchemy.exc import OperationalError

from db import db

from models.immunediscoverdata import ImmuneDiscoverDataModel
from resources.immunediscoverdata import blp as ImmuneDiscoverDataBlueprint
from utils import load_tsv_to_db

def create_app(db_url=None):
    load_dotenv(override=True)

    app = Flask(__name__)
    # CORS(app, origins=[os.getenv("FRONTEND_URL")])

    app.config["PROPAGATE_EXCEPTIONS"] = True
    app.config["API_TITLE"] = "Precision Medicine Portal REST API"
    app.config["API_VERSION"] = "v1"
    app.config["OPENAPI_VERSION"] = "3.0.3"
    app.config["OPENAPI_URL_PREFIX"] = "/"
    app.config["OPENAPI_SWAGGER_UI_PATH"] = "/swagger-ui"
    app.config["OPENAPI_SWAGGER_UI_URL"] = "https://cdn.jsdelivr.net/npm/swagger-ui-dist/"
    app.config["SQLALCHEMY_DATABASE_URI"] = db_url or os.getenv("DATABASE_URL", "sqlite:///data.db")
    app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False

    db.init_app(app)
    migrate = Migrate(app, db)
    api = Api(app)

    api.register_blueprint(ImmuneDiscoverDataBlueprint)

    with app.app_context():
        try:
            data = ImmuneDiscoverDataModel.query.all()
            if len(data) == 0:
                load_tsv_to_db()
        except OperationalError:
            print("---DB not initialized---")

    return app