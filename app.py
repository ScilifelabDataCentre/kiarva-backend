import os
from dotenv import load_dotenv

from flask import Flask, request
from flask_smorest import Api
from flask_migrate import Migrate
from flask_cors import CORS
import pandas as pd
from sqlalchemy import column
from sqlalchemy.exc import OperationalError, IntegrityError

from constants import ROOT_DIR
from db import db

from models.immunediscoverdata import ImmuneDiscoverDataModel
from resources.immunediscoverdata import blp as ImmuneDiscoverDataBlueprint
from utils import load_tsv_to_db, write_fastas

def create_app(db_url=None):
    load_dotenv(override=True)

    app = Flask(__name__)
    # CORS(app, origins=[os.getenv("FRONTEND_URL")])
    CORS(app, origins=["*"])

    app.config["PROPAGATE_EXCEPTIONS"] = True
    app.config["API_TITLE"] = "Precision Medicine Portal REST API"
    app.config["API_VERSION"] = "v1"
    app.config["OPENAPI_VERSION"] = "3.0.3"
    app.config["OPENAPI_URL_PREFIX"] = "/"
    app.config["OPENAPI_SWAGGER_UI_PATH"] = "/swagger-ui"
    app.config["OPENAPI_SWAGGER_UI_URL"] = "https://cdn.jsdelivr.net/npm/swagger-ui-dist/"
    app.config["SQLALCHEMY_DATABASE_URI"] = db_url or os.getenv("DATABASE_URL", "sqlite:///data.db")
    app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
    app.config["FASTA_DIR"] = ROOT_DIR + "/data/out/"

    db.init_app(app)
    migrate = Migrate(app, db)
    api = Api(app)

    api.register_blueprint(ImmuneDiscoverDataBlueprint)

    # loads tsv data to db only if db is initialized and upgraded and tsv data has not been loaded yet
    # writes fasta files from db data to /data/out/ once db is loaded
    with app.app_context():
        data_dir = ROOT_DIR + '/data/'
        tsv_files = [file for file in os.listdir(data_dir+'in/') if file.endswith('.tsv')]

        try:
            files_in_db = ImmuneDiscoverDataModel.query.with_entities(ImmuneDiscoverDataModel.loaded_from_tsv).distinct().all()

            loaded_files = [loaded_file[0] for loaded_file in files_in_db]
                    
            for file in tsv_files:
                if file not in loaded_files:
                    load_tsv_to_db(file)
            
            df = pd.read_sql_table('immunediscoverdata', db.engine)
            write_fastas(df, data_dir +'out/')
        except OperationalError:
            print("---DB not initialized---")

    return app