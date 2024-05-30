import json
from flask.views import MethodView
from flask_smorest import Blueprint, abort
from sqlalchemy.exc import SQLAlchemyError, IntegrityError

from db import db
from models import ImmuneDiscoverDataModel
from schemas import ImmuneDiscoverDataGetSchema, ImmuneDiscoverDataUploadSchema

from security import api_key_required
from utils import load_tsv_to_db


blp = Blueprint("ImmuneDiscoverData", __name__, description="Operations on ImmuneDiscover Data")

@blp.route("/data")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataUploadSchema(many=True))
    def get(self):
        # load_tsv_to_db()
        data = ImmuneDiscoverDataModel.query.all()
        # data_out = []
        # for row in data:
        #     data_out.append({'case': row.case, 'db_name': row.db_name, 'sequence': row.sequence})
        return data