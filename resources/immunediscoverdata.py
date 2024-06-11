import json
from flask.views import MethodView
from flask_smorest import Blueprint, abort
from sqlalchemy.exc import SQLAlchemyError, IntegrityError
from flask import send_from_directory

from constants import ROOT_DIR
from models import ImmuneDiscoverDataModel
from schemas import ImmuneDiscoverDataUploadSchema

from security import api_key_required


blp = Blueprint("ImmuneDiscoverData", __name__, description="Operations on ImmuneDiscover Data")

@blp.route("/data")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataUploadSchema(many=True))
    def get(self):
        data = ImmuneDiscoverDataModel.query.all()
        # data_out = []
        # for row in data:
        #     data_out.append({'case': row.case, 'db_name': row.db_name, 'sequence': row.sequence})
        return data
    
@blp.route("/fastas/<path:name>")
def send_fasta(name):
    data_out_path = ROOT_DIR + '/data/out/'
    return send_from_directory(data_out_path, name)
        