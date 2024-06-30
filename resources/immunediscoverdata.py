from flask.views import MethodView
from flask_smorest import Blueprint, abort
# from sqlalchemy.exc import SQLAlchemyError, IntegrityError
from flask import send_from_directory
from flask import send_file

from constants import ROOT_DIR
from models import ImmuneDiscoverDataModel
from schemas import ImmuneDiscoverDataFrequencySchema, ImmuneDiscoverDataGetAllSchema

# from security import api_key_required
from utils import calculate_allele_frequencies


blp = Blueprint("ImmuneDiscoverData", __name__, description="Operations on ImmuneDiscover Data")

@blp.route("/data")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataGetAllSchema(many=True))
    def get(self):
        data = ImmuneDiscoverDataModel.query.all()
        # data_out = []
        # for row in data:
        #     data_out.append({'case': row.case, 'db_name': row.db_name, 'sequence': row.sequence})
        return data

@blp.route("/data/frequencies/superpopulations/<allele_name>")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataFrequencySchema(many=True))
    def get(self, allele_name):        
        return calculate_allele_frequencies(allele_name, True)
    
@blp.route("/data/frequencies/populations/<allele_name>")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataFrequencySchema(many=True))
    def get(self, allele_name):
        return calculate_allele_frequencies(allele_name)
    
@blp.route("/fasta/<file_name>")
def send_fasta(file_name):
    data_out_path = ROOT_DIR + '/data/out/'
    return send_file(data_out_path + file_name, as_attachment=True)
        