from io import BytesIO
from flask.views import MethodView
from flask_smorest import Blueprint, abort
# from sqlalchemy.exc import SQLAlchemyError, IntegrityError
from flask import current_app, send_file
from werkzeug.utils import safe_join

from constants import ROOT_DIR
from models import ImmuneDiscoverDataModel
from schemas import ImmuneDiscoverDataFrequencySchema, ImmuneDiscoverDataGetAllSchema, ImmuneDiscoverPopulationRegionSchema

# from security import api_key_required
from utils import calculate_allele_frequencies
from utils.generate_fasta import generate_fasta


blp = Blueprint("ImmuneDiscoverData", __name__, description="Operations on ImmuneDiscover Data")

@blp.route("/data")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataGetAllSchema(many=True))
    def get(self):
        data = ImmuneDiscoverDataModel.query.all()
        return data
    
@blp.route("/data/<gene_segment>")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataGetAllSchema(many=True))
    def get(self, gene_segment):
        data = ImmuneDiscoverDataModel.query.distinct().filter(ImmuneDiscoverDataModel.db_name.like(gene_segment+'%')).all()
        return data

@blp.route("/data/frequencies/superpopulations/<allele_name>")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataFrequencySchema(many=True))
    def get(self, allele_name):        
        return calculate_allele_frequencies(allele_name, "superpopulation")
    
@blp.route("/data/frequencies/populations/<allele_name>")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataFrequencySchema(many=True))
    def get(self, allele_name):
        return calculate_allele_frequencies(allele_name, "population")
    
@blp.route("/data/populationregions")
class ImmuneDiscoverPopulationDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverPopulationRegionSchema(many=True))
    def get(self):
        data = ImmuneDiscoverDataModel.query.with_entities(
                ImmuneDiscoverDataModel.superpopulation,
                ImmuneDiscoverDataModel.population,
                ).distinct().all()
        return data

@blp.route("/fasta/<file_name>")
def send_fasta(file_name):
    # send_file expects bytes rather than str
    buffer = BytesIO()
    buffer.write(str.encode(generate_fasta(file_name)))
    buffer.seek(0)
    return send_file(
        buffer,
        as_attachment=True,
        download_name=file_name + '.fasta'
    )

@blp.route("/fasta/genomic/<file_name>")
def send_fasta(file_name):
    # send_file expects bytes rather than str
    buffer = BytesIO()
    buffer.write(str.encode(generate_fasta(file_name, genomic=True)))
    buffer.seek(0)
    return send_file(
        buffer,
        as_attachment=True,
        download_name=file_name + '_genomic.fasta'
    )   