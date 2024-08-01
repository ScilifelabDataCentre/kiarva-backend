from io import BytesIO
import time
from flask.views import MethodView
from flask_smorest import Blueprint, abort
# from sqlalchemy.exc import SQLAlchemyError, IntegrityError
from flask import current_app, send_file
from werkzeug.utils import safe_join

from constants import ROOT_DIR
from models import ImmuneDiscoverDataModel
from schemas import ImmuneDiscoverDataFrequencySchema, ImmuneDiscoverDataGetAllSchema, ImmuneDiscoverIgSNPerDataSchema, ImmuneDiscoverPopulationRegionSchema

# from security import api_key_required
from utils import calculate_allele_frequencies
from utils.generate_fasta import generate_fasta
from utils.igSNPer_data import get_igSNPer_data


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
        print("---")
        time_start = time.time()
        data_out = calculate_allele_frequencies(allele_name, "superpopulation")
        time_end = time.time()
        print("superpopulations frequency calc delta: " + str((time_end-time_start)*1000) + " ms")
        print("---")
        return data_out
    
@blp.route("/data/frequencies/populations/<allele_name>")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverDataFrequencySchema(many=True))
    def get(self, allele_name):
        print("---")
        time_start = time.time()
        data_out = calculate_allele_frequencies(allele_name, "population")
        time_end = time.time()
        print("subpopulations frequency calc delta: " + str((time_end-time_start)*1000) + " ms")
        print("---")
        return data_out
    
@blp.route("/data/igsnperdata/<allele_name>")
class ImmuneDiscoverDataList(MethodView):
    # @api_key_required
    @blp.response(200, ImmuneDiscoverIgSNPerDataSchema)
    def get(self, allele_name):
        data_out = get_igSNPer_data(allele_name)
        return data_out
    
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