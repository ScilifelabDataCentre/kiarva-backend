# Flask resources file containing all endpoints related to ImmuneDiscoverData

from io import BytesIO
from flask_smorest import Blueprint, abort
from flask import current_app, request, send_file

from security import api_key_required
from repositories import *
from services import *

from constants import allele_superpopulation_frequencies, allele_population_frequencies, aminoacid_allele_superpopulation_frequencies, aminoacid_allele_population_frequencies
from services.frequencies import create_frequencies_table


blp = Blueprint("ImmuneDiscoverData", __name__, description="Operations on ImmuneDiscover Data")

@blp.route("/health")
def health():
    return {"status": "ok"}, 200

@blp.route("/data/db_name")
@api_key_required
# API to fetch the corresponding db_name (true allele name) if corresponding values from the
# "gene" and "allele" columns are supplied in the URL.
# For example: 
# Rows with db_name: IGHV3-30*02/IGHV3-30-5*02 have gene: IGHV3-30-5 and allele: IGHV3-30*02/IGHV3-30-5*02.
# If the request:
#  /data/db_name?selection=IGHV3-30-5,IGHV3-30*02/IGHV3-30-5*02
# is sent to the server, the requester then gets the response
# {db_name: IGHV3-30*02/IGHV3-30-5*02}.
def get_db_name():
    gene, allele = request.args.get("selection").split(",")
    db_name = get_db_name_from_options(gene, allele)
    return {"db_name": db_name}

# The four frequency endpoints below all have the same shape: read an allele name
# from the query string, then either serve it from the dict pre-calculated at
# startup (prod) or calculate it on the spot (debug and testing, where the
# pre-load is skipped).
def frequency_data(allele_name, param_name, precalculated, population_type, plot_type):
    if not allele_name:
        abort(400, message="Missing required query parameter '" + param_name + "'.")

    if not current_app.debug and not current_app.config.get("TESTING"):
        # A miss means the allele is not plottable: a flanking ('_F') variant, an
        # IGHD/IGHJ allele, or a name that is not in the data at all. None of those
        # are pre-calculated, and the resulting KeyError used to surface as a 500.
        if allele_name not in precalculated:
            abort(404, message="No plot data for allele '" + allele_name + "'.")
        return precalculated[allele_name]

    return calculate_frequencies(allele_name, population_type, plot_type)

@blp.route("/data/frequencies/superpopulations")
@api_key_required
def get_superpopulation_allele_frequencies():
    return frequency_data(request.args.get("allele_name"), "allele_name",
                          allele_superpopulation_frequencies, "superpopulation", "genomic")

@blp.route("/data/aminoacidfrequencies/superpopulations")
@api_key_required
def get_superpopulation_aminoacid_frequencies():
    return frequency_data(request.args.get("aa_allele_name"), "aa_allele_name",
                          aminoacid_allele_superpopulation_frequencies, "superpopulation", "aminoacid")

@blp.route("/data/frequencies/populations")
@api_key_required
def get_subpopulation_allele_frequencies():
    return frequency_data(request.args.get("allele_name"), "allele_name",
                          allele_population_frequencies, "population", "genomic")

@blp.route("/data/aminoacidfrequencies/populations")
@api_key_required
def get_subpopulation_aminoacid_frequencies():
    return frequency_data(request.args.get("aa_allele_name"), "aa_allele_name",
                          aminoacid_allele_population_frequencies, "population", "aminoacid")

@blp.route("/data/frequencies/table/allele")
@api_key_required
def get_subpopulation_allele_frequencies():
    allele_name = request.args.get("allele_name")
    buffer = BytesIO()
    buffer.write(str.encode(create_frequencies_table(allele_name, "genomic")))
    buffer.seek(0)
    return send_file(
        buffer,
        as_attachment=True,
        download_name=allele_name + '_frequencies_genomic.tsv'
    )

@blp.route("/data/aminoacidfrequencies/table/allele")
@api_key_required
def get_subpopulation_allele_frequencies():
    aa_allele_name = request.args.get("aa_allele_name")
    buffer = BytesIO()
    buffer.write(str.encode(create_frequencies_table(aa_allele_name, "aminoacid")))
    buffer.seek(0)
    return send_file(
        buffer,
        as_attachment=True,
        download_name=aa_allele_name + '_frequencies_aminoacid.tsv'
    )

@blp.route("/data/frequencies/table/gene")
@api_key_required
def get_subpopulation_allele_frequencies():
    gene_name = request.args.get("gene_name")
    buffer = BytesIO()
    buffer.write(str.encode(create_frequencies_table(gene_name, "genomic", full_gene=True)))
    buffer.seek(0)
    return send_file(
        buffer,
        as_attachment=True,
        download_name=gene_name + '_frequencies_genomic.tsv'
    )

@blp.route("/data/aminoacidfrequencies/table/gene")
@api_key_required
def get_subpopulation_allele_frequencies():
    aa_gene_name = request.args.get("aa_gene_name")
    buffer = BytesIO()
    buffer.write(str.encode(create_frequencies_table(aa_gene_name, "aminoacid", full_gene=True)))
    buffer.seek(0)
    return send_file(
        buffer,
        as_attachment=True,
        download_name=aa_gene_name + '_frequencies_aminoacid.tsv'
    )
    
@blp.route("/data/igsnperdata")
@api_key_required
def get_igsnper_data():
    allele_name = request.args.get("allele_name")
    data_out = get_igSNPer_data(allele_name)
    return data_out

@blp.route("/data/aminoacidalleles")
@api_key_required
def get_aa_top_allele():
    aa_allele_name = request.args.get("aa_allele_name")
    data_out = get_aminoacid_top_allele(aa_allele_name)
    return data_out
    
@blp.route("/data/aminoacidlist")
@api_key_required
def get_aminoacid_list():
    aa_allele_name = request.args.get("aa_allele_name")
    data_out = get_aminoacid_allele_list(aa_allele_name)
    return data_out
    
@blp.route("/data/populationregions")
@api_key_required
def get_population_regions():
    data_out = get_populations()
    return data_out
    

@blp.route("/data/plotoptions", methods=["GET"])
@api_key_required
def get_next_selection_option():
    gene = request.args.get("current_selection")
    data_out = get_plot_options(gene)
    return data_out

@blp.route("/data/sequences/alignedsequences")
@api_key_required
def get_aligned_sequences():
    gene = request.args.get("gene_name")
    aligned_seqs = align_sequences(gene)
    return aligned_seqs

@blp.route("/data/sequences")
@api_key_required
def get_sequence_search():
    sequence_str = request.args.get("sequence_str")
    data_out = sequence_search(sequence_str)
    return data_out

@blp.route("/fasta/genomic")
@api_key_required
def send_fasta():
    file_name = request.args.get("file_name")
    # send_file expects bytes rather than str
    buffer = BytesIO()
    buffer.write(str.encode(generate_fasta(file_name, type="genomic")))
    buffer.seek(0)
    return send_file(
        buffer,
        as_attachment=True,
        download_name=file_name + '_genomic.fasta'
    )

@blp.route("/fasta/genomic_fl")
@api_key_required
def send_fasta():
    file_name = request.args.get("file_name")
    # send_file expects bytes rather than str
    buffer = BytesIO()
    buffer.write(str.encode(generate_fasta(file_name, type="genomic_fl")))
    buffer.seek(0)
    return send_file(
        buffer,
        as_attachment=True,
        download_name=file_name + '_genomic_fl.fasta'
    )

@blp.route("/fasta/translated")
@api_key_required
def send_fasta():
    file_name = request.args.get("file_name")
    # send_file expects bytes rather than str
    buffer = BytesIO()
    buffer.write(str.encode(generate_fasta(file_name, type="translated")))
    buffer.seek(0)
    return send_file(
        buffer,
        as_attachment=True,
        download_name=file_name + '_translated.fasta'
    )

@blp.route("/checkapikey")
@api_key_required
def check_api_key():
    return "Correct key!"