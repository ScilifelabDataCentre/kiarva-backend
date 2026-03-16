# Flask resources file containing all endpoints related to ImmuneDiscoverData

from io import BytesIO
from flask_smorest import Blueprint
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

@blp.route("/data/frequencies/superpopulations")
@api_key_required
def get_superpopulation_allele_frequencies():
    allele_name = request.args.get("allele_name")
    if not current_app.debug and not current_app.config.get("TESTING"):
        data_out = allele_superpopulation_frequencies[allele_name]
    else:
        data_out = calculate_frequencies(allele_name, "superpopulation", "genomic")
    return data_out

@blp.route("/data/aminoacidfrequencies/superpopulations")
@api_key_required
def get_superpopulation_aminoacid_frequencies():
    aa_allele_name = request.args.get("aa_allele_name")
    if not current_app.debug and not current_app.config.get("TESTING"):
        data_out = aminoacid_allele_superpopulation_frequencies[aa_allele_name]
    else:
        data_out = calculate_frequencies(aa_allele_name, "superpopulation", "aminoacid")
    return data_out
    
@blp.route("/data/frequencies/populations")
@api_key_required
def get_subpopulation_allele_frequencies():
    allele_name = request.args.get("allele_name")
    if not current_app.debug and not current_app.config.get("TESTING"):
        data_out = allele_population_frequencies[allele_name]
    else:
        data_out = calculate_frequencies(allele_name, "population", "genomic")
    return data_out
    
@blp.route("/data/aminoacidfrequencies/populations")
@api_key_required
def get_subpopulation_aminoacid_frequencies():
    aa_allele_name = request.args.get("aa_allele_name")
    if not current_app.debug and not current_app.config.get("TESTING"):
        data_out = aminoacid_allele_population_frequencies[aa_allele_name]
    else:
        data_out = calculate_frequencies(aa_allele_name, "population", "aminoacid")
    return data_out

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
    print(data_out)
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