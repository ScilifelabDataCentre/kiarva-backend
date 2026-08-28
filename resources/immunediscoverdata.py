# Flask resources file containing all endpoints related to ImmuneDiscoverData
#
# Every endpoint declares its query arguments and its response through the marshmallow
# schemas in schemas.py. flask_smorest then rejects a malformed request with a 422 before
# the view runs, serialises the response down to the declared fields, and documents both
# in the OpenAPI spec at /swagger-ui.
#
# Decorator order matters: api_key_required sits outside blp.arguments so an unauthorised
# caller is turned away before any request parsing happens.

from io import BytesIO
from flask_smorest import Blueprint, abort
from flask import current_app, send_file

from security import api_key_required
from repositories import *
from services import *

from constants import allele_superpopulation_frequencies, allele_population_frequencies, aminoacid_allele_superpopulation_frequencies, aminoacid_allele_population_frequencies
from schemas import (
    AlleleNameArgs,
    AlignedSequenceSchema,
    AminoAcidAlleleListSchema,
    AminoAcidAlleleNameArgs,
    AminoAcidGeneNameArgs,
    AminoAcidTopAlleleSchema,
    DbNameSchema,
    FastaFileNameArgs,
    FrequencySchema,
    GeneNameArgs,
    HealthSchema,
    IgSNPerSchema,
    PlotOptionArgs,
    PopulationRegionSchema,
    SelectionArgs,
    SequenceSearchArgs,
    SequenceSearchSchema,
)
from services.frequencies import create_frequencies_table

TSV_CONTENT_TYPE = "text/tab-separated-values"
FASTA_CONTENT_TYPE = "text/x-fasta"
# blp.response(content_type=...) only labels a body when there is a schema to label, so
# the downloads describe theirs through blp.doc instead. The header itself is set by
# send_file - see file_attachment below.
FILE_BODY = {"schema": {"type": "string", "format": "binary"}}

blp = Blueprint("ImmuneDiscoverData", __name__, description="Operations on ImmuneDiscover Data")

@blp.route("/health")
@blp.response(200, HealthSchema)
def health():
    return {"status": "ok"}

@blp.route("/data/db_name")
@api_key_required
@blp.arguments(SelectionArgs, location="query")
@blp.response(200, DbNameSchema)
# API to fetch the corresponding db_name (true allele name) if corresponding values from the
# "gene" and "allele" columns are supplied in the URL.
# For example:
# Rows with db_name: IGHV3-30*02/IGHV3-30-5*02 have gene: IGHV3-30-5 and allele: IGHV3-30*02/IGHV3-30-5*02.
# If the request:
#  /data/db_name?selection=IGHV3-30-5,IGHV3-30*02/IGHV3-30-5*02
# is sent to the server, the requester then gets the response
# {db_name: IGHV3-30*02/IGHV3-30-5*02}.
def get_db_name(args):
    selection = args["selection"].split(",")
    # Unpacking a selection that is not exactly gene,allele used to raise and surface as
    # a 500. The schema cannot express this, since a gene name may itself contain a comma.
    if len(selection) != 2:
        abort(422, message="'selection' must be a gene and an allele, comma separated.")

    gene, allele = selection
    return {"db_name": get_db_name_from_options(gene, allele)}

# The four frequency endpoints below all have the same shape: take the allele name from
# the validated query arguments, then either serve it from the dict pre-calculated at
# startup (prod) or calculate it on the spot (debug and testing, where the pre-load is
# skipped).
def frequency_data(allele_name, precalculated, population_type, plot_type):
    if not current_app.debug and not current_app.config.get("TESTING"):
        # A miss means the allele is not plottable: a flanking ('_F') variant, an
        # IGHD/IGHJ allele, or a name that is not in the data at all. None of those
        # are pre-calculated, and the resulting KeyError used to surface as a 500.
        if allele_name not in precalculated:
            abort(404, message="No plot data for the requested allele.")
        return precalculated[allele_name]

    # Nothing is pre-calculated under debug or pytest, so plottability is asked of the data
    # instead. Without this the same request answered 404 in prod and 200 with an all-zero
    # plot here, and the 404 could not be covered by a test.
    if not is_plottable_allele(allele_name, plot_type):
        abort(404, message="No plot data for the requested allele.")

    return calculate_frequencies(allele_name, population_type, plot_type)

@blp.route("/data/frequencies/superpopulations")
@api_key_required
@blp.arguments(AlleleNameArgs, location="query")
@blp.response(200, FrequencySchema(many=True))
@blp.alt_response(404, description="No plot data for the requested allele")
def get_superpopulation_allele_frequencies(args):
    return frequency_data(args["allele_name"], allele_superpopulation_frequencies,
                          "superpopulation", "genomic")

@blp.route("/data/aminoacidfrequencies/superpopulations")
@api_key_required
@blp.arguments(AminoAcidAlleleNameArgs, location="query")
@blp.response(200, FrequencySchema(many=True))
@blp.alt_response(404, description="No plot data for the requested allele")
def get_superpopulation_aminoacid_frequencies(args):
    return frequency_data(args["aa_allele_name"], aminoacid_allele_superpopulation_frequencies,
                          "superpopulation", "aminoacid")

@blp.route("/data/frequencies/populations")
@api_key_required
@blp.arguments(AlleleNameArgs, location="query")
@blp.response(200, FrequencySchema(many=True))
@blp.alt_response(404, description="No plot data for the requested allele")
def get_subpopulation_allele_frequencies(args):
    return frequency_data(args["allele_name"], allele_population_frequencies,
                          "population", "genomic")

@blp.route("/data/aminoacidfrequencies/populations")
@api_key_required
@blp.arguments(AminoAcidAlleleNameArgs, location="query")
@blp.response(200, FrequencySchema(many=True))
@blp.alt_response(404, description="No plot data for the requested allele")
def get_subpopulation_aminoacid_frequencies(args):
    return frequency_data(args["aa_allele_name"], aminoacid_allele_population_frequencies,
                          "population", "aminoacid")

# The four table downloads and the three FASTA downloads return a file rather than JSON, so
# they have no response schema. The content type is passed to send_file because that is what
# actually sets the header - left to mimetypes, '.fasta' is unrecognised and the downloads
# went out as application/octet-stream.
def file_attachment(body, download_name, content_type):
    # send_file expects bytes rather than str
    buffer = BytesIO()
    buffer.write(str.encode(body))
    buffer.seek(0)
    return send_file(buffer, as_attachment=True, download_name=download_name,
                     mimetype=content_type)

@blp.route("/data/frequencies/table/allele")
@api_key_required
@blp.arguments(AlleleNameArgs, location="query")
@blp.response(200, description="Frequency table as a tab separated file")
@blp.doc(responses={200: {"content": {TSV_CONTENT_TYPE: FILE_BODY}}})
def get_allele_frequencies_table(args):
    allele_name = args["allele_name"]
    return file_attachment(create_frequencies_table(allele_name, "genomic"),
                           allele_name + '_frequencies_genomic.tsv', TSV_CONTENT_TYPE)

@blp.route("/data/aminoacidfrequencies/table/allele")
@api_key_required
@blp.arguments(AminoAcidAlleleNameArgs, location="query")
@blp.response(200, description="Frequency table as a tab separated file")
@blp.doc(responses={200: {"content": {TSV_CONTENT_TYPE: FILE_BODY}}})
def get_aminoacid_allele_frequencies_table(args):
    aa_allele_name = args["aa_allele_name"]
    return file_attachment(create_frequencies_table(aa_allele_name, "aminoacid"),
                           aa_allele_name + '_frequencies_aminoacid.tsv', TSV_CONTENT_TYPE)

@blp.route("/data/frequencies/table/gene")
@api_key_required
@blp.arguments(GeneNameArgs, location="query")
@blp.response(200, description="Frequency table as a tab separated file")
@blp.doc(responses={200: {"content": {TSV_CONTENT_TYPE: FILE_BODY}}})
def get_gene_frequencies_table(args):
    gene_name = args["gene_name"]
    return file_attachment(create_frequencies_table(gene_name, "genomic", full_gene=True),
                           gene_name + '_frequencies_genomic.tsv', TSV_CONTENT_TYPE)

@blp.route("/data/aminoacidfrequencies/table/gene")
@api_key_required
@blp.arguments(AminoAcidGeneNameArgs, location="query")
@blp.response(200, description="Frequency table as a tab separated file")
@blp.doc(responses={200: {"content": {TSV_CONTENT_TYPE: FILE_BODY}}})
def get_aminoacid_gene_frequencies_table(args):
    aa_gene_name = args["aa_gene_name"]
    return file_attachment(create_frequencies_table(aa_gene_name, "aminoacid", full_gene=True),
                           aa_gene_name + '_frequencies_aminoacid.tsv', TSV_CONTENT_TYPE)

@blp.route("/data/igsnperdata")
@api_key_required
@blp.arguments(AlleleNameArgs, location="query")
@blp.response(200, IgSNPerSchema)
@blp.alt_response(404, description="No IgSNPer data for the requested allele")
def get_igsnper_data(args):
    data = get_igSNPer_data(args["allele_name"])
    if not data:
        abort(404, message="No IgSNPer data for the requested allele.")
    return data

@blp.route("/data/aminoacidalleles")
@api_key_required
@blp.arguments(AminoAcidAlleleNameArgs, location="query")
@blp.response(200, AminoAcidTopAlleleSchema)
def get_aa_top_allele(args):
    return get_aminoacid_top_allele(args["aa_allele_name"])

@blp.route("/data/aminoacidlist")
@api_key_required
@blp.arguments(AminoAcidAlleleNameArgs, location="query")
@blp.response(200, AminoAcidAlleleListSchema)
def get_aminoacid_list(args):
    return get_aminoacid_allele_list(args["aa_allele_name"])

@blp.route("/data/populationregions")
@api_key_required
@blp.response(200, PopulationRegionSchema(many=True))
def get_population_regions():
    return get_populations()

@blp.route("/data/plotoptions", methods=["GET"])
@api_key_required
@blp.arguments(PlotOptionArgs, location="query")
# Responds with a bare list of names rather than a list of objects, which a marshmallow
# schema cannot describe, so the response is documented without one.
@blp.response(200, description="Selectable names for the next part of the selection")
def get_next_selection_option(args):
    return get_plot_options(args["current_selection"])

@blp.route("/data/sequences/alignedsequences")
@api_key_required
@blp.arguments(GeneNameArgs, location="query")
@blp.response(200, AlignedSequenceSchema(many=True))
def get_aligned_sequences(args):
    return align_sequences(args["gene_name"])

@blp.route("/data/sequences")
@api_key_required
@blp.arguments(SequenceSearchArgs, location="query")
@blp.response(200, SequenceSearchSchema(many=True))
def get_sequence_search(args):
    return sequence_search(args["sequence_str"])

@blp.route("/fasta/genomic")
@api_key_required
@blp.arguments(FastaFileNameArgs, location="query")
@blp.response(200, description="Sequences in FASTA format")
@blp.doc(responses={200: {"content": {FASTA_CONTENT_TYPE: FILE_BODY}}})
def send_genomic_fasta(args):
    file_name = args["file_name"]
    return file_attachment(generate_fasta(file_name, type="genomic"),
                           file_name + '_genomic.fasta', FASTA_CONTENT_TYPE)

@blp.route("/fasta/genomic_fl")
@api_key_required
@blp.arguments(FastaFileNameArgs, location="query")
@blp.response(200, description="Sequences in FASTA format")
@blp.doc(responses={200: {"content": {FASTA_CONTENT_TYPE: FILE_BODY}}})
def send_flanking_genomic_fasta(args):
    file_name = args["file_name"]
    return file_attachment(generate_fasta(file_name, type="genomic_fl"),
                           file_name + '_genomic_fl.fasta', FASTA_CONTENT_TYPE)

@blp.route("/fasta/translated")
@api_key_required
@blp.arguments(FastaFileNameArgs, location="query")
@blp.response(200, description="Sequences in FASTA format")
@blp.doc(responses={200: {"content": {FASTA_CONTENT_TYPE: FILE_BODY}}})
def send_translated_fasta(args):
    file_name = args["file_name"]
    return file_attachment(generate_fasta(file_name, type="translated"),
                           file_name + '_translated.fasta', FASTA_CONTENT_TYPE)

@blp.route("/checkapikey")
@api_key_required
# Deliberately undecorated: this returns a bare string, and blp.response would JSON-encode
# it into '"Correct key!"' with quotes, which is a breaking change for the caller.
def check_api_key():
    return "Correct key!"
