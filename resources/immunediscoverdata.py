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
from flask import send_file

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
# The api key requirement is declared at document root, so every operation inherits it.
# This one has no api_key_required, and a Kubernetes probe wired from the spec would
# otherwise be built to send a header it does not need.
@blp.doc(security=[])
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

# The four frequency endpoints below all have the same shape: take the allele name from the
# validated query arguments and serve it from the dictionary pre-calculated at startup.
#
# There is no separate on-demand path any more. The pre-load runs in every mode, so a
# request is answered the same way under pytest, under 'flask run --debug' and in
# production - which is what makes the miss below reachable by a test at all. It used to be
# prod-only, so the same request 404d there and returned an all-zero plot everywhere else.
def frequency_data(allele_name, precalculated):
    # A miss means the allele is not plottable: a flanking ('_F') variant, an IGHD/IGHJ
    # allele, or a name that is not in the data at all.
    if allele_name not in precalculated:
        abort(404, message="No plot data for the requested allele.")

    return precalculated[allele_name]

@blp.route("/data/frequencies/superpopulations")
@api_key_required
@blp.arguments(AlleleNameArgs, location="query")
@blp.response(200, FrequencySchema(many=True))
@blp.alt_response(404, description="No plot data for the requested allele")
def get_superpopulation_allele_frequencies(args):
    return frequency_data(args["allele_name"], allele_superpopulation_frequencies)

@blp.route("/data/aminoacidfrequencies/superpopulations")
@api_key_required
@blp.arguments(AminoAcidAlleleNameArgs, location="query")
@blp.response(200, FrequencySchema(many=True))
@blp.alt_response(404, description="No plot data for the requested allele")
def get_superpopulation_aminoacid_frequencies(args):
    return frequency_data(args["aa_allele_name"], aminoacid_allele_superpopulation_frequencies)

@blp.route("/data/frequencies/populations")
@api_key_required
@blp.arguments(AlleleNameArgs, location="query")
@blp.response(200, FrequencySchema(many=True))
@blp.alt_response(404, description="No plot data for the requested allele")
def get_subpopulation_allele_frequencies(args):
    return frequency_data(args["allele_name"], allele_population_frequencies)

@blp.route("/data/aminoacidfrequencies/populations")
@api_key_required
@blp.arguments(AminoAcidAlleleNameArgs, location="query")
@blp.response(200, FrequencySchema(many=True))
@blp.alt_response(404, description="No plot data for the requested allele")
def get_subpopulation_aminoacid_frequencies(args):
    return frequency_data(args["aa_allele_name"], aminoacid_allele_population_frequencies)

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

# The four frequency table downloads. create_frequencies_table returns None when the
# requested allele or gene is not one the plots cover, which is the same set the download is
# offered for in the frontend - the option only appears once an allele has been selected to
# plot.
def frequency_table_attachment(name, plot_type, download_suffix, full_gene = False):
    table = create_frequencies_table(name, plot_type, full_gene=full_gene)
    if table is None:
        abort(404, message="No frequency table for the requested "
                           + ("gene" if full_gene else "allele") + ".")

    return file_attachment(table, name + download_suffix, TSV_CONTENT_TYPE)

@blp.route("/data/frequencies/table/allele")
@api_key_required
@blp.arguments(AlleleNameArgs, location="query")
@blp.response(200, description="Frequency table as a tab separated file")
@blp.doc(responses={200: {"content": {TSV_CONTENT_TYPE: FILE_BODY}}})
@blp.alt_response(404, description="No frequency table for the requested allele or gene")
def get_allele_frequencies_table(args):
    return frequency_table_attachment(args["allele_name"], "genomic",
                                      '_frequencies_genomic.tsv')

@blp.route("/data/aminoacidfrequencies/table/allele")
@api_key_required
@blp.arguments(AminoAcidAlleleNameArgs, location="query")
@blp.response(200, description="Frequency table as a tab separated file")
@blp.doc(responses={200: {"content": {TSV_CONTENT_TYPE: FILE_BODY}}})
@blp.alt_response(404, description="No frequency table for the requested allele or gene")
def get_aminoacid_allele_frequencies_table(args):
    return frequency_table_attachment(args["aa_allele_name"], "aminoacid",
                                      '_frequencies_aminoacid.tsv')

@blp.route("/data/frequencies/table/gene")
@api_key_required
@blp.arguments(GeneNameArgs, location="query")
@blp.response(200, description="Frequency table as a tab separated file")
@blp.doc(responses={200: {"content": {TSV_CONTENT_TYPE: FILE_BODY}}})
@blp.alt_response(404, description="No frequency table for the requested allele or gene")
def get_gene_frequencies_table(args):
    return frequency_table_attachment(args["gene_name"], "genomic",
                                      '_frequencies_genomic.tsv', full_gene=True)

@blp.route("/data/aminoacidfrequencies/table/gene")
@api_key_required
@blp.arguments(AminoAcidGeneNameArgs, location="query")
@blp.response(200, description="Frequency table as a tab separated file")
@blp.doc(responses={200: {"content": {TSV_CONTENT_TYPE: FILE_BODY}}})
@blp.alt_response(404, description="No frequency table for the requested allele or gene")
def get_aminoacid_gene_frequencies_table(args):
    return frequency_table_attachment(args["aa_gene_name"], "aminoacid",
                                      '_frequencies_aminoacid.tsv', full_gene=True)

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
@blp.alt_response(404, description="No sequences to align for the requested gene")
def get_aligned_sequences(args):
    aligned = align_sequences(args["gene_name"])
    if aligned is None:
        abort(404, message="No sequences to align for the requested gene.")

    return aligned

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
