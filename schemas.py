# Marshmallow schemas for the request arguments and responses of the endpoints in
# resources/. flask_smorest uses them for three things at once:
#
# - rejecting malformed requests with a 422 before a view function runs, so no untrusted
#   value reaches a database query,
# - serialising responses to exactly the fields declared here and nothing else,
# - generating the OpenAPI document served at /swagger-ui, which previously described
#   none of the request or response shapes.

from marshmallow import EXCLUDE, Schema, fields, validate

# Allele, gene and plot-selection names in the source data are built from a narrow
# character set - letters, digits and * / , _ - - and the longest value in the dataset is
# 28 characters. Rejecting anything else keeps regex metacharacters and SQL LIKE
# wildcards out of the queries these names are used in.
NAME_PATTERN = r"\A[A-Za-z0-9*/,_-]{1,64}\Z"

# Sequences are nucleotide (ACGT) or amino acid single-letter codes, '*' for a stop codon,
# with the longest stored sequence at 330 characters. Digits are allowed because a search
# term is free text supplied by the user rather than a stored sequence.
SEQUENCE_PATTERN = r"\A[A-Za-z0-9*-]{1,400}\Z"

def name_field(description):
    return fields.Str(
        required = True,
        validate = validate.Regexp(NAME_PATTERN),
        metadata = {"description": description},
    )

class ArgsSchema(Schema):
    """Base for query argument schemas.

    Unknown query parameters are ignored rather than rejected. Marshmallow 4 defaults to
    raising on them, which would turn a previously working request that carries a stray
    parameter into a 422.
    """
    class Meta:
        unknown = EXCLUDE

# ---------------------------------------------------------------- request arguments

class AlleleNameArgs(ArgsSchema):
    allele_name = name_field("Genomic allele name (db_name), e.g. IGHV1-2*02")

class AminoAcidAlleleNameArgs(ArgsSchema):
    aa_allele_name = name_field("Amino acid allele name (db_name_AA)")

class GeneNameArgs(ArgsSchema):
    gene_name = name_field("Gene name, e.g. IGHV1-2")

class AminoAcidGeneNameArgs(ArgsSchema):
    aa_gene_name = name_field("Gene name to report amino acid frequencies for")

class FastaFileNameArgs(ArgsSchema):
    file_name = name_field("Gene or gene segment to build the FASTA file from")

class PlotOptionArgs(ArgsSchema):
    current_selection = name_field(
        "Selection made so far. A trailing '*' asks for allele names, otherwise gene names"
    )

class SelectionArgs(ArgsSchema):
    selection = name_field("Gene and allele as displayed in the plot options, comma separated")

class SequenceSearchArgs(ArgsSchema):
    sequence_str = fields.Str(
        required = True,
        validate = validate.Regexp(SEQUENCE_PATTERN),
        metadata = {"description": "Nucleotide or amino acid sequence to search for"},
    )

# ---------------------------------------------------------------- responses

class HealthSchema(Schema):
    status = fields.Str(required = True)

class DbNameSchema(Schema):
    db_name = fields.Str(required = True, metadata = {"description": "'Not found' if the selection does not resolve"})

class FrequencySchema(Schema):
    population = fields.Str(required = True, metadata = {"description": "Population, superpopulation, or 'ALL'"})
    n = fields.Int(required = True, metadata = {"description": "Cases carrying the allele"})
    frequency = fields.Float(required = True, metadata = {"description": "n divided by the cases in the population"})

class PopulationRegionSchema(Schema):
    population = fields.Str(required = True)
    superpopulation = fields.Str(required = True)

class IgSNPerSchema(Schema):
    # Null for rows with no IgSNPer columns at all, such as the homozygous deletions.
    igSNPer_score = fields.Float(required = True, allow_none = True)
    igSNPer_SNPs = fields.List(fields.Str(), required = True)

class AminoAcidTopAlleleSchema(Schema):
    # Both omitted when the allele is not found, which is why neither is required.
    allele = fields.Str()
    allele_aa = fields.Str()

class AminoAcidAlleleListSchema(Schema):
    aa_allele_list = fields.List(fields.Str(), required = True, allow_none = True)

class AlignedSequenceSchema(Schema):
    allele = fields.Str(required = True)
    sequence_nt = fields.Str(required = True, metadata = {"description": "Aligned nucleotide sequence, '-' for gaps"})
    sequence_aa = fields.Str(required = True, metadata = {"description": "Translation, 'X' where a frameshift gap cuts it short"})

class SequenceSearchSchema(Schema):
    allele = fields.Str(required = True)
    sequence = fields.Str(required = True)
    positions = fields.List(fields.Int(), required = True,
                            metadata = {"description": "Every offset the search term occurs at"})
