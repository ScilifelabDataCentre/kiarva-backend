# Script for generating fasta files from db data

from models.immunediscoverdata import ImmuneDiscoverDataModel
from repositories.filters import is_deletion
from utils import dict_to_fasta_str

# Three different types of fasta files can be generated:
# "Standard" genomic, using the main sequence of a gene/allele
# "Flanking" genomic, using sequences with flanking segments included
# "Translated", using the translated amino acids
#
# All three match on the requested gene segment with startswith(autoescape=True) rather
# than like(gene_segment + '%'). The segment comes from a request, and LIKE reads '_' as a
# single-character wildcard, so a file_name of "_" returned most of the table in one file.
# '%' is already rejected by the schema, but '_' cannot be - it is a legitimate character
# in allele names.
def generate_fasta(gene_segment, type = "genomic"):
    if type not in ("genomic", "genomic_fl", "translated"):
        # Raise rather than fall through the branches below leaving distinct_sequences
        # unassigned, which failed with UnboundLocalError - an error about a local variable
        # rather than about the argument that was wrong.
        raise ValueError("unknown fasta type: " + repr(type))

    if type == "genomic":
        distinct_sequences = ImmuneDiscoverDataModel.query.with_entities(
                ImmuneDiscoverDataModel.db_name,
                ImmuneDiscoverDataModel.sequence,
                ).distinct().filter(ImmuneDiscoverDataModel.db_name.startswith(gene_segment, autoescape=True)).filter(~ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True)).filter(~is_deletion()).all()
    elif type == "genomic_fl":
        distinct_sequences = ImmuneDiscoverDataModel.query.with_entities(
                ImmuneDiscoverDataModel.db_name,
                ImmuneDiscoverDataModel.sequence,
                ).distinct().filter(ImmuneDiscoverDataModel.db_name.startswith(gene_segment, autoescape=True)).filter(ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True)).all()
    elif type == "translated":
                distinct_sequences = ImmuneDiscoverDataModel.query.with_entities(
                ImmuneDiscoverDataModel.db_name_AA,
                ImmuneDiscoverDataModel.sequence_AA,
                ).distinct().filter(ImmuneDiscoverDataModel.db_name_AA.startswith(gene_segment, autoescape=True)).all()

    seq_data = {}
    for row in distinct_sequences:
        seq_data[row[0]] = row[1]

    seq_data_sorted = dict(sorted(seq_data.items()))

    fasta_out = dict_to_fasta_str(seq_data_sorted)

    return fasta_out