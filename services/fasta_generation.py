# Script for generating fasta files from db data

from models.immunediscoverdata import ImmuneDiscoverDataModel
from utils import dict_to_fasta_str

# Three different types of fasta files can be generated: 
# "Standard" genomic, using the main sequence of a gene/allele
# "Flanking" genomic, using sequences with flanking segments included
# "Translated", using the translated amino acids
def generate_fasta(gene_segment, type = "genomic"):
    if type == "genomic":
        distinct_sequences = ImmuneDiscoverDataModel.query.with_entities(
                ImmuneDiscoverDataModel.db_name,
                ImmuneDiscoverDataModel.sequence,
                ).distinct().filter(ImmuneDiscoverDataModel.db_name.like(gene_segment+'%')).filter(ImmuneDiscoverDataModel.db_name.notlike('%_F%')).filter(ImmuneDiscoverDataModel.db_name.notlike('%*DEL')).all()
    elif type == "genomic_fl":
        distinct_sequences = ImmuneDiscoverDataModel.query.with_entities(
                ImmuneDiscoverDataModel.db_name,
                ImmuneDiscoverDataModel.sequence,
                ).distinct().filter(ImmuneDiscoverDataModel.db_name.like(gene_segment+'%_F%')).all()
    elif type == "translated":
                distinct_sequences = ImmuneDiscoverDataModel.query.with_entities(
                ImmuneDiscoverDataModel.db_name_AA,
                ImmuneDiscoverDataModel.sequence_AA,
                ).distinct().filter(ImmuneDiscoverDataModel.db_name_AA.like(gene_segment+'%')).all()

    seq_data = {}
    for row in distinct_sequences:
        seq_data[row[0]] = row[1]

    seq_data_sorted = dict(sorted(seq_data.items()))

    fasta_out = dict_to_fasta_str(seq_data_sorted)

    return fasta_out