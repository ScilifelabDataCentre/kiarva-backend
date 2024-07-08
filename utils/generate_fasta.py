from models.immunediscoverdata import ImmuneDiscoverDataModel

def table_to_fasta(gene_segment):
    distinct_sequences = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.db_name,
            ImmuneDiscoverDataModel.sequence,
            ).distinct().filter(ImmuneDiscoverDataModel.db_name.like(gene_segment+'%')).all()
    fasta_out = ''
    for row in distinct_sequences:
        fasta_out += '>' + row[0] + '\n' +row[1] + '\n'
    return fasta_out