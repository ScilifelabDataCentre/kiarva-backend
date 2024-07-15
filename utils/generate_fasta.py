from models.immunediscoverdata import ImmuneDiscoverDataModel

def generate_fasta(gene_segment, genomic = False):
    prefix = ''
    suffix = ''
    match gene_segment:
        case 'IGHV':
            prefix = 'prefix'
            suffix = 'heptamer'
        case 'IGHD':
            prefix = 'pre_heptamer'
            suffix = 'post_heptamer'
        case 'IGHJ':
            prefix = 'heptamer'
            suffix = 'suffix'

    if genomic:
        distinct_sequences = ImmuneDiscoverDataModel.query.with_entities(
                ImmuneDiscoverDataModel.db_name,
                ImmuneDiscoverDataModel.sequence,
                getattr(ImmuneDiscoverDataModel, prefix),
                getattr(ImmuneDiscoverDataModel, suffix),
                ImmuneDiscoverDataModel.flank_index
                ).distinct().filter(ImmuneDiscoverDataModel.db_name.like(gene_segment+'%')).all()
    else:
        distinct_sequences = ImmuneDiscoverDataModel.query.with_entities(
                ImmuneDiscoverDataModel.db_name,
                ImmuneDiscoverDataModel.sequence,
                ).distinct().filter(ImmuneDiscoverDataModel.db_name.like(gene_segment+'%')).all()
    

    fasta_out = ''
    for row in distinct_sequences:
        if genomic:
            fasta_out += '>' + row[0] + '_F' + str(row[4]) + '\n' + row[2] + row[1] + row[3] + '\n'
        else:
            fasta_out += '>' + row[0] + '\n' +row[1] + '\n'
    return fasta_out