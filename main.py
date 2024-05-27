import pandas as pd
from constants import ROOT_DIR

def df_to_fasta(df, db_name_col, seq_col, seq_type):
    df = df[df[db_name_col].str.startswith(seq_type)]
    df = df[[db_name_col, seq_col]].drop_duplicates()
    fasta_out = ''
    for index, row in df.iterrows():
        fasta_out += '>' + row[db_name_col] + '\n' +row[seq_col] + '\n'
    return fasta_out

def write_fastas(df, output_dir):
    seq_types = ['IGHV', 'IGHD', 'IGHJ']
    for seq_type in seq_types: 
        f = open(output_dir + seq_type + '.fasta', 'w')
        f.write(df_to_fasta(df, 'db_name', 'sequence', seq_type))
        f.close()


def main():
    data_in_dir = ROOT_DIR + "/data/in/"
    data_out_dir = ROOT_DIR + "/data/out/"

    file_name = '1KGP_ImmuneDiscover_IGHV1-2.tsv'

    df = pd.read_csv(data_in_dir + file_name,sep='\t')

    print(df.dtypes)
    print(df[df.columns[:5]].head())
    print()
    print(df[df.columns[5:10]].head())
    print()
    print(df[df.columns[10:15]].head())
    # write_fastas(df, data_out_dir)


main()