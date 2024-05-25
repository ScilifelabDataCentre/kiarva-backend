import os
import pandas as pd

def df_to_fasta(df, desc_col, seq_col):
    df = df[[desc_col, seq_col]].drop_duplicates()
    fasta_out = ''
    for index, row in df.iterrows():
        fasta_out += '>' + row[desc_col] + '\n' +row[seq_col] + '\n'
    return fasta_out


def main():
    file_name = '1KGP_ImmuneDiscover_IGHV1-2.tsv'
    file_dir = os.path.dirname(os.path.realpath('./data/'+file_name))
    df = pd.read_csv(file_dir + '/' + file_name,sep='\t')
    print(df.dtypes)
    pd.set_option('display.max_colwidth', None)
    # print(df[['db_name', 'sequence']].drop_duplicates())
    print(df_to_fasta(df, 'db_name', 'sequence'))


main()