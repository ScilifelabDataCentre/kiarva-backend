# import pandas as pd

# from constants import ROOT_DIR

# def main():
#     data_dir = ROOT_DIR + '/data/in/'
#     file_name = '1KGP_ImmuneDiscover_IGHV1-2.tsv'
#     df = pd.read_csv(data_dir + file_name, sep='\t')
#     # print(len(df))
#     # for column in df.columns:
#     #     unique_values = df[column].unique()
#     #     print(f"'{column}': {len(unique_values)}")
#     subset = ['case', 'db_name']
#     print(f"{subset}: '{len(df.drop_duplicates(subset=subset))}'")
#     pd.set_option('display.max_columns', None)
#     pd.set_option('display.max_rows', None)
#     print(df[df.duplicated(subset, keep=False)].head())

# main()
