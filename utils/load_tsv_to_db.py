import csv
import pandas as pd
from db import db
import datetime
from constants import ROOT_DIR
from models import ImmuneDiscoverDataModel



def load_tsv_to_db():
    data_in_dir = ROOT_DIR + "/data/in/"

    file_name = '1KGP_ImmuneDiscover_IGHV1-2.tsv'

    # df = pd.read_csv(data_in_dir + file_name,sep='\t')
    # df['load_tsv_metadata']=file_name + '_' + str(datetime.datetime.now())

    # # Insert to DB
    # df.to_sql('immunediscoverdata',
    #         con=db.engine,
    #         if_exists='append',
    #         index=True,
    #         index_label='id')


    with open(data_in_dir+file_name, encoding='utf-8', newline='') as tsv_file:
        tsvreader = csv.DictReader(tsv_file, delimiter='\t')

        immune_discover_data = [ImmuneDiscoverDataModel(**row) for row in tsvreader]

        db.session.add(immune_discover_data)
        db.session.commit()