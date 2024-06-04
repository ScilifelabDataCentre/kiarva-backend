import csv
from db import db
from constants import ROOT_DIR
from models import ImmuneDiscoverDataModel



def load_tsv_to_db():
    data_in_dir = ROOT_DIR + "/data/in/"

    file_name = '1KGP_ImmuneDiscover_IGHV1-2.tsv'

    # Open file 
    with open(data_in_dir+file_name, encoding='utf-8', newline='') as tsv_file:
        tsvreader = csv.DictReader(tsv_file, delimiter='\t')

        immune_discover_data = [ImmuneDiscoverDataModel(**row) for row in tsvreader]

        for row in immune_discover_data:
            db.session.add(row)
        db.session.commit()