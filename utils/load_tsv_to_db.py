import csv
from datetime import datetime
from db import db
from constants import ROOT_DIR
from models import ImmuneDiscoverDataModel



def load_tsv_to_db(file_name):
    data_in_dir = ROOT_DIR + "/data/in/"

    # Open file 
    with open(data_in_dir+file_name, encoding='utf-8', newline='') as tsv_file:
        tsvreader = csv.DictReader(tsv_file, delimiter='\t')

        immune_discover_data = []
        for row in tsvreader:
            row["loaded_from_tsv"] = file_name
            row["loaded_at"] = str(datetime.now())
            immune_discover_data.append(ImmuneDiscoverDataModel(**row))

        for row in immune_discover_data:
            db.session.add(row)
        db.session.commit()