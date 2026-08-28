# Research data for this project is provided by the Hedestam group in the format
# of .tsv files. The scripts in this file relates to loading the tsv data from file, 
# validating it, some minor fixes to column data types, and storing in project db.

import csv
from datetime import datetime
import os
from db import db
from constants import ROOT_DIR
from models import ImmuneDiscoverDataModel
from loaders.validation import SourceDataError, validate_loaded_data
from sqlalchemy.exc import IntegrityError
from flask import current_app
import pyzipper

validation_schema = {
    "cohort": str,
    "case": str,
    "db_name": str,
    "gene": str,
    "allele": str,
    "sequence": (str, type(None)),
    "flank_index": int,
    "IgSNPer_uncommon": (str, type(None)),
    "IgSNPer_SNPs": (str, type(None)),
    "db_name_AA": (str, type(None)),
    "db_name_AA_list": (str, type(None)),
    "sequence_AA": (str, type(None)),
    "superpopulation": str,
    "population": str,
    "loaded_from_tsv": str,
    "loaded_at": str,
}

# The following two dictionaries + 1 list contain temporary changes to data, which will later be changed 
# at a source data level and removed from here. For now they are used as a proof of concept that these changes
# will fix issues we currently have with how plot options are displayed.
gene_name_change_dict = {
                            'IGHV1-69/1-69D': 'IGHV1-69',
                            'IGHV3-23/3-23D': 'IGHV3-23',
                            'IGHV3-30+': 'IGHV3-30',
                        }

db_name_gene_name_change_dict = {
                            'IGHV3-30*02/IGHV3-30-5*02': 'IGHV3-30,IGHV3-30-5',
                            'IGHV3-30*18/IGHV3-30-5*01': 'IGHV3-30,IGHV3-30-5',
                            'IGHV3-30*04/IGHV3-30-3*03': 'IGHV3-30,IGHV3-30-3',
                            'IGHV3-30-5*03': 'IGHV3-30-5',
                            'IGHV3-30-5*03_S1123': 'IGHV3-30-5',
                            'IGHV3-30-5*03_S5145': 'IGHV3-30-5',
                            'IGHV3-30-5*03_S8990': 'IGHV3-30-5',
                                }

alleles_to_change = ['IGHV3-30*02/IGHV3-30-5*02']

# genes for which *DEL should not be shown in the data. Not sure if these entries will be removed from source
# data at some point, or if they want to keep them but not show them to the users yet.
genes_without_del = ['IGHV1-69/1-69D']

def validate_row(row_num, data):
    for k, expected in validation_schema.items():
        v = data.get(k)

        if isinstance(expected, tuple):  # allow multiple types
            if not isinstance(v, expected):
                print(f"[Row {row_num}] {k} has invalid type: {type(v)} → {repr(v)}")
        else:
            if not isinstance(v, expected):
                print(f"[Row {row_num}] {k} has invalid type: {type(v)} → {repr(v)}")

        # Check for weird values
        if isinstance(v, str):
            val_clean = v.strip().lower()
            if val_clean in {"none", "nan", "null"}:
                print(f"[Row {row_num}] {k} contains suspicious string value: {repr(v)}")

        # Optional: log empty strings
        if v == "":
            print(f"[Row {row_num}] {k} is an empty string")

        # Optional: log None where not expected
        if v is None and not isinstance(expected, tuple):
            print(f"[Row {row_num}] {k} is None, expected {expected.__name__}")

# tsv_files are compressed in repo, unpack them before loading.
#
# Every archive found is extracted on a best-effort basis, and an archive that will not open
# is skipped rather than raised on. That is deliberate, not an oversight:
#
# Prepub data - results the research group has not published yet - is password protected and
# unlocked in k8s with a sealed secret. It is served by a second deployment sitting behind
# basic auth, so the researchers can use the site's own functionality on their unpublished
# data without the app needing logins or sessions. The public deployment has no password, so
# the prepub archive fails to open there and its data is simply left out, which is the
# intended behaviour. The same holds locally: anyone can set the app up and run it without
# access to the password, and gets exactly the public deployment's dataset.
#
# If no archive at all extracts, load_tsv_to_db() below raises on finding no .tsv files.
def unpack_compressed_tsv_files(data_dir):
    paths_to_compressed_files = [data_dir + 'compressed/' + file for file in os.listdir(data_dir + 'compressed/') if file.endswith('.zip')]
    zip_pass = os.getenv("ZIP_PASSWORD") or None
    
    for path in paths_to_compressed_files:
        # using pyzipper instead of the native zipfile lib, due to the latter not supporting
        # handling AES256-encrypted files
        with pyzipper.AESZipFile(path, 'r') as filezip:
            # if zip_pass was set using env variable, set it as password for current zipfile.
            # Must be converted to byte-type using .encode('utf-8')
            if zip_pass:
                filezip.setpassword(zip_pass.encode('utf-8'))
            try:
                filezip.extractall(data_dir + 'in/')
            except RuntimeError:
                print("Password incorrect/not found for " + path)

# Loading one file is its own function so the whole of it sits inside a single try in the
# caller, which is where the one commit for the file lives. Nothing here commits, so nothing
# it writes is visible until that commit succeeds - a file that fails partway is absent from
# loaded_from_tsv and the next start reads it again instead of skipping it as already done.
#
# The memory bound comes from the batch size, which caps how many objects are built at once.
# bulk_save_objects() emits its own INSERTs into the open transaction, so there is nothing
# pending for a flush to write - an explicit flush here measured identically, 78 MB either
# way, and the suite passed with it deleted.
def load_one_tsv(data_in_dir, file):
    with open(data_in_dir+file, encoding='utf-8', newline='') as tsv_file:
        tsvreader = csv.DictReader(tsv_file, delimiter='\t')

        # load data in batches, to avoid spiking memory usage
        batch_size = 1000
        batch = []
        for row_index, row in enumerate(tsvreader):
            # omit *DEL data for specified genes by skipping an iteration
            if row["gene"].strip() in genes_without_del and row["allele"] == "DEL":
                continue

            data_to_add = {}
            # population names are part of "case" column string,
            # split them out to add them to their own columns
            population_data_split = row["case"].strip().split("_")

            db_name = row["db_name"].strip().split("_F")[0]
            gene = row["gene"].strip()
            # add row content to dict, stripping whitespaces,
            # setting null values to fitting values if necessary
            # (if col is nullable, set to None)
            data_to_add["cohort"] = row["cohort"].strip()
            data_to_add["case"] = row["case"].strip()
            data_to_add["db_name"] = row["db_name"].strip()

            if db_name in db_name_gene_name_change_dict.keys():
                data_to_add["gene"] = db_name_gene_name_change_dict[db_name]
            elif gene in gene_name_change_dict.keys():
                data_to_add["gene"] = gene_name_change_dict[gene]
            else:
                data_to_add["gene"] = row["gene"].strip()

            if db_name in alleles_to_change:
                data_to_add["allele"] = db_name
            else:
                data_to_add["allele"] = row["allele"].strip()

            if not row["sequence"]:
                data_to_add["sequence"] = None
            else:
                data_to_add["sequence"] = row["sequence"].strip()
            if not row["flank_index"]:
                # non-nullable but contains positive integer null values.
                # set null values to -1
                data_to_add["flank_index"] = -1
            else:
                data_to_add["flank_index"] = int(float(row["flank_index"].strip()))
            if not row["IgSNPer_uncommon"]:
                data_to_add["IgSNPer_uncommon"] = None
            else:
                data_to_add["IgSNPer_uncommon"] = row["IgSNPer_uncommon"].strip()
            if not row["IgSNPer_SNPs"]:
                data_to_add["IgSNPer_SNPs"] = None
            else:
                data_to_add["IgSNPer_SNPs"] = row["IgSNPer_SNPs"].strip()
            if not row["db_name_AA"]:
                data_to_add["db_name_AA"] = None
            else:
                data_to_add["db_name_AA"] = row["db_name_AA"].strip()
            if not row["db_name_AA_list"] or row["db_name_AA_list"] == "nan":
                data_to_add["db_name_AA_list"] = None
            else:
                data_to_add["db_name_AA_list"] = row["db_name_AA_list"].strip()
            if not row["sequence_AA"]:
                data_to_add["sequence_AA"] = None
            else:
                data_to_add["sequence_AA"] = row["sequence_AA"].strip()
            data_to_add["superpopulation"] = population_data_split[2].strip()
            data_to_add["population"] = population_data_split[1].strip()
            data_to_add["loaded_from_tsv"] = file.strip()
            data_to_add["loaded_at"] = str(datetime.now())

            validate_row(row_index, data_to_add)
            batch.append(ImmuneDiscoverDataModel(**data_to_add))

            if len(batch) >= batch_size:
                db.session.bulk_save_objects(batch)
                batch = []

        if batch:
            db.session.bulk_save_objects(batch)

# if db is set up - load all tsv files which have not yet been loaded
def load_tsv_to_db():
    data_dir = current_app.config.get("DATA_DIR")
    unpack_compressed_tsv_files(data_dir)
        
    data_in_dir = data_dir + 'in/'
    # get all tsv files from the current data/in dir
    tsv_files = [file for file in os.listdir(data_in_dir) if file.endswith('.tsv')]
    if not tsv_files:
        raise SourceDataError(
            "No .tsv files found in " + data_in_dir + ", so there is no data to serve. "
            "Check that data/compressed/ holds the archives and that ZIP_PASSWORD is set "
            "for the password protected ones.")

    print("Loading files: " + str(tsv_files))

    # check db for which files have already been loaded ones, to not process
    # them again
    files_in_db = ImmuneDiscoverDataModel.query.with_entities(ImmuneDiscoverDataModel.loaded_from_tsv).distinct().all()
    loaded_files = [loaded_file[0] for loaded_file in files_in_db]
            
    for file in tsv_files:
        if file not in loaded_files:
            # One transaction per file. Committing each batch as it filled meant an
            # IntegrityError partway through left the earlier batches committed, which
            # recorded the file in loaded_from_tsv - so every later startup skipped it and
            # the data stayed silently truncated.
            #
            # Rolled back on anything, not only IntegrityError: a missing column, a
            # non-numeric flank_index or a case value with too few parts raises from inside
            # load_one_tsv, and leaving those rows in an open transaction hands the next
            # user of the session a failed one.
            try:
                load_one_tsv(data_in_dir, file)
                db.session.commit()
            except IntegrityError as e:
                db.session.rollback()
                raise SourceDataError(
                    file + " duplicates rows already in the database and was not loaded. "
                    "Every row must be a unique combination of case, db_name, gene and "
                    "flank_index." ) from e
            except Exception:
                db.session.rollback()
                raise

    # Validated here rather than by the caller so that no code path can load data without
    # checking it - the pytest fixtures call this function directly rather than going
    # through create_app(), and a broken TSV should fail the tests too.
    validate_loaded_data()