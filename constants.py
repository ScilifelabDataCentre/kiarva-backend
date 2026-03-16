# Constants that are reused in the rest of the app

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# Dicts that will contain pre-loaded frequency plot data, if not running
# 'TESTING' or 'DEBUG'
allele_superpopulation_frequencies = {}
allele_population_frequencies = {}
aminoacid_allele_superpopulation_frequencies = {}
aminoacid_allele_population_frequencies = {}
