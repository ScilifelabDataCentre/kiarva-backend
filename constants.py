# Constants that are reused in the rest of the app

import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# Dicts holding the pre-calculated frequency plot data, filled by
# loaders/plot_loader.py at startup in every mode. Module-level, and imported by name
# elsewhere, so plot_loader fills them with update() rather than rebinding - and the pytest
# fixture clears them between tests, since they outlive the app that filled them.
allele_superpopulation_frequencies = {}
allele_population_frequencies = {}
aminoacid_allele_superpopulation_frequencies = {}
aminoacid_allele_population_frequencies = {}
