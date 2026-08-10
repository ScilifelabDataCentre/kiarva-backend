# Pytest tests for all endpoints.
#
# The dataset these run against is built by tests/generate_mock_data.py, which
# documents why each row exists. Assertions here are written against exact output
# rather than shapes wherever the dataset makes that possible, so both files have to
# be read together when either changes.

import os

from urllib.parse import quote

api_key_header = {"X-api-key": os.getenv("API_KEY")}

# The allele carried by (almost) every case in the mock data, and the gene it belongs
# to. The frontend E2E fixtures are recorded for the same allele.
TEST_ALLELE = "IGHV1-8*01"
TEST_GENE = "IGHV1-8"

# db_name that belongs to two genes at once, stored with allele == db_name.
TWO_GENE_ALLELE = "IGHV3-30*02/IGHV3-30-5*02"

# Populations in the order the frequency endpoints report them, which is the order the
# research group asked for rather than alphabetical.
EXPECTED_POPULATIONS = [
    "ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI",
    "FIN", "GBR", "IBS", "TSI",
    "CDX", "CHB", "CHS", "JPT", "KHV",
    "BEB", "GIH", "ITU", "PJL", "STU",
    "CLM", "MXL", "PEL", "PUR",
    "ALL",
]
EXPECTED_SUPERPOPULATIONS = ["AFR", "EUR", "EAS", "SAS", "AMR", "ALL"]

POPULATION_REGIONS = {
    "ACB": "AFR", "ASW": "AFR", "ESN": "AFR", "GWD": "AFR", "LWK": "AFR",
    "MSL": "AFR", "YRI": "AFR",
    "FIN": "EUR", "GBR": "EUR", "IBS": "EUR", "TSI": "EUR",
    "CDX": "EAS", "CHB": "EAS", "CHS": "EAS", "JPT": "EAS", "KHV": "EAS",
    "BEB": "SAS", "GIH": "SAS", "ITU": "SAS", "PJL": "SAS", "STU": "SAS",
    "CLM": "AMR", "MXL": "AMR", "PEL": "AMR", "PUR": "AMR",
    "ALL": "ALL",
}

# The three IGHV1-8 sequences, aligned against each other by MAFFT, and their
# translations: no gap, a whole-codon gap, and a frameshift gap that cuts the
# translation short and pads it with X.
ALIGNED_IGHV1_8 = [
    {
        "allele": "IGHV1-8*01",
        "sequence_nt": "CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC",
        "sequence_aa": "LDSPLLALLCSGCDRLVDNA",
    },
    {
        "allele": "IGHV1-8*02",
        "sequence_nt": "CTGGATTCACCTTTAC---GTATGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC",
        "sequence_aa": "LDSPL-RMLCSGCDRLVDNA",
    },
    {
        "allele": "IGHV1-8*04",
        "sequence_nt": "CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGA-AGGCTCGTGGACAACGCC",
        "sequence_aa": "LDSPLLALLCSGCXXXXXXX",
    },
]

SEARCH_SEQUENCE = "THISISASEQUENCESEARCHTEST123"
# Sequence of the flanking (_F) row only, so searching for it must not return anything.
FLANKING_ONLY_SEQUENCE = "AGGTGCCCACTCC"


def get(client, path):
    return client.get(path, headers=api_key_header)


def url(path, **params):
    query = "&".join(f"{k}={quote(str(v), safe='')}" for k, v in params.items())
    return f"{path}?{query}" if query else path


def frequency_dict(data):
    return {entry["population"]: (entry["n"], entry["frequency"]) for entry in data}


def test_health(client):
    res = client.get("/health")
    assert res.status_code == 200
    assert res.json["status"] == "ok"


def test_checkapikey(client):
    res = get(client, "/checkapikey")
    assert res.status_code == 200
    assert res.get_data(as_text=True) == "Correct key!"


def test_populationregions(client):
    res = get(client, "/data/populationregions")
    assert res.status_code == 200
    data = res.get_json()

    # All 25 populations, their superpopulation, and the aggregated "ALL" entry, with
    # nothing extra. Order is left to the database, so compare as a mapping.
    assert {entry["population"]: entry["superpopulation"] for entry in data} == POPULATION_REGIONS
    assert len(data) == len(POPULATION_REGIONS)


def test_superpopfrequencies(client):
    res = get(client, url("/data/frequencies/superpopulations", allele_name=TEST_ALLELE))
    assert res.status_code == 200
    data = res.get_json()

    # AFR is missing IGHV1-8*01 in one of its 7 cases (MSL carries *DEL), EUR in one of
    # 5 (one of the two TSI cases carries *02) and EAS in one of 5 (CDX carries *04).
    assert data == [
        {"population": "AFR", "n": 6, "frequency": 0.85714},
        {"population": "EUR", "n": 4, "frequency": 0.8},
        {"population": "EAS", "n": 4, "frequency": 0.8},
        {"population": "SAS", "n": 5, "frequency": 1.0},
        {"population": "AMR", "n": 4, "frequency": 1.0},
        {"population": "ALL", "n": 23, "frequency": 0.88462},
    ]


def test_aminoacid_superpopfrequencies(client):
    res = get(
        client,
        url("/data/aminoacidfrequencies/superpopulations", aa_allele_name=TEST_ALLELE),
    )
    assert res.status_code == 200
    data = res.get_json()

    # Same allele, but counted through db_name_AA: *02 and *04 collapse into
    # IGHV1-8*01, so only the population without any IGHV1-8 at all (MSL, in AFR) is
    # missing it. A bug serving genomic data here would show up as the EUR/EAS values
    # from the test above.
    assert data == [
        {"population": "AFR", "n": 6, "frequency": 0.85714},
        {"population": "EUR", "n": 5, "frequency": 1.0},
        {"population": "EAS", "n": 5, "frequency": 1.0},
        {"population": "SAS", "n": 5, "frequency": 1.0},
        {"population": "AMR", "n": 4, "frequency": 1.0},
        {"population": "ALL", "n": 25, "frequency": 0.96154},
    ]


def test_populationfrequencies(client):
    res = get(client, url("/data/frequencies/populations", allele_name=TEST_ALLELE))
    assert res.status_code == 200
    data = res.get_json()

    assert [entry["population"] for entry in data] == EXPECTED_POPULATIONS

    frequencies = frequency_dict(data)
    # MSL carries only IGHV1-8*DEL and CDX only *04, so neither has the allele.
    assert frequencies["MSL"] == (0, 0.0)
    assert frequencies["CDX"] == (0, 0.0)
    # TSI is the only population with two cases, one of which carries *02 instead.
    assert frequencies["TSI"] == (1, 0.5)
    # 23 of the 26 cases in the dataset.
    assert frequencies["ALL"] == (23, 0.88462)
    # Every other population has a single case, carrying the allele.
    for population in EXPECTED_POPULATIONS:
        if population not in {"MSL", "CDX", "TSI", "ALL"}:
            assert frequencies[population] == (1, 1.0)


def test_aminoacid_populationfrequencies(client):
    res = get(
        client, url("/data/aminoacidfrequencies/populations", aa_allele_name=TEST_ALLELE)
    )
    assert res.status_code == 200
    data = res.get_json()

    assert [entry["population"] for entry in data] == EXPECTED_POPULATIONS

    frequencies = frequency_dict(data)
    # Only MSL, which has no IGHV1-8 sequence at all, is missing the amino acid allele.
    assert frequencies["MSL"] == (0, 0.0)
    # CDX (*04) and TSI (*02 plus *01) both collapse into IGHV1-8*01, so unlike the
    # genomic frequencies above they come out at 1.0.
    assert frequencies["CDX"] == (1, 1.0)
    assert frequencies["TSI"] == (2, 1.0)
    assert frequencies["ALL"] == (25, 0.96154)


def test_igsnperdata(client):
    res = get(client, url("/data/igsnperdata", allele_name=TEST_ALLELE))
    assert res.status_code == 200
    assert res.get_json() == {
        "igSNPer_score": 2.0,
        "igSNPer_SNPs": [
            "rs12345678(A:87,106539163)",
            "rs87654321(G:204,106539280)",
        ],
    }


def test_igsnperdata_without_rs_ids(client):
    # An allele with no rs ids still has one position-only entry, not an empty list.
    res = get(client, url("/data/igsnperdata", allele_name="IGHV3-30*01"))
    assert res.status_code == 200
    assert res.get_json() == {
        "igSNPer_score": 0.0,
        "igSNPer_SNPs": ["(:0,106791243)"],
    }


def test_igsnperdata_missing(client):
    # *DEL rows have no IgSNPer columns at all.
    res = get(client, url("/data/igsnperdata", allele_name="IGHV1-8*DEL"))
    assert res.status_code == 200
    assert res.get_json() == {"igSNPer_score": None, "igSNPer_SNPs": []}


def test_aminoacidalleles(client):
    # *02 is one of the alleles collapsed into the IGHV1-8*01 amino acid allele, so
    # asking for it has to return the collapsed name rather than the one asked for.
    res = get(client, url("/data/aminoacidalleles", aa_allele_name="IGHV1-8*02"))
    assert res.status_code == 200
    assert res.get_json() == {"allele": "IGHV1-8*02", "allele_aa": TEST_ALLELE}


def test_aminoacidlist(client):
    res = get(client, url("/data/aminoacidlist", aa_allele_name=TEST_ALLELE))
    assert res.status_code == 200
    assert res.get_json() == {
        "aa_allele_list": ["IGHV1-8*01", "IGHV1-8*02", "IGHV1-8*04"]
    }


def test_aminoacidlist_for_two_gene_allele(client):
    # The '*' and '/' in this name have to be escaped before it is matched against the
    # comma separated list it appears in.
    res = get(client, url("/data/aminoacidlist", aa_allele_name=TWO_GENE_ALLELE))
    assert res.status_code == 200
    assert res.get_json() == {
        "aa_allele_list": [TWO_GENE_ALLELE, "IGHV3-30*02_S9143"]
    }


def test_plotoptions_genes(client):
    # Gene lists per gene type. IGHV3-30 and IGHV3-30-5 both come out of the single
    # "IGHV3-30,IGHV3-30-5" gene value, and the renamed genes appear under the name
    # they were rewritten to ("IGHV1-69/1-69D" -> IGHV1-69, "IGHV3-23/3-23D" ->
    # IGHV3-23, "IGHV3-30+" -> IGHV3-30).
    assert get(client, url("/data/plotoptions", current_selection="IGHV")).get_json() == [
        "1-69", "1-8", "3-23", "3-30", "3-30-5",
    ]
    assert get(client, url("/data/plotoptions", current_selection="IGHD")).get_json() == ["1-1"]
    assert get(client, url("/data/plotoptions", current_selection="IGHJ")).get_json() == ["6"]


def test_plotoptions_alleles(client):
    # Allele lists per gene. The flanking IGHV1-8*01_F1 row repeats allele "01" and so
    # adds nothing here, while the homozygous deletion is offered as "DEL".
    res = get(client, url("/data/plotoptions", current_selection=TEST_GENE + "*"))
    assert res.status_code == 200
    assert res.get_json() == ["01", "02", "04", "DEL"]

    # For IGHV1-69 the loader drops *DEL rows instead of storing them, so the same
    # gene must not offer a DEL option.
    res = get(client, url("/data/plotoptions", current_selection="IGHV1-69*"))
    assert res.get_json() == ["04_S7754"]

    # Both genes of the two-gene db_name offer it as an allele, alongside their own.
    assert get(client, url("/data/plotoptions", current_selection="IGHV3-30*")).get_json() == [
        "01", TWO_GENE_ALLELE,
    ]
    assert get(client, url("/data/plotoptions", current_selection="IGHV3-30-5*")).get_json() == [
        "03_S1123", TWO_GENE_ALLELE,
    ]

    # D and J flanking rows, unlike V ones, do suffix the allele column.
    assert get(client, url("/data/plotoptions", current_selection="IGHD1-1*")).get_json() == [
        "01", "01_F1",
    ]


def test_db_name(client):
    # Plain gene + allele. The flanking row shares gene and allele with the row below
    # it, so this also asserts that the _F name is not what comes back.
    res = get(client, url("/data/db_name", selection=f"{TEST_GENE},01"))
    assert res.status_code == 200
    assert res.get_json() == {"db_name": TEST_ALLELE}

    # A gene whose alleles are only reachable through the two-gene db_name.
    res = get(client, url("/data/db_name", selection=f"IGHV3-30-5,{TWO_GENE_ALLELE}"))
    assert res.get_json() == {"db_name": TWO_GENE_ALLELE}

    # A gene the allele does not belong to.
    res = get(client, url("/data/db_name", selection=f"{TEST_GENE},03"))
    assert res.get_json() == {"db_name": "Not found"}


def test_alignedsequences(client):
    res = get(client, url("/data/sequences/alignedsequences", gene_name=TEST_GENE))
    assert res.status_code == 200
    # The *DEL row (no sequence) and the flanking row (a longer sequence that would
    # skew the whole alignment) are both left out.
    assert res.get_json() == ALIGNED_IGHV1_8


def test_sequencesearch(client):
    res = get(client, url("/data/sequences", sequence_str="ESEARCHTES"))
    assert res.status_code == 200
    assert res.get_json() == [
        {"allele": "IGHV3-23*01", "positions": [14], "sequence": SEARCH_SEQUENCE}
    ]


def test_sequencesearch_repeated_match(client):
    # Every position of the searched-for sequence is reported, not just the first.
    res = get(client, url("/data/sequences", sequence_str="SE"))
    assert res.status_code == 200
    assert res.get_json() == [
        {"allele": "IGHV3-23*01", "positions": [7, 15], "sequence": SEARCH_SEQUENCE}
    ]


def test_sequencesearch_skips_flanking_rows(client):
    # Sequences with flanking regions are a separate set of rows and are not searched,
    # so a match inside a flanking region is not a match.
    res = get(client, url("/data/sequences", sequence_str=FLANKING_ONLY_SEQUENCE))
    assert res.status_code == 200
    assert res.get_json() == [{"allele": "", "sequence": "", "position": []}]


def test_sequencesearch_no_match(client):
    res = get(client, url("/data/sequences", sequence_str="NOMATCHFORTHIS"))
    assert res.status_code == 200
    assert res.get_json() == [{"allele": "", "sequence": "", "position": []}]


def test_send_fasta_genomic(client):
    res = get(client, url("/fasta/genomic", file_name="IGHV1-8"))
    assert res.status_code == 200
    assert "text/plain" in res.content_type or "application/octet-stream" in res.content_type
    # Alleles in name order, with the *DEL and flanking rows left out.
    assert res.get_data(as_text=True) == (
        ">IGHV1-8*01\nCTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC\n"
        ">IGHV1-8*02\nCTGGATTCACCTTTACGTATGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC\n"
        ">IGHV1-8*04\nCTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGAAGGCTCGTGGACAACGCC\n"
    )


def test_send_fasta_genomicwithflanking(client):
    # The frontend asks for a whole gene type at a time, as here.
    res = get(client, url("/fasta/genomic_fl", file_name="IGHV"))
    assert res.status_code == 200
    assert "text/plain" in res.content_type or "application/octet-stream" in res.content_type
    # Only the flanking rows, and their sequence includes the prefix and suffix.
    assert res.get_data(as_text=True) == (
        ">IGHV1-8*01_F1\nAGGTGCCCACTCC"
        "CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC"
        "CACAGTGTGAAA\n"
    )


def test_send_fasta_translated(client):
    res = get(client, url("/fasta/translated", file_name="IGHV1-8"))
    assert res.status_code == 200
    assert "text/plain" in res.content_type or "application/octet-stream" in res.content_type
    # One entry per amino acid allele: the three IGHV1-8 alleles share one.
    assert res.get_data(as_text=True) == ">IGHV1-8*01\nLDSPLLALLCSGCDRLVDNA\n"


def test_send_fasta_only_covers_requested_gene_type(client):
    # D genes are not translated, so they have no amino acid rows to return.
    assert get(client, url("/fasta/genomic", file_name="IGHD")).get_data(as_text=True) == (
        ">IGHD1-1*01\nGGTACAACTGGAACGAC\n"
    )
    assert get(client, url("/fasta/translated", file_name="IGHD")).get_data(as_text=True) == ""


def test_frequencies_table_allele(client):
    res = get(client, url("/data/frequencies/table/allele", allele_name=TEST_ALLELE))
    assert res.status_code == 200
    rows = res.get_data(as_text=True).split("\n")

    assert rows[0] == "allele\tpopulation\tsuperpopulation\tfrequency\tn"
    # Superpopulations first, then populations, with the aggregated "ALL" of each
    # dropped: 5 + 25 rows.
    assert len(rows) == 31
    assert rows[1] == f"{TEST_ALLELE}\tAFR\tAFR\t0.85714\t6"
    assert f"{TEST_ALLELE}\tTSI\tEUR\t0.5\t1" in rows
    assert not [row for row in rows if "\tALL\t" in row]


def test_aminoacidfrequencies_table_gene(client):
    # Requesting a whole gene returns one block of rows per amino acid allele of that
    # gene, each carrying the list of alleles collapsed into it.
    res = get(client, url("/data/aminoacidfrequencies/table/gene", aa_gene_name=TEST_GENE))
    assert res.status_code == 200
    rows = res.get_data(as_text=True).split("\n")

    assert rows[0] == "collapsed_translated_sequence\tallele\tpopulation\tsuperpopulation\tfrequency\tn"
    assert len(rows) == 31
    assert rows[1] == (
        "['IGHV1-8*01', 'IGHV1-8*02', 'IGHV1-8*04']"
        f"\t{TEST_ALLELE}\tAFR\tAFR\t0.85714\t6"
    )
