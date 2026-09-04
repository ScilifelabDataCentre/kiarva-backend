# Pytest tests for all endpoints.
#
# The dataset these run against is built by tests/generate_mock_data.py, which
# documents why each row exists. Assertions here are written against exact output
# rather than shapes wherever the dataset makes that possible, so both files have to
# be read together when either changes.

import copy
import os

import pytest

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
    # MSL carries only IGHV1-8*DEL and CDX only *04, so neither has the allele. Asserting
    # the full population list above and these zeroes here is what pins the zero fill:
    # GROUP BY returns no row for a population with no carriers, so they have to be filled
    # back in rather than left out.
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
    # These exact figures are what pin the amino acid denominator. case_MSL_AFR carries only
    # IGHV1-8*DEL, so it has no amino acid row at all - and it still counts as an individual,
    # which is why ALL below is 25 of 26 rather than 25 of 25. Restricting the denominator to
    # cases that have amino acid data is a silent mistake on the real dataset, where every
    # case happens to have some; here it turns these numbers into 1.0.
    assert frequencies["MSL"] == (0, 0.0)
    # CDX (*04) and TSI (*02 plus *01) both collapse into IGHV1-8*01, so unlike the
    # genomic frequencies above they come out at 1.0.
    assert frequencies["CDX"] == (1, 1.0)
    assert frequencies["TSI"] == (2, 1.0)
    assert frequencies["ALL"] == (25, 0.96154)


FREQUENCY_ENDPOINTS = (
    ("/data/frequencies/superpopulations", "allele_name"),
    ("/data/frequencies/populations", "allele_name"),
    ("/data/aminoacidfrequencies/superpopulations", "aa_allele_name"),
    ("/data/aminoacidfrequencies/populations", "aa_allele_name"),
)


def test_frequencies_missing_allele_name(client):
    # A request with no allele name is a client error, not a plot of zeroes. In prod the
    # amino acid dicts used to hold a None key - db_name_AA is null on every row that is
    # not plottable, and SELECT DISTINCT returns that null as an allele - so the missing
    # parameter looked up that key and returned an all-zero plot with a 200. The schemas
    # make it a 422, the code webargs uses for a request that fails validation.
    for endpoint, _ in FREQUENCY_ENDPOINTS:
        assert get(client, endpoint).status_code == 422


def test_frequencies_malformed_allele_name(client):
    # Allele names use a narrow character set, so the schema rejects anything outside it
    # before the view runs. The payloads below are the shapes a scanner worries about:
    # markup, a quote break-out, an injection attempt and an over-long value.
    for endpoint, param in FREQUENCY_ENDPOINTS:
        for value in ("<script>alert(1)</script>", "IGHV1-8*01'\"", "IGHV1-8*01 OR 1=1", "A" * 65):
            res = get(client, url(endpoint, **{param: value}))
            assert res.status_code == 422, f"{endpoint} accepted {value!r}"

            # Whatever the response says, it must not contain the rejected value: that
            # reflection is the actual vulnerability, independent of content type.
            assert value not in res.get_data(as_text=True)

    # A legitimate name with every special character the dataset uses still gets through.
    res = get(client, url("/data/frequencies/superpopulations", allele_name=TWO_GENE_ALLELE))
    assert res.status_code == 200


def test_frequency_table_does_not_mutate_the_precalculated_dicts(app):
    # Building the tsv writes 'allele' and 'superpopulation' into each entry. Done in place,
    # those keys leaked into every later response from the plot endpoints for that allele.
    # No flag flipping needed any more: the dictionaries are the only source in every mode.
    from constants import allele_superpopulation_frequencies
    from services.frequencies import create_frequencies_table

    cached_before = copy.deepcopy(allele_superpopulation_frequencies[TEST_ALLELE])

    table = create_frequencies_table(TEST_ALLELE, "genomic")

    assert allele_superpopulation_frequencies[TEST_ALLELE] == cached_before
    assert table.splitlines()[0].split("\t") == [
        "allele", "population", "superpopulation", "frequency", "n",
    ]


def test_frequency_tables_are_offered_only_for_plottable_alleles(client):
    # The download option only appears in the frontend once an allele has been selected to
    # plot, so the endpoint covers the same set. Everything reported comes from the
    # pre-calculated dictionaries, so a non-plottable name has no table rather than a table
    # of zeroes naming an allele that does not exist.
    for endpoint, param in (
        ("/data/frequencies/table/allele", "allele_name"),
        ("/data/aminoacidfrequencies/table/allele", "aa_allele_name"),
    ):
        assert get(client, url(endpoint, **{param: TEST_ALLELE})).status_code == 200
        for absent in ("IGHV9-99*99", TEST_ALLELE + "_F1", "IGHD1-1*01"):
            res = get(client, url(endpoint, **{param: absent}))
            assert res.status_code == 404, f"{endpoint} served a table for {absent}"
            assert absent not in res.get_data(as_text=True)

    for endpoint, param in (
        ("/data/frequencies/table/gene", "gene_name"),
        ("/data/aminoacidfrequencies/table/gene", "aa_gene_name"),
    ):
        assert get(client, url(endpoint, **{param: TEST_GENE})).status_code == 200
        for absent in ("IGHV9-99", "IGHD1-1"):
            res = get(client, url(endpoint, **{param: absent}))
            assert res.status_code == 404, f"{endpoint} served a table for {absent}"


def test_underscore_is_not_a_like_wildcard(client):
    # LIKE reads '_' as a single-character wildcard. These endpoints matched on
    # like(value + '%') with the raw parameter, so a value of "_" matched everything:
    # /fasta/genomic returned 11 of the 14 sequences in the mock data instead of none.
    # '_' cannot be rejected by the schema the way '%' is - it is legitimate in allele
    # names like IGHV1-58*02_S3393 - so the query has to escape it.
    res = get(client, url("/data/plotoptions", current_selection="_"))
    assert res.status_code == 200
    assert res.get_json() == []

    for endpoint in ("/fasta/genomic", "/fasta/genomic_fl", "/fasta/translated"):
        res = get(client, url(endpoint, file_name="_"))
        assert res.status_code == 200
        assert res.get_data(as_text=True) == "", f"{endpoint} treated '_' as a wildcard"

    # A real prefix still matches, so the escaping has not broken ordinary lookups.
    res = get(client, url("/fasta/genomic", file_name=TEST_GENE))
    assert res.get_data(as_text=True).count(">") == 3


def test_igsnperdata_score_without_snps(client):
    # A score with no SNPs is the majority shape in the real data - 209,867 rows - and not a
    # contradiction: the score counts uncommon SNPs, so 0.0 means there were none to list and
    # the column is empty, which the loader stores as NULL. Taking len() of that null raised
    # TypeError and answered 500 for 29 of the 732 plottable alleles. IGHV1-69*04_S7754 is
    # the mock allele carrying this shape; asserting on an allele that has SNPs would pass
    # either way and guard nothing.
    res = get(client, url("/data/igsnperdata", allele_name="IGHV1-69*04_S7754"))
    assert res.status_code == 200
    assert res.get_json() == {"igSNPer_score": 0.0, "igSNPer_SNPs": []}


def test_openapi_documents_the_404s_and_the_api_key(client):
    # The spec is what /swagger-ui renders. Without the security scheme the page offers no
    # way to send X-api-key, so every "Try it out" comes back 400.
    spec = get(client, "/openapi.json").get_json()

    assert spec["components"]["securitySchemes"]["ApiKeyAuth"]["name"] == "X-api-key"
    assert spec["security"] == [{"ApiKeyAuth": []}]

    for path in (
        "/data/frequencies/superpopulations",
        "/data/frequencies/populations",
        "/data/aminoacidfrequencies/superpopulations",
        "/data/aminoacidfrequencies/populations",
        "/data/igsnperdata",
    ):
        assert "404" in spec["paths"][path]["get"]["responses"], f"{path} does not document its 404"

    # The downloads have no schema, so blp.response(content_type=...) documented nothing for
    # them; the body is declared through blp.doc instead.
    fasta = spec["paths"]["/fasta/genomic"]["get"]["responses"]["200"]
    assert "text/x-fasta" in fasta["content"]
    table = spec["paths"]["/data/frequencies/table/allele"]["get"]["responses"]["200"]
    assert "text/tab-separated-values" in table["content"]


def test_igsnperdata_unknown_allele(client):
    # An allele that is not in the data at all yields no rows. Indexing the empty result
    # raised IndexError and surfaced as a 500; it is a 404. This is distinct from an allele
    # that exists with no IgSNPer columns, which is still a 200 - see the test below.
    res = get(client, url("/data/igsnperdata", allele_name="IGHV9-99*99"))
    assert res.status_code == 404


def test_frequencies_unknown_allele_is_404_in_every_mode(client):
    # Plottability used to be read off the dictionaries pre-calculated at startup, which are
    # only populated in prod, so this request 404d there and returned a 200 all-zero plot
    # under pytest - leaving the 404 impossible to cover. Both modes now agree.
    for endpoint, param in FREQUENCY_ENDPOINTS:
        assert get(client, url(endpoint, **{param: TEST_ALLELE})).status_code == 200
        for absent in ("IGHV9-99*99", TEST_ALLELE + "_F1", "IGHD1-1*01"):
            res = get(client, url(endpoint, **{param: absent}))
            assert res.status_code == 404, f"{endpoint} returned {res.status_code} for {absent}"


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
    # Only IGHV and TRGV are plotted. The research group has specified that IGHD and
    # IGHJ are not to be offered as plot or MSA selections, so they yield no options
    # even though their rows are loaded for the FASTA downloads.
    assert get(client, url("/data/plotoptions", current_selection="IGHD")).get_json() == []
    assert get(client, url("/data/plotoptions", current_selection="IGHJ")).get_json() == []

    # A '*' with no gene in front of it selects nothing. plot_options_regex("") builds
    # ',{0,1}(,|$)', which matches every gene, so this answered with every plottable allele
    # in the data. '*' is legal in an allele name, so NAME_PATTERN cannot reject it.
    assert get(client, url("/data/plotoptions", current_selection="*")).get_json() == []


def test_plotoptions_alleles(client):
    # Allele lists per gene. The flanking IGHV1-8*01_F1 row is excluded outright, and
    # for V genes it repeats allele "01" anyway so it would add nothing here. The
    # homozygous deletion is offered as "DEL".
    res = get(client, url("/data/plotoptions", current_selection=TEST_GENE + "*"))
    assert res.status_code == 200
    assert res.get_json() == ["01", "02", "04", "DEL"]

    # For IGHV1-69 the loader drops *DEL rows instead of storing them, so the same
    # gene must not offer a DEL option.
    res = get(client, url("/data/plotoptions", current_selection="IGHV1-69*"))
    assert res.get_json() == ["04_S7754"]

    # Both genes of the two-gene db_name offer it as an allele, alongside their own.
    assert get(client, url("/data/plotoptions", current_selection="IGHV3-30*")).get_json() == [
        "01", "02_S9143", TWO_GENE_ALLELE,
    ]
    assert get(client, url("/data/plotoptions", current_selection="IGHV3-30-5*")).get_json() == [
        "03_S1123", TWO_GENE_ALLELE,
    ]

    # D and J flanking rows, unlike V ones, do suffix the allele column - so before
    # IGHD was excluded this gene offered a selectable "01_F1". It now offers nothing,
    # both because the locus is not plotted and because _F rows are filtered out.
    assert get(client, url("/data/plotoptions", current_selection="IGHD1-1*")).get_json() == []


def test_db_name(client):
    # Plain gene + allele. The flanking row shares gene and allele with the row below
    # it, so this also asserts that the _F name is not what comes back.
    res = get(client, url("/data/db_name", selection=f"{TEST_GENE},01"))
    assert res.status_code == 200
    assert res.get_json() == {"db_name": TEST_ALLELE}

    # A gene whose alleles are only reachable through the two-gene db_name.
    res = get(client, url("/data/db_name", selection=f"IGHV3-30-5,{TWO_GENE_ALLELE}"))
    assert res.get_json() == {"db_name": TWO_GENE_ALLELE}
    res = get(client, url("/data/db_name", selection=f"IGHV3-30,{TWO_GENE_ALLELE}"))
    assert res.get_json() == {"db_name": TWO_GENE_ALLELE}

    # The same row named by the composite gene the loader actually stored. The selection is
    # split on its last comma, because a gene value can contain one itself; splitting on
    # every comma made this a 422. The frontend sends one component at a time, which is why
    # that never showed up in production.
    res = get(client, url("/data/db_name", selection=f"IGHV3-30,IGHV3-30-5,{TWO_GENE_ALLELE}"))
    assert res.get_json() == {"db_name": TWO_GENE_ALLELE}

    # A selection with no comma at all still cannot be split into a gene and an allele.
    assert get(client, url("/data/db_name", selection=TEST_GENE)).status_code == 422

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
    # The no-match entry used to key its empty list as "position" while every matching
    # entry used "positions", so the frontend - whose type declares positions - read
    # undefined here. The response schema does not allow the two to disagree.
    assert res.get_json() == [{"allele": "", "sequence": "", "positions": []}]


def test_sequencesearch_no_match(client):
    res = get(client, url("/data/sequences", sequence_str="NOMATCHFORTHIS"))
    assert res.status_code == 200
    assert res.get_json() == [{"allele": "", "sequence": "", "positions": []}]


def test_sequencesearch_rejects_wildcards_and_metacharacters(client):
    # Sequences are letters and digits, so the schema rejects both the SQL LIKE wildcards
    # ('%', '_') and the regex metacharacters that used to reach re.finditer.
    for value in ("%", "_", "__________", "(", "[", "(a+)+", "<script>"):
        res = get(client, url("/data/sequences", sequence_str=value))
        assert res.status_code == 422, f"{value!r} was accepted"


def test_sequence_search_escapes_like_wildcards(app):
    # The schema now stops a wildcard reaching the query, so this covers the layer below
    # it: .contains() without autoescape treated '%' and '_' as LIKE wildcards, matching
    # every row, so a search for a single '%' returned the whole sequence table.
    from services.sequence_search import sequence_search

    empty = [{"allele": "", "sequence": "", "positions": []}]
    assert sequence_search("%") == empty
    assert sequence_search("_") == empty
    # A real substring still matches, so the escaping has not broken ordinary searching.
    assert sequence_search("ESEARCHTES")[0]["allele"] == "IGHV3-23*01"


def test_send_fasta_genomic(client):
    res = get(client, url("/fasta/genomic", file_name="IGHV1-8"))
    assert res.status_code == 200
    # Served as text/x-fasta rather than the application/octet-stream mimetypes falls
    # back to for an unknown ".fasta" extension.
    assert "text/x-fasta" in res.content_type
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
    # Served as text/x-fasta rather than the application/octet-stream mimetypes falls
    # back to for an unknown ".fasta" extension.
    assert "text/x-fasta" in res.content_type
    # Only the flanking rows, and their sequence includes the prefix and suffix.
    assert res.get_data(as_text=True) == (
        ">IGHV1-8*01_F1\nAGGTGCCCACTCC"
        "CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC"
        "CACAGTGTGAAA\n"
    )


def test_send_fasta_translated(client):
    res = get(client, url("/fasta/translated", file_name="IGHV1-8"))
    assert res.status_code == 200
    # Served as text/x-fasta rather than the application/octet-stream mimetypes falls
    # back to for an unknown ".fasta" extension.
    assert "text/x-fasta" in res.content_type
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


def test_aminoacid_table_for_a_name_with_no_master_is_a_404(client):
    # get_aminoacid_top_allele returns {} when the name appears in no db_name_AA_list, and
    # subscripting that was a KeyError and a 500. Reachable with any schema-valid name that
    # has no translation: an allele not in the data, a flanking variant, an IGHD allele.
    for name in ("IGHV9-99*99", TEST_ALLELE + "_F1", "IGHD1-1*01"):
        res = get(client, url("/data/aminoacidfrequencies/table/allele", aa_allele_name=name))
        assert res.status_code == 404, name

    # The genomic table of the same allele is unaffected - it needs no translation.
    assert get(client, url("/data/frequencies/table/allele",
                           allele_name=TEST_ALLELE)).status_code == 200


def test_a_master_amino_acid_allele_of_null_has_no_table(app):
    # get_aminoacid_top_allele returns {'allele': ..., 'allele_aa': None} for a row that
    # carries db_name_AA_list without db_name_AA, and that dict is truthy - so guarding on it
    # rather than on the master let allele_name become None and produced a table of zeroes
    # naming an allele called "None". validate_loaded_data() rejects that row shape at
    # startup, which is why reaching it means writing the row in behind the loader.
    from db import db
    from models.immunediscoverdata import ImmuneDiscoverDataModel
    from services.frequencies import create_frequencies_table

    row = ImmuneDiscoverDataModel.query.filter_by(db_name=TEST_ALLELE).first()
    orphan = ImmuneDiscoverDataModel(**{
        column.name: getattr(row, column.name)
        for column in ImmuneDiscoverDataModel.__table__.columns if column.name != "id"
    })
    orphan.db_name, orphan.allele, orphan.case = "IGHV9-99*DEL", "DEL", "case_XTRA_EUR"
    orphan.db_name_AA, orphan.db_name_AA_list = None, "IGHV9-99*01"
    db.session.add(orphan)
    db.session.flush()

    assert create_frequencies_table("IGHV9-99*01", "aminoacid") is None


def test_the_amino_acid_lookups_404_a_name_they_cannot_resolve(client):
    # Both used to answer 200 with an empty body - {} and {"aa_allele_list": null} - while
    # every sibling on the blueprint 404s. The frontend reads them through
    # getJson(url, fallback), which turns a rejected request into the same fallback the empty
    # body produced, so the answer is unchanged for the only caller there is.
    for endpoint, param in (("/data/aminoacidalleles", "aa_allele_name"),
                            ("/data/aminoacidlist", "aa_allele_name")):
        assert get(client, url(endpoint, **{param: TEST_ALLELE})).status_code == 200
        for absent in ("IGHV9-99*99", TEST_ALLELE + "_F1", "IGHD1-1*01"):
            res = get(client, url(endpoint, **{param: absent}))
            assert res.status_code == 404, f"{endpoint} answered for {absent}"


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


def rows_per_allele(table, allele_column):
    counts = {}
    for row in table.split("\n")[1:]:
        allele = row.split("\t")[allele_column]
        counts[allele] = counts.get(allele, 0) + 1
    return counts


def test_gene_tables_list_each_allele_once(client):
    # IGHV3-30 is the gene whose alleles are stored under more than one gene value:
    # IGHV3-30*02/IGHV3-30-5*02 under "IGHV3-30,IGHV3-30-5" and IGHV3-30*02_S9143 under
    # "IGHV3-30", both collapsing into the same amino acid allele. Selecting the gene
    # column alongside them would make DISTINCT report that allele twice, once per gene
    # value, and repeat its 30 rows in the downloaded table.
    genomic = get(client, url("/data/frequencies/table/gene", gene_name="IGHV3-30"))
    assert genomic.status_code == 200
    assert rows_per_allele(genomic.get_data(as_text=True), 0) == {
        "IGHV3-30*01": 30,
        "IGHV3-30*02_S9143": 30,
        TWO_GENE_ALLELE: 30,
    }

    aminoacid = get(client, url("/data/aminoacidfrequencies/table/gene", aa_gene_name="IGHV3-30"))
    assert aminoacid.status_code == 200
    assert rows_per_allele(aminoacid.get_data(as_text=True), 1) == {
        "IGHV3-30*01": 30,
        TWO_GENE_ALLELE: 30,
    }


def test_openapi_spec_options_are_not_shared_between_apps(monkeypatch):
    # from_object copies the config class attribute by reference, and apispec's to_dict()
    # ends by deep-merging the spec it just built into it. Shared, the first app's schemas
    # accumulated there and won for every app built afterwards in the process.
    from app import create_app
    from config import TestConfig

    # Neither app is the server and neither has a schema; only the generated spec is read.
    monkeypatch.setenv("SKIP_DATA_LOAD", "1")
    first, second = create_app(TestConfig), create_app(TestConfig)
    assert first.config["API_SPEC_OPTIONS"] is not second.config["API_SPEC_OPTIONS"]

    with first.test_client() as client:
        client.get("/openapi.json")

    assert sorted(TestConfig.API_SPEC_OPTIONS["components"]) == ["securitySchemes"], \
        "the config class attribute accumulated a built spec"


def test_health_is_documented_as_needing_no_api_key(client):
    # The requirement is declared at document root so every operation inherits it, but this
    # route has no api_key_required. Anyone wiring a Kubernetes probe from the spec would
    # build it to send a header it does not need.
    spec = get(client, "/openapi.json").get_json()
    assert spec["paths"]["/health"]["get"]["security"] == []
    assert client.get("/health").status_code == 200


def test_unknown_fasta_type_and_plot_type_raise(app):
    # Both chose on the argument with if/elif and no else: the fasta branches left
    # distinct_sequences unassigned and failed with UnboundLocalError, and the frequency
    # table silently treated anything unrecognised as amino acid, answering 200 with a
    # header-only tsv for a full gene.
    from services.fasta_generation import generate_fasta
    from services.frequencies import create_frequencies_table

    with pytest.raises(ValueError, match="unknown fasta type"):
        generate_fasta(TEST_GENE, type="genomicc")
    with pytest.raises(ValueError, match="unknown plot_type"):
        create_frequencies_table(TEST_GENE, "genomicc", full_gene=True)


def test_empty_plot_selection_is_not_an_error(client):
    # Nothing chosen yet is a legitimate state, not a malformed request. The frontend sends
    # an empty selection on first load and whenever a selection is reset, so requiring at
    # least one character turned both into a 422 that getJson quietly swallowed.
    for u in ("/data/plotoptions?current_selection=", "/data/plotoptions"):
        res = get(client, u)
        assert res.status_code == 200, f"{u} returned {res.status_code}"
        assert res.get_json() == []

    # A real selection is unaffected.
    assert get(client, url("/data/plotoptions", current_selection="IGHV")).status_code == 200


def test_aminoacid_table_resolves_the_master_allele(app):
    # The frontend asks for this table by *genomic* allele name - that is what it resolved
    # from the plot selection - so it has to be translated to the master amino acid allele it
    # collapses into. 243 of the 732 plottable alleles are not their own master, and checking
    # the cache before resolving made every one of them a 404.
    from services.frequencies import create_frequencies_table

    # IGHV1-8*02 and *04 both collapse into IGHV1-8*01 in the mock data.
    for genomic in (TEST_ALLELE, "IGHV1-8*02", "IGHV1-8*04"):
        table = create_frequencies_table(genomic, "aminoacid")
        assert table is not None, f"no amino acid table for {genomic}"
        assert TEST_ALLELE in table, f"{genomic} did not resolve to its master"

    # A *DEL has no amino acid master, and a non-plottable allele none either.
    assert create_frequencies_table("IGHV1-8*DEL", "aminoacid") is None
    assert create_frequencies_table("IGHD1-1*01", "aminoacid") is None


def test_alignment_covers_only_the_plotted_loci(client):
    # repositories/filters.py says IGHD and IGHJ must never be offered as a plot or
    # alignment selection, and this was the one query in that family without the locus
    # restriction. Handed no sequences, MAFFT also produced one blank record rather than
    # nothing, so the answer was 200 with a row of empty strings.
    assert get(client, url("/data/sequences/alignedsequences", gene_name=TEST_GENE)).status_code == 200

    for absent in ("IGHD1-1", "IGHJ6", "IGHV9-99"):
        res = get(client, url("/data/sequences/alignedsequences", gene_name=absent))
        assert res.status_code == 404, f"aligned {absent}, which no plot offers"



def test_the_preload_replaces_rather_than_merges(app):
    # target.update() alone left alleles from an earlier load sitting alongside the new ones.
    # The dictionaries are module-level and shared by name, so they are cleared in place
    # rather than rebound - and the invariant belongs here, not in the test fixture.
    from constants import allele_superpopulation_frequencies
    from loaders import load_plot_data_to_dict

    before = len(allele_superpopulation_frequencies)
    allele_superpopulation_frequencies["STALE*99"] = [{"population": "AFR", "n": 0, "frequency": 0.0}]

    load_plot_data_to_dict()

    assert "STALE*99" not in allele_superpopulation_frequencies
    assert len(allele_superpopulation_frequencies) == before
