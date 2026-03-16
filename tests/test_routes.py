# Pytest tests for all endpoints

import os

from urllib.parse import quote

api_key_header = {"X-api-key": os.getenv("API_KEY")}
expected_superpops = {"AFR", "EUR", "EAS", "SAS", "AMR", "ALL"}
expected_pops = {
    "ACB",
    "ASW",
    "PJL",
    "STU",
    "ITU",
    "BEB",
    "GIH",
    "CHS",
    "CDX",
    "KHV",
    "CHB",
    "JPT",
    "CLM",
    "GWD",
    "ESN",
    "MSL",
    "YRI",
    "FIN",
    "GBR",
    "IBS",
    "LWK",
    "TSI",
    "MXL",
    "PEL",
    "PUR",
    "ALL"
    }
test_allele_name = "TEST1-8*01"


def test_health(client):
    res = client.get("/health")
    assert res.status_code == 200
    assert res.json["status"] == "ok"

def test_superpopfrequencies(client):
    res = client.get(
        f"/data/frequencies/superpopulations?allele_name={quote(test_allele_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    data = res.get_json()

    # Invariant 1: Expected keys exist
    for entry in data:
        assert "population" in entry
        assert "n" in entry
        assert "frequency" in entry

    # Invariant 2: All known superpopulations + 'ALL' exist in output
    returned_pops = {entry["population"] for entry in data}
    assert expected_superpops.issubset(returned_pops)

    # Optional: Compare to manual recalculation if feasible
    # expected = calculate_frequencies("TEST1-8*01", "superpopulation", "genomic")
    # assert data == expected

    # Optional: Snapshot testing
    # from approvaltests import verify
    # verify(data)

def test_aminoacid_superpopfrequencies(client):
    res = client.get(
        f"/data/aminoacidfrequencies/superpopulations?aa_allele_name={quote(test_allele_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    data = res.get_json()
    # Invariant 1: Expected keys exist
    for entry in data:
        assert "population" in entry
        assert "n" in entry
        assert "frequency" in entry

    # Invariant 2: All known superpopulations + 'ALL' exist in output
    returned_pops = {entry["population"] for entry in data}
    assert expected_superpops.issubset(returned_pops)

def test_populationfrequencies(client):
    res = client.get(
        f"/data/frequencies/populations?allele_name={quote(test_allele_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    data = res.get_json()

    # Invariant 1: Expected keys exist
    for entry in data:
        assert "population" in entry
        assert "n" in entry
        assert "frequency" in entry

    # Invariant 2: All known superpopulations + 'ALL' exist in output
    returned_pops = {entry["population"] for entry in data}
    assert expected_pops.issubset(returned_pops)

def test_aminoacid_populationfrequencies(client):
    res = client.get(
        f"/data/aminoacidfrequencies/populations?allele_name={quote(test_allele_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    data = res.get_json()

    # Invariant 1: Expected keys exist
    for entry in data:
        assert "population" in entry
        assert "n" in entry
        assert "frequency" in entry

    # Invariant 2: All known superpopulations + 'ALL' exist in output
    returned_pops = {entry["population"] for entry in data}
    assert expected_pops.issubset(returned_pops)

def test_igsnperdata(client):
    res = client.get(
        f"/data/igsnperdata?allele_name={quote(test_allele_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    data = res.get_json()
    assert "igSNPer_SNPs" in data
    assert "igSNPer_score" in data

def test_aminoacidalleles(client):
    res = client.get(
        f"/data/aminoacidalleles?aa_allele_name={quote(test_allele_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    data = res.get_json()
    assert "allele" in data
    assert "allele_aa" in data

def test_aminoacidlist(client):
    res = client.get(
        f"/data/aminoacidlist?aa_allele_name={quote(test_allele_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    data = res.get_json()
    assert "aa_allele_list" in data

def test_populationregions(client):
    expected_output_content = [
        {
            "population": "ACB",
            "superpopulation": "AFR"
        },
        {
            "population": "ASW",
            "superpopulation": "AFR"
        },
        {
            "population": "PJL",
            "superpopulation": "SAS"
        },
        {
            "population": "STU",
            "superpopulation": "SAS"
        },
        {
            "population": "ITU",
            "superpopulation": "SAS"
        },
        {
            "population": "BEB",
            "superpopulation": "SAS"
        },
        {
            "population": "GIH",
            "superpopulation": "SAS"
        },
        {
            "population": "CHS",
            "superpopulation": "EAS"
        },
        {
            "population": "CDX",
            "superpopulation": "EAS"
        },
        {
            "population": "KHV",
            "superpopulation": "EAS"
        },
        {
            "population": "CHB",
            "superpopulation": "EAS"
        },
        {
            "population": "JPT",
            "superpopulation": "EAS"
        },
        {
            "population": "CLM",
            "superpopulation": "AMR"
        },
        {
            "population": "GWD",
            "superpopulation": "AFR"
        },
        {
            "population": "ESN",
            "superpopulation": "AFR"
        },
        {
            "population": "MSL",
            "superpopulation": "AFR"
        },
        {
            "population": "YRI",
            "superpopulation": "AFR"
        },
        {
            "population": "FIN",
            "superpopulation": "EUR"
        },
        {
            "population": "GBR",
            "superpopulation": "EUR"
        },
        {
            "population": "IBS",
            "superpopulation": "EUR"
        },
        {
            "population": "LWK",
            "superpopulation": "AFR"
        },
        {
            "population": "TSI",
            "superpopulation": "EUR"
        },
        {
            "population": "MXL",
            "superpopulation": "AMR"
        },
        {
            "population": "PEL",
            "superpopulation": "AMR"
        },
        {
            "population": "PUR",
            "superpopulation": "AMR"
        },
        {
            "population": "ALL",
            "superpopulation": "ALL"
        }
    ]
    res = client.get(
        "/data/populationregions",
        headers=api_key_header
        )
    assert res.status_code == 200
    data = res.get_json()
    for expected_item in expected_output_content:
        assert(expected_item in data)

def test_plotoptions(client):
    current_selection = "TEST1-8*"
    res = client.get(
        f"/data/plotoptions?current_selection={quote(current_selection, safe='')}",
        headers=api_key_header
    )
    current_selection2 = "TEST1-69/1-69D*"
    res2 = client.get(
        f"/data/plotoptions?current_selection={quote(current_selection2, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    assert res2.status_code == 200
    data = res.get_json()
    data2 = res2.get_json()
    assert data == ['01', '04']
    assert data2 == ['01_S01']

def test_alignedsequences(client):
    alignment_gene_test = "ALIGNMENTTEST1-2"
    res = client.get(
        f"/data/sequences/alignedsequences?gene_name={quote(alignment_gene_test, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    data = res.get_json()
    assert data == [{'allele': 'ALIGNMENTTEST1-2*01', 'sequence_aa': 'LDSPLLALLCSGCDRLVDNA', 'sequence_nt': 'CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC'}, {'allele': 'ALIGNMENTTEST1-2*02', 'sequence_aa': 'LDSPL-RMLCSGCDRLVDNA', 'sequence_nt': 'CTGGATTCACCTTTAC---GTATGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC'}, {'allele': 'ALIGNMENTTEST1-2*03', 'sequence_aa': 'LDSPLLALLCSGCXXXXXXX', 'sequence_nt': 'CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGA-AGGCTCGTGGACAACGCC'}]

def test_sequencesearch(client):
    sequence_str = "ESEARCHTES"
    res = client.get(
        f"/data/sequences?sequence_str={quote(sequence_str, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    data = res.get_json()
    assert data == [{'allele': 'SEQTEST1-2*01', 'positions': [14], 'sequence': 'THISISASEQUENCESEARCHTEST123'}]

def test_send_fasta_genomic(client):
    file_name = "testfile.tsv.gz"
    res = client.get(
        f"/fasta/genomic?file_name={quote(file_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    assert "text/plain" in res.content_type or "application/octet-stream" in res.content_type

def test_send_fasta_genomicwithflanking(client):
    file_name = "testfile.tsv.gz"
    res = client.get(
        f"/fasta/genomic_fl?file_name={quote(file_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    assert "text/plain" in res.content_type or "application/octet-stream" in res.content_type

def test_send_fasta_translated(client):
    file_name = "testfile.tsv.gz"
    res = client.get(
        f"/fasta/translated?file_name={quote(file_name, safe='')}",
        headers=api_key_header
    )
    assert res.status_code == 200
    assert "text/plain" in res.content_type or "application/octet-stream" in res.content_type

def test_checkapikey(client):
    res = client.get(
        "/checkapikey",
        headers=api_key_header
    )
    assert res.status_code == 200
    assert res.get_data(as_text=True) == "Correct key!"