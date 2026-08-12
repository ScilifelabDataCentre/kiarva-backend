# Tests for the standard genetic code translation in utils/translation.py.
#
# These previously came for free from biopython. Since that dependency was
# dropped, the expected values below were generated from biopython's
# Bio.Seq.Seq.translate() and checked to match it exactly for every codon over
# the full IUPAC alphabet, so they pin the behaviour that used to be inherited.

import pytest

from utils.translation import CODON_TABLE, translate, translate_codon

# The standard genetic code, in the conventional NCBI base ordering.
EXPECTED_CODE = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)


def test_codon_table_is_the_standard_genetic_code():
    bases = "TCAG"
    codons = [b1 + b2 + b3 for b1 in bases for b2 in bases for b3 in bases]
    assert len(CODON_TABLE) == 64
    assert "".join(CODON_TABLE[codon] for codon in codons) == EXPECTED_CODE


@pytest.mark.parametrize(
    "codon,expected",
    [
        ("ATG", "M"),   # start
        ("TAA", "*"),   # stop
        ("TAG", "*"),
        ("TGA", "*"),
        ("---", "-"),   # a fully gapped codon stays a gap
    ],
)
def test_unambiguous_codons(codon, expected):
    assert translate_codon(codon) == expected


@pytest.mark.parametrize(
    "codon,expected",
    [
        ("GGN", "G"),   # every GGx is glycine, so this resolves
        ("CGN", "R"),
        ("ACN", "T"),
        ("TTY", "F"),
        ("RAY", "B"),   # Asn or Asp
        ("SAR", "Z"),   # Gln or Glu
        ("MTA", "J"),   # Leu or Ile
        ("TTN", "X"),   # Phe or Leu, no single letter for that pair
        ("NNN", "X"),
        ("ATN", "X"),
    ],
)
def test_ambiguity_codes(codon, expected):
    assert translate_codon(codon) == expected


def test_unknown_characters_translate_to_x():
    assert translate_codon("AZG") == "X"
    assert translate_codon("A-G") == "X"


@pytest.mark.parametrize(
    "sequence,expected",
    [
        ("", ""),
        ("ATGGCCATTGTAATGGGCCGC", "MAIVMGR"),
        ("ATG---GCC", "M-A"),
        ("TAAATG", "*M"),
    ],
)
def test_sequences(sequence, expected):
    assert translate(sequence) == expected


@pytest.mark.parametrize("sequence", ["A", "AT", "ATGA", "ATGAT"])
def test_trailing_partial_codon_is_ignored(sequence):
    assert translate(sequence) == translate(sequence[: len(sequence) // 3 * 3])
