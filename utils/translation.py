# Translation of nucleotide triplets to amino acids using the standard genetic
# code (NCBI translation table 1).
#
# This replaces biopython's Bio.Seq.Seq.translate(), which was the only thing
# biopython was used for. Dropping it also drops numpy, which biopython
# requires. The genetic code is a fixed biological constant, so there is nothing
# here that needs to track upstream releases.
#
# Behaviour is matched to biopython's defaults, so results are unchanged:
#   - stop codons translate to '*'
#   - a codon of three gaps ('---') translates to '-'
#   - IUPAC ambiguity codes are resolved when every base they could stand for
#     gives the same amino acid (so 'GGN' -> 'G'), and give 'X' otherwise
#     (so 'TTN' -> 'X')
#   - a trailing partial codon is ignored

# The standard genetic code, in the conventional NCBI base ordering: the amino
# acids below line up with codons generated as T/C/A/G in the first position,
# then the second, then the third.
_BASES = "TCAG"
_AMINO_ACIDS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)

CODON_TABLE = {
    b1 + b2 + b3: _AMINO_ACIDS[i]
    for i, (b1, b2, b3) in enumerate(
        (x, y, z) for x in _BASES for y in _BASES for z in _BASES
    )
}

# Which concrete bases each IUPAC ambiguity code can stand for.
_AMBIGUITY_CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

# Where an ambiguous codon narrows down to exactly one of these amino acid
# pairs, there is a standard letter for that pair rather than a plain 'X'.
_AMBIGUOUS_AMINO_ACIDS = {
    frozenset("ND"): "B",   # Asx
    frozenset("QE"): "Z",   # Glx
    frozenset("LI"): "J",   # Xle
}


def translate_codon(codon):
    """Translate a single three-base codon to its amino acid letter."""
    if codon in CODON_TABLE:
        return CODON_TABLE[codon]

    if codon == "---":
        return "-"

    # An ambiguity code is only translatable if every concrete codon it could
    # stand for yields the same amino acid.
    options = [_AMBIGUITY_CODES.get(base) for base in codon]
    if any(option is None for option in options):
        return "X"

    amino_acids = {
        CODON_TABLE[b1 + b2 + b3]
        for b1 in options[0]
        for b2 in options[1]
        for b3 in options[2]
    }
    if len(amino_acids) == 1:
        return amino_acids.pop()
    return _AMBIGUOUS_AMINO_ACIDS.get(frozenset(amino_acids), "X")


def translate(sequence):
    """Translate a nucleotide sequence, ignoring any trailing partial codon."""
    return "".join(
        translate_codon(sequence[i:i + 3])
        for i in range(0, len(sequence) - len(sequence) % 3, 3)
    )
