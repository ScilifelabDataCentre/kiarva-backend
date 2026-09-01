# Shared query criteria for the plot and MSA data paths.
#
# Two categories of row are loaded into the db but must never be offered as a
# plot or alignment selection:
#
# - Flanking-region variants, which carry a '_F' suffix on db_name. The research
#   group has specified these are not to be shown in plots or alignments. They
#   stay in the db because the "Flanking genomic" FASTA download is built from
#   exactly these rows (services/fasta_generation.py, type='genomic_fl').
# - IGHD and IGHJ genes. Out of the IGH loci only IGHV is plotted, per the
#   research group; TRGV is the other plotted locus. IGHD/IGHJ rows are still
#   used by the FASTA downloads.
#
# Excluding '_F' has a second effect worth knowing about: it makes (gene, allele)
# resolve to exactly one db_name. With flanking variants included, 553 of 732
# IGHV/TRGV (gene, allele) pairs are ambiguous, because a flanking variant shares
# the 'allele' value of its parent for these loci - IGHV1-2*02_F1 has allele '02',
# the same as IGHV1-2*02. See get_db_name_from_options in repositories/allele.py.

from sqlalchemy import or_

from models.immunediscoverdata import ImmuneDiscoverDataModel

# Loci offered as plot / MSA selections. These are the V genes, which is also exactly the
# set whose alleles have an amino acid sequence - V genes are translated, D and J genes are
# not - so loaders/validation.py derives "never translated" from this rather than listing
# the other loci separately. Keeping one definition means a locus that turns up in future
# source data cannot be plotted by one rule and expected to be translated by another.
PLOT_LOCI = ("IGHV", "TRGV")

# Every locus the source data is known to contain. PLOT_LOCI is the subset that is plotted;
# the rest are loaded for the FASTA downloads. A locus outside this list means the source
# data has gained something nobody has decided how to present, which loaders/validation.py
# reports rather than quietly leaving out of every plot.
KNOWN_LOCI = ("IGHV", "IGHD", "IGHJ", "TRGV")

# Gene names are prefixed with their locus, e.g. IGHV1-2, IGHD5-18/5-5, TRGV9.
def locus_of(gene):
    """The locus a gene name belongs to, or None if it is not one we know.

    Takes one gene name, not a gene column value: the column can hold a comma-separated
    list, and prefix-matching the whole thing answers with the first component's locus and
    says nothing about the rest. loaders/validation.py calls this per component and refuses
    to boot on a value whose components do not all resolve to one locus, which is what makes
    the whole-value matching in in_plot_loci() below correct rather than lucky.

    Matched against KNOWN_LOCI rather than sliced at a fixed width. Every locus here happens
    to be four characters, but slicing means a shorter gene name yields a truncated string
    that is in no list - so the service would refuse to boot citing a locus that does not
    exist, which is the opposite of the clear diagnosis this check is for.
    """
    for locus in KNOWN_LOCI:
        if gene.startswith(locus):
            return locus
    return None

def allele_column(plot_type):
    """The column holding allele names for a plot type.

    Raises rather than falling off the end of an if/elif. Three call sites chose this
    column that way, leaving the variable unassigned for an unrecognised plot_type - so
    the next line raised UnboundLocalError and surfaced as a 500 complaining about a
    local variable rather than about the argument.
    """
    if plot_type == "genomic":
        return ImmuneDiscoverDataModel.db_name
    if plot_type == "aminoacid":
        return ImmuneDiscoverDataModel.db_name_AA
    raise ValueError("unknown plot_type: " + repr(plot_type))

def is_deletion():
    """True for the homozygous deletion rows.

    One definition, because there were two: loaders/validation.py matched allele == 'DEL'
    while the FASTA and MSA queries matched db_name NOT LIKE '%*DEL'. A *DEL row whose allele
    column was spelled differently would have been reported as "should be translated" and
    stopped the service booting, while the FASTA queries excluded it correctly.
    """
    return or_(ImmuneDiscoverDataModel.allele == 'DEL',
               ImmuneDiscoverDataModel.db_name.like('%*DEL'))

def in_plot_loci():
    """True for rows belonging to a locus that is plotted, and therefore translated.

    Matched against the whole gene column value rather than per comma-separated component,
    which is safe only because validate_loaded_data() rejects a value whose components span
    more than one locus. Without that, a composite like "IGHV1-8,TRGJ1" would be plotted on
    the strength of its first component and offer the second's allele under the wrong gene.
    Kept as one LIKE per locus because this has to be SQL the database can index, not Python.
    """
    return or_(*(ImmuneDiscoverDataModel.gene.like(locus + '%') for locus in PLOT_LOCI))

def plot_selection_criteria():
    """Criteria restricting a query to the rows eligible for plots and MSA.

    Splat into an existing filter chain: .filter(*plot_selection_criteria())
    """
    return (
        ~ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True),
        in_plot_loci(),
    )
