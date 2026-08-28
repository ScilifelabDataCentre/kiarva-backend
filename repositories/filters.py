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
LOCUS_PREFIX_LENGTH = 4

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

def in_plot_loci():
    """True for rows belonging to a locus that is plotted, and therefore translated."""
    return or_(*(ImmuneDiscoverDataModel.gene.like(locus + '%') for locus in PLOT_LOCI))

def plot_selection_criteria():
    """Criteria restricting a query to the rows eligible for plots and MSA.

    Splat into an existing filter chain: .filter(*plot_selection_criteria())
    """
    return (
        ~ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True),
        in_plot_loci(),
    )
