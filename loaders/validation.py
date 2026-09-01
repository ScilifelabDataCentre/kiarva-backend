# Validation of the source data once it has been loaded into the db.
#
# These checks raise rather than warn. The research group asked to be told about data
# problems immediately: a wrong frequency in a published resource is worse than a service
# that will not start, and a warning in a log nobody reads is not a report. They run in
# every mode - pytest, local dev and production - so a broken TSV is caught while it is
# still being worked on rather than after it ships.
#
# Raising during app creation means the Gunicorn worker never boots and the pod never
# passes its readiness check. With a Deployment and maxUnavailable: 0 the previous pod goes
# on serving, so bad data cannot reach users and CrashLoopBackOff is the alarm. The data is
# baked into the image, so it can only change when a new image ships.

from sqlalchemy import and_, func, or_

from db import db
from models.immunediscoverdata import ImmuneDiscoverDataModel
from repositories.filters import in_plot_loci, is_deletion, locus_of, plot_selection_criteria
from services.frequencies import subpopulation_order, superpopulation_order

class SourceDataError(Exception):
    """Raised when loaded data breaks an assumption the rest of the app relies on."""

# A row carries no amino acid data exactly when it is one of the kinds that is never
# translated: a flanking region variant, a homozygous deletion, or a gene outside the V loci.
#
# The last of those is derived from PLOT_LOCI rather than listing IGHD and IGHJ, which is
# what the D and J genes in the current data happen to be. Spelling them out made the two
# definitions hand-maintained complements of each other: a TRGJ row appearing in future
# source data would have been outside PLOT_LOCI, so never plotted, and simultaneously
# outside the IGHD/IGHJ list, so expected to be translated - and the service would have
# refused to boot reporting "should be translated" about a gene that never is.
def never_translated():
    return or_(
        ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True),
        is_deletion(),
        ~in_plot_loci(),
    )

def sample_names(rows, limit = 10):
    names = sorted(row[0] for row in rows)
    if len(names) <= limit:
        return ", ".join(names)
    return ", ".join(names[:limit]) + ", ... (" + str(len(names) - limit) + " more)"

def amino_acid_coverage_problems(unknown_loci = ()):
    """db_name_AA must be null on exactly the rows that are never translated.

    Both directions matter. A translated allele missing db_name_AA drops out of the amino
    acid plots and the translated FASTA download without any error; a flanking or *DEL row
    that has db_name_AA set would appear in them, and would also break the translated
    FASTA branch in services/fasta_generation.py, which has no exclusion of its own and
    relies on db_name_AA being null to leave those rows out.

    Rows of an unrecognised locus are excluded, because locus_problems() already reports
    them and this check would say something true but unhelpful about the same rows. Only
    those rows: suppressing the whole check whenever any unknown locus exists would hide a
    genuine missing db_name_AA on a known allele until the unrelated locus was dealt with,
    which is the extra round trip reporting everything at once exists to avoid.

    db_name_AA_list is checked alongside it, in both directions too. The list is a
    property of db_name_AA - the master name plus every genomic allele collapsing into it
    - so the two are set together or not at all, which is exactly how the data reads
    today. Three callers take that for granted: get_aminoacid_allele_list() measures the
    list, the full-gene amino acid download in services/frequencies.py splits it, and
    get_aminoacid_top_allele() matches a regex against it to find the master. A list
    missing from a translated row is a 500 in the first two; a list on a row with no
    master resolves that master to None and 404s instead. Asserted here rather than
    guarded at each caller, so one broken row is reported once as a data problem, with its
    allele name, instead of as three different symptoms.
    """
    problems = []
    known_locus_rows = (ImmuneDiscoverDataModel.gene.notin_(unknown_loci)
                        if unknown_loci else True)

    missing = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name
        ).filter(ImmuneDiscoverDataModel.db_name_AA == None
        ).filter(known_locus_rows
        ).filter(~never_translated()).distinct().all()
    if missing:
        problems.append(
            str(len(missing)) + " allele(s) have no db_name_AA but should be translated: "
            + sample_names(missing))

    unexpected = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name
        ).filter(ImmuneDiscoverDataModel.db_name_AA != None
        ).filter(known_locus_rows
        ).filter(never_translated()).distinct().all()
    if unexpected:
        problems.append(
            str(len(unexpected)) + " allele(s) have db_name_AA set but are never translated"
            " (a flanking variant, a *DEL, or a gene outside the plotted V loci): "
            + sample_names(unexpected))

    # Restricted the same way masterless below is, and for the same reason: a row that is
    # never translated yet carries db_name_AA is already reported by the check above, and
    # adding "and no list" names the same allele twice under two headings.
    listless = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name
        ).filter(ImmuneDiscoverDataModel.db_name_AA != None
        ).filter(ImmuneDiscoverDataModel.db_name_AA_list == None
        ).filter(known_locus_rows
        ).filter(~never_translated()).distinct().all()
    if listless:
        problems.append(
            str(len(listless)) + " allele(s) have db_name_AA but no db_name_AA_list, which"
            " three callers read as the list of alleles collapsing into it: "
            + sample_names(listless))

    # Restricted to the rows that should carry neither. A translated row missing db_name_AA
    # while keeping its list is the same row the first check above already reports, and
    # reporting it twice under two headings is the overlap this module avoids elsewhere -
    # so what is left to say here is about a flanking, *DEL or non-V row, where the null
    # db_name_AA is correct and the stray list is the whole problem.
    masterless = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name
        ).filter(ImmuneDiscoverDataModel.db_name_AA_list != None
        ).filter(ImmuneDiscoverDataModel.db_name_AA == None
        ).filter(known_locus_rows
        ).filter(never_translated()).distinct().all()
    if masterless:
        problems.append(
            str(len(masterless)) + " allele(s) are never translated but have a"
            " db_name_AA_list, which get_aminoacid_top_allele() matches on and would resolve"
            " to a master of None: " + sample_names(masterless))

    return problems

def gene_column_values():
    return {row[0] for row in ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.gene).distinct().all()}

def unrecognised_loci():
    """Gene column values naming a gene whose locus is not one the app knows to present.

    Tested per comma-separated component. The column can hold two gene names at once - the
    loader builds "IGHV3-30,IGHV3-30-5" itself - and locus_of() prefix-matches whatever it is
    given, so asking it about the whole value returns the first component's locus and never
    looks at the rest. A composite whose second component belonged to an unknown locus was
    therefore reported by nothing, and in_plot_loci() matched it on the first component's
    prefix - so its allele was offered under the wrong gene and resolved to the wrong
    db_name, which is the silent wrong answer this check exists to prevent.

    The whole column value is returned rather than the offending component, because that is
    what amino_acid_coverage_problems() filters those rows out by.
    """
    return sorted(gene for gene in gene_column_values()
                  if any(locus_of(component) is None for component in gene.split(",")))

def mixed_locus_genes():
    """Gene column values naming genes from more than one locus.

    Everything that matches a locus by prefix - in_plot_loci(), never_translated(),
    locus_of() itself - reads the whole column value, so a value spanning two loci is plotted
    or translated according to whichever comes first. The comma-separated values in the data
    are all one locus (they come from one rewrite dictionary, every entry IGHV), and this is
    what keeps that true, so the prefix matching elsewhere does not have to be per component.
    """
    mixed = []
    for gene in gene_column_values():
        loci = {locus_of(component) for component in gene.split(",")}
        if len(loci) > 1 and None not in loci:
            mixed.append(gene)
    return sorted(mixed)

def locus_problems(unknown):
    """Every locus in the data must be one the app knows how to present.

    Deriving "never translated" from PLOT_LOCI means an unrecognised locus is treated as
    untranslated and left out of the plots without complaint, which is consistent but
    silent. A locus nobody has decided about is a change the research group should have
    announced, so it is reported here - with the right diagnosis, rather than as a
    translation problem.
    """
    if not unknown:
        return []

    return ["Gene(s) belonging to no known locus: " + sample_names([(g,) for g in unknown])
            + ". Add the locus to KNOWN_LOCI in repositories/filters.py once it is decided "
            "whether it should be plotted, which also decides whether its alleles are "
            "expected to be translated"]

def mixed_locus_problems():
    """A gene column value must not name genes from more than one locus."""
    mixed = mixed_locus_genes()
    if not mixed:
        return []

    return [str(len(mixed)) + " gene value(s) name genes from more than one locus, so every "
            "prefix match on that column - whether the row is plotted, whether it is expected "
            "to be translated, which locus it is reported as - answers on whichever component "
            "comes first: " + sample_names([(g,) for g in mixed])]

def igsnper_consistency_problems():
    """An allele's IgSNPer values must not differ between the rows carrying it.

    The researchers look these up per allele name, so cohort cannot change them - which is
    why get_igSNPer_data() can answer from the first row it gets back. It queries with
    .distinct() and no ORDER BY, so if one allele ever carried two different values the row
    that won would be whichever the query plan returned first, and a real score could
    disappear silently. Zero alleles diverge today; this is what keeps that true.
    """
    # A group diverges if it holds more than one non-null value, or one non-null value on
    # some rows and nothing on others. Counted rather than coalesced to a sentinel:
    # IgSNPer_uncommon is a Float, and coalescing it to "~" only works because SQLite types
    # values rather than columns. SQLALCHEMY_DATABASE_URI is a DATABASE_URL override, so
    # against PostgreSQL that sentinel would make this startup check itself the crash
    # ("invalid input syntax for type double precision"). COUNT(DISTINCT col) ignores nulls
    # on both engines, which is what makes the two comparisons enough.
    def diverges(column):
        values = func.count(func.distinct(column))
        return or_(values > 1, and_(values == 1, func.count(column) < func.count()))

    divergent = db.session.query(ImmuneDiscoverDataModel.db_name).group_by(
        ImmuneDiscoverDataModel.db_name
        ).having(or_(diverges(ImmuneDiscoverDataModel.IgSNPer_uncommon),
                     diverges(ImmuneDiscoverDataModel.IgSNPer_SNPs))).all()

    if not divergent:
        return []

    return [str(len(divergent)) + " allele(s) carry more than one set of IgSNPer values, so "
            "which one is reported depends on the query plan: " + sample_names(divergent)]

def amino_acid_list_problems():
    """An amino acid master must carry exactly one db_name_AA_list.

    The same shape as igsnper_consistency_problems() above, for the other column whose value
    is a property of a name rather than of a row. db_name_AA_list is the master plus every
    genomic allele collapsing into it, so every row sharing a db_name_AA has to repeat the
    same list - and the full-gene amino acid download selects the two columns together and
    relies on a plain DISTINCT returning one row per master. Two lists for one master make
    DISTINCT return two, and the allele is written into the table twice: the mock data goes
    from 30 rows for IGHV1-8*01 to 60 with one divergent row added.

    No coalesce or null handling needed here, unlike the IgSNPer check: the pair check in
    amino_acid_coverage_problems() already reports a master without a list, so anything
    reaching this one has a list on every row.
    """
    divergent = db.session.query(ImmuneDiscoverDataModel.db_name_AA).filter(
        ImmuneDiscoverDataModel.db_name_AA != None
        ).group_by(ImmuneDiscoverDataModel.db_name_AA
        ).having(func.count(
            func.distinct(ImmuneDiscoverDataModel.db_name_AA_list)) > 1).all()

    if not divergent:
        return []

    return [str(len(divergent)) + " amino acid allele(s) carry more than one db_name_AA_list, "
            "so the full-gene download repeats them once per list: " + sample_names(divergent)]

def population_problems():
    """Every population in the data must appear in the hardcoded display order.

    The frequency endpoints report populations in the order the research group asked for,
    and a value missing from those lists is silently left out of every plot rather than
    reported - so an unrecognised population is most likely a typo in the source data.
    """
    problems = []

    for column, order, label in (
        (ImmuneDiscoverDataModel.superpopulation, superpopulation_order, "superpopulation"),
        (ImmuneDiscoverDataModel.population, subpopulation_order, "population"),
    ):
        in_data = {row[0] for row in
                   ImmuneDiscoverDataModel.query.with_entities(column).distinct().all()}
        unknown = sorted(in_data - set(order))
        if unknown:
            problems.append(
                "Unknown " + label + "(s) " + sample_names([(u,) for u in unknown])
                + ", which are not in the "
                "display order in services/frequencies.py and would be left out of every plot")

    return problems

def allele_resolution_problems():
    """(gene, allele) must identify exactly one db_name once flanking rows are excluded.

    get_db_name_from_options resolves a plot selection by taking the first row it gets back,
    which is only correct because this holds. It does not hold with flanking variants
    included - they share their parent's allele value - which is why they are excluded from
    the plot and MSA paths.

    Grouped over the comma-separated components of the gene column rather than its literal
    value, because that is how the resolver matches: plot_options_regex treats a composite
    value like "IGHV3-30,IGHV3-30-5" as naming both genes, and /data/plotoptions offers them
    as two separate options. So the rows a selection can reach span every literal value that
    names that gene, and grouping literally split them into separate groups that each looked
    unambiguous - a new composite row colliding with an existing plain one would have passed
    while making the resolver's answer arbitrary.

    Splitting on commas covers what the resolver matches for whole components, but not quite
    everything it matches: plot_options_regex makes the leading comma optional, so a
    selection that is only the tail of a component - "30-5" against "IGHV3-30-5" - resolves
    too, and no group here holds it. That is left alone deliberately. Every gene name in the
    data is prefixed with its locus, so none is a suffix of another (93 components, 0 such
    pairs), and /data/plotoptions only ever offers whole components - so a selection of that
    shape cannot come from the frontend and cannot collide with a real gene.
    """
    rows = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.gene,
        ImmuneDiscoverDataModel.allele,
        ImmuneDiscoverDataModel.db_name,
        ).filter(*plot_selection_criteria()).distinct().all()

    by_selection = {}
    for gene, allele, db_name in rows:
        for component in gene.split(","):
            by_selection.setdefault((component, allele), set()).add(db_name)

    ambiguous = sorted(pair for pair, names in by_selection.items() if len(names) > 1)

    if not ambiguous:
        return []

    return [str(len(ambiguous)) + " gene/allele selection(s) resolve to more than one allele "
            "name, so which one a plot shows is arbitrary: "
            + sample_names([(g + "," + a,) for g, a in ambiguous])]

def validate_loaded_data():
    """Check the loaded data, raising SourceDataError if anything is wrong.

    Every check runs before anything is raised, so one crash reports everything that needs
    fixing rather than one problem per restart - but each problem is reported once. An
    unrecognised locus makes its rows look untranslated too, so reporting the locus and
    stopping there says one true thing rather than two overlapping ones about the same row.
    """
    unknown_loci = unrecognised_loci()

    problems = (locus_problems(unknown_loci)
                + mixed_locus_problems()
                + amino_acid_coverage_problems(unknown_loci)
                + igsnper_consistency_problems()
                + amino_acid_list_problems()
                + population_problems()
                + allele_resolution_problems())

    if problems:
        raise SourceDataError(
            "Loaded data breaks " + str(len(problems)) + " assumption(s) this app relies "
            "on. These need reporting to the research group so the source data can be "
            "corrected:\n  - " + "\n  - ".join(problems))

    print("Source data validated: " + str(
        ImmuneDiscoverDataModel.query.with_entities(ImmuneDiscoverDataModel.id).count()
    ) + " rows.", flush = True)
