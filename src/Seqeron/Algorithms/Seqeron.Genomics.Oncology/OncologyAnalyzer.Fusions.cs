namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region Fusion Gene Detection (ONCO-FUSION-001)

    /// <summary>
    /// Default minimum number of junction-spanning (split) reads required to support a fusion.
    /// Source: STAR-Fusion source (Haas et al. 2017), default <c>my $MIN_JUNCTION_READS = 1;</c>
    /// https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
    /// </summary>
    public const int DefaultMinJunctionReads = 1;

    /// <summary>
    /// Default minimum total fusion support = (junction reads + spanning/discordant fragments), applied
    /// when at least one junction read is present.
    /// Source: STAR-Fusion source, default <c>my $MIN_SUM_FRAGS = 2;</c> ("requires at least one junction
    /// read"). https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
    /// </summary>
    public const int DefaultMinSumFrags = 2;

    /// <summary>
    /// Default minimum number of spanning (discordant) fragments required when there are NO junction reads.
    /// Source: STAR-Fusion source, default <c>my $MIN_SPANNING_FRAGS_ONLY = 5;</c>
    /// https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
    /// </summary>
    public const int DefaultMinSpanningFragsOnly = 5;

    /// <summary>Number of nucleotides per codon (reading frame is read in triplets).
    /// Source: Wikipedia "Reading frame" (Badger &amp; Olsen 1999; Lodish 6th ed.):
    /// a reading frame reads nucleotides "as a sequence of triplets". https://en.wikipedia.org/wiki/Reading_frame</summary>
    private const int CodonLength = 3;

    /// <summary>
    /// A candidate gene-fusion breakpoint with its per-class supporting-read counts.
    /// Read classes follow the Arriba output schema (Uhrig et al. 2021): split reads anchored in each
    /// partner and discordant (spanning) mate-pairs.
    /// </summary>
    /// <param name="Gene5Prime">5' (upstream) fusion partner gene symbol.</param>
    /// <param name="Gene3Prime">3' (downstream) fusion partner gene symbol.</param>
    /// <param name="SplitReads5Prime">Split (junction) reads whose longer segment (anchor) maps to the 5' gene (Arriba split_reads1).</param>
    /// <param name="SplitReads3Prime">Split (junction) reads whose longer segment (anchor) maps to the 3' gene (Arriba split_reads2).</param>
    /// <param name="DiscordantMates">Discordant (spanning/bridge) mate-pairs supporting the fusion (Arriba discordant_mates).</param>
    /// <param name="FivePrimeCodingBases">Coding bases the 5' partner contributes upstream of the breakpoint (for the reading-frame test); -1 if unknown.</param>
    /// <param name="ThreePrimeStartPhase">Coding-start phase (0, 1, or 2) of the 3' partner at the breakpoint; -1 if unknown.</param>
    public readonly record struct FusionCandidate(
        string Gene5Prime,
        string Gene3Prime,
        int SplitReads5Prime,
        int SplitReads3Prime,
        int DiscordantMates,
        int FivePrimeCodingBases = -1,
        int ThreePrimeStartPhase = -1);

    /// <summary>
    /// Minimum-support thresholds for fusion calling. Defaults mirror STAR-Fusion (Haas et al.).
    /// </summary>
    /// <param name="MinJunctionReads">Minimum junction (split) reads required (STAR-Fusion min_junction_reads).</param>
    /// <param name="MinSumFrags">Minimum total support when ≥1 junction read present (STAR-Fusion min_sum_frags).</param>
    /// <param name="MinSpanningFragsOnly">Minimum discordant fragments required when there are no junction reads (STAR-Fusion min_spanning_frags_only).</param>
    // S3427: the explicit parameterless ctor deliberately overlaps the positional ctor's all-defaults
    // call — that overlap is the fix (it restores the STAR-Fusion defaults that `new()`/default(T) would
    // otherwise zero out), not an accident. Suppressed rather than removed.
#pragma warning disable S3427
    public readonly record struct FusionDetectionThresholds(
        int MinJunctionReads = DefaultMinJunctionReads,
        int MinSumFrags = DefaultMinSumFrags,
        int MinSpanningFragsOnly = DefaultMinSpanningFragsOnly)
    {
        // A record struct's *implicit* parameterless constructor (e.g. default(T) / new()) bypasses the
        // positional defaults above and zero-initialises every field, which would silently disable the
        // STAR-Fusion thresholds. Declaring an explicit parameterless constructor restores the defaults.
        public FusionDetectionThresholds()
            : this(DefaultMinJunctionReads, DefaultMinSumFrags, DefaultMinSpanningFragsOnly)
        {
        }
    }
#pragma warning restore S3427

    /// <summary>Reading-frame status of the 3' partner across the fusion junction.</summary>
    public enum FusionReadingFrame
    {
        /// <summary>Reading frame is preserved across the junction (3' partner stays in phase).</summary>
        InFrame,

        /// <summary>Reading frame is shifted across the junction (3' partner translated out of phase).</summary>
        OutOfFrame,

        /// <summary>Coding-phase information was not supplied, so frame could not be determined.</summary>
        Unknown
    }

    /// <summary>A detected (passing) gene fusion with its supporting evidence.</summary>
    /// <param name="Gene5Prime">5' fusion partner.</param>
    /// <param name="Gene3Prime">3' fusion partner.</param>
    /// <param name="JunctionReads">Total junction (split) reads = SplitReads5Prime + SplitReads3Prime.</param>
    /// <param name="DiscordantMates">Discordant (spanning) mate-pairs.</param>
    /// <param name="TotalSupport">Total supporting reads = SplitReads5Prime + SplitReads3Prime + DiscordantMates (Arriba).</param>
    /// <param name="ReadingFrame">In-frame / out-of-frame / unknown status of the junction.</param>
    public readonly record struct FusionCall(
        string Gene5Prime,
        string Gene3Prime,
        int JunctionReads,
        int DiscordantMates,
        int TotalSupport,
        FusionReadingFrame ReadingFrame);

    /// <summary>
    /// Total supporting reads for a candidate = split_reads1 + split_reads2 + discordant_mates.
    /// Source: Arriba output spec — "The total number of supporting reads can be obtained by summing up the
    /// reads given in the columns split_reads1, split_reads2, discordant_mates".
    /// https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    public static int ComputeTotalSupport(FusionCandidate candidate)
        => candidate.SplitReads5Prime + candidate.SplitReads3Prime + candidate.DiscordantMates;

    /// <summary>
    /// Determines whether the 3' partner is fused in-frame using codon phase.
    /// A fusion is in-frame iff the coding bases the 5' partner contributes upstream of the breakpoint,
    /// taken modulo 3 relative to the 3' partner's coding-start phase, keep the downstream codons in phase:
    /// <c>(fivePrimeCodingBases - threePrimeStartPhase) mod 3 == 0</c>.
    /// Source: Genomics England exon-phase rule ("if one exon finishes after the second letter of a triplet
    /// (end phase 2), the next one should start with the third letter") +
    /// reading frames read in triplets / modulo 3 (Wikipedia "Reading frame", Badger &amp; Olsen 1999).
    /// https://www.genomicsengland.co.uk/blog/gene-fusion-reporting , https://en.wikipedia.org/wiki/Reading_frame
    /// </summary>
    /// <param name="fivePrimeCodingBases">Coding bases contributed by the 5' partner upstream of the breakpoint (≥ 0).</param>
    /// <param name="threePrimeStartPhase">Coding-start phase of the 3' partner at the breakpoint (0, 1, or 2).</param>
    /// <returns><see langword="true"/> if the junction preserves the reading frame.</returns>
    /// <exception cref="ArgumentOutOfRangeException">A negative base count, or a phase outside {0,1,2}.</exception>
    public static bool IsInFrame(int fivePrimeCodingBases, int threePrimeStartPhase)
    {
        if (fivePrimeCodingBases < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(fivePrimeCodingBases),
                "Coding-base count cannot be negative.");
        }

        if (threePrimeStartPhase < 0 || threePrimeStartPhase >= CodonLength)
        {
            throw new ArgumentOutOfRangeException(nameof(threePrimeStartPhase),
                "Coding-start phase must be 0, 1, or 2.");
        }

        // Modulo arithmetic over a non-negative dividend; result is in [0, CodonLength).
        return (fivePrimeCodingBases - threePrimeStartPhase) % CodonLength == 0;
    }

    /// <summary>
    /// Detects candidate gene fusions from breakpoint supporting-read counts using the STAR-Fusion
    /// minimum-support rule, and reports each passing fusion with its total support and reading-frame status.
    /// </summary>
    /// <remarks>
    /// A candidate is reported as a fusion iff:
    /// <list type="bullet">
    /// <item><description>gene5p ≠ gene3p (a gene is not fused with itself — Registry invariant); and</description></item>
    /// <item><description>when junction (split) reads ≥ 1: junctionReads ≥ MinJunctionReads AND totalSupport ≥ MinSumFrags
    /// (STAR-Fusion min_junction_reads / min_sum_frags); OR</description></item>
    /// <item><description>when there are no junction reads: discordantMates ≥ MinSpanningFragsOnly
    /// (STAR-Fusion min_spanning_frags_only).</description></item>
    /// </list>
    /// Results are returned ordered by descending total support (STAR-Fusion scores fusions by the
    /// abundance of supporting reads), with the gene pair as a deterministic tie-breaker.
    /// </remarks>
    /// <param name="candidates">Candidate breakpoints with per-class supporting-read counts.</param>
    /// <param name="thresholds">Optional support thresholds; defaults to STAR-Fusion defaults.</param>
    /// <returns>Detected fusions ordered by descending total support.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="candidates"/> is null.</exception>
    /// <exception cref="ArgumentException">A candidate has a negative supporting-read count.</exception>
    public static IReadOnlyList<FusionCall> DetectFusions(
        IEnumerable<FusionCandidate> candidates,
        FusionDetectionThresholds? thresholds = null)
    {
        ArgumentNullException.ThrowIfNull(candidates);

        FusionDetectionThresholds t = thresholds ?? new FusionDetectionThresholds();
        var calls = new List<FusionCall>();

        foreach (FusionCandidate candidate in candidates)
        {
            if (candidate.SplitReads5Prime < 0 || candidate.SplitReads3Prime < 0 || candidate.DiscordantMates < 0)
            {
                throw new ArgumentException(
                    "Supporting-read counts cannot be negative.", nameof(candidates));
            }

            // INV-1: a gene is not fused with itself (case-insensitive symbol comparison).
            if (string.Equals(candidate.Gene5Prime, candidate.Gene3Prime, StringComparison.OrdinalIgnoreCase))
            {
                continue;
            }

            int junctionReads = candidate.SplitReads5Prime + candidate.SplitReads3Prime;
            int totalSupport = ComputeTotalSupport(candidate);

            bool passes = junctionReads >= t.MinJunctionReads
                ? totalSupport >= t.MinSumFrags
                // No junction reads: fall back to the spanning-fragments-only rule.
                : candidate.DiscordantMates >= t.MinSpanningFragsOnly;

            if (!passes)
            {
                continue;
            }

            FusionReadingFrame frame = ResolveReadingFrame(candidate);
            calls.Add(new FusionCall(
                candidate.Gene5Prime,
                candidate.Gene3Prime,
                junctionReads,
                candidate.DiscordantMates,
                totalSupport,
                frame));
        }

        // INV-4: order by descending total support (abundance of supporting reads), deterministic ties.
        return calls
            .OrderByDescending(c => c.TotalSupport)
            .ThenBy(c => c.Gene5Prime, StringComparer.Ordinal)
            .ThenBy(c => c.Gene3Prime, StringComparer.Ordinal)
            .ToArray();
    }

    private static FusionReadingFrame ResolveReadingFrame(FusionCandidate candidate)
    {
        if (candidate.FivePrimeCodingBases < 0 || candidate.ThreePrimeStartPhase < 0
            || candidate.ThreePrimeStartPhase >= CodonLength)
        {
            return FusionReadingFrame.Unknown;
        }

        return IsInFrame(candidate.FivePrimeCodingBases, candidate.ThreePrimeStartPhase)
            ? FusionReadingFrame.InFrame
            : FusionReadingFrame.OutOfFrame;
    }

    #endregion


    #region Known Fusion Database Lookup (ONCO-FUSION-002)

    /// <summary>
    /// HGNC fusion-designation separator: a double colon between the 5' and 3' partner symbols.
    /// Source: Bruford et al. (2021), HGNC recommendations for the designation of gene fusions —
    /// "HGNC recommends that a new separator—a double colon (::)—be used in describing gene fusions,
    /// e.g., BCR::ABL1." https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/
    /// </summary>
    public const string FusionDesignationSeparator = "::";

    /// <summary>
    /// The result of looking a fusion up against a known-fusion set.
    /// </summary>
    /// <param name="Designation">The HGNC <c>5'::3'</c> designation of the queried fusion (e.g. <c>BCR::ABL1</c>).</param>
    /// <param name="IsKnown"><see langword="true"/> if the directional designation was present in the supplied set.</param>
    /// <param name="Annotation">The caller-supplied annotation for the matched designation, or <see langword="null"/> if not known.</param>
    public readonly record struct KnownFusionMatch(string Designation, bool IsKnown, string? Annotation);

    /// <summary>
    /// Formats the HGNC designation of a gene fusion as <c>gene5p::gene3p</c>.
    /// The 5' partner is always written first, before the double colon, irrespective of chromosomal
    /// location or gene orientation; the designation is therefore directional (A::B ≠ B::A).
    /// Source: Bruford et al. (2021), HGNC recommendations for the designation of gene fusions —
    /// "the 5′ partner gene should always be listed first in the description of a fusion gene, i.e.,
    /// before the double colon" and "a double colon (::) … e.g., BCR::ABL1".
    /// https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/
    /// </summary>
    /// <param name="gene5p">5' (upstream) partner gene symbol; must be non-empty.</param>
    /// <param name="gene3p">3' (downstream) partner gene symbol; must be non-empty.</param>
    /// <returns>The designation string <c>gene5p + "::" + gene3p</c>.</returns>
    /// <exception cref="ArgumentException">Either symbol is null, empty, or whitespace.</exception>
    public static string GetFusionAnnotation(string gene5p, string gene3p)
    {
        if (string.IsNullOrWhiteSpace(gene5p))
        {
            throw new ArgumentException("5' partner gene symbol must be non-empty.", nameof(gene5p));
        }

        if (string.IsNullOrWhiteSpace(gene3p))
        {
            throw new ArgumentException("3' partner gene symbol must be non-empty.", nameof(gene3p));
        }

        return gene5p + FusionDesignationSeparator + gene3p;
    }

    /// <summary>
    /// Looks a detected fusion up against a caller-supplied set of known fusions, keyed by their
    /// HGNC <c>5'::3'</c> designation.
    /// </summary>
    /// <remarks>
    /// The lookup is <b>directional</b>: the key is built with the 5' partner first (per Bruford et al. 2021),
    /// so a reciprocal fusion (partners swapped) is a different designation and does NOT match.
    /// Symbol comparison is case-insensitive (ordinal-ignore-case); the known-fusion set membership and the
    /// annotation text are entirely caller-supplied — this library bundles no curated fusion database
    /// (Mitelman / COSMIC / ChimerDB content is the caller's responsibility).
    /// Source (designation format and directional keying): Bruford et al. (2021),
    /// https://pmc.ncbi.nlm.nih.gov/articles/PMC8550944/
    /// </remarks>
    /// <param name="fusion">The fusion to look up (its 5'/3' partners define the key).</param>
    /// <param name="knownFusions">
    /// Caller-supplied map from <c>5'::3'</c> designation to its annotation. For case-insensitive matching,
    /// supply a dictionary built with <see cref="StringComparer.OrdinalIgnoreCase"/> (the method also probes
    /// case-insensitively when the dictionary is not already case-insensitive).
    /// </param>
    /// <returns>A <see cref="KnownFusionMatch"/> reporting the designation and, if present, its annotation.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="knownFusions"/> is null.</exception>
    /// <exception cref="ArgumentException">A fusion partner symbol is null, empty, or whitespace.</exception>
    public static KnownFusionMatch MatchKnownFusions(
        FusionCall fusion,
        IReadOnlyDictionary<string, string> knownFusions)
    {
        ArgumentNullException.ThrowIfNull(knownFusions);

        string designation = GetFusionAnnotation(fusion.Gene5Prime, fusion.Gene3Prime);

        // Directional key (5'::3'). Try the supplied dictionary's own comparer first; if it is not
        // case-insensitive, fall back to an explicit case-insensitive scan so that e.g. "eml4::alk"
        // matches a stored "EML4::ALK" (HGNC symbols are case-defined, but inputs vary in case).
        if (knownFusions.TryGetValue(designation, out string? annotation))
        {
            return new KnownFusionMatch(designation, IsKnown: true, annotation);
        }

        foreach (KeyValuePair<string, string> entry in knownFusions)
        {
            if (string.Equals(entry.Key, designation, StringComparison.OrdinalIgnoreCase))
            {
                return new KnownFusionMatch(designation, IsKnown: true, entry.Value);
            }
        }

        return new KnownFusionMatch(designation, IsKnown: false, Annotation: null);
    }

    #endregion


    #region Fusion Breakpoint Analysis (ONCO-FUSION-003)

    /// <summary>
    /// Location category of a fusion breakpoint within a partner transcript, mirroring the Arriba
    /// <c>site1</c>/<c>site2</c> output column. Source: Arriba output spec — "Possible values are: 5'UTR,
    /// 3'UTR, UTR, CDS, exon, intron, and intergenic". https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    public enum BreakpointSite
    {
        /// <summary>Coding sequence (the only site at which two reading frames are joined).</summary>
        Cds,

        /// <summary>5' untranslated region.</summary>
        FivePrimeUtr,

        /// <summary>3' untranslated region.</summary>
        ThreePrimeUtr,

        /// <summary>Untranslated region (strand/UTR side not resolved).</summary>
        Utr,

        /// <summary>Exon (non-coding part).</summary>
        Exon,

        /// <summary>Intron.</summary>
        Intron,

        /// <summary>Intergenic region.</summary>
        Intergenic
    }

    /// <summary>
    /// Reading-frame consequence reported by <see cref="AnalyzeBreakpoint"/>, mirroring the Arriba
    /// <c>reading_frame</c> column. Source: Arriba output spec — reading_frame ∈ {in-frame, out-of-frame,
    /// stop-codon, .}; the dot ("not predicted") is used when the peptide cannot be predicted because a
    /// breakpoint is not in coding context. https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    public enum BreakpointFrameStatus
    {
        /// <summary>Both breakpoints are in CDS and the junction preserves the reading frame (in-frame).</summary>
        InFrame,

        /// <summary>Both breakpoints are in CDS but the junction shifts the reading frame (out-of-frame).</summary>
        OutOfFrame,

        /// <summary>The chimeric ORF reaches a stop codon at/after the junction (Arriba <c>stop-codon</c>).</summary>
        StopCodon,

        /// <summary>Frame cannot be called because a breakpoint is not in coding context (Arriba <c>.</c>).</summary>
        NotPredicted
    }

    /// <summary>
    /// A fusion breakpoint described by its per-partner site categories and the coding-frame quantities
    /// needed to call the junction reading frame. The 5'/3' partner symbols carry over from a
    /// <see cref="FusionCall"/>; <see cref="FivePrimeCodingBases"/> and <see cref="ThreePrimeStartPhase"/>
    /// are the same coding-frame inputs used by <see cref="IsInFrame"/> (ONCO-FUSION-001). Read class
    /// schema follows Arriba (Uhrig et al. 2021). https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    /// <param name="Gene5Prime">5' (upstream) fusion partner gene symbol.</param>
    /// <param name="Gene3Prime">3' (downstream) fusion partner gene symbol.</param>
    /// <param name="Site5Prime">Site category of the breakpoint in the 5' partner.</param>
    /// <param name="Site3Prime">Site category of the breakpoint in the 3' partner.</param>
    /// <param name="FivePrimeCodingBases">Coding bases the 5' partner contributes upstream of the breakpoint (≥ 0).</param>
    /// <param name="ThreePrimeStartPhase">Coding-start phase (0, 1, or 2) of the 3' partner at the breakpoint.</param>
    public readonly record struct FusionBreakpoint(
        string Gene5Prime,
        string Gene3Prime,
        BreakpointSite Site5Prime,
        BreakpointSite Site3Prime,
        int FivePrimeCodingBases,
        int ThreePrimeStartPhase);

    /// <summary>The breakpoint-analysis result for a fusion junction.</summary>
    /// <param name="Gene5Prime">5' partner (carried through unchanged).</param>
    /// <param name="Gene3Prime">3' partner (carried through unchanged).</param>
    /// <param name="Site5Prime">Site category of the 5' breakpoint.</param>
    /// <param name="Site3Prime">Site category of the 3' breakpoint.</param>
    /// <param name="BreakpointInCoding">True iff both breakpoints lie in CDS (a coding-to-coding junction).</param>
    /// <param name="FrameStatus">Reading-frame consequence of the junction.</param>
    public readonly record struct BreakpointAnalysis(
        string Gene5Prime,
        string Gene3Prime,
        BreakpointSite Site5Prime,
        BreakpointSite Site3Prime,
        bool BreakpointInCoding,
        BreakpointFrameStatus FrameStatus);

    /// <summary>The predicted protein product of a fusion junction.</summary>
    /// <param name="ChimericCds">The chimeric coding sequence = 5' CDS prefix ++ 3' CDS suffix (uppercased DNA).</param>
    /// <param name="Peptide">The translated fusion peptide, truncated at the first stop codon.</param>
    /// <param name="Effect">Reading-frame effect of the junction (in-frame / out-of-frame).</param>
    /// <param name="HasPrematureStop">True iff a stop codon was reached before the end of the chimeric ORF.</param>
    public readonly record struct FusionProteinPrediction(
        string ChimericCds,
        string Peptide,
        BreakpointFrameStatus Effect,
        bool HasPrematureStop);

    /// <summary>
    /// Analyzes a fusion breakpoint: reports the per-partner site categories and calls the junction
    /// reading frame. A reading-frame call (<see cref="BreakpointFrameStatus.InFrame"/> /
    /// <see cref="BreakpointFrameStatus.OutOfFrame"/>) is made ONLY when both breakpoints lie in CDS; if
    /// either breakpoint is in a UTR, intron, exon (non-coding part), or intergenic region the frame is
    /// <see cref="BreakpointFrameStatus.NotPredicted"/> (Arriba <c>reading_frame = .</c>). The in-frame
    /// test reuses the codon-phase rule of <see cref="IsInFrame"/>:
    /// <c>(fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0</c>.
    /// Source: Arriba output spec (site / reading_frame), AGFusion frame rule (Murphy &amp; Elemento 2016).
    /// https://github.com/suhrig/arriba/wiki/05-Output-files
    /// </summary>
    /// <param name="fusion">The breakpoint to analyze.</param>
    /// <returns>The site categories and reading-frame consequence of the junction.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Both breakpoints are CDS but <see cref="FusionBreakpoint.FivePrimeCodingBases"/> is negative or
    /// <see cref="FusionBreakpoint.ThreePrimeStartPhase"/> is outside {0, 1, 2}.
    /// </exception>
    public static BreakpointAnalysis AnalyzeBreakpoint(FusionBreakpoint fusion)
    {
        bool bothCoding = fusion.Site5Prime == BreakpointSite.Cds
                       && fusion.Site3Prime == BreakpointSite.Cds;

        // A frame call is only defined for a coding-to-coding junction (Arriba reading_frame = '.' otherwise).
        BreakpointFrameStatus frame;
        if (!bothCoding)
            frame = BreakpointFrameStatus.NotPredicted;
        else if (IsInFrame(fusion.FivePrimeCodingBases, fusion.ThreePrimeStartPhase))
            frame = BreakpointFrameStatus.InFrame;
        else
            frame = BreakpointFrameStatus.OutOfFrame;

        return new BreakpointAnalysis(
            fusion.Gene5Prime,
            fusion.Gene3Prime,
            fusion.Site5Prime,
            fusion.Site3Prime,
            bothCoding,
            frame);
    }

    /// <summary>
    /// Predicts the protein product of a fusion from the two partners' coding sequences. The chimeric CDS
    /// is the 5' partner's CDS taken up to the breakpoint (a prefix) concatenated with the 3' partner's CDS
    /// taken from the breakpoint onward (a suffix); it is then translated with the standard genetic code
    /// (NCBI Table 1) and truncated at the first stop codon.
    /// Source (verbatim from AGFusion model.py, Murphy &amp; Elemento 2016):
    /// <c>cds_5prime = transcript1.coding_sequence[0:junction5]</c>,
    /// <c>cds_3prime = transcript2.coding_sequence[junction3:]</c>,
    /// <c>seq = cds_5prime + cds_3prime</c>,
    /// <c>protein_seq = cds.seq.translate(); protein_seq = protein_seq[0:protein_seq.find("*")]</c>.
    /// When the junction is out-of-frame the chimeric CDS is first trimmed to a whole number of codons
    /// (<c>cds[0:3*(len//3)]</c>) so the 3' partner is read in its (shifted) frame.
    /// https://raw.githubusercontent.com/murphycj/AGFusion/master/agfusion/model.py
    /// </summary>
    /// <param name="fusion">The breakpoint (its site categories and codon-frame quantities drive the effect call).</param>
    /// <param name="transcripts">
    /// The two partner CDS sequences as (fivePrimeCds, threePrimeCds): the 5' partner's full coding sequence
    /// and the 3' partner's full coding sequence (DNA, A/C/G/T). The breakpoint offsets are taken from
    /// <paramref name="fusion"/>: the 5' prefix length is <see cref="FusionBreakpoint.FivePrimeCodingBases"/>
    /// and the 3' suffix starts at <see cref="FusionBreakpoint.ThreePrimeStartPhase"/>.
    /// </param>
    /// <returns>The chimeric CDS, the translated (first-stop-truncated) peptide, the frame effect, and a premature-stop flag.</returns>
    /// <exception cref="ArgumentNullException">A CDS sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// An offset is out of range for its CDS, or the 3' start phase is outside {0, 1, 2}.
    /// </exception>
    public static FusionProteinPrediction PredictFusionProtein(
        FusionBreakpoint fusion,
        (string FivePrimeCds, string ThreePrimeCds) transcripts)
    {
        ArgumentNullException.ThrowIfNull(transcripts.FivePrimeCds);
        ArgumentNullException.ThrowIfNull(transcripts.ThreePrimeCds);

        string fivePrimeCds = transcripts.FivePrimeCds.ToUpperInvariant();
        string threePrimeCds = transcripts.ThreePrimeCds.ToUpperInvariant();

        int junction5 = fusion.FivePrimeCodingBases;
        int junction3 = fusion.ThreePrimeStartPhase;

        if (junction5 < 0 || junction5 > fivePrimeCds.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(fusion),
                "FivePrimeCodingBases (5' prefix length) is out of range for the 5' CDS.");
        }

        if (junction3 < 0 || junction3 > threePrimeCds.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(fusion),
                "ThreePrimeStartPhase (3' suffix start) is out of range for the 3' CDS.");
        }

        // AGFusion: cds_5prime = transcript1.coding_sequence[0:junction5]; cds_3prime = transcript2[junction3:].
        string cds5 = fivePrimeCds.Substring(0, junction5);
        string cds3 = threePrimeCds.Substring(junction3);
        string chimericCds = cds5 + cds3;

        // Frame effect by the AGFusion / IsInFrame codon-phase rule (only meaningful junction3 in {0,1,2}).
        bool inFrame = junction3 < CodonLength
            && IsInFrame(junction5, junction3);
        BreakpointFrameStatus effect = inFrame
            ? BreakpointFrameStatus.InFrame
            : BreakpointFrameStatus.OutOfFrame;

        // AGFusion: an out-of-frame CDS is trimmed to whole codons before translation; an in-frame CDS is
        // translated as-is (a trailing partial codon, if any, is not translatable and is dropped).
        int translatableLength = chimericCds.Length - (chimericCds.Length % CodonLength);
        var peptide = new System.Text.StringBuilder(translatableLength / CodonLength);
        bool hasStop = false;

        for (int i = 0; i < translatableLength; i += CodonLength)
        {
            string codon = chimericCds.Substring(i, CodonLength);
            char aminoAcid = GeneticCode.Standard.Translate(codon);
            if (aminoAcid == StopCodonSymbol)
            {
                // AGFusion truncates the peptide at the first stop codon (protein_seq[0:find("*")]).
                hasStop = true;
                break;
            }

            peptide.Append(aminoAcid);
        }

        return new FusionProteinPrediction(chimericCds, peptide.ToString(), effect, hasStop);
    }

    /// <summary>
    /// Stop-codon marker returned by <see cref="GeneticCode.Translate"/> for UAA/UAG/UGA (NCBI Table 1).
    /// Source: AGFusion truncation at the first <c>'*'</c>; <see cref="GeneticCode"/> emits '*' for stops.
    /// </summary>
    private const char StopCodonSymbol = '*';

    #endregion

}
