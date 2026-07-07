namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region Driver Mutation Detection (20/20 rule)

    /// <summary>
    /// Fraction-of-mutations threshold of the Vogelstein 20/20 rule: a gene is classified as a driver
    /// only when more than this fraction of its mutations meet the oncogene or tumor-suppressor criterion.
    /// Source: Vogelstein B et al. (2013), Science 339(6127):1546–1558 — "more than 20%"; restated verbatim
    /// by Tokheim &amp; Karchin (2020), Bioinformatics 36(6):1712–1719 ("OGs have &gt;20% ... TSGs have &gt;20% ...").
    /// The comparison is strict (&gt;), so an exact 20% fraction is not sufficient. Value = 0.20.
    /// </summary>
    public const double DriverGeneFractionThreshold = 0.20;

    /// <summary>
    /// Minimum number of mutations at the same protein position for that position to count as
    /// <i>recurrent</i> (a hotspot). Source: Miller ML et al. (2017), Oncotarget 8(20):33321–33333 — a
    /// recurrent position requires "at least two mutations of the same class" at an identical location.
    /// Value = 2.
    /// </summary>
    public const int RecurrentPositionMinCount = 2;

    /// <summary>
    /// The functional consequence of a coding mutation, restricted to the categories the 20/20 rule
    /// distinguishes. Truncating categories are those listed by Schroeder MP et al. (2014),
    /// Bioinformatics 30(17):i549–i555 and Miller ML et al. (2017): nonsense (stop gain/loss), frameshift
    /// indels, and splice donor/acceptor mutations. Missense at recurrent positions drives the oncogene call.
    /// </summary>
    public enum MutationConsequence
    {
        /// <summary>Amino-acid-changing substitution (drives the oncogene criterion when recurrent).</summary>
        Missense,

        /// <summary>Premature/lost stop codon (nonsense) — truncating/inactivating.</summary>
        Nonsense,

        /// <summary>Insertion/deletion shifting the reading frame — truncating/inactivating.</summary>
        Frameshift,

        /// <summary>Mutation at a splice donor/acceptor site — truncating/inactivating.</summary>
        SpliceSite,

        /// <summary>Synonymous or other non-truncating, non-missense change (counts toward the denominator only).</summary>
        Other
    }

    /// <summary>The 20/20-rule role assigned to a gene from its mutation spectrum.</summary>
    public enum DriverGeneRole
    {
        /// <summary>&gt;20% of mutations are missense at recurrent positions (Vogelstein 2013).</summary>
        Oncogene,

        /// <summary>&gt;20% of mutations are truncating/inactivating (Vogelstein 2013).</summary>
        TumorSuppressor,

        /// <summary>Neither criterion exceeds 20% (or an exact tie): not classified as a driver gene.</summary>
        Ambiguous
    }

    /// <summary>
    /// A single coding mutation observed in a gene, reduced to the features the 20/20 rule needs.
    /// </summary>
    /// <param name="Gene">Gene symbol the mutation falls in.</param>
    /// <param name="ProteinPosition">1-based codon / amino-acid position used to detect recurrence.</param>
    /// <param name="Consequence">Functional consequence category.</param>
    public readonly record struct GeneMutation(string Gene, int ProteinPosition, MutationConsequence Consequence);

    /// <summary>
    /// The 20/20-rule classification of one gene, with the two criterion fractions that produced it.
    /// </summary>
    /// <param name="Gene">Gene symbol.</param>
    /// <param name="Role">Assigned <see cref="DriverGeneRole"/>.</param>
    /// <param name="TruncatingFraction">Fraction of mutations that are truncating/inactivating, in [0, 1].</param>
    /// <param name="RecurrentMissenseFraction">Fraction of mutations that are missense at a recurrent position, in [0, 1].</param>
    /// <param name="MutationCount">Total number of mutations considered (the denominator).</param>
    public readonly record struct DriverGeneClassification(
        string Gene,
        DriverGeneRole Role,
        double TruncatingFraction,
        double RecurrentMissenseFraction,
        int MutationCount);

    /// <summary>
    /// Classifies a single gene by the Vogelstein 20/20 rule from its observed coding mutations. The gene is
    /// an <see cref="DriverGeneRole.Oncogene"/> when more than 20% of its mutations are missense at recurrent
    /// positions (a position observed ≥ <see cref="RecurrentPositionMinCount"/> times), and a
    /// <see cref="DriverGeneRole.TumorSuppressor"/> when more than 20% of its mutations are truncating
    /// (nonsense, frameshift, or splice-site). If both criteria are met the dominant fraction decides; an exact
    /// tie or neither criterion met yields <see cref="DriverGeneRole.Ambiguous"/>.
    /// Source: Vogelstein et al. (2013); Tokheim &amp; Karchin (2020); Schroeder et al. (2014); Miller et al. (2017).
    /// </summary>
    /// <param name="mutations">All coding mutations observed in one gene.</param>
    /// <returns>The gene's 20/20-rule classification with its criterion fractions.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="mutations"/> is null.</exception>
    public static DriverGeneClassification ClassifyGene(IEnumerable<GeneMutation> mutations)
    {
        ArgumentNullException.ThrowIfNull(mutations);

        var list = mutations as IReadOnlyList<GeneMutation> ?? mutations.ToList();
        int total = list.Count;
        string gene = total > 0 ? list[0].Gene : string.Empty;

        if (total == 0)
        {
            return new DriverGeneClassification(gene, DriverGeneRole.Ambiguous, 0.0, 0.0, 0);
        }

        int truncating = list.Count(IsTruncating);
        int recurrentMissense = CountRecurrentMissense(list);

        double truncatingFraction = (double)truncating / total;
        double recurrentMissenseFraction = (double)recurrentMissense / total;

        bool isTsg = truncatingFraction > DriverGeneFractionThreshold;
        bool isOg = recurrentMissenseFraction > DriverGeneFractionThreshold;

        DriverGeneRole role;
        if (isTsg && isOg)
        {
            // Both criteria pass (atypical per Vogelstein 2013 — well-documented genes far surpass one
            // criterion). Resolve by the dominant signal; an exact tie is genuinely ambiguous.
            role = truncatingFraction > recurrentMissenseFraction ? DriverGeneRole.TumorSuppressor
                 : recurrentMissenseFraction > truncatingFraction ? DriverGeneRole.Oncogene
                 : DriverGeneRole.Ambiguous;
        }
        else if (isTsg)
        {
            role = DriverGeneRole.TumorSuppressor;
        }
        else if (isOg)
        {
            role = DriverGeneRole.Oncogene;
        }
        else
        {
            role = DriverGeneRole.Ambiguous;
        }

        return new DriverGeneClassification(gene, role, truncatingFraction, recurrentMissenseFraction, total);
    }

    /// <summary>
    /// Computes the 20/20-rule driver-signal score for a gene: the larger of its truncating fraction and its
    /// recurrent-missense fraction, in [0, 1]. This is the transparent, source-derived strength of the driver
    /// signal underlying <see cref="ClassifyGene"/>; it is NOT an external pathogenicity model (CADD/SIFT/
    /// PolyPhen), which are caller-supplied / not implemented. Source: Vogelstein et al. (2013) 20/20 rule.
    /// </summary>
    /// <param name="mutations">All coding mutations observed in one gene.</param>
    /// <returns>max(truncating fraction, recurrent-missense fraction), in [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="mutations"/> is null.</exception>
    public static double ScoreDriverPotential(IEnumerable<GeneMutation> mutations)
    {
        var c = ClassifyGene(mutations);
        return Math.Max(c.TruncatingFraction, c.RecurrentMissenseFraction);
    }

    /// <summary>
    /// Tests whether a mutation's (gene, protein position) is present in a caller-supplied set of known
    /// cancer hotspot positions. The 20/20 rule treats recurrent positions as the activating signal of
    /// oncogenes (Miller et al. 2017); curated hotspot catalogs (COSMIC, Cancer Hotspots, OncoKB) are
    /// supplied by the caller rather than hardcoded, since they cannot be reproduced authoritatively here.
    /// </summary>
    /// <param name="mutation">The mutation to test.</param>
    /// <param name="knownHotspots">Set of known hotspots as (gene, protein position) pairs.</param>
    /// <returns><c>true</c> when (gene, position) is in <paramref name="knownHotspots"/>; otherwise <c>false</c>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="knownHotspots"/> is null.</exception>
    public static bool MatchCancerHotspots(
        GeneMutation mutation,
        IReadOnlySet<(string Gene, int ProteinPosition)> knownHotspots)
    {
        ArgumentNullException.ThrowIfNull(knownHotspots);
        return knownHotspots.Contains((mutation.Gene, mutation.ProteinPosition));
    }

    /// <summary>
    /// Identifies driver mutations across a set of somatic coding mutations by applying the 20/20 rule
    /// per gene: a mutation is a driver if its gene classifies as an <see cref="DriverGeneRole.Oncogene"/>
    /// or <see cref="DriverGeneRole.TumorSuppressor"/>, OR if its (gene, position) is a known hotspot in
    /// <paramref name="knownHotspots"/>. The returned mutations are always a subset of the input, in input
    /// order (invariant: driver_mutations ⊆ somatic_mutations). Source: Vogelstein et al. (2013); Miller
    /// et al. (2017).
    /// </summary>
    /// <param name="mutations">Somatic coding mutations (one entry per observed mutation).</param>
    /// <param name="knownHotspots">Optional caller-supplied hotspot set; null is treated as empty.</param>
    /// <returns>The subset of input mutations that fall in a driver gene or at a known hotspot.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="mutations"/> is null.</exception>
    public static IReadOnlyList<GeneMutation> IdentifyDriverMutations(
        IEnumerable<GeneMutation> mutations,
        IReadOnlySet<(string Gene, int ProteinPosition)>? knownHotspots = null)
    {
        ArgumentNullException.ThrowIfNull(mutations);

        var list = mutations as IReadOnlyList<GeneMutation> ?? mutations.ToList();
        var hotspots = knownHotspots ?? EmptyHotspots;

        // Classify each gene once from its full mutation spectrum.
        var driverGenes = new HashSet<string>();
        foreach (var byGene in list.GroupBy(m => m.Gene, StringComparer.Ordinal))
        {
            if (ClassifyGene(byGene).Role != DriverGeneRole.Ambiguous)
            {
                driverGenes.Add(byGene.Key);
            }
        }

        var drivers = new List<GeneMutation>();
        foreach (var mutation in list)
        {
            if (driverGenes.Contains(mutation.Gene) || hotspots.Contains((mutation.Gene, mutation.ProteinPosition)))
            {
                drivers.Add(mutation);
            }
        }

        return drivers;
    }

    private static readonly IReadOnlySet<(string Gene, int ProteinPosition)> EmptyHotspots =
        new HashSet<(string, int)>();

    private static bool IsTruncating(GeneMutation mutation) =>
        mutation.Consequence is MutationConsequence.Nonsense
            or MutationConsequence.Frameshift
            or MutationConsequence.SpliceSite;

    /// <summary>
    /// Counts mutations that are missense AND located at a recurrent position (a protein position carrying
    /// ≥ <see cref="RecurrentPositionMinCount"/> missense mutations), per Miller et al. (2017).
    /// </summary>
    private static int CountRecurrentMissense(IReadOnlyList<GeneMutation> mutations)
    {
        var missenseByPosition = new Dictionary<int, int>();
        foreach (var mutation in mutations)
        {
            if (mutation.Consequence == MutationConsequence.Missense)
            {
                missenseByPosition.TryGetValue(mutation.ProteinPosition, out int count);
                missenseByPosition[mutation.ProteinPosition] = count + 1;
            }
        }

        int recurrent = 0;
        foreach (var positionCount in missenseByPosition.Values)
        {
            if (positionCount >= RecurrentPositionMinCount)
            {
                recurrent += positionCount;
            }
        }

        return recurrent;
    }

    #endregion


    #region Sequencing Artifact Detection

    /// <summary>
    /// GIV (Global Imbalance Value) threshold above which a library is declared damaged for a given
    /// substitution type. Source: Chen L. et al. (2017), Science 355(6326):752–756, as summarized in
    /// Nature Methods (2017) 14:330 ("DNA variants or DNA damage?"): "A GIV score of 1 indicates there
    /// is no DNA damage and a GIV score above 1.5 is defined as damaged DNA." Value = 1.5.
    /// </summary>
    public const double DamagedGivThreshold = 1.5;

    /// <summary>
    /// GIV value of a perfectly balanced (undamaged) library: equal G&gt;T counts in read 1 and read 2.
    /// Source: Chen et al. (2017) / Nature Methods (2017) — GIV = 1 means no DNA damage. Value = 1.0.
    /// </summary>
    public const double UndamagedGivScore = 1.0;

    /// <summary>
    /// Minimum two-sided Fisher exact p-value used before Phred-scaling the strand-bias score, mirroring
    /// GATK's <c>FisherStrand.MIN_PVALUE</c>. Source: Broad Institute GATK, FisherStrand.java
    /// (<c>static final double MIN_PVALUE = 1E-320;</c>). Caps FS so a p-value of 0 does not produce
    /// an infinite Phred score.
    /// </summary>
    private const double MinFisherPValue = 1E-320;

    /// <summary>
    /// Classes of sequencing artifact distinguished by substitution type (and, for OxoG, read-orientation
    /// imbalance). The two artifact classes are disjoint by substitution: deamination is C:G&gt;T:A
    /// (C&gt;T / G&gt;A), oxidation is G:C&gt;T:A read as G&gt;T / C&gt;A
    /// (Do &amp; Dobrovic 2015; Chen et al. 2017).
    /// </summary>
    public enum ArtifactType
    {
        /// <summary>Not a recognized substitution-class artifact (a candidate true variant).</summary>
        None,

        /// <summary>
        /// FFPE cytosine-deamination artifact: C&gt;T or G&gt;A (collectively C:G&gt;T:A). Deaminated
        /// cytosine becomes uracil, which pairs with adenine. Source: Do &amp; Dobrovic (2015).
        /// </summary>
        FfpeDeamination,

        /// <summary>
        /// OxoG (8-oxoguanine) oxidative artifact: G&gt;T (read 1) or C&gt;A (read 2, reverse complement).
        /// Source: Chen et al. (2017).
        /// </summary>
        OxoG
    }

    /// <summary>
    /// One observed candidate variant together with the strand- and read-orientation read evidence needed
    /// for artifact classification. The strand counts feed the GATK FisherStrand 2×2 contingency table; the
    /// read-mate counts (<paramref name="AltReadsR1"/> / <paramref name="AltReadsR2"/>) feed the OxoG GIV
    /// imbalance. (The repository has no BAM reader; these counts are supplied directly rather than parsed
    /// from a BAM file — an API-shape decision that does not change the classification rules.)
    /// </summary>
    /// <param name="ReferenceAllele">Single reference base (A/C/G/T).</param>
    /// <param name="AlternateAllele">Single alternate base (A/C/G/T).</param>
    /// <param name="RefForward">Reference-supporting reads on the forward strand.</param>
    /// <param name="RefReverse">Reference-supporting reads on the reverse strand.</param>
    /// <param name="AltForward">Alternate-supporting reads on the forward strand.</param>
    /// <param name="AltReverse">Alternate-supporting reads on the reverse strand.</param>
    /// <param name="AltReadsR1">Alternate-supporting reads from read 1 of the pair (for GIV).</param>
    /// <param name="AltReadsR2">Alternate-supporting reads from read 2 of the pair (for GIV).</param>
    public readonly record struct ArtifactObservation(
        char ReferenceAllele,
        char AlternateAllele,
        int RefForward,
        int RefReverse,
        int AltForward,
        int AltReverse,
        int AltReadsR1,
        int AltReadsR2);

    /// <summary>
    /// Classification of a single candidate variant as a sequencing artifact (or not), with the supporting
    /// substitution class, GIV score and Phred-scaled strand-bias score.
    /// </summary>
    /// <param name="Type">Artifact class by substitution (and read orientation for OxoG).</param>
    /// <param name="GivScore">GIV imbalance for this variant's substitution (R1/R2); 1.0 means balanced.</param>
    /// <param name="StrandBiasPhred">Phred-scaled two-sided Fisher strand-bias score (FS); 0 means none.</param>
    /// <param name="IsArtifact">True when the variant is flagged as a likely artifact.</param>
    public readonly record struct ArtifactCall(
        ArtifactType Type,
        double GivScore,
        double StrandBiasPhred,
        bool IsArtifact);

    /// <summary>
    /// Computes the GIV (Global Imbalance Value) for one substitution type as the ratio of read-1 to read-2
    /// alternate-supporting read counts: GIV = r1Count / r2Count. For OxoG the canonical substitution is
    /// G&gt;T, which appears in excess in read 1 (its reverse complement C&gt;A appearing in read 2), so an
    /// elevated GIV signals oxidative damage. GIV = 1 means a balanced, undamaged library; GIV &gt; 1.5 is
    /// defined as damaged. When both counts are 0 there is no imbalance evidence and GIV = 1; when only the
    /// read-2 count is 0 the imbalance is maximal and GIV = <see cref="double.PositiveInfinity"/>.
    /// Source: Chen et al. (2017); Nature Methods (2017); Ettwiller Damage-estimator.
    /// </summary>
    /// <param name="r1Count">Alternate-supporting reads from read 1 (≥ 0).</param>
    /// <param name="r2Count">Alternate-supporting reads from read 2 (≥ 0).</param>
    /// <returns>The GIV ratio r1Count / r2Count (≥ 0); 1.0 when both are 0; +∞ when only r2Count is 0.</returns>
    /// <exception cref="ArgumentOutOfRangeException">A count is negative.</exception>
    public static double CalculateGivScore(int r1Count, int r2Count)
    {
        if (r1Count < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(r1Count), "Read-1 count cannot be negative.");
        }

        if (r2Count < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(r2Count), "Read-2 count cannot be negative.");
        }

        if (r2Count == 0)
        {
            // No read-2 support: balanced when read 1 is also empty (no imbalance evidence), otherwise a
            // maximal one-sided imbalance (Chen et al. 2017 — GIV is an R1/R2 ratio).
            return r1Count == 0 ? UndamagedGivScore : double.PositiveInfinity;
        }

        return (double)r1Count / r2Count;
    }

    /// <summary>
    /// Computes the GATK FisherStrand score FS: the Phred-scaled p-value of a two-sided Fisher exact test on
    /// the 2×2 strand contingency table [refForward, refReverse, altForward, altReverse], testing whether the
    /// reference and alternate alleles are distributed differently across forward/reverse strands (strand
    /// bias). FS = −10·log₁₀(max(p, MIN_PVALUE)); FS = 0 when there is no bias (p = 1) and grows as the
    /// alleles segregate by strand. Source: Broad Institute GATK FisherStrand / StrandBiasTest
    /// (table cell ordering ref-fwd, ref-rev, alt-fwd, alt-rev; FS = phredScaleErrorRate(p)).
    /// </summary>
    /// <param name="refForward">Reference reads on the forward strand (≥ 0).</param>
    /// <param name="refReverse">Reference reads on the reverse strand (≥ 0).</param>
    /// <param name="altForward">Alternate reads on the forward strand (≥ 0).</param>
    /// <param name="altReverse">Alternate reads on the reverse strand (≥ 0).</param>
    /// <returns>The Phred-scaled FisherStrand score FS (≥ 0).</returns>
    /// <exception cref="ArgumentOutOfRangeException">A count is negative.</exception>
    public static double CalculateStrandBias(int refForward, int refReverse, int altForward, int altReverse)
    {
        RequireNonNegative(refForward, nameof(refForward));
        RequireNonNegative(refReverse, nameof(refReverse));
        RequireNonNegative(altForward, nameof(altForward));
        RequireNonNegative(altReverse, nameof(altReverse));

        double pValue = FisherExactTwoSided(refForward, refReverse, altForward, altReverse);
        double floored = Math.Max(pValue, MinFisherPValue);

        // Phred-scaled error rate: FS = -10 * log10(p) (GATK QualityUtils.phredScaleErrorRate).
        return -10.0 * Math.Log10(floored);
    }

    /// <summary>
    /// Classifies one candidate variant as an artifact by substitution class. C&gt;T / G&gt;A are FFPE
    /// cytosine-deamination artifacts (Do &amp; Dobrovic 2015); G&gt;T / C&gt;A are OxoG oxidative artifacts
    /// (Chen et al. 2017). Any other substitution is <see cref="ArtifactType.None"/>. The returned call also
    /// carries the GIV score (from the read-mate alt counts) and the Phred-scaled strand-bias FS, and flags
    /// the variant as an artifact when it is a deamination class, OR an OxoG class whose GIV exceeds the
    /// damaged threshold (1.5).
    /// </summary>
    /// <param name="observation">The candidate variant with its strand / read-mate read evidence.</param>
    /// <returns>The artifact classification for this variant.</returns>
    public static ArtifactCall ClassifyArtifact(ArtifactObservation observation)
    {
        char reference = char.ToUpperInvariant(observation.ReferenceAllele);
        char alternate = char.ToUpperInvariant(observation.AlternateAllele);

        ArtifactType type = ClassifySubstitution(reference, alternate);
        double giv = CalculateGivScore(observation.AltReadsR1, observation.AltReadsR2);
        double strandBias = CalculateStrandBias(
            observation.RefForward, observation.RefReverse, observation.AltForward, observation.AltReverse);

        // Deamination (C>T/G>A) is flagged by substitution class alone; OxoG (G>T/C>A) is confirmed by the
        // read-orientation imbalance (GIV > 1.5 = damaged, Chen et al. 2017).
        bool isArtifact = type switch
        {
            ArtifactType.FfpeDeamination => true,
            ArtifactType.OxoG => giv > DamagedGivThreshold,
            _ => false
        };

        return new ArtifactCall(type, giv, strandBias, isArtifact);
    }

    /// <summary>
    /// Detects OxoG (8-oxoguanine) artifacts among candidate variants: returns the calls for variants whose
    /// substitution is the OxoG class (G&gt;T / C&gt;A) AND whose GIV read-orientation imbalance exceeds the
    /// damaged threshold (GIV &gt; 1.5). Source: Chen et al. (2017).
    /// </summary>
    /// <param name="variants">Candidate variant observations.</param>
    /// <returns>The OxoG-artifact calls, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    public static IReadOnlyList<ArtifactCall> DetectOxoGArtifacts(IEnumerable<ArtifactObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var oxoG = new List<ArtifactCall>();
        foreach (var variant in variants)
        {
            ArtifactCall call = ClassifyArtifact(variant);
            if (call.Type == ArtifactType.OxoG && call.IsArtifact)
            {
                oxoG.Add(call);
            }
        }

        return oxoG;
    }

    /// <summary>
    /// Filters sequencing artifacts out of a candidate variant set, returning only the variants that are NOT
    /// flagged as artifacts (per <see cref="ClassifyArtifact"/>). The result is always a subset of the input,
    /// in input order. Source: composition of the FFPE-deamination and OxoG substitution-class rules
    /// (Do &amp; Dobrovic 2015; Chen et al. 2017).
    /// </summary>
    /// <param name="variants">Candidate variant observations.</param>
    /// <returns>The subset of <paramref name="variants"/> not classified as artifacts, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    public static IReadOnlyList<ArtifactObservation> FilterArtifacts(IEnumerable<ArtifactObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var kept = new List<ArtifactObservation>();
        foreach (var variant in variants)
        {
            if (!ClassifyArtifact(variant).IsArtifact)
            {
                kept.Add(variant);
            }
        }

        return kept;
    }

    /// <summary>
    /// Maps a single-base substitution to its artifact class. C&gt;T / G&gt;A = FFPE deamination;
    /// G&gt;T / C&gt;A = OxoG oxidation; everything else = none (Do &amp; Dobrovic 2015; Chen et al. 2017).
    /// </summary>
    private static ArtifactType ClassifySubstitution(char reference, char alternate)
    {
        return (reference, alternate) switch
        {
            ('C', 'T') => ArtifactType.FfpeDeamination,
            ('G', 'A') => ArtifactType.FfpeDeamination,
            ('G', 'T') => ArtifactType.OxoG,
            ('C', 'A') => ArtifactType.OxoG,
            _ => ArtifactType.None
        };
    }

    private static void RequireNonNegative(int value, string paramName)
    {
        if (value < 0)
        {
            throw new ArgumentOutOfRangeException(paramName, value, "Read count cannot be negative.");
        }
    }

    /// <summary>
    /// Two-sided Fisher exact test p-value for the 2×2 table
    /// <code>[[a, b], [c, d]]</code> (a = refForward, b = refReverse, c = altForward, d = altReverse).
    /// Sums the hypergeometric probabilities of all tables with the same margins whose probability is
    /// ≤ that of the observed table (the conventional two-sided definition). Source: GATK FisherStrand uses
    /// <c>FisherExactTest.twoSidedPValue</c> on this 2×2 strand table.
    /// </summary>
    private static double FisherExactTwoSided(int a, int b, int c, int d)
    {
        int rowOne = a + b;
        int rowTwo = c + d;
        int colOne = a + c;
        int total = a + b + c + d;

        if (total == 0)
        {
            // An empty table provides no evidence of strand bias: p = 1.
            return 1.0;
        }

        double observedLogProb = HypergeometricLogProbability(a, rowOne, rowTwo, colOne, total);

        // The cell 'a' ranges over the values compatible with the fixed margins.
        int minA = Math.Max(0, colOne - rowTwo);
        int maxA = Math.Min(rowOne, colOne);

        double pValue = 0.0;
        for (int x = minA; x <= maxA; x++)
        {
            double logProb = HypergeometricLogProbability(x, rowOne, rowTwo, colOne, total);
            // Include tables at least as extreme (≤ observed probability), with a small tolerance for
            // floating-point comparison of equal-probability tables.
            if (logProb <= observedLogProb + 1e-7)
            {
                pValue += Math.Exp(logProb);
            }
        }

        return Math.Min(1.0, pValue);
    }

    /// <summary>
    /// log P(table) under the hypergeometric distribution for a 2×2 table with the given margins, where the
    /// top-left cell equals <paramref name="a"/>:
    /// P = C(rowOne, a)·C(rowTwo, colOne−a) / C(total, colOne).
    /// </summary>
    private static double HypergeometricLogProbability(int a, int rowOne, int rowTwo, int colOne, int total)
    {
        return LogChoose(rowOne, a)
             + LogChoose(rowTwo, colOne - a)
             - LogChoose(total, colOne);
    }

    /// <summary>log of the binomial coefficient C(n, k) via log-gamma (numerically stable for read counts).</summary>
    private static double LogChoose(int n, int k)
    {
        if (k < 0 || k > n)
        {
            return double.NegativeInfinity;
        }

        return LogFactorial(n) - LogFactorial(k) - LogFactorial(n - k);
    }

    /// <summary>log(n!) = logΓ(n+1).</summary>
    private static double LogFactorial(int n) => LogGamma(n + 1.0);

    /// <summary>
    /// Lanczos approximation of the natural log of the gamma function. Coefficients g = 7, n = 9 per the
    /// standard Lanczos series (Numerical Recipes / Lanczos 1964); accurate to ~1e-13 for the positive
    /// arguments used here, so the resulting binomial coefficients are exact to floating precision for the
    /// read-count magnitudes encountered in strand-bias tables.
    /// </summary>
    private static double LogGamma(double x)
    {
        // Lanczos coefficients (g = 7).
        double[] coefficients =
        {
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        };

        const double g = 7.0;
        x -= 1.0;
        double sum = coefficients[0];
        for (int i = 1; i < coefficients.Length; i++)
        {
            sum += coefficients[i] / (x + i);
        }

        double t = x + g + 0.5;
        return 0.5 * Math.Log(2.0 * Math.PI) + (x + 0.5) * Math.Log(t) - t + Math.Log(sum);
    }

    #endregion


    #region Cancer Variant Annotation (AMP/ASCO/CAP 2017 tiers)

    /// <summary>
    /// Minor-allele-frequency (MAF) cutoff at or above which a variant is treated as a common
    /// polymorphism and classified Tier IV (benign / likely benign). Source: Li MM et al. (2017),
    /// J Mol Diagn 19(1):4–23 — "the work group recommends using 1% (0.01) as a primary cutoff" for
    /// eliminating polymorphic or benign variants (Population Databases section), and Table 7 (Tier IV)
    /// lists "MAF ≥ 1% in the general population" as the population-database criterion. Value = 0.01.
    /// </summary>
    public const double BenignPopulationMafThreshold = 0.01;

    /// <summary>
    /// Strength of clinical/experimental evidence supporting a variant as a biomarker, per the
    /// four evidence levels of Li MM et al. (2017), Table 3 / Figure 2. Levels A and B map to Tier I
    /// (strong clinical significance); Levels C and D map to Tier II (potential clinical significance).
    /// </summary>
    public enum ClinicalEvidenceLevel
    {
        /// <summary>No biomarker evidence level assigned (the variant is not a known biomarker).</summary>
        None,

        /// <summary>
        /// Level A: biomarkers that predict response/resistance to FDA-approved therapies for a specific
        /// tumor type, or are included in professional guidelines (therapeutic/diagnostic/prognostic).
        /// Maps to Tier I. Source: Li et al. (2017), Table 3.
        /// </summary>
        A,

        /// <summary>
        /// Level B: biomarkers based on well-powered studies with expert consensus. Maps to Tier I.
        /// Source: Li et al. (2017), Table 3.
        /// </summary>
        B,

        /// <summary>
        /// Level C: FDA-approved/guideline therapies for a different tumor type (off-label), clinical-trial
        /// inclusion criteria, or diagnostic/prognostic significance from multiple small studies. Maps to
        /// Tier II. Source: Li et al. (2017), Table 3.
        /// </summary>
        C,

        /// <summary>
        /// Level D: plausible therapeutic significance from preclinical studies, or diagnostic/prognostic
        /// support from small studies / case reports without consensus. Maps to Tier II.
        /// Source: Li et al. (2017), Table 3.
        /// </summary>
        D
    }

    /// <summary>
    /// AMP/ASCO/CAP 2017 four-tier clinical-significance classification of a somatic sequence variant.
    /// Source: Li MM et al. (2017), J Mol Diagn 19(1):4–23, Figure 2.
    /// </summary>
    public enum VariantTier
    {
        /// <summary>Tier I: variants of strong clinical significance (Level A or B evidence).</summary>
        TierI_StrongClinicalSignificance,

        /// <summary>Tier II: variants of potential clinical significance (Level C or D evidence).</summary>
        TierII_PotentialClinicalSignificance,

        /// <summary>
        /// Tier III: variants of unknown clinical significance — not common in population databases and
        /// with no convincing published evidence of cancer association.
        /// </summary>
        TierIII_UnknownClinicalSignificance,

        /// <summary>
        /// Tier IV: benign or likely benign variants — observed at a significant allele frequency
        /// (MAF ≥ 1%) in population databases, or with no evidence of cancer association.
        /// </summary>
        TierIV_BenignOrLikelyBenign
    }

    /// <summary>
    /// Caller-supplied evidence for one somatic variant, reduced to the features the AMP/ASCO/CAP 2017
    /// tiering rule consumes (Li et al. 2017, Figure 2 / Tables 4–7). The guideline classifies variants
    /// from external knowledge (professional guidelines, population databases, somatic databases,
    /// literature); this library does not reproduce those curated resources — the relevant facts are
    /// supplied by the caller, who has performed the database lookups.
    /// </summary>
    /// <param name="Gene">Gene symbol the variant falls in.</param>
    /// <param name="ProteinChange">HGVS protein change (e.g. p.V600E); informational.</param>
    /// <param name="EvidenceLevel">Strongest assigned clinical evidence level (A–D, or None).</param>
    /// <param name="PopulationMaf">
    /// Minor allele frequency in a population database (e.g. gnomAD/ExAC/1000 Genomes), in [0, 1].
    /// A value ≥ <see cref="BenignPopulationMafThreshold"/> indicates a common polymorphism.
    /// </param>
    /// <param name="HasCancerAssociation">
    /// True when there is published evidence associating the variant with cancer (somatic database
    /// presence, functional/population study). Distinguishes Tier III from Tier IV when MAF is low.
    /// </param>
    public readonly record struct CancerVariantAnnotationInput(
        string Gene,
        string ProteinChange,
        ClinicalEvidenceLevel EvidenceLevel,
        double PopulationMaf,
        bool HasCancerAssociation);

    /// <summary>The tier classification of one variant, with the input evidence that produced it.</summary>
    /// <param name="Variant">The variant evidence that was classified.</param>
    /// <param name="Tier">Assigned AMP/ASCO/CAP 2017 tier.</param>
    public readonly record struct CancerVariantAnnotation(
        CancerVariantAnnotationInput Variant,
        VariantTier Tier);

    /// <summary>
    /// Classifies a single somatic variant into the AMP/ASCO/CAP 2017 four-tier system from caller-supplied
    /// evidence, applying the decision criteria of Li MM et al. (2017), Figure 2 in priority order:
    /// <list type="number">
    /// <item><description>Level A or B evidence ⇒ <see cref="VariantTier.TierI_StrongClinicalSignificance"/>.</description></item>
    /// <item><description>Level C or D evidence ⇒ <see cref="VariantTier.TierII_PotentialClinicalSignificance"/>.</description></item>
    /// <item><description>Otherwise, MAF ≥ 1% (common polymorphism) OR no cancer association ⇒
    /// <see cref="VariantTier.TierIV_BenignOrLikelyBenign"/> (Table 7).</description></item>
    /// <item><description>Otherwise (rare, no clinical evidence, but a cancer association exists) ⇒
    /// <see cref="VariantTier.TierIII_UnknownClinicalSignificance"/> (Table 6).</description></item>
    /// </list>
    /// Clinical evidence (Tier I/II) is evaluated before the benign-frequency rule because a Level A/B
    /// biomarker remains strongly significant even if it also appears in population databases; Table 4
    /// (Tier I) and Table 5 (Tier II) note such variants are "absent or extremely low MAF" but the
    /// guideline assigns them by evidence level, not frequency.
    /// </summary>
    /// <param name="variant">Caller-supplied evidence for the variant.</param>
    /// <returns>The variant's tier classification.</returns>
    /// <exception cref="ArgumentOutOfRangeException">PopulationMaf is NaN or outside [0, 1].</exception>
    public static VariantTier ClassifyVariantTier(CancerVariantAnnotationInput variant)
    {
        if (double.IsNaN(variant.PopulationMaf) || variant.PopulationMaf < 0.0 || variant.PopulationMaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.PopulationMaf, "Population MAF must be in the range [0, 1].");
        }

        // Tier I — strong clinical significance: Level A or B evidence (Li et al. 2017, Figure 2).
        if (variant.EvidenceLevel is ClinicalEvidenceLevel.A or ClinicalEvidenceLevel.B)
        {
            return VariantTier.TierI_StrongClinicalSignificance;
        }

        // Tier II — potential clinical significance: Level C or D evidence (Li et al. 2017, Figure 2).
        if (variant.EvidenceLevel is ClinicalEvidenceLevel.C or ClinicalEvidenceLevel.D)
        {
            return VariantTier.TierII_PotentialClinicalSignificance;
        }

        // No clinical evidence level. Tier IV — benign/likely benign: observed at a significant allele
        // frequency (MAF ≥ 1%, Table 7), OR no published evidence of cancer association (Figure 2).
        if (variant.PopulationMaf >= BenignPopulationMafThreshold || !variant.HasCancerAssociation)
        {
            return VariantTier.TierIV_BenignOrLikelyBenign;
        }

        // Tier III — unknown clinical significance: rare (low MAF), no clinical evidence, but a cancer
        // association exists so it cannot be called benign (Li et al. 2017, Table 6).
        return VariantTier.TierIII_UnknownClinicalSignificance;
    }

    /// <summary>
    /// Annotates a set of somatic variants with their AMP/ASCO/CAP 2017 clinical-significance tiers by
    /// applying <see cref="ClassifyVariantTier"/> to each variant. The output preserves input order and
    /// has one entry per input variant. Source: Li MM et al. (2017), J Mol Diagn 19(1):4–23.
    /// </summary>
    /// <param name="variants">Caller-supplied variant evidence records.</param>
    /// <returns>One <see cref="CancerVariantAnnotation"/> per input variant, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A variant's PopulationMaf is outside [0, 1].</exception>
    public static IReadOnlyList<CancerVariantAnnotation> AnnotateCancerVariants(
        IEnumerable<CancerVariantAnnotationInput> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var annotations = new List<CancerVariantAnnotation>();
        foreach (var variant in variants)
        {
            annotations.Add(new CancerVariantAnnotation(variant, ClassifyVariantTier(variant)));
        }

        return annotations;
    }

    /// <summary>
    /// Looks up a variant's COSMIC (Catalogue Of Somatic Mutations In Cancer) annotation in a
    /// caller-supplied catalog keyed by (gene, protein change). COSMIC is a large, expert-curated
    /// somatic-mutation database (Tate JG et al. 2019, Nucleic Acids Res 47:D941–D947) that cannot be
    /// reproduced or hardcoded here; the caller passes the relevant records (e.g. a COSMIC export),
    /// and this method performs the exact-match lookup the AMP/ASCO/CAP workflow uses to flag a variant
    /// as present in a somatic database (Li et al. 2017, Tables 4–6, "Somatic database: COSMIC...").
    /// </summary>
    /// <param name="variant">The variant to look up.</param>
    /// <param name="cosmicCatalog">
    /// Caller-supplied COSMIC records keyed by (gene, protein change); e.g. COSMIC identifier strings.
    /// </param>
    /// <returns>
    /// The catalog value (e.g. a COSMIC ID) for the variant's (gene, protein change), or <c>null</c>
    /// when the variant is not present in the supplied catalog.
    /// </returns>
    /// <exception cref="ArgumentNullException"><paramref name="cosmicCatalog"/> is null.</exception>
    public static string? GetCOSMICAnnotation(
        CancerVariantAnnotationInput variant,
        IReadOnlyDictionary<(string Gene, string ProteinChange), string> cosmicCatalog)
    {
        ArgumentNullException.ThrowIfNull(cosmicCatalog);
        return cosmicCatalog.TryGetValue((variant.Gene, variant.ProteinChange), out string? id) ? id : null;
    }

    #endregion

}
