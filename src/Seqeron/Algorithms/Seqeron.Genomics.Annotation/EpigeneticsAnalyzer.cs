using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Provides algorithms for epigenetic analysis including methylation and chromatin state.
/// </summary>
public static class EpigeneticsAnalyzer
{
    #region Records and Types

    /// <summary>
    /// Represents a methylation site.
    /// </summary>
    public readonly record struct MethylationSite(
        int Position,
        MethylationType Type,
        string Context,
        double MethylationLevel,
        int Coverage);

    /// <summary>
    /// Type of DNA methylation.
    /// </summary>
    public enum MethylationType
    {
        CpG,      // CG context
        CHG,      // CHG context (H = A, C, or T)
        CHH,      // CHH context
        N6A,      // 6-methyladenine (bacterial)
        N4C       // 4-methylcytosine (bacterial)
    }

    /// <summary>
    /// Represents a differentially methylated region.
    /// </summary>
    public readonly record struct DifferentiallyMethylatedRegion(
        int Start,
        int End,
        double MeanDifference,
        double PValue,
        int CpGCount,
        string Annotation);

    /// <summary>
    /// Represents a methylation profile.
    /// </summary>
    public readonly record struct MethylationProfile(
        double GlobalMethylation,
        double CpGMethylation,
        double CHGMethylation,
        double CHHMethylation,
        int TotalCpGSites,
        int MethylatedCpGSites,
        IReadOnlyList<(int Position, double Level)> MethylationByPosition);

    /// <summary>
    /// Represents a chromatin state.
    /// </summary>
    /// <summary>
    /// Chromatin states inferred from combinatorial histone-mark signatures.
    /// Names follow the Roadmap Epigenomics core/expanded model mnemonics
    /// (TssA, Enh, Tx, ReprPC, Het, TssBiv, EnhBiv, Quies).
    /// </summary>
    /// <remarks>Roadmap Epigenomics chromatin state learning (15/18-state models).</remarks>
    public enum ChromatinState
    {
        /// <summary>Active TSS / promoter (Roadmap TssA) — H3K4me3 present.</summary>
        ActivePromoter,

        /// <summary>Active enhancer (Roadmap active Enh) — H3K4me1 + H3K27ac present.</summary>
        ActiveEnhancer,

        /// <summary>Weak / poised enhancer (Roadmap Enh) — H3K4me1 present without H3K27ac.</summary>
        WeakEnhancer,

        /// <summary>Transcribed gene body (Roadmap Tx) — H3K36me3 present.</summary>
        Transcribed,

        /// <summary>Polycomb-repressed (Roadmap ReprPC) — H3K27me3 present.</summary>
        Repressed,

        /// <summary>Constitutive heterochromatin (Roadmap Het) — H3K9me3 present.</summary>
        Heterochromatin,

        /// <summary>Bivalent / poised promoter (Roadmap TssBiv) — H3K4me3 + H3K27me3.</summary>
        BivalentPromoter,

        /// <summary>Bivalent enhancer (Roadmap EnhBiv) — H3K4me1 + H3K27me3.</summary>
        BivalentEnhancer,

        /// <summary>Quiescent / low signal (Roadmap Quies) — no mark present.</summary>
        LowSignal
    }

    /// <summary>
    /// Represents a histone modification.
    /// </summary>
    public readonly record struct HistoneModification(
        int Start,
        int End,
        string Mark,
        double Signal,
        ChromatinState PredictedState);

    /// <summary>
    /// Represents a chromatin accessibility region.
    /// </summary>
    public readonly record struct AccessibilityRegion(
        int Start,
        int End,
        double AccessibilityScore,
        string PeakType,
        IReadOnlyList<string> NearbyGenes);

    /// <summary>
    /// Represents an imprinted gene prediction.
    /// </summary>
    public readonly record struct ImprintedGene(
        string GeneId,
        int Start,
        int End,
        double ImprintingScore,
        string ParentalOrigin,
        bool HasDMR);

    /// <summary>
    /// Represents a genomic feature annotation (e.g. a gene or regulatory element)
    /// used to label differentially methylated regions by genomic context.
    /// Coordinates are 0-based, half-open in the same coordinate space as the
    /// methylation site positions.
    /// </summary>
    public readonly record struct GeneAnnotation(
        string Feature,
        int Start,
        int End);

    #endregion

    #region CpG Site Detection

    /// <summary>
    /// Finds all CpG dinucleotides in a sequence.
    /// </summary>
    public static IEnumerable<int> FindCpGSites(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        sequence = sequence.ToUpperInvariant();

        for (int i = 0; i < sequence.Length - 1; i++)
        {
            if (sequence[i] == 'C' && sequence[i + 1] == 'G')
            {
                yield return i;
            }
        }
    }

    /// <summary>
    /// Classifies the cytosine at <paramref name="index"/> into its DNA-methylation
    /// sequence context: CpG (C followed by G), CHG (C, H, G), or CHH (C, H, H),
    /// where H is the IUPAC ambiguity code for A, C, or T (i.e. "not G").
    /// </summary>
    /// <param name="sequence">DNA sequence (case-insensitive).</param>
    /// <param name="index">0-based position; must reference a cytosine.</param>
    /// <returns>
    /// The methylation context, or <c>null</c> when the base is not a cytosine,
    /// the downstream context is incomplete, or a context base is not in {A,C,G,T}.
    /// </returns>
    /// <remarks>
    /// Context definitions: Krueger &amp; Andrews (2011) Bismark; Lister et al. (2009).
    /// H = A, C, or T per IUPAC nucleotide nomenclature (Cornish-Bowden 1985).
    /// </remarks>
    public static MethylationType? GetMethylationContext(string sequence, int index)
    {
        if (string.IsNullOrEmpty(sequence) || index < 0 || index >= sequence.Length)
            return null;

        if (ToUpperBase(sequence[index]) != 'C')
            return null;

        // Need at least one downstream base to classify any context.
        if (index + 1 >= sequence.Length)
            return null;

        char next = ToUpperBase(sequence[index + 1]);

        // CpG: cytosine immediately followed by guanine (Bismark; Lister 2009).
        if (next == 'G')
            return MethylationType.CpG;

        // Non-CpG contexts require the second base to be H (A, C, or T).
        if (!IsHBase(next))
            return null;

        // CHG / CHH require a third base in the window.
        if (index + 2 >= sequence.Length)
            return null;

        char third = ToUpperBase(sequence[index + 2]);

        // CHG: C, H, G (symmetric context). H already validated above.
        if (third == 'G')
            return MethylationType.CHG;

        // CHH: C, H, H (asymmetric context); third base must also be H.
        if (IsHBase(third))
            return MethylationType.CHH;

        return null;
    }

    /// <summary>
    /// Finds all potential cytosine methylation sites with their sequence context.
    /// </summary>
    /// <remarks>
    /// Each classifiable cytosine (CpG, CHG, or CHH per
    /// <see cref="GetMethylationContext(string,int)"/>) is reported with its 0-based
    /// position. <see cref="MethylationSite.MethylationLevel"/> and
    /// <see cref="MethylationSite.Coverage"/> are 0 here because sequence-only input
    /// carries no bisulfite read evidence; measured levels come from
    /// <see cref="CalculateMethylationFromBisulfite"/> or caller-supplied sites.
    /// </remarks>
    public static IEnumerable<MethylationSite> FindMethylationSites(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        for (int i = 0; i < sequence.Length; i++)
        {
            MethylationType? context = GetMethylationContext(sequence, i);
            if (context is null)
                continue;

            // Context window: 3 bases when available, otherwise the terminal CpG pair.
            int windowSize = Math.Min(MethylationContextWindow, sequence.Length - i);
            string contextString = sequence.Substring(i, windowSize).ToUpperInvariant();

            yield return new MethylationSite(
                Position: i,
                Type: context.Value,
                Context: contextString,
                MethylationLevel: 0, // Unknown without bisulfite data
                Coverage: 0);
        }
    }

    // CpG/CHG/CHH context spans the cytosine plus up to two downstream bases.
    private const int MethylationContextWindow = 3;

    // IUPAC code H = A, C, or T ("not G") — Cornish-Bowden (1985), Nucleic Acids Res. 13:3021–3030.
    private static bool IsHBase(char upperBase) =>
        upperBase is 'A' or 'C' or 'T';

    private static char ToUpperBase(char c) =>
        (char)(c >= 'a' && c <= 'z' ? c - ('a' - 'A') : c);

    #endregion

    #region CpG Island Analysis

    /// <summary>
    /// Calculates CpG observed/expected ratio.
    /// </summary>
    public static double CalculateCpGObservedExpected(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return 0;

        sequence = sequence.ToUpperInvariant();

        int c = sequence.Count(ch => ch == 'C');
        int g = sequence.Count(ch => ch == 'G');
        int cpg = 0;

        for (int i = 0; i < sequence.Length - 1; i++)
        {
            if (sequence[i] == 'C' && sequence[i + 1] == 'G')
                cpg++;
        }

        // Compute C×G in double to avoid int overflow on long sequences:
        // for c, g near 10^5 the product exceeds Int32.MaxValue and would wrap,
        // corrupting the Gardiner-Garden & Frommer O/E ratio.
        double expected = ((double)c * g) / sequence.Length;
        return expected > 0 ? cpg / expected : 0;
    }

    /// <summary>
    /// Identifies CpG islands using Gardiner-Garden and Frommer criteria.
    /// </summary>
    public static IEnumerable<(int Start, int End, double GcContent, double CpGRatio)>
        FindCpGIslands(
            string sequence,
            int minLength = 200,
            double minGc = 0.5,
            double minCpGRatio = 0.6)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minLength)
            yield break;

        sequence = sequence.ToUpperInvariant();

        int? islandStart = null;
        int islandEnd = 0;

        for (int i = 0; i <= sequence.Length - minLength; i++)
        {
            int windowSize = Math.Min(minLength, sequence.Length - i);
            string window = sequence.Substring(i, windowSize);

            double gc = CalculateGcContent(window);
            double cpgRatio = CalculateCpGObservedExpected(window);

            bool isCpGIsland = gc >= minGc && cpgRatio >= minCpGRatio;

            if (isCpGIsland)
            {
                if (islandStart == null)
                    islandStart = i;
                islandEnd = i + windowSize;
            }
            else if (islandStart != null)
            {
                if (islandEnd - islandStart.Value >= minLength)
                {
                    string island = sequence.Substring(islandStart.Value, islandEnd - islandStart.Value);
                    double islandGc = CalculateGcContent(island);
                    double islandCpG = CalculateCpGObservedExpected(island);
                    if (islandGc >= minGc && islandCpG >= minCpGRatio)
                    {
                        yield return (islandStart.Value, islandEnd, islandGc, islandCpG);
                    }
                }
                islandStart = null;
            }
        }

        // Handle island at end
        if (islandStart != null && islandEnd - islandStart.Value >= minLength)
        {
            string island = sequence.Substring(islandStart.Value, islandEnd - islandStart.Value);
            double islandGc = CalculateGcContent(island);
            double islandCpG = CalculateCpGObservedExpected(island);
            if (islandGc >= minGc && islandCpG >= minCpGRatio)
            {
                yield return (islandStart.Value, islandEnd, islandGc, islandCpG);
            }
        }
    }

    private static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcFractionFast();

    #endregion

    #region Methylation Analysis

    /// <summary>
    /// Simulates in-silico sodium-bisulfite conversion of a single DNA strand.
    /// </summary>
    /// <remarks>
    /// Per Frommer et al. (1992): bisulfite treatment converts unmethylated cytosine to
    /// uracil (which reads/amplifies as thymine), while 5-methylcytosine remains
    /// nonreactive and stays a cytosine. This method substitutes each unprotected
    /// cytosine with thymine (case preserved: <c>C</c>→<c>T</c>, <c>c</c>→<c>t</c>),
    /// leaves cytosines whose 0-based index is in <paramref name="methylatedPositions"/>
    /// as cytosines, and returns every non-cytosine base unchanged. Only the supplied
    /// strand is converted; the complementary strand is a separate molecule and is not
    /// synthesized or merged (Frommer's protocol is strand-specific).
    /// </remarks>
    /// <param name="sequence">DNA sequence to convert. Null or empty yields an empty string.</param>
    /// <param name="methylatedPositions">
    /// 0-based indices of protected (methylated) cytosines. Null means none are protected.
    /// </param>
    /// <returns>The bisulfite-converted strand; same length as <paramref name="sequence"/>.</returns>
    public static string SimulateBisulfiteConversion(
        string sequence,
        IReadOnlySet<int>? methylatedPositions = null)
    {
        if (string.IsNullOrEmpty(sequence))
            return "";

        methylatedPositions ??= new HashSet<int>();
        var result = new StringBuilder(sequence.Length);

        for (int i = 0; i < sequence.Length; i++)
        {
            char c = sequence[i];

            if (c == 'C' || c == 'c')
            {
                // Methylated cytosines are protected from conversion
                if (methylatedPositions.Contains(i))
                {
                    result.Append(c);
                }
                else
                {
                    // Unmethylated C converts to T (U in bisulfite chemistry)
                    result.Append(c == 'C' ? 'T' : 't');
                }
            }
            else
            {
                result.Append(c);
            }
        }

        return result.ToString();
    }

    /// <summary>
    /// Calls per-CpG methylation levels from aligned bisulfite reads against a reference.
    /// </summary>
    /// <remarks>
    /// Methylation call rule (Krueger &amp; Andrews 2011, Bismark): at a reference CpG
    /// cytosine, a read base of <c>C</c> indicates a methylated cytosine (protected from
    /// conversion) and a read base of <c>T</c> indicates an unmethylated cytosine
    /// (converted). Any other read base is not a valid bisulfite call and is ignored.
    /// The reported methylation level is the Bismark fraction
    /// <c>methylated / (methylated + unmethylated)</c> and <c>Coverage</c> is the number of
    /// valid C/T calls. CpG sites with zero coverage have an undefined percentage and are
    /// excluded from the result. Read bases falling outside the reference are ignored.
    /// </remarks>
    /// <param name="referenceSequence">Reference DNA sequence (case-insensitive).</param>
    /// <param name="bisulfiteReads">Aligned reads as (read sequence, 0-based start position) pairs.</param>
    /// <returns>One <see cref="MethylationSite"/> per covered CpG site, in reference order.</returns>
    public static IEnumerable<MethylationSite> CalculateMethylationFromBisulfite(
        string referenceSequence,
        IEnumerable<(string ReadSequence, int StartPosition)> bisulfiteReads)
    {
        if (string.IsNullOrEmpty(referenceSequence))
            yield break;

        referenceSequence = referenceSequence.ToUpperInvariant();

        // Find all CpG sites
        var cpgSites = FindCpGSites(referenceSequence).ToList();
        var siteData = cpgSites.ToDictionary(
            site => site,
            site => (Methylated: 0, Total: 0));

        foreach (var (readSeq, startPos) in bisulfiteReads)
        {
            string read = readSeq.ToUpperInvariant();

            for (int i = 0; i < read.Length && startPos + i < referenceSequence.Length - 1; i++)
            {
                int refPos = startPos + i;

                if (siteData.ContainsKey(refPos))
                {
                    var (meth, total) = siteData[refPos];

                    // C in read at CpG site means methylated (protected from conversion)
                    // T in read means unmethylated (C converted to T)
                    if (read[i] == 'C')
                    {
                        siteData[refPos] = (meth + 1, total + 1);
                    }
                    else if (read[i] == 'T')
                    {
                        siteData[refPos] = (meth, total + 1);
                    }
                }
            }
        }

        foreach (var site in cpgSites)
        {
            var (meth, total) = siteData[site];
            if (total == 0)
                continue;

            string context = site + 2 < referenceSequence.Length
                ? referenceSequence.Substring(site, 3)
                : referenceSequence.Substring(site);

            yield return new MethylationSite(
                Position: site,
                Type: MethylationType.CpG,
                Context: context,
                MethylationLevel: (double)meth / total,
                Coverage: total);
        }
    }

    /// <summary>
    /// Generates a methylation profile (per-context weighted methylation levels and
    /// CpG counts) from a set of measured methylation sites.
    /// </summary>
    /// <remarks>
    /// Per-context levels are <b>weighted methylation levels</b> per Schultz et al. (2012):
    /// the sum of methylated reads divided by the sum of total reads over all cytosines
    /// of that context, i.e. Σ(level·coverage) / Σ(coverage). When every site has equal
    /// coverage this equals the unweighted mean of per-site fractions; under unequal
    /// coverage the weighted level is the correct aggregate. Sites with zero coverage
    /// (e.g. sequence-only sites from <see cref="FindMethylationSites"/>) fall back to an
    /// unweighted mean of per-site levels so they are not silently dropped.
    /// </remarks>
    public static MethylationProfile GenerateMethylationProfile(IEnumerable<MethylationSite> sites)
    {
        var siteList = sites.ToList();

        if (siteList.Count == 0)
        {
            return new MethylationProfile(0, 0, 0, 0, 0, 0, new List<(int, double)>());
        }

        var cpgSites = siteList.Where(s => s.Type == MethylationType.CpG).ToList();
        var chgSites = siteList.Where(s => s.Type == MethylationType.CHG).ToList();
        var chhSites = siteList.Where(s => s.Type == MethylationType.CHH).ToList();

        double globalMeth = WeightedMethylationLevel(siteList);
        double cpgMeth = WeightedMethylationLevel(cpgSites);
        double chgMeth = WeightedMethylationLevel(chgSites);
        double chhMeth = WeightedMethylationLevel(chhSites);

        int totalCpG = cpgSites.Count;
        // Descriptive count of "methylated" CpG sites using a fractional cutoff.
        int methylatedCpG = cpgSites.Count(s => s.MethylationLevel >= MethylatedSiteThreshold);

        var byPosition = siteList
            .Select(s => (s.Position, s.MethylationLevel))
            .OrderBy(x => x.Position)
            .ToList();

        return new MethylationProfile(
            GlobalMethylation: globalMeth,
            CpGMethylation: cpgMeth,
            CHGMethylation: chgMeth,
            CHHMethylation: chhMeth,
            TotalCpGSites: totalCpG,
            MethylatedCpGSites: methylatedCpG,
            MethylationByPosition: byPosition);
    }

    // Fractional-methylation cutoff used only for the descriptive MethylatedCpGSites count.
    // Schultz et al. (2012) recommend a binomial test; 0.5 is a convenience threshold and
    // does not affect the continuous (weighted) methylation-level outputs.
    private const double MethylatedSiteThreshold = 0.5;

    /// <summary>
    /// Weighted methylation level per Schultz et al. (2012): Σ(methylated reads) / Σ(total reads).
    /// With per-site fraction = methylated/coverage, this is Σ(level·coverage) / Σ(coverage).
    /// Falls back to the unweighted mean of per-site levels when total coverage is zero.
    /// </summary>
    private static double WeightedMethylationLevel(IReadOnlyCollection<MethylationSite> sites)
    {
        if (sites.Count == 0)
            return 0;

        double totalCoverage = sites.Sum(s => (double)s.Coverage);
        if (totalCoverage <= 0)
            return sites.Average(s => s.MethylationLevel);

        double methylatedReads = sites.Sum(s => s.MethylationLevel * s.Coverage);
        return methylatedReads / totalCoverage;
    }

    #endregion

    #region Differentially Methylated Regions

    // methylKit default tiling window size in bp (tileMethylCounts win.size=1000).
    // Akalin et al. (2012) Genome Biology 13:R87; tileMethylCounts man page.
    private const int DefaultWindowSize = 1000;

    // methylKit default %methylation-difference cutoff = 25% = 0.25 as a fraction.
    // getMethylDiff difference=25; Akalin et al. (2012). The repository expresses
    // methylation as a fraction in [0,1], so 25 percentage points == 0.25.
    private const double DefaultMinDifference = 0.25;

    // A DMR is a region of *adjacent* CpG sites (Akalin et al. 2012); a single
    // site is not a region. Minimum number of covered cytosines per reported region.
    private const int DefaultMinCpGCount = 3;

    /// <summary>
    /// Identifies differentially methylated regions (DMRs) between two single-sample
    /// methylation profiles using the methylKit tiling-window model: positions are
    /// grouped into fixed-size genomic windows; within each window the per-site
    /// methylation differences (sample2 − sample1, in fraction units) are averaged,
    /// and significance is assessed with a two-sided Fisher's exact test on the pooled
    /// methylated/unmethylated read counts of the window.
    /// </summary>
    /// <param name="sample1">Control sample methylation sites (per-site level + coverage).</param>
    /// <param name="sample2">Treatment sample methylation sites.</param>
    /// <param name="windowSize">Tiling window width in bp. Default 1000 (methylKit win.size).</param>
    /// <param name="minDifference">
    /// Minimum absolute mean methylation difference to report, as a fraction in [0,1].
    /// Default 0.25 (methylKit difference=25%). The cutoff is strict: |meanDiff| must be
    /// greater than this value (methylKit getMethylDiff uses <c>meth.diff &gt; difference</c>).
    /// </param>
    /// <param name="minCpGCount">Minimum number of covered cytosines for a window to be a region. Default 3.</param>
    /// <returns>
    /// DMRs ordered by start position. Each region is annotated "Hypermethylated" when the
    /// treatment has higher methylation than control (positive difference) and
    /// "Hypomethylated" when lower (negative difference), per methylKit getMethylDiff.
    /// </returns>
    /// <exception cref="ArgumentNullException">A sample enumerable is null.</exception>
    public static IEnumerable<DifferentiallyMethylatedRegion> FindDMRs(
        IEnumerable<MethylationSite> sample1,
        IEnumerable<MethylationSite> sample2,
        int windowSize = DefaultWindowSize,
        double minDifference = DefaultMinDifference,
        int minCpGCount = DefaultMinCpGCount)
    {
        ArgumentNullException.ThrowIfNull(sample1);
        ArgumentNullException.ThrowIfNull(sample2);

        return FindDMRsIterator(sample1, sample2, windowSize, minDifference, minCpGCount);
    }

    private static IEnumerable<DifferentiallyMethylatedRegion> FindDMRsIterator(
        IEnumerable<MethylationSite> sample1,
        IEnumerable<MethylationSite> sample2,
        int windowSize,
        double minDifference,
        int minCpGCount)
    {
        var sites1 = new Dictionary<int, MethylationSite>();
        foreach (var s in sample1)
            sites1[s.Position] = s;

        var sites2 = new Dictionary<int, MethylationSite>();
        foreach (var s in sample2)
            sites2[s.Position] = s;

        var allPositions = sites1.Keys.Union(sites2.Keys).OrderBy(p => p).ToList();
        if (allPositions.Count == 0)
            yield break;

        int start = allPositions[0];
        var window = new List<int>();

        foreach (var pos in allPositions)
        {
            // methylKit tiling: a window spans [start, start+windowSize). A position
            // at or beyond windowSize from the window start opens a new window.
            if (pos - start >= windowSize)
            {
                if (TryBuildRegion(window, sites1, sites2, start, minDifference, minCpGCount, out var dmr))
                    yield return dmr;

                start = pos;
                window.Clear();
            }

            window.Add(pos);
        }

        if (TryBuildRegion(window, sites1, sites2, start, minDifference, minCpGCount, out var lastDmr))
            yield return lastDmr;
    }

    private static bool TryBuildRegion(
        List<int> positions,
        IReadOnlyDictionary<int, MethylationSite> sites1,
        IReadOnlyDictionary<int, MethylationSite> sites2,
        int start,
        double minDifference,
        int minCpGCount,
        out DifferentiallyMethylatedRegion region)
    {
        region = default;

        if (positions.Count < minCpGCount)
            return false;

        double sumDiff = 0;
        // Pooled 2x2 contingency table for the window's Fisher's exact test:
        //              methylated      unmethylated
        //   sample1    numC1           numT1
        //   sample2    numC2           numT2
        long numC1 = 0, numT1 = 0, numC2 = 0, numT2 = 0;

        foreach (var pos in positions)
        {
            double level1 = sites1.TryGetValue(pos, out var s1) ? s1.MethylationLevel : 0;
            double level2 = sites2.TryGetValue(pos, out var s2) ? s2.MethylationLevel : 0;
            sumDiff += level2 - level1;

            // Reconstruct integer C/T read counts from fractional level x coverage.
            // methylKit operates on numCs/numTs directly (C/(C+T)); MethylationSite
            // stores a fraction + coverage, so counts are recovered by rounding.
            if (sites1.TryGetValue(pos, out var c1) && c1.Coverage > 0)
            {
                long meth = (long)Math.Round(c1.MethylationLevel * c1.Coverage);
                numC1 += meth;
                numT1 += c1.Coverage - meth;
            }

            if (sites2.TryGetValue(pos, out var c2) && c2.Coverage > 0)
            {
                long meth = (long)Math.Round(c2.MethylationLevel * c2.Coverage);
                numC2 += meth;
                numT2 += c2.Coverage - meth;
            }
        }

        double meanDiff = sumDiff / positions.Count;

        // methylKit getMethylDiff cutoff is STRICT: |meth.diff| > difference.
        if (Math.Abs(meanDiff) <= minDifference)
            return false;

        region = new DifferentiallyMethylatedRegion(
            Start: start,
            End: positions[^1],
            MeanDifference: meanDiff,
            PValue: FisherExactTwoSided(numC1, numT1, numC2, numT2),
            CpGCount: positions.Count,
            // hyper = treatment higher than control; hypo = lower (Akalin et al. 2012).
            Annotation: meanDiff > 0 ? "Hypermethylated" : "Hypomethylated");

        return true;
    }

    /// <summary>
    /// Annotates differentially methylated regions with overlapping genomic features.
    /// Each DMR is labelled with the <c>Feature</c> of the first annotation whose
    /// [Start,End) interval overlaps the DMR; DMRs with no overlap keep their methylation
    /// annotation ("Hypermethylated"/"Hypomethylated"). This is the genomic-context
    /// labelling step that follows DMR calling (Akalin et al. 2012, gene-region annotation).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="dmrs"/> or <paramref name="annotations"/> is null.</exception>
    public static IEnumerable<DifferentiallyMethylatedRegion> AnnotateDMRs(
        IEnumerable<DifferentiallyMethylatedRegion> dmrs,
        IEnumerable<GeneAnnotation> annotations)
    {
        ArgumentNullException.ThrowIfNull(dmrs);
        ArgumentNullException.ThrowIfNull(annotations);

        var features = annotations.ToList();
        return AnnotateDMRsIterator(dmrs, features);
    }

    private static IEnumerable<DifferentiallyMethylatedRegion> AnnotateDMRsIterator(
        IEnumerable<DifferentiallyMethylatedRegion> dmrs,
        List<GeneAnnotation> features)
    {
        foreach (var dmr in dmrs)
        {
            string? label = null;
            foreach (var f in features)
            {
                // Half-open interval overlap: [dmr.Start, dmr.End] vs [f.Start, f.End).
                if (dmr.Start < f.End && f.Start <= dmr.End)
                {
                    label = f.Feature;
                    break;
                }
            }

            yield return label is null
                ? dmr
                : dmr with { Annotation = label };
        }
    }

    /// <summary>
    /// Probability of a single 2×2 contingency table under fixed margins
    /// (the hypergeometric probability that underlies Fisher's exact test):
    /// p = (a+b)!(c+d)!(a+c)!(b+d)! / (a! b! c! d! n!), with n = a+b+c+d.
    /// </summary>
    /// <remarks>Fisher's exact test, hypergeometric form (Fisher 1922, 1935).</remarks>
    public static double FisherExactProbability(long a, long b, long c, long d)
    {
        if (a < 0 || b < 0 || c < 0 || d < 0)
            throw new ArgumentOutOfRangeException(nameof(a), "Contingency-table cells must be non-negative.");

        long n = a + b + c + d;
        if (n == 0)
            return 1.0;

        // Compute in log space to avoid factorial overflow, then exponentiate.
        double logP = LogFactorial(a + b) + LogFactorial(c + d)
                    + LogFactorial(a + c) + LogFactorial(b + d)
                    - LogFactorial(a) - LogFactorial(b)
                    - LogFactorial(c) - LogFactorial(d) - LogFactorial(n);

        return Math.Exp(logP);
    }

    /// <summary>
    /// Two-sided Fisher's exact test p-value for a 2×2 table: the sum of the
    /// probabilities of every table with the same margins whose hypergeometric
    /// probability is ≤ that of the observed table (Fisher's exact test, two-sided rule).
    /// </summary>
    private static double FisherExactTwoSided(long a, long b, long c, long d)
    {
        long n = a + b + c + d;
        if (n == 0)
            return 1.0;

        long row1 = a + b;          // sample1 total reads
        long colMeth = a + c;       // total methylated reads
        // Degenerate margins (a row or column total is 0) admit only one table → p = 1.
        if (row1 == 0 || (c + d) == 0 || colMeth == 0 || (b + d) == 0)
            return 1.0;

        double pObserved = FisherExactProbability(a, b, c, d);

        // Enumerate all tables sharing the fixed margins by varying the top-left cell x.
        // Constraints: x in [max(0, colMeth - (c+d)), min(row1, colMeth)].
        long xMin = Math.Max(0, colMeth - (c + d));
        long xMax = Math.Min(row1, colMeth);

        // Numerical tolerance so a table tied with the observed one is included.
        const double Tolerance = 1e-7;

        double pValue = 0;
        for (long x = xMin; x <= xMax; x++)
        {
            long xa = x;
            long xb = row1 - x;
            long xc = colMeth - x;
            long xd = (c + d) - (colMeth - x);

            double p = FisherExactProbability(xa, xb, xc, xd);
            if (p <= pObserved * (1 + Tolerance))
                pValue += p;
        }

        return Math.Min(1.0, pValue);
    }

    private static double LogFactorial(long k)
    {
        // log(k!) = lgamma(k+1). Math has no Lgamma; use the standard series via
        // sum of logs for the modest counts arising from read coverage.
        if (k < 2)
            return 0.0;

        double sum = 0;
        for (long i = 2; i <= k; i++)
            sum += Math.Log(i);

        return sum;
    }

    #endregion

    #region Histone Modification Analysis

    // ChromHMM models each chromatin mark as present/absent (binarized) before
    // assigning a state (Ernst & Kellis 2012, Nat Methods 9:215-216; ChromHMM
    // BinarizeBed/BinarizeBam). A mark is "present" when its (normalized [0,1])
    // enrichment signal meets the presence call. The exact numeric call is not
    // fixed by the sources; 0.5 is the default midpoint for a normalized signal.
    private const double DefaultMarkPresenceThreshold = 0.5;

    /// <summary>
    /// Predicts the chromatin state at a locus from six histone-modification signals,
    /// using the canonical Roadmap Epigenomics mark-combination signatures.
    /// </summary>
    /// <remarks>
    /// Each mark is binarized to present/absent at <paramref name="presenceThreshold"/>
    /// (ChromHMM binary model; Ernst &amp; Kellis 2012). The state is then assigned from the
    /// present-mark set per the Roadmap 15/18-state definitions:
    /// H3K4me3→active promoter (TssA; Liang 2004); H3K4me1+H3K27ac→active enhancer
    /// (active Enh; Creyghton 2010); H3K4me1 alone→weak/poised enhancer (Enh;
    /// Rada-Iglesias 2018); H3K36me3→transcribed (Tx); H3K27me3→Polycomb-repressed
    /// (ReprPC; Ferrari 2014); H3K9me3→heterochromatin (Het; Nicetto 2019);
    /// H3K4me3+H3K27me3→bivalent promoter (TssBiv); H3K4me1+H3K27me3→bivalent enhancer
    /// (EnhBiv); no mark present→quiescent/low (Quies). Active promoter signature takes
    /// precedence over the enhancer signature when both occur (TSS ranks above Enh).
    /// </remarks>
    /// <param name="h3k4me3">H3K4me3 signal (active-promoter mark).</param>
    /// <param name="h3k4me1">H3K4me1 signal (enhancer mark).</param>
    /// <param name="h3k27ac">H3K27ac signal (active-enhancer/active-regulatory mark).</param>
    /// <param name="h3k36me3">H3K36me3 signal (transcribed gene-body mark).</param>
    /// <param name="h3k27me3">H3K27me3 signal (Polycomb-repression mark).</param>
    /// <param name="h3k9me3">H3K9me3 signal (heterochromatin mark).</param>
    /// <param name="presenceThreshold">
    /// Inclusive presence call: a mark counts as present when its signal is
    /// &gt;= this value. Defaults to <see cref="DefaultMarkPresenceThreshold"/>.
    /// </param>
    /// <returns>The predicted <see cref="ChromatinState"/>.</returns>
    public static ChromatinState PredictChromatinState(
        double h3k4me3,
        double h3k4me1,
        double h3k27ac,
        double h3k36me3,
        double h3k27me3,
        double h3k9me3,
        double presenceThreshold = DefaultMarkPresenceThreshold)
    {
        // Binarize each mark (ChromHMM present/absent model; Ernst & Kellis 2012).
        bool hasK4me3 = h3k4me3 >= presenceThreshold;
        bool hasK4me1 = h3k4me1 >= presenceThreshold;
        bool hasK27ac = h3k27ac >= presenceThreshold;
        bool hasK36me3 = h3k36me3 >= presenceThreshold;
        bool hasK27me3 = h3k27me3 >= presenceThreshold;
        bool hasK9me3 = h3k9me3 >= presenceThreshold;

        // Bivalent signatures: an active promoter/enhancer mark co-occurring with the
        // Polycomb mark H3K27me3 (Roadmap TssBiv / EnhBiv). Checked first because the
        // combination, not the active mark alone, defines the state.
        if (hasK4me3 && hasK27me3)
            return ChromatinState.BivalentPromoter;

        if (hasK4me1 && hasK27me3)
            return ChromatinState.BivalentEnhancer;

        // Active promoter (Roadmap TssA; Liang 2004). H3K4me3 dominates the enhancer
        // mark when both are present (TSS ranks above Enh in the Roadmap model).
        if (hasK4me3)
            return ChromatinState.ActivePromoter;

        // Enhancer: active when accompanied by H3K27ac (Creyghton 2010), otherwise
        // weak/poised (Roadmap Enh; Rada-Iglesias 2018).
        if (hasK4me1)
            return hasK27ac ? ChromatinState.ActiveEnhancer : ChromatinState.WeakEnhancer;

        // Transcribed gene body (Roadmap Tx).
        if (hasK36me3)
            return ChromatinState.Transcribed;

        // Polycomb-repressed (Roadmap ReprPC; Ferrari 2014).
        if (hasK27me3)
            return ChromatinState.Repressed;

        // Constitutive heterochromatin (Roadmap Het; Nicetto 2019).
        if (hasK9me3)
            return ChromatinState.Heterochromatin;

        // No mark present → quiescent/low (Roadmap Quies).
        return ChromatinState.LowSignal;
    }

    /// <summary>
    /// Annotates genomic regions with the chromatin state implied by their single
    /// histone mark, labelling each region by that mark's canonical Roadmap state.
    /// </summary>
    /// <remarks>
    /// A region is labelled by the present-mark rule: a mark below
    /// <paramref name="presenceThreshold"/> yields <see cref="ChromatinState.LowSignal"/>;
    /// otherwise the mark maps to its canonical state per the Roadmap definitions
    /// (see <see cref="PredictChromatinState"/>). H3K27ac is an active-enhancer mark
    /// (Creyghton 2010); H3K9 acetylation is linked to transcriptional activation
    /// (Wang et al. 2008, Nat Genet 40:897-903).
    /// </remarks>
    public static IEnumerable<HistoneModification> AnnotateHistoneModifications(
        IEnumerable<(int Start, int End, string Mark, double Signal)> modifications,
        double presenceThreshold = DefaultMarkPresenceThreshold)
    {
        foreach (var (start, end, mark, signal) in modifications)
        {
            var state = InferStateFromMark(mark, signal, presenceThreshold);

            yield return new HistoneModification(
                Start: start,
                End: end,
                Mark: mark,
                Signal: signal,
                PredictedState: state);
        }
    }

    private static ChromatinState InferStateFromMark(string mark, double signal, double presenceThreshold)
    {
        // Mark not present (below the ChromHMM-style binarization call) → no enrichment.
        if (signal < presenceThreshold)
            return ChromatinState.LowSignal;

        // Single-mark Roadmap signature mapping (Roadmap Epigenomics 15/18-state models).
        return mark.ToUpperInvariant() switch
        {
            "H3K4ME3" => ChromatinState.ActivePromoter,       // TssA (Liang 2004)
            "H3K4ME1" => ChromatinState.WeakEnhancer,         // Enh, weak without H3K27ac (Rada-Iglesias 2018)
            "H3K27AC" => ChromatinState.ActiveEnhancer,       // active enhancer mark (Creyghton 2010)
            "H3K36ME3" => ChromatinState.Transcribed,         // Tx (transcribed gene body)
            "H3K27ME3" => ChromatinState.Repressed,           // ReprPC (Ferrari 2014)
            "H3K9ME3" => ChromatinState.Heterochromatin,      // Het (Nicetto 2019)
            "H3K9AC" => ChromatinState.ActivePromoter,        // H3K9 acetylation → activation (Wang 2008)
            _ => ChromatinState.LowSignal
        };
    }

    #endregion

    #region Chromatin Accessibility

    /// <summary>
    /// Identifies accessible chromatin regions (ATAC-seq like analysis).
    /// </summary>
    public static IEnumerable<AccessibilityRegion> FindAccessibleRegions(
        IEnumerable<(int Position, double Signal)> accessibilitySignal,
        double threshold = 0.5,
        int minWidth = 100,
        int maxGap = 50)
    {
        var signalList = accessibilitySignal.OrderBy(s => s.Position).ToList();

        if (signalList.Count == 0)
            yield break;

        int? regionStart = null;
        int lastPosition = 0;
        double maxSignal = 0;

        foreach (var (pos, signal) in signalList)
        {
            if (signal >= threshold)
            {
                if (regionStart == null)
                {
                    regionStart = pos;
                    maxSignal = signal;
                }
                else if (pos - lastPosition > maxGap)
                {
                    // End previous region, start new one
                    if (lastPosition - regionStart.Value >= minWidth)
                    {
                        yield return new AccessibilityRegion(
                            Start: regionStart.Value,
                            End: lastPosition,
                            AccessibilityScore: maxSignal,
                            PeakType: ClassifyPeakType(maxSignal),
                            NearbyGenes: new List<string>());
                    }
                    regionStart = pos;
                    maxSignal = signal;
                }
                else
                {
                    maxSignal = Math.Max(maxSignal, signal);
                }

                lastPosition = pos;
            }
            else if (regionStart != null)
            {
                if (pos - lastPosition > maxGap)
                {
                    if (lastPosition - regionStart.Value >= minWidth)
                    {
                        yield return new AccessibilityRegion(
                            Start: regionStart.Value,
                            End: lastPosition,
                            AccessibilityScore: maxSignal,
                            PeakType: ClassifyPeakType(maxSignal),
                            NearbyGenes: new List<string>());
                    }
                    regionStart = null;
                    maxSignal = 0;
                }
            }
        }

        // Handle final region
        if (regionStart != null && lastPosition - regionStart.Value >= minWidth)
        {
            yield return new AccessibilityRegion(
                Start: regionStart.Value,
                End: lastPosition,
                AccessibilityScore: maxSignal,
                PeakType: ClassifyPeakType(maxSignal),
                NearbyGenes: new List<string>());
        }
    }

    // Descriptive strength labels for an accessibility peak's normalized [0,1] score.
    // These cutoffs set only the cosmetic PeakType label; they do not affect which
    // regions are detected or returned (that is governed by the caller's threshold).
    private const double StrongPeakScoreCutoff = 0.8;
    private const double ModeratePeakScoreCutoff = 0.5;

    private static string ClassifyPeakType(double score)
    {
        return score switch
        {
            > StrongPeakScoreCutoff => "Strong",
            > ModeratePeakScoreCutoff => "Moderate",
            _ => "Weak"
        };
    }

    #endregion

    #region Imprinting Analysis

    /// <summary>
    /// Predicts imprinted genes based on allele-specific methylation.
    /// </summary>
    public static IEnumerable<ImprintedGene> PredictImprintedGenes(
        IEnumerable<(string GeneId, int Start, int End, double MaternalMethylation, double PaternalMethylation)> genes,
        double minDifference = 0.4)
    {
        foreach (var (geneId, start, end, maternal, paternal) in genes)
        {
            double diff = Math.Abs(maternal - paternal);

            if (diff >= minDifference)
            {
                string origin = maternal > paternal ? "Maternal" : "Paternal";
                double score = diff / (maternal + paternal + 0.01);

                yield return new ImprintedGene(
                    GeneId: geneId,
                    Start: start,
                    End: end,
                    ImprintingScore: Math.Min(score, 1.0),
                    ParentalOrigin: origin,
                    HasDMR: diff > 0.5);
            }
        }
    }

    #endregion

    #region DNA Methylation Age (Epigenetic Clock)

    // Adult-age constant of the Horvath (2013) calibration function F. The relationship
    // between chronological age and DNAm is logarithmic up to this age and linear after.
    // Source: Horvath S (2013), Genome Biology 14:R115, "A transformed version of age"
    // (Additional file 2); reference implementations use adult.age = 20.
    private const double HorvathAdultAge = 20.0;

    /// <summary>
    /// Estimates DNA methylation (epigenetic) age from methylation values at clock CpGs using a
    /// caller-supplied linear predictor (generalised Horvath-style epigenetic clock).
    /// </summary>
    /// <remarks>
    /// Computes the linear predictor Y = intercept + Σ(coefficient_i · β_i) over the CpGs present in
    /// both <paramref name="methylationAtClockCpGs"/> and <paramref name="coefficients"/>, then maps it
    /// to age with the Horvath (2013) inverse calibration F⁻¹ (<see cref="HorvathAntiTransform"/>).
    /// Clock coefficients (e.g. the 353-CpG Horvath set) are large published tables and are NOT bundled;
    /// the caller MUST supply the coefficient table and intercept for the clock they intend to use.
    /// </remarks>
    /// <param name="methylationAtClockCpGs">Methylation β-values keyed by CpG identifier (0..1).</param>
    /// <param name="coefficients">Clock coefficients keyed by CpG identifier. Required.</param>
    /// <param name="intercept">Model intercept added to the weighted sum before the inverse transform.</param>
    /// <returns>Estimated DNAm age in years.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="methylationAtClockCpGs"/> or <paramref name="coefficients"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="coefficients"/> is empty.</exception>
    public static double CalculateEpigeneticAge(
        IReadOnlyDictionary<string, double> methylationAtClockCpGs,
        IReadOnlyDictionary<string, double> coefficients,
        double intercept = 0.0)
    {
        if (methylationAtClockCpGs == null)
            throw new ArgumentNullException(nameof(methylationAtClockCpGs));
        if (coefficients == null)
            throw new ArgumentNullException(nameof(coefficients));
        if (coefficients.Count == 0)
            throw new ArgumentException("Clock coefficient table cannot be empty.", nameof(coefficients));

        // Linear predictor: Y = intercept + Σ coef_i · β_i over the clock CpGs.
        // Source: Horvath (2013) reference R code — predictedAge = anti.trafo(
        //   CoefficientTraining[1] + datMethClock %*% CoefficientTraining[-1]).
        double linearPredictor = intercept;

        foreach (var (cpg, methylation) in methylationAtClockCpGs)
        {
            if (coefficients.TryGetValue(cpg, out double coef))
            {
                linearPredictor += coef * methylation;
            }
        }

        return HorvathAntiTransform(linearPredictor);
    }

    /// <summary>
    /// Horvath (2013) inverse calibration F⁻¹ mapping a transformed-age linear predictor to years.
    /// </summary>
    /// <remarks>
    /// anti.trafo(x) = (1 + adult.age)·exp(x) − 1 for x &lt; 0; (1 + adult.age)·x + adult.age otherwise.
    /// Source: Horvath (2013) reference R code (anti.trafo), with adult.age = 20.
    /// </remarks>
    public static double HorvathAntiTransform(double transformedAge)
    {
        return transformedAge < 0
            ? (1 + HorvathAdultAge) * Math.Exp(transformedAge) - 1
            : (1 + HorvathAdultAge) * transformedAge + HorvathAdultAge;
    }

    #endregion
}
