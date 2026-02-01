using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics;

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
    public enum ChromatinState
    {
        ActivePromoter,
        WeakPromoter,
        ActiveEnhancer,
        WeakEnhancer,
        Transcribed,
        Repressed,
        Heterochromatin,
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
    /// Finds all potential methylation sites with context.
    /// </summary>
    public static IEnumerable<MethylationSite> FindMethylationSites(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        sequence = sequence.ToUpperInvariant();

        for (int i = 0; i < sequence.Length; i++)
        {
            if (sequence[i] != 'C')
                continue;

            if (i + 1 >= sequence.Length)
                continue;

            char next = sequence[i + 1];
            string context;
            MethylationType type;

            if (next == 'G')
            {
                type = MethylationType.CpG;
                context = i + 2 < sequence.Length ? sequence.Substring(i, 3) : sequence.Substring(i);
            }
            else if (i + 2 < sequence.Length)
            {
                char third = sequence[i + 2];
                context = sequence.Substring(i, 3);

                if (third == 'G')
                {
                    type = MethylationType.CHG;
                }
                else
                {
                    type = MethylationType.CHH;
                }
            }
            else
            {
                continue;
            }

            yield return new MethylationSite(
                Position: i,
                Type: type,
                Context: context,
                MethylationLevel: 0, // Unknown without bisulfite data
                Coverage: 0);
        }
    }

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

        double expected = (c * g) / (double)sequence.Length;
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

        for (int i = 0; i <= sequence.Length - minLength; i += 10)
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
                    yield return (
                        islandStart.Value,
                        islandEnd,
                        CalculateGcContent(island),
                        CalculateCpGObservedExpected(island));
                }
                islandStart = null;
            }
        }

        // Handle island at end
        if (islandStart != null && islandEnd - islandStart.Value >= minLength)
        {
            string island = sequence.Substring(islandStart.Value, islandEnd - islandStart.Value);
            yield return (
                islandStart.Value,
                islandEnd,
                CalculateGcContent(island),
                CalculateCpGObservedExpected(island));
        }
    }

    private static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcFractionFast();

    #endregion

    #region Methylation Analysis

    /// <summary>
    /// Simulates bisulfite conversion of a DNA sequence.
    /// </summary>
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
    /// Calculates methylation level from bisulfite sequencing data.
    /// </summary>
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
    /// Generates a methylation profile for a sequence.
    /// </summary>
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

        double globalMeth = siteList.Average(s => s.MethylationLevel);
        double cpgMeth = cpgSites.Count > 0 ? cpgSites.Average(s => s.MethylationLevel) : 0;
        double chgMeth = chgSites.Count > 0 ? chgSites.Average(s => s.MethylationLevel) : 0;
        double chhMeth = chhSites.Count > 0 ? chhSites.Average(s => s.MethylationLevel) : 0;

        int totalCpG = cpgSites.Count;
        int methylatedCpG = cpgSites.Count(s => s.MethylationLevel >= 0.5);

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

    #endregion

    #region Differentially Methylated Regions

    /// <summary>
    /// Identifies differentially methylated regions between two samples.
    /// </summary>
    public static IEnumerable<DifferentiallyMethylatedRegion> FindDMRs(
        IEnumerable<MethylationSite> sample1,
        IEnumerable<MethylationSite> sample2,
        int windowSize = 1000,
        double minDifference = 0.25,
        int minCpGCount = 3)
    {
        var sites1 = sample1.ToDictionary(s => s.Position, s => s);
        var sites2 = sample2.ToDictionary(s => s.Position, s => s);

        var allPositions = sites1.Keys.Union(sites2.Keys).OrderBy(p => p).ToList();

        if (allPositions.Count == 0)
            yield break;

        int start = allPositions[0];
        var windowSites = new List<(int Pos, double Diff)>();

        foreach (var pos in allPositions)
        {
            // Start new window if needed
            if (pos - start >= windowSize && windowSites.Count >= minCpGCount)
            {
                double meanDiff = windowSites.Average(s => s.Diff);

                if (Math.Abs(meanDiff) >= minDifference)
                {
                    yield return new DifferentiallyMethylatedRegion(
                        Start: start,
                        End: windowSites.Last().Pos,
                        MeanDifference: meanDiff,
                        PValue: CalculateDMRPValue(windowSites.Select(s => s.Diff)),
                        CpGCount: windowSites.Count,
                        Annotation: meanDiff > 0 ? "Hypermethylated" : "Hypomethylated");
                }

                start = pos;
                windowSites.Clear();
            }

            double level1 = sites1.TryGetValue(pos, out var s1) ? s1.MethylationLevel : 0;
            double level2 = sites2.TryGetValue(pos, out var s2) ? s2.MethylationLevel : 0;
            double diff = level2 - level1;

            windowSites.Add((pos, diff));
        }

        // Handle last window
        if (windowSites.Count >= minCpGCount)
        {
            double meanDiff = windowSites.Average(s => s.Diff);

            if (Math.Abs(meanDiff) >= minDifference)
            {
                yield return new DifferentiallyMethylatedRegion(
                    Start: start,
                    End: windowSites.Last().Pos,
                    MeanDifference: meanDiff,
                    PValue: CalculateDMRPValue(windowSites.Select(s => s.Diff)),
                    CpGCount: windowSites.Count,
                    Annotation: meanDiff > 0 ? "Hypermethylated" : "Hypomethylated");
            }
        }
    }

    private static double CalculateDMRPValue(IEnumerable<double> differences)
    {
        var diffs = differences.ToList();
        if (diffs.Count < 2)
            return 1.0;

        double mean = diffs.Average();
        double variance = diffs.Sum(d => (d - mean) * (d - mean)) / (diffs.Count - 1);
        double se = Math.Sqrt(variance / diffs.Count);

        if (se == 0)
            return mean == 0 ? 1.0 : 0.0;

        double t = Math.Abs(mean) / se;

        // Approximate p-value
        return 2 * (1 - StatisticsHelper.NormalCDF(t));
    }

    #endregion

    #region Histone Modification Analysis

    /// <summary>
    /// Predicts chromatin state from histone modification signals.
    /// </summary>
    public static ChromatinState PredictChromatinState(
        double h3k4me3,  // Active promoter
        double h3k4me1,  // Enhancer
        double h3k27ac,  // Active enhancer
        double h3k36me3, // Transcription
        double h3k27me3, // Polycomb repression
        double h3k9me3)  // Heterochromatin
    {
        // Simple decision rules based on histone marks
        if (h3k4me3 > 0.5 && h3k27ac > 0.5)
            return ChromatinState.ActivePromoter;

        if (h3k4me3 > 0.3 && h3k27ac < 0.3)
            return ChromatinState.WeakPromoter;

        if (h3k4me1 > 0.5 && h3k27ac > 0.5)
            return ChromatinState.ActiveEnhancer;

        if (h3k4me1 > 0.3 && h3k27ac < 0.3)
            return ChromatinState.WeakEnhancer;

        if (h3k36me3 > 0.5)
            return ChromatinState.Transcribed;

        if (h3k27me3 > 0.5)
            return ChromatinState.Repressed;

        if (h3k9me3 > 0.5)
            return ChromatinState.Heterochromatin;

        return ChromatinState.LowSignal;
    }

    /// <summary>
    /// Annotates regions with histone modifications.
    /// </summary>
    public static IEnumerable<HistoneModification> AnnotateHistoneModifications(
        IEnumerable<(int Start, int End, string Mark, double Signal)> modifications)
    {
        foreach (var (start, end, mark, signal) in modifications)
        {
            var state = InferStateFromMark(mark, signal);

            yield return new HistoneModification(
                Start: start,
                End: end,
                Mark: mark,
                Signal: signal,
                PredictedState: state);
        }
    }

    private static ChromatinState InferStateFromMark(string mark, double signal)
    {
        if (signal < 0.3)
            return ChromatinState.LowSignal;

        return mark.ToUpperInvariant() switch
        {
            "H3K4ME3" => ChromatinState.ActivePromoter,
            "H3K4ME1" => ChromatinState.WeakEnhancer,
            "H3K27AC" => signal > 0.5 ? ChromatinState.ActiveEnhancer : ChromatinState.WeakEnhancer,
            "H3K36ME3" => ChromatinState.Transcribed,
            "H3K27ME3" => ChromatinState.Repressed,
            "H3K9ME3" => ChromatinState.Heterochromatin,
            "H3K9AC" => ChromatinState.ActivePromoter,
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

    private static string ClassifyPeakType(double score)
    {
        return score switch
        {
            > 0.8 => "Strong",
            > 0.5 => "Moderate",
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

    /// <summary>
    /// Calculates epigenetic age using simplified Horvath clock-like model.
    /// </summary>
    public static double CalculateEpigeneticAge(
        IReadOnlyDictionary<string, double> methylationAtClockCpGs,
        IReadOnlyDictionary<string, double>? coefficients = null)
    {
        if (methylationAtClockCpGs == null || methylationAtClockCpGs.Count == 0)
            return 0;

        coefficients ??= GetDefaultClockCoefficients();

        double age = 0;

        foreach (var (cpg, methylation) in methylationAtClockCpGs)
        {
            if (coefficients.TryGetValue(cpg, out double coef))
            {
                age += coef * methylation;
            }
        }

        // Anti-log transformation (Horvath uses transformed age)
        return age > 0 ? Math.Exp(age) - 1 : 0;
    }

    private static Dictionary<string, double> GetDefaultClockCoefficients()
    {
        // Simplified example coefficients
        return new Dictionary<string, double>
        {
            { "cg00000029", 0.0127 },
            { "cg00000165", -0.0312 },
            { "cg00000236", 0.0089 },
            { "cg00000289", -0.0156 },
            { "cg00000363", 0.0245 }
        };
    }

    #endregion
}
