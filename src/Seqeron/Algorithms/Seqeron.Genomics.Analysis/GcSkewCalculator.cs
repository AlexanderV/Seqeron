namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Calculates GC skew and related metrics for identifying replication origins and termini.
/// GC skew = (G - C) / (G + C), useful for finding origin of replication in bacterial genomes.
/// </summary>
public static class GcSkewCalculator
{
    #region GC Skew Calculation

    /// <summary>
    /// Calculates GC skew for a single sequence or window.
    /// GC skew = (G - C) / (G + C).
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <returns>GC skew value (-1 to 1).</returns>
    public static double CalculateGcSkew(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateGcSkewCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates GC skew from a raw sequence string.
    /// </summary>
    public static double CalculateGcSkew(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateGcSkewCore(sequence.ToUpperInvariant());
    }

    private static double CalculateGcSkewCore(string seq)
    {
        int gCount = seq.Count(c => c == 'G');
        int cCount = seq.Count(c => c == 'C');
        int total = gCount + cCount;

        return total > 0 ? (double)(gCount - cCount) / total : 0;
    }

    #endregion

    #region Sliding Window GC Skew

    /// <summary>
    /// Calculates GC skew using a sliding window across the sequence.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the sliding window (default: 1000).</param>
    /// <param name="stepSize">Step size for window movement (default: 100).</param>
    /// <returns>Collection of GC skew values with positions.</returns>
    public static IEnumerable<GcSkewPoint> CalculateWindowedGcSkew(
        DnaSequence sequence,
        int windowSize = 1000,
        int stepSize = 100)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        return CalculateWindowedGcSkewCore(sequence.Sequence, windowSize, stepSize);
    }

    /// <summary>
    /// Calculates windowed GC skew from a raw sequence string.
    /// </summary>
    public static IEnumerable<GcSkewPoint> CalculateWindowedGcSkew(
        string sequence,
        int windowSize = 1000,
        int stepSize = 100)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var point in CalculateWindowedGcSkewCore(sequence.ToUpperInvariant(), windowSize, stepSize))
            yield return point;
    }

    private static IEnumerable<GcSkewPoint> CalculateWindowedGcSkewCore(
        string seq,
        int windowSize,
        int stepSize)
    {
        for (int i = 0; i + windowSize <= seq.Length; i += stepSize)
        {
            string window = seq.Substring(i, windowSize);
            double skew = CalculateGcSkewCore(window);

            yield return new GcSkewPoint(
                Position: i + windowSize / 2,
                GcSkew: skew,
                WindowStart: i,
                WindowEnd: i + windowSize - 1);
        }
    }

    #endregion

    #region Cumulative GC Skew

    /// <summary>
    /// Calculates cumulative GC skew across the sequence.
    /// Useful for identifying origin and terminus of replication.
    /// Minimum = origin of replication, Maximum = terminus.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the window for cumulative calculation (default: 1000).</param>
    /// <returns>Collection of cumulative GC skew values.</returns>
    public static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkew(
        DnaSequence sequence,
        int windowSize = 1000)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);

        return CalculateCumulativeGcSkewCore(sequence.Sequence, windowSize);
    }

    /// <summary>
    /// Calculates cumulative GC skew from a raw sequence string.
    /// </summary>
    public static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkew(
        string sequence,
        int windowSize = 1000)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var point in CalculateCumulativeGcSkewCore(sequence.ToUpperInvariant(), windowSize))
            yield return point;
    }

    private static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkewCore(
        string seq,
        int windowSize)
    {
        double cumulative = 0;
        int stepSize = windowSize;

        for (int i = 0; i + windowSize <= seq.Length; i += stepSize)
        {
            string window = seq.Substring(i, windowSize);
            double skew = CalculateGcSkewCore(window);
            cumulative += skew;

            yield return new CumulativeGcSkewPoint(
                Position: i + windowSize / 2,
                GcSkew: skew,
                CumulativeGcSkew: cumulative);
        }
    }

    #endregion

    #region AT Skew Calculation

    /// <summary>
    /// Calculates AT skew for a sequence: (A - T) / (A + T).
    /// The result lies in [-1, 1]; +1 when no T, -1 when no A. Returns 0 when the
    /// sequence contains no A and no T (A + T = 0).
    /// </summary>
    /// <remarks>
    /// AT skew = (A - T) / (A + T) per Charneski et al. (2011) PLoS Genet 7(9):e1002283
    /// and Lobry (1996) Mol Biol Evol 13(5):660-665. Only A and T are counted; all other
    /// symbols are ignored (cf. Biopython Bio.SeqUtils.GC_skew, which ignores non-G/C bases).
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static double CalculateAtSkew(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateAtSkewCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates AT skew from a raw sequence string: (A - T) / (A + T).
    /// Counting is case-insensitive; symbols other than A/T are ignored. Returns 0 for
    /// null/empty input or when A + T = 0.
    /// </summary>
    /// <remarks>
    /// AT skew = (A - T) / (A + T) per Charneski et al. (2011) PLoS Genet 7(9):e1002283
    /// and Lobry (1996) Mol Biol Evol 13(5):660-665.
    /// </remarks>
    public static double CalculateAtSkew(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateAtSkewCore(sequence.ToUpperInvariant());
    }

    private static double CalculateAtSkewCore(string seq)
    {
        int aCount = seq.Count(c => c == 'A');
        int tCount = seq.Count(c => c == 'T');
        int total = aCount + tCount;

        // (A - T) / (A + T); zero denominator (no A and no T) -> 0
        // per Biopython GC_skew ZeroDivisionError -> 0.0 convention.
        return total > 0 ? (double)(aCount - tCount) / total : 0;
    }

    #endregion

    #region Origin/Terminus Prediction

    // Per-nucleotide skew increments for the cumulative skew diagram:
    // G contributes +1, C contributes -1, A/T contribute 0.
    // Grigoriev A (1998) Nucleic Acids Res 26(10):2286-2290; Rosalind BA1F "Minimum Skew Problem".
    private const int GuanineSkewIncrement = +1;
    private const int CytosineSkewIncrement = -1;

    /// <summary>
    /// Predicts the origin and terminus of replication from the cumulative GC-skew diagram.
    /// </summary>
    /// <remarks>
    /// The cumulative skew Skew_i is the running difference (#G − #C) over the prefix
    /// Genome[0..i): Skew_0 = 0 and each base updates the running total by +1 for G, −1 for C,
    /// and 0 for A/T (Grigoriev 1998; Rosalind BA1F). The global <b>minimum</b> of this
    /// diagram marks the replication <b>origin</b> and the global <b>maximum</b> marks the
    /// <b>terminus</b> (Lobry 1996; Grigoriev 1998; GC-skew Wikipedia citing both). Positions
    /// are 0-based prefix indices i ∈ [0, n], so position i refers to the boundary <i>before</i>
    /// base i, matching the Rosalind BA1F convention (its sample returns 53 and 97). When
    /// several positions tie for the extreme value, the first (smallest index) is reported.
    /// </remarks>
    /// <param name="sequence">DNA sequence (typically a complete bacterial chromosome).</param>
    /// <returns>Predicted origin and terminus positions and their cumulative skew values.
    /// <see cref="ReplicationOriginPrediction.IsSignificant"/> is true when the diagram has a
    /// non-zero amplitude (max &gt; min), i.e. a detectable strand-composition asymmetry.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static ReplicationOriginPrediction PredictReplicationOrigin(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return PredictReplicationOriginCore(sequence.Sequence);
    }

    /// <summary>
    /// Predicts the origin and terminus of replication from the cumulative GC-skew diagram of a
    /// raw sequence string. Counting is case-insensitive; only G and C affect the skew.
    /// Returns a zero prediction with <c>IsSignificant = false</c> for null/empty input.
    /// </summary>
    public static ReplicationOriginPrediction PredictReplicationOrigin(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return new ReplicationOriginPrediction(0, 0, 0, 0, false);

        return PredictReplicationOriginCore(sequence.ToUpperInvariant());
    }

    private static ReplicationOriginPrediction PredictReplicationOriginCore(string seq)
    {
        if (seq.Length == 0)
            return new ReplicationOriginPrediction(0, 0, 0, 0, false);

        // Build the cumulative skew diagram and track its first global min/max prefix index.
        int cumulative = 0;          // Skew_0 = 0
        int minSkew = 0, maxSkew = 0;
        int minPos = 0, maxPos = 0;

        for (int i = 0; i < seq.Length; i++)
        {
            char c = seq[i];
            if (c == 'G') cumulative += GuanineSkewIncrement;
            else if (c == 'C') cumulative += CytosineSkewIncrement;
            // A, T and any other symbol leave the cumulative skew unchanged.

            int prefixIndex = i + 1; // Skew_{i+1} is defined after consuming base i.
            if (cumulative < minSkew) { minSkew = cumulative; minPos = prefixIndex; }
            if (cumulative > maxSkew) { maxSkew = cumulative; maxPos = prefixIndex; }
        }

        // Amplitude > 0 means the strands differ in G/C composition (a detectable origin signal).
        bool isSignificant = maxSkew > minSkew;

        return new ReplicationOriginPrediction(
            PredictedOrigin: minPos,
            PredictedTerminus: maxPos,
            OriginSkew: minSkew,
            TerminusSkew: maxSkew,
            IsSignificant: isSignificant);
    }

    #endregion

    #region Comprehensive GC Analysis

    /// <summary>
    /// Gets comprehensive GC analysis including overall GC content, GC skew, AT skew, sliding-window
    /// GC-skew/GC-content profiles, and the compositional variability of those windows.
    /// </summary>
    /// <remarks>
    /// Combines the per-metric definitions used elsewhere in this class:
    /// GC content = (G+C)/(A+T+G+C)·100 (Madigan &amp; Martinko, <i>Brock Biology of Microorganisms</i>,
    /// via Wikipedia "GC-content"); GC skew = (G−C)/(G+C) and AT skew = (A−T)/(A+T) (Lobry 1996;
    /// Charneski et al. 2011). "Variability" is the <b>population</b> variance σ² = Σ(xᵢ−μ)²/N of the
    /// per-window values (the windows form the complete population for this sequence; cf. the
    /// population-variance definition Σ(x−μ)²/N). When the sequence is shorter than the window no full
    /// window exists, so the windowed lists are empty and both window-derived variances are 0; the
    /// overall scalar metrics are still computed over the whole sequence.
    /// </remarks>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Sliding-window length for the profiles (default: 1000).</param>
    /// <param name="stepSize">Step between window starts (default: 100).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static GcAnalysisResult AnalyzeGcContent(
        DnaSequence sequence,
        int windowSize = 1000,
        int stepSize = 100,
        bool fraction = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return AnalyzeGcContentCore(sequence.Sequence, windowSize, stepSize, fraction);
    }

    /// <summary>
    /// Gets comprehensive GC analysis from a raw sequence string. Counting is case-insensitive; only
    /// A/T/G/C contribute to the metrics, other symbols are ignored. Returns a zero result with empty
    /// windowed profiles for null/empty input.
    /// </summary>
    /// <remarks>See <see cref="AnalyzeGcContent(DnaSequence,int,int)"/> for the formulas and conventions.</remarks>
    public static GcAnalysisResult AnalyzeGcContent(
        string sequence,
        int windowSize = 1000,
        int stepSize = 100,
        bool fraction = false)
    {
        if (string.IsNullOrEmpty(sequence))
            return new GcAnalysisResult(0, 0, 0, 0, 0, Array.Empty<GcSkewPoint>(), Array.Empty<GcContentPoint>(), 0);

        return AnalyzeGcContentCore(sequence.ToUpperInvariant(), windowSize, stepSize, fraction);
    }

    private static GcAnalysisResult AnalyzeGcContentCore(string seq, int windowSize, int stepSize, bool fraction)
    {
        var windowedSkew = CalculateWindowedGcSkewCore(seq, windowSize, stepSize).ToList();
        var windowedContent = CalculateWindowedGcContentCore(seq, windowSize, stepSize, fraction).ToList();

        // Opt-in Biopython convention: fraction == true reports GC content in [0,1]
        // (Bio.SeqUtils.gc_fraction); the default (false) keeps the existing percentage [0,100].
        double overallGcContent = CalculateGcContent(seq, fraction);
        double overallGcSkew = CalculateGcSkewCore(seq);
        double overallAtSkew = CalculateAtSkewCore(seq);

        double gcContentVariance = windowedContent.Count > 0
            ? CalculateVariance(windowedContent.Select(w => w.GcContent).ToList())
            : 0;

        double gcSkewVariance = windowedSkew.Count > 0
            ? CalculateVariance(windowedSkew.Select(w => w.GcSkew).ToList())
            : 0;

        return new GcAnalysisResult(
            OverallGcContent: overallGcContent,
            OverallGcSkew: overallGcSkew,
            OverallAtSkew: overallAtSkew,
            GcContentVariance: gcContentVariance,
            GcSkewVariance: gcSkewVariance,
            WindowedGcSkew: windowedSkew,
            WindowedGcContent: windowedContent,
            SequenceLength: seq.Length);
    }

    private static IEnumerable<GcContentPoint> CalculateWindowedGcContentCore(
        string seq,
        int windowSize,
        int stepSize,
        bool fraction = false)
    {
        for (int i = 0; i + windowSize <= seq.Length; i += stepSize)
        {
            string window = seq.Substring(i, windowSize);
            double gcContent = CalculateGcContent(window, fraction);

            yield return new GcContentPoint(
                Position: i + windowSize / 2,
                GcContent: gcContent,
                WindowStart: i,
                WindowEnd: i + windowSize - 1);
        }
    }

    // GC content as a percentage of all bases: GC% = (G+C)/(A+T+G+C)·100
    // per Madigan & Martinko, Brock Biology of Microorganisms (via Wikipedia "GC-content").
    private const double PercentScale = 100.0;

    private static double CalculateGcContent(string seq, bool fraction = false)
    {
        if (string.IsNullOrEmpty(seq)) return 0;
        // Only A/C/G/T are counted; ambiguous/other symbols are ignored in BOTH the
        // numerator and the denominator (Biopython gc_fraction "remove" / Comprehensive
        // GC Analysis §2.2/§3.3). The denominator is A+T+G+C, NOT the raw seq length.
        int gcCount = 0;
        int atgcCount = 0;
        foreach (char c in seq)
        {
            switch (c)
            {
                case 'G' or 'C':
                    gcCount++;
                    atgcCount++;
                    break;
                case 'A' or 'T':
                    atgcCount++;
                    break;
            }
        }

        // Opt-in Biopython convention: fraction == true reports [0,1] (Bio.SeqUtils.gc_fraction);
        // default (false) keeps the percentage GC% = (G+C)/(A+T+G+C)·100.
        double scale = fraction ? 1.0 : PercentScale;
        return atgcCount > 0 ? (double)gcCount / atgcCount * scale : 0;
    }

    // Population variance σ² = Σ(xᵢ−μ)²/N (division by N, not Bessel-corrected N−1):
    // the windows are the complete population for this sequence. Population-variance definition
    // Σ(x−μ)²/N (Cuemath "Population Variance"; worked example {12,13,12,14,19} -> 6.8).
    private static double CalculateVariance(IList<double> values)
    {
        if (values.Count == 0) return 0;
        double mean = values.Average();
        return values.Sum(v => (v - mean) * (v - mean)) / values.Count;
    }

    #endregion
}

/// <summary>
/// A point in GC skew analysis with position and value.
/// </summary>
public readonly record struct GcSkewPoint(
    int Position,
    double GcSkew,
    int WindowStart,
    int WindowEnd);

/// <summary>
/// A point in cumulative GC skew analysis.
/// </summary>
public readonly record struct CumulativeGcSkewPoint(
    int Position,
    double GcSkew,
    double CumulativeGcSkew);

/// <summary>
/// A point in GC content analysis.
/// </summary>
public readonly record struct GcContentPoint(
    int Position,
    double GcContent,
    int WindowStart,
    int WindowEnd);

/// <summary>
/// Predicted origin and terminus of replication.
/// </summary>
public readonly record struct ReplicationOriginPrediction(
    int PredictedOrigin,
    int PredictedTerminus,
    double OriginSkew,
    double TerminusSkew,
    bool IsSignificant);

/// <summary>
/// Comprehensive GC analysis results.
/// </summary>
public sealed record GcAnalysisResult(
    double OverallGcContent,
    double OverallGcSkew,
    double OverallAtSkew,
    double GcContentVariance,
    double GcSkewVariance,
    IReadOnlyList<GcSkewPoint> WindowedGcSkew,
    IReadOnlyList<GcContentPoint> WindowedGcContent,
    int SequenceLength);
