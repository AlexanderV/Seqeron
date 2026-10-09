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

    private static double CalculateGcSkewCore(ReadOnlySpan<char> seq) => CalculateSkewCore(seq, 'G', 'C');

    /// <summary>
    /// Canonical single-pass nucleotide-skew kernel shared by GC skew and AT skew:
    /// (X − Y) / (X + Y) over the (already upper-cased) sequence, counting only
    /// <paramref name="plus"/> (X) and <paramref name="minus"/> (Y); every other symbol is ignored.
    /// Zero denominator (no X and no Y) ⇒ 0, per Biopython <c>GC_skew</c>'s ZeroDivisionError → 0.0.
    /// GC skew = (G−C)/(G+C), AT skew = (A−T)/(A+T) (Lobry 1996; Charneski et al. 2011).
    /// </summary>
    private static double CalculateSkewCore(ReadOnlySpan<char> seq, char plus, char minus)
    {
        int plusCount = 0, minusCount = 0;
        foreach (char c in seq)
        {
            if (c == plus) plusCount++;
            else if (c == minus) minusCount++;
        }

        int total = plusCount + minusCount;
        return total > 0 ? (double)(plusCount - minusCount) / total : 0;
    }

    #endregion

    #region Sliding Window GC Skew

    /// <summary>
    /// Calculates GC skew using a sliding window across the sequence.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the sliding window (default: 1000).</param>
    /// <param name="stepSize">Step size for window movement (default: 100).</param>
    /// <returns>Collection of GC skew values with positions. Only complete windows are emitted
    /// (window starts 0, step, 2·step, … while start + windowSize ≤ length); <c>Position</c> is the
    /// 0-based window centre <c>WindowStart + windowSize / 2</c> (integer division).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1.</exception>
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
    /// Calculates windowed GC skew from a raw sequence string (case-insensitive; only G and C are
    /// counted). Only complete windows are emitted: a trailing partial window is not reported
    /// (cf. SkewIT gcskew.py, which also skips it; Biopython <c>GC_skew</c> instead appends it — use the
    /// <c>includePartialWindow</c> overload for that behaviour).
    /// Returns an empty sequence for null/empty input.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (validated eagerly, as in the
    /// <see cref="DnaSequence"/> overload; a zero step would otherwise never terminate).</exception>
    public static IEnumerable<GcSkewPoint> CalculateWindowedGcSkew(
        string sequence,
        int windowSize = 1000,
        int stepSize = 100)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        if (string.IsNullOrEmpty(sequence))
            return Enumerable.Empty<GcSkewPoint>();

        return CalculateWindowedGcSkewCore(sequence.ToUpperInvariant(), windowSize, stepSize);
    }

    /// <summary>
    /// Calculates windowed GC skew, optionally emitting the trailing partial window(s).
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the sliding window (≥ 1).</param>
    /// <param name="stepSize">Step size for window movement (≥ 1).</param>
    /// <param name="includePartialWindow">When false, identical to
    /// <see cref="CalculateWindowedGcSkew(DnaSequence,int,int)"/> (complete windows only). When true,
    /// every window start <c>i = 0, step, 2·step, …</c> with <c>i &lt; length</c> is emitted; a window
    /// running past the end is truncated to <c>[i, length−1]</c>. With <c>stepSize == windowSize</c> the
    /// skew values equal Biopython 1.88 <c>Bio.SeqUtils.GC_skew(seq, window)</c> exactly
    /// (<c>for i in range(0, len(seq), window): s = seq[i:i+window]</c>; no G/C ⇒ 0.0).
    /// <c>Position</c> is <c>WindowStart + (actual window length) / 2</c> (integer division), which is
    /// <c>WindowStart + windowSize / 2</c> for every complete window.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1.</exception>
    public static IEnumerable<GcSkewPoint> CalculateWindowedGcSkew(
        DnaSequence sequence,
        int windowSize,
        int stepSize,
        bool includePartialWindow)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        return CalculateWindowedGcSkewCore(sequence.Sequence, windowSize, stepSize, includePartialWindow);
    }

    /// <summary>
    /// Calculates windowed GC skew from a raw sequence string (case-insensitive), optionally emitting
    /// the trailing partial window(s); see
    /// <see cref="CalculateWindowedGcSkew(DnaSequence,int,int,bool)"/> for the partial-window contract
    /// (Biopython 1.88 <c>GC_skew</c> parity when <c>stepSize == windowSize</c>).
    /// Returns an empty sequence for null/empty input.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1.</exception>
    public static IEnumerable<GcSkewPoint> CalculateWindowedGcSkew(
        string sequence,
        int windowSize,
        int stepSize,
        bool includePartialWindow)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        if (string.IsNullOrEmpty(sequence))
            return Enumerable.Empty<GcSkewPoint>();

        return CalculateWindowedGcSkewCore(sequence.ToUpperInvariant(), windowSize, stepSize, includePartialWindow);
    }

    private static IEnumerable<GcSkewPoint> CalculateWindowedGcSkewCore(
        string seq,
        int windowSize,
        int stepSize,
        bool includePartialWindow = false) =>
        EnumerateSkewWindows(seq, windowSize, stepSize, includePartialWindow, 'G', 'C')
            .Select(w => new GcSkewPoint(w.Position, w.Skew, w.Start, w.End));

    /// <summary>
    /// Shared window enumerator for windowed/cumulative GC and AT skew. Window starts are
    /// 0, step, 2·step, …; complete windows only (start + windowSize ≤ length) unless
    /// <paramref name="includePartialWindow"/>, in which case every start &lt; length is emitted and
    /// the last window(s) are truncated at the sequence end (Biopython <c>GC_skew</c> slicing
    /// <c>seq[i:i+window]</c>). Position = start + actualLength / 2.
    /// </summary>
    private static IEnumerable<(int Position, double Skew, int Start, int End)> EnumerateSkewWindows(
        string seq,
        int windowSize,
        int stepSize,
        bool includePartialWindow,
        char plus,
        char minus)
    {
        for (int i = 0; includePartialWindow ? i < seq.Length : i + windowSize <= seq.Length; i += stepSize)
        {
            int length = Math.Min(windowSize, seq.Length - i);
            double skew = CalculateSkewCore(seq.AsSpan(i, length), plus, minus);
            yield return (i + length / 2, skew, i, i + length - 1);

            // Guard against int overflow of i + stepSize on huge steps.
            if (stepSize > seq.Length - i)
                yield break;
        }
    }

    #endregion

    #region Cumulative GC Skew

    /// <summary>
    /// Calculates cumulative GC skew across the sequence: the running sum of (G−C)/(G+C) over
    /// adjacent, non-overlapping windows from the sequence start (Grigoriev 1998, NAR 26:2286).
    /// Useful for identifying origin and terminus of replication.
    /// Minimum = origin of replication, Maximum = terminus. Only complete windows are used;
    /// a trailing partial window is not reported (Biopython <c>GC_skew</c> would append it — use the
    /// <c>includePartialWindow</c> overload for that behaviour).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1.</exception>
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
    /// Calculates cumulative GC skew from a raw sequence string (case-insensitive): the running sum
    /// of (G−C)/(G+C) over adjacent, non-overlapping windows (Grigoriev 1998). Only complete windows
    /// are used; a trailing partial window is not reported. Returns an empty sequence for null/empty input.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1
    /// (validated eagerly, as in the <see cref="DnaSequence"/> overload; a zero window would
    /// otherwise never terminate).</exception>
    public static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkew(
        string sequence,
        int windowSize = 1000)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);

        if (string.IsNullOrEmpty(sequence))
            return Enumerable.Empty<CumulativeGcSkewPoint>();

        return CalculateCumulativeGcSkewCore(sequence.ToUpperInvariant(), windowSize);
    }

    /// <summary>
    /// Calculates cumulative GC skew, optionally including the trailing partial window.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the adjacent, non-overlapping windows (≥ 1).</param>
    /// <param name="includePartialWindow">When false, identical to
    /// <see cref="CalculateCumulativeGcSkew(DnaSequence,int)"/>. When true, a trailing window shorter
    /// than <paramref name="windowSize"/> is also emitted and added to the running sum, so the
    /// per-window <c>GcSkew</c> values equal Biopython 1.88 <c>Bio.SeqUtils.GC_skew(seq, window)</c>
    /// exactly and <c>CumulativeGcSkew</c> equals <c>itertools.accumulate</c> of that list.
    /// <c>Position</c> of the partial window is <c>start + (actual length) / 2</c>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1.</exception>
    public static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkew(
        DnaSequence sequence,
        int windowSize,
        bool includePartialWindow)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);

        return CalculateCumulativeGcSkewCore(sequence.Sequence, windowSize, includePartialWindow);
    }

    /// <summary>
    /// Calculates cumulative GC skew from a raw sequence string (case-insensitive), optionally
    /// including the trailing partial window; see
    /// <see cref="CalculateCumulativeGcSkew(DnaSequence,int,bool)"/>. Returns an empty sequence for
    /// null/empty input.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1.</exception>
    public static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkew(
        string sequence,
        int windowSize,
        bool includePartialWindow)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);

        if (string.IsNullOrEmpty(sequence))
            return Enumerable.Empty<CumulativeGcSkewPoint>();

        return CalculateCumulativeGcSkewCore(sequence.ToUpperInvariant(), windowSize, includePartialWindow);
    }

    private static IEnumerable<CumulativeGcSkewPoint> CalculateCumulativeGcSkewCore(
        string seq,
        int windowSize,
        bool includePartialWindow = false)
    {
        double cumulative = 0;

        // Adjacent, non-overlapping windows: step == windowSize (Grigoriev 1998).
        foreach (var w in EnumerateSkewWindows(seq, windowSize, windowSize, includePartialWindow, 'G', 'C'))
        {
            cumulative += w.Skew;
            yield return new CumulativeGcSkewPoint(
                Position: w.Position,
                GcSkew: w.Skew,
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

    // (A - T) / (A + T) via the shared skew kernel; zero denominator (no A and no T) -> 0
    // per Biopython GC_skew ZeroDivisionError -> 0.0 convention.
    private static double CalculateAtSkewCore(ReadOnlySpan<char> seq) => CalculateSkewCore(seq, 'A', 'T');

    /// <summary>
    /// Calculates AT skew (A−T)/(A+T) in sliding windows (Charneski et al. 2011; Lobry 1996),
    /// mirroring <see cref="CalculateWindowedGcSkew(DnaSequence,int,int,bool)"/>: window starts
    /// 0, step, 2·step, …; complete windows only unless <paramref name="includePartialWindow"/>
    /// (then every start &lt; length is emitted, truncated at the sequence end). Only A and T are
    /// counted; a window with no A and no T has skew 0. <c>Position</c> = start + (actual length) / 2.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1.</exception>
    public static IEnumerable<AtSkewPoint> CalculateWindowedAtSkew(
        DnaSequence sequence,
        int windowSize = 1000,
        int stepSize = 100,
        bool includePartialWindow = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        return CalculateWindowedAtSkewCore(sequence.Sequence, windowSize, stepSize, includePartialWindow);
    }

    /// <summary>
    /// Calculates windowed AT skew from a raw sequence string (case-insensitive); see
    /// <see cref="CalculateWindowedAtSkew(DnaSequence,int,int,bool)"/>. Returns an empty sequence
    /// for null/empty input.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (validated eagerly).</exception>
    public static IEnumerable<AtSkewPoint> CalculateWindowedAtSkew(
        string sequence,
        int windowSize = 1000,
        int stepSize = 100,
        bool includePartialWindow = false)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        if (string.IsNullOrEmpty(sequence))
            return Enumerable.Empty<AtSkewPoint>();

        return CalculateWindowedAtSkewCore(sequence.ToUpperInvariant(), windowSize, stepSize, includePartialWindow);
    }

    private static IEnumerable<AtSkewPoint> CalculateWindowedAtSkewCore(
        string seq,
        int windowSize,
        int stepSize,
        bool includePartialWindow) =>
        EnumerateSkewWindows(seq, windowSize, stepSize, includePartialWindow, 'A', 'T')
            .Select(w => new AtSkewPoint(w.Position, w.Skew, w.Start, w.End));

    /// <summary>
    /// Calculates cumulative AT skew: the running sum of (A−T)/(A+T) over adjacent, non-overlapping
    /// windows from the sequence start (the AT analogue of Grigoriev's 1998 cumulative GC skew),
    /// mirroring <see cref="CalculateCumulativeGcSkew(DnaSequence,int,bool)"/>. Only complete
    /// windows are used unless <paramref name="includePartialWindow"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1.</exception>
    public static IEnumerable<CumulativeAtSkewPoint> CalculateCumulativeAtSkew(
        DnaSequence sequence,
        int windowSize = 1000,
        bool includePartialWindow = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);

        return CalculateCumulativeAtSkewCore(sequence.Sequence, windowSize, includePartialWindow);
    }

    /// <summary>
    /// Calculates cumulative AT skew from a raw sequence string (case-insensitive); see
    /// <see cref="CalculateCumulativeAtSkew(DnaSequence,int,bool)"/>. Returns an empty sequence for
    /// null/empty input.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1
    /// (validated eagerly).</exception>
    public static IEnumerable<CumulativeAtSkewPoint> CalculateCumulativeAtSkew(
        string sequence,
        int windowSize = 1000,
        bool includePartialWindow = false)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);

        if (string.IsNullOrEmpty(sequence))
            return Enumerable.Empty<CumulativeAtSkewPoint>();

        return CalculateCumulativeAtSkewCore(sequence.ToUpperInvariant(), windowSize, includePartialWindow);
    }

    private static IEnumerable<CumulativeAtSkewPoint> CalculateCumulativeAtSkewCore(
        string seq,
        int windowSize,
        bool includePartialWindow)
    {
        double cumulative = 0;
        foreach (var w in EnumerateSkewWindows(seq, windowSize, windowSize, includePartialWindow, 'A', 'T'))
        {
            cumulative += w.Skew;
            yield return new CumulativeAtSkewPoint(w.Position, w.Skew, cumulative);
        }
    }

    #endregion

    #region Origin/Terminus Prediction

    // Grigoriev (1998) cumulative skew = running sum of (G−C)/(G+C) over adjacent windows. With a
    // one-base window each window's skew is +1 (G), −1 (C) or 0 (A/T/other), so the diagram is the
    // per-nucleotide running #G − #C of Rosalind BA1F ("Minimum Skew Problem"). The prediction
    // therefore folds over the canonical CalculateCumulativeGcSkewCore with this window size.
    private const int PerNucleotideWindow = 1;

    /// <summary>
    /// Predicts the origin and terminus of replication from the cumulative GC-skew diagram.
    /// </summary>
    /// <remarks>
    /// The cumulative skew Skew_i is the running difference (#G − #C) over the prefix
    /// Genome[0..i): Skew_0 = 0 and each base updates the running total by +1 for G, −1 for C,
    /// and 0 for A/T (Rosalind BA1F). This is Grigoriev's (1998) cumulative skew — the running sum
    /// of (G−C)/(G+C) over adjacent windows — at a one-base window, and is computed by
    /// <see cref="CalculateCumulativeGcSkew(string,int)"/>'s canonical kernel with windowSize = 1
    /// (the per-base resolution needed to reproduce BA1F). The global <b>minimum</b> of this
    /// diagram marks the replication <b>origin</b> and the global <b>maximum</b> marks the
    /// <b>terminus</b> (Lobry 1996; Grigoriev 1998; GC-skew Wikipedia citing both). Positions
    /// are 0-based prefix indices i ∈ [0, n], so position i refers to the boundary <i>before</i>
    /// base i, matching the Rosalind BA1F convention (its sample returns 53 and 97). When
    /// several positions tie for the extreme value, the first (smallest index) is reported.
    /// The input is treated as a linear string read from index 0 (Grigoriev's "arbitrary start");
    /// for a circular chromosome prefix index n denotes the same junction as index 0.
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

        // Skew_0 = 0 (empty prefix) is part of the diagram; the canonical cumulative point for the
        // one-base window starting at i carries Skew_{i+1}. Strict comparisons keep the first
        // (smallest prefix index) extremum on ties.
        double minSkew = 0, maxSkew = 0;
        int minPos = 0, maxPos = 0;
        int prefixIndex = 0;

        foreach (var point in CalculateCumulativeGcSkewCore(seq, PerNucleotideWindow))
        {
            prefixIndex++;
            double cumulative = point.CumulativeGcSkew;
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
    /// GC content (overall and windowed) is computed by the canonical
    /// <see cref="SequenceExtensions.CalculateGcFraction(ReadOnlySpan{char})"/>: G+C over A+C+G+T+U
    /// (U is the RNA counterpart of T — "adenine and uracil in RNA", Wikipedia "GC-content"; Biopython
    /// <c>gc_fraction(seq, "remove")</c> likewise counts U), all other symbols excluded from both counts.
    /// </remarks>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Sliding-window length for the profiles (default: 1000); must be ≥ 1.</param>
    /// <param name="stepSize">Step between window starts (default: 100); must be ≥ 1.</param>
    /// <param name="fraction">When true, GC content is reported in [0,1] (Biopython <c>gc_fraction</c>);
    /// default false reports a percentage in [0,100].</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (a zero step would otherwise never terminate;
    /// Biopython <c>GC_skew</c> likewise rejects window 0).</exception>
    public static GcAnalysisResult AnalyzeGcContent(
        DnaSequence sequence,
        int windowSize = 1000,
        int stepSize = 100,
        bool fraction = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);
        return AnalyzeGcContentCore(sequence.Sequence, windowSize, stepSize, fraction);
    }

    /// <summary>
    /// Gets comprehensive GC analysis from a raw sequence string. Counting is case-insensitive.
    /// GC content counts G+C over A+C+G+T+U (RNA U included in the denominator, as in Biopython
    /// <c>gc_fraction</c>); the GC skew counts only G/C and the AT skew only A/T; every other symbol
    /// is ignored. Returns a zero result with empty windowed profiles for null/empty input.
    /// </summary>
    /// <remarks>See <see cref="AnalyzeGcContent(DnaSequence,int,int,bool)"/> for the formulas and conventions.</remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (validated eagerly, before the null/empty check).</exception>
    public static GcAnalysisResult AnalyzeGcContent(
        string sequence,
        int windowSize = 1000,
        int stepSize = 100,
        bool fraction = false)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

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

        // Population variance Σ(xᵢ−μ)²/N via the canonical StatisticsHelper (0 when no windows).
        double gcContentVariance = StatisticsHelper.PopulationVariance(
            windowedContent.Select(w => w.GcContent).ToList());
        double gcSkewVariance = StatisticsHelper.PopulationVariance(
            windowedSkew.Select(w => w.GcSkew).ToList());

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
            double gcContent = CalculateGcContent(seq.AsSpan(i, windowSize), fraction);

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

    // Delegates to the canonical SequenceExtensions.CalculateGcFraction (case-insensitive;
    // G/C over A/C/G/T/U — U is the RNA counterpart of T, as in Biopython gc_fraction
    // "remove", whose denominator counts ATWU; every other symbol is excluded from both counts).
    private static double CalculateGcContent(ReadOnlySpan<char> seq, bool fraction = false)
    {
        // Opt-in Biopython convention: fraction == true reports [0,1] (Bio.SeqUtils.gc_fraction);
        // default (false) keeps the percentage GC% = fraction·100.
        double gcFraction = seq.CalculateGcFraction();
        return fraction ? gcFraction : gcFraction * PercentScale;
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
/// A point in windowed AT skew analysis: (A−T)/(A+T) of window [WindowStart, WindowEnd].
/// </summary>
public readonly record struct AtSkewPoint(
    int Position,
    double AtSkew,
    int WindowStart,
    int WindowEnd);

/// <summary>
/// A point in cumulative AT skew analysis.
/// </summary>
public readonly record struct CumulativeAtSkewPoint(
    int Position,
    double AtSkew,
    double CumulativeAtSkew);

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
