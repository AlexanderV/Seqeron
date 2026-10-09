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
    /// several positions tie for the extreme value, the first (smallest index) is reported; the
    /// full BA1F answer (all minimizers / maximizers) is given by
    /// <see cref="FindMinimumSkewPositions(DnaSequence,bool)"/> /
    /// <see cref="FindMaximumSkewPositions(DnaSequence,bool)"/>.
    /// The input is treated as a linear string read from index 0 (Grigoriev's "arbitrary start");
    /// for a circular chromosome prefix index n denotes the same junction as index 0 (use the
    /// <see cref="PredictReplicationOrigin(DnaSequence,bool)"/> overload to report positions mod n).
    /// </remarks>
    /// <param name="sequence">DNA sequence (typically a complete bacterial chromosome).</param>
    /// <returns>Predicted origin and terminus positions and their cumulative skew values.
    /// <see cref="ReplicationOriginPrediction.IsSignificant"/> is true when the diagram has a
    /// non-zero amplitude (max &gt; min), i.e. a detectable strand-composition asymmetry (no statistical
    /// cutoff; for SkewIT's sourced per-genus test see <see cref="IsSkewIBelowGenusThreshold"/>).</returns>
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

    private static ReplicationOriginPrediction PredictReplicationOriginCore(string seq, bool circular = false)
    {
        if (seq.Length == 0)
            return new ReplicationOriginPrediction(0, 0, 0, 0, false);

        var e = ScanSkewExtrema(seq, circular, collectPositions: false);

        return new ReplicationOriginPrediction(
            PredictedOrigin: e.FirstMin,
            PredictedTerminus: e.FirstMax,
            OriginSkew: e.Min,
            TerminusSkew: e.Max,
            // Amplitude > 0 means the strands differ in G/C composition (a detectable origin signal).
            IsSignificant: e.Max > e.Min);
    }

    /// <summary>
    /// Single pass over the per-nucleotide cumulative skew Skew_0 … Skew_n (Rosalind BA1F), folded
    /// from the canonical <see cref="CalculateCumulativeGcSkewCore"/> at window 1. Skew_0 = 0 (empty
    /// prefix) is part of the diagram; the cumulative point of the one-base window starting at i
    /// carries Skew_{i+1}. Strict comparisons keep the first (smallest prefix index) extremum.
    /// When <paramref name="circular"/>, prefix index n is identified with index 0 (the same
    /// junction of a circular chromosome), so Skew_n is not visited and positions lie in [0, n−1].
    /// </summary>
    private static (double Min, double Max, int FirstMin, int FirstMax, List<int>? MinPositions, List<int>? MaxPositions)
        ScanSkewExtrema(string seq, bool circular, bool collectPositions)
    {
        double minSkew = 0, maxSkew = 0;
        int minPos = 0, maxPos = 0;
        List<int>? minPositions = collectPositions ? new List<int> { 0 } : null;
        List<int>? maxPositions = collectPositions ? new List<int> { 0 } : null;
        int lastPrefixIndex = circular ? seq.Length - 1 : seq.Length;
        int prefixIndex = 0;

        foreach (var point in CalculateCumulativeGcSkewCore(seq, PerNucleotideWindow))
        {
            prefixIndex++;
            if (prefixIndex > lastPrefixIndex)
                break;

            double cumulative = point.CumulativeGcSkew;
            if (cumulative < minSkew)
            {
                minSkew = cumulative; minPos = prefixIndex;
                if (minPositions is not null) { minPositions.Clear(); minPositions.Add(prefixIndex); }
            }
            else if (cumulative == minSkew)
            {
                minPositions?.Add(prefixIndex);
            }

            if (cumulative > maxSkew)
            {
                maxSkew = cumulative; maxPos = prefixIndex;
                if (maxPositions is not null) { maxPositions.Clear(); maxPositions.Add(prefixIndex); }
            }
            else if (cumulative == maxSkew)
            {
                maxPositions?.Add(prefixIndex);
            }
        }

        return (minSkew, maxSkew, minPos, maxPos, minPositions, maxPositions);
    }

    /// <summary>
    /// Minimum Skew Problem (Rosalind BA1F): returns <b>all</b> prefix indices i minimizing the
    /// cumulative skew Skew_i = #G − #C over Genome[0..i) (Skew_0 = 0), in ascending order —
    /// the candidate replication origins (Lobry 1996; Grigoriev 1998). The first element equals
    /// <see cref="ReplicationOriginPrediction.PredictedOrigin"/> of
    /// <see cref="PredictReplicationOrigin(DnaSequence)"/> (linear) or
    /// <see cref="PredictReplicationOrigin(DnaSequence,bool)"/> (circular).
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="circular">When false (default) positions range over [0, n] as in BA1F. When true the
    /// input is a circular chromosome read from index 0: prefix index n denotes the same junction as 0
    /// and is not reported separately, so positions range over [0, n−1]. See
    /// <see cref="PredictReplicationOrigin(DnaSequence,bool)"/> for the rotation behaviour.</param>
    /// <returns>Ascending minimizing prefix indices (never empty: Skew_0 is always a candidate value).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static IReadOnlyList<int> FindMinimumSkewPositions(DnaSequence sequence, bool circular = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindSkewExtremumPositionsCore(sequence.Sequence, circular, minimum: true);
    }

    /// <summary>
    /// Minimum Skew Problem (Rosalind BA1F) on a raw string (case-insensitive); see
    /// <see cref="FindMinimumSkewPositions(DnaSequence,bool)"/>. Null/empty input → <c>[0]</c>
    /// (the diagram is the single value Skew_0 = 0).
    /// </summary>
    public static IReadOnlyList<int> FindMinimumSkewPositions(string sequence, bool circular = false) =>
        FindSkewExtremumPositionsCore(string.IsNullOrEmpty(sequence) ? string.Empty : sequence.ToUpperInvariant(), circular, minimum: true);

    /// <summary>
    /// Returns <b>all</b> prefix indices i maximizing the cumulative skew Skew_i (the candidate
    /// replication termini; Grigoriev 1998), in ascending order — the BA1F definition with the
    /// maximum in place of the minimum. Indexing and <paramref name="circular"/> as in
    /// <see cref="FindMinimumSkewPositions(DnaSequence,bool)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static IReadOnlyList<int> FindMaximumSkewPositions(DnaSequence sequence, bool circular = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindSkewExtremumPositionsCore(sequence.Sequence, circular, minimum: false);
    }

    /// <summary>
    /// All maximizing prefix indices of the cumulative skew of a raw string (case-insensitive); see
    /// <see cref="FindMaximumSkewPositions(DnaSequence,bool)"/>. Null/empty input → <c>[0]</c>.
    /// </summary>
    public static IReadOnlyList<int> FindMaximumSkewPositions(string sequence, bool circular = false) =>
        FindSkewExtremumPositionsCore(string.IsNullOrEmpty(sequence) ? string.Empty : sequence.ToUpperInvariant(), circular, minimum: false);

    private static IReadOnlyList<int> FindSkewExtremumPositionsCore(string seq, bool circular, bool minimum)
    {
        if (seq.Length == 0)
            return new[] { 0 };

        var e = ScanSkewExtrema(seq, circular, collectPositions: true);
        return minimum ? e.MinPositions! : e.MaxPositions!;
    }

    /// <summary>
    /// Predicts the origin and terminus of replication, optionally treating the input as a
    /// circular chromosome.
    /// </summary>
    /// <remarks>
    /// With <paramref name="circular"/> = false this is identical to
    /// <see cref="PredictReplicationOrigin(DnaSequence)"/>. With <paramref name="circular"/> = true
    /// the same prefix walk Skew_0 … Skew_{n−1} is computed from index 0 (Grigoriev's 1998
    /// "arbitrary start" on a circular chromosome), but prefix index n is identified with 0, so
    /// positions are reported modulo n in [0, n−1] and Skew_n is not visited (the junction carries
    /// the value Skew_0 = 0). Let D = Skew_n = total #G − #C:
    /// <list type="bullet">
    /// <item>D = 0: the walk closes on the circle and rotating the input left by r gives
    /// Skew'_j = Skew_{(j+r) mod n} − Skew_r, so every minimizing/maximizing position shifts to
    /// (p − r) mod n — the prediction is rotation-equivariant (ties aside, see
    /// <see cref="FindMinimumSkewPositions(DnaSequence,bool)"/> for the full sets).</item>
    /// <item>D ≠ 0: the circle's cumulative skew is not single-valued — the rotated walk equals the
    /// original one shifted by −Skew_r, plus a step of D for every position past the junction
    /// (Skew'_j = Skew_{j+r−n} + D − Skew_r when j + r ≥ n), so the extrema can move with the start.
    /// No detrending is applied: the result is the walk from the supplied start, so supply the
    /// chromosome from its annotated coordinate 0 (ASM-02).</item>
    /// </list>
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static ReplicationOriginPrediction PredictReplicationOrigin(DnaSequence sequence, bool circular)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return PredictReplicationOriginCore(sequence.Sequence, circular);
    }

    /// <summary>
    /// Raw-string (case-insensitive) form of <see cref="PredictReplicationOrigin(DnaSequence,bool)"/>.
    /// Null/empty input → zero prediction, <c>IsSignificant = false</c>.
    /// </summary>
    public static ReplicationOriginPrediction PredictReplicationOrigin(string sequence, bool circular)
    {
        if (string.IsNullOrEmpty(sequence))
            return new ReplicationOriginPrediction(0, 0, 0, 0, false);

        return PredictReplicationOriginCore(sequence.ToUpperInvariant(), circular);
    }

    /// <summary>
    /// Predicts the origin and terminus from Grigoriev's (1998) <b>windowed</b> cumulative skew
    /// diagram: the running sum of (G−C)/(G+C) over adjacent, non-overlapping windows of
    /// <paramref name="windowSize"/> bases (exactly the points of
    /// <see cref="CalculateCumulativeGcSkew(DnaSequence,int)"/>; complete windows only).
    /// </summary>
    /// <remarks>
    /// Origin = <c>Position</c> (window centre, start + windowSize/2) of the first point with the
    /// minimum cumulative value; terminus = that of the first point with the maximum (Grigoriev
    /// 1998: minimum = origin, maximum = terminus). <c>OriginSkew</c>/<c>TerminusSkew</c> are those
    /// cumulative values; <c>IsSignificant</c> = max &gt; min. Only emitted window points are
    /// candidates (there is no Skew_0 = 0 baseline point), so unlike the per-nucleotide method
    /// <c>OriginSkew</c> may be positive. With windowSize = 1 the values equal the per-nucleotide
    /// diagram without Skew_0. A sequence shorter than one window yields no point → zero prediction,
    /// not significant.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1.</exception>
    public static ReplicationOriginPrediction PredictReplicationOrigin(DnaSequence sequence, int windowSize)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        return PredictReplicationOriginWindowedCore(sequence.Sequence, windowSize);
    }

    /// <summary>
    /// Raw-string (case-insensitive) form of <see cref="PredictReplicationOrigin(DnaSequence,int)"/>.
    /// Null/empty input → zero prediction, <c>IsSignificant = false</c>.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1.</exception>
    public static ReplicationOriginPrediction PredictReplicationOrigin(string sequence, int windowSize)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        if (string.IsNullOrEmpty(sequence))
            return new ReplicationOriginPrediction(0, 0, 0, 0, false);

        return PredictReplicationOriginWindowedCore(sequence.ToUpperInvariant(), windowSize);
    }

    private static ReplicationOriginPrediction PredictReplicationOriginWindowedCore(string seq, int windowSize)
    {
        bool any = false;
        double minSkew = 0, maxSkew = 0;
        int minPos = 0, maxPos = 0;

        foreach (var point in CalculateCumulativeGcSkewCore(seq, windowSize))
        {
            double cumulative = point.CumulativeGcSkew;
            if (!any || cumulative < minSkew) { minSkew = cumulative; minPos = point.Position; }
            if (!any || cumulative > maxSkew) { maxSkew = cumulative; maxPos = point.Position; }
            any = true;
        }

        return new ReplicationOriginPrediction(minPos, maxPos, minSkew, maxSkew, IsSignificant: maxSkew > minSkew);
    }

    #endregion

    #region Skew Index (SkewIT)

    /// <summary>Default SkewIT window (skewi.py <c>-k</c>, 20 kb).</summary>
    public const int DefaultSkewIndexWindow = 20000;

    /// <summary>
    /// Skew Index (SkewI) of Lu &amp; Salzberg (2020, PLoS Comput Biol 16:e1008439, "SkewIT"),
    /// computed exactly as the authors' <c>src/skewi.py</c>: a single [0, 1] measure of how strongly a
    /// complete chromosome shows the two-strand GC-skew pattern (higher = stronger).
    /// </summary>
    /// <remarks>
    /// Steps (verbatim from <c>skewi.py</c>, github.com/jenniferlu717/SkewIT):
    /// <list type="number">
    /// <item>Windows start at 0, k, 2k, … (the last one may be partial); each gets sign(#G − #C) ∈ {+1, 0, −1}.</item>
    /// <item>With L windows, h = round(L/2) and r = round(0.04·L) (Python 3 round, half-to-even), the
    /// window list is doubled (circular) and, for each start i ∈ [0, L) and split t ∈ [i+h−r, i+h+r),
    /// D = |Σ skew[i:t] − Σ skew[t:i+L]|; maxDiff = max D.</item>
    /// <item>SkewI = maxDiff / n · k, capped at 1.0.</item>
    /// </list>
    /// skewi.py writes no SkewI for a sequence when maxDiff ≤ 0 (all windows sign 0, or fewer than 13
    /// windows so that r = 0); this method then returns <c>null</c>. Differences from the CLI, which
    /// are input filters rather than part of the metric: skewi.py by default skips sequences shorter
    /// than 500 kb, without "complete" in the FASTA header, or with "plasmid" in it; its <c>-f</c>
    /// option is parsed but not used by the computation. skewi.py counts only upper-case G/C; this
    /// method upper-cases the input (identical on upper-case FASTA). There is no universal SkewI
    /// cutoff: SkewIT publishes per-genus thresholds (genus mean − 2 SD, genera with ≥ 10 RefSeq-97
    /// genomes; data/RefSeq97_Bacteria_GenusSkewIThresholds.txt, e.g. Escherichia 0.7110); a SkewI
    /// below its genus threshold flags a possible mis-assembly — see <see cref="ParseSkewIGenusThresholds"/>,
    /// <see cref="IsSkewIBelowGenusThreshold"/> and <see cref="IsSkewIBelowThreshold"/>. <see cref="ReplicationOriginPrediction.IsSignificant"/>
    /// is unrelated and unchanged.
    /// </remarks>
    /// <param name="sequence">Complete chromosome sequence.</param>
    /// <param name="windowSize">SkewIT window k (default 20 000).</param>
    /// <returns>SkewI in (0, 1], or null when skewi.py would report none.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1.</exception>
    public static double? CalculateSkewIndex(DnaSequence sequence, int windowSize = DefaultSkewIndexWindow)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        return CalculateSkewIndexCore(sequence.Sequence, windowSize);
    }

    /// <summary>
    /// Raw-string form of <see cref="CalculateSkewIndex(DnaSequence,int)"/> (case-insensitive).
    /// Null/empty input → null.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1.</exception>
    public static double? CalculateSkewIndex(string sequence, int windowSize = DefaultSkewIndexWindow)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        if (string.IsNullOrEmpty(sequence))
            return null;

        return CalculateSkewIndexCore(sequence.ToUpperInvariant(), windowSize);
    }

    private static double? CalculateSkewIndexCore(string seq, int windowSize)
    {
        // Window signs: skewi.py `for i in range(0, len(seq), k)` over seq[i:i+k] (partial tail kept);
        // sign((G−C)/(G+C)) == sign(G−C), and the shared kernel returns 0 for G+C = 0.
        var signs = EnumerateSkewWindows(seq, windowSize, windowSize, includePartialWindow: true, 'G', 'C')
            .Select(w => Math.Sign(w.Skew))
            .ToArray();

        int fullLen = signs.Length;
        int halfLen = (int)Math.Round(fullLen / 2.0, MidpointRounding.ToEven);
        int currRange = (int)Math.Round(fullLen * 0.04, MidpointRounding.ToEven);

        // Prefix sums over the doubled list (skew += skew[:full_len]); Sum(a, b) = Σ skew[a:b]
        // with Python slice semantics for non-negative indices.
        var prefix = new long[2 * fullLen + 1];
        for (int j = 0; j < 2 * fullLen; j++)
            prefix[j + 1] = prefix[j] + signs[j % fullLen];
        long Sum(int a, int b)
        {
            a = Math.Min(a, 2 * fullLen);
            b = Math.Min(b, 2 * fullLen);
            return b > a ? prefix[b] - prefix[a] : 0;
        }

        long maxDiff = -1;
        for (int i = 0; i < fullLen; i++)
        {
            for (int t = i + halfLen - currRange; t < i + halfLen + currRange; t++)
            {
                long diff = Math.Abs(Sum(i, t) - Sum(t, i + fullLen));
                if (diff > maxDiff)
                    maxDiff = diff;
            }
        }

        if (maxDiff <= 0)
            return null;

        double skewI = (double)maxDiff / seq.Length * windowSize;
        return skewI > 1 ? 1.0 : skewI;
    }

    /// <summary>
    /// Parses SkewIT's per-genus SkewI threshold table (Lu &amp; Salzberg 2020,
    /// <c>data/RefSeq97_Bacteria_GenusSkewIThresholds.txt</c>, github.com/jenniferlu717/SkewIT).
    /// </summary>
    /// <remarks>
    /// Format: tab-separated, header row <c>Genus  Num_Genomes  Mean  STDEV  Threshold</c>, one row per
    /// genus named <c>g__&lt;Genus&gt;</c>, CRLF or LF line endings. Threshold = genus mean − 2 SD and is
    /// published only for genera with ≥ 10 RefSeq-97 genomes; rows with an empty threshold column are
    /// skipped. The <c>g__</c> prefix is removed; the returned map is case-insensitive (ordinal).
    /// The table is <b>not bundled</b> with Seqeron: the SkewIT repository is licensed GPL-3.0, Seqeron
    /// MIT, so callers download the file themselves and pass it here.
    /// </remarks>
    /// <param name="reader">Reader positioned at the start of the table.</param>
    /// <returns>Genus → threshold (the RefSeq-97 file gives 160 genera, out of 1 147 rows).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="reader"/> is null.</exception>
    /// <exception cref="FormatException">A row has a malformed threshold or a duplicate genus.</exception>
    public static IReadOnlyDictionary<string, double> ParseSkewIGenusThresholds(TextReader reader)
    {
        ArgumentNullException.ThrowIfNull(reader);
        var map = new Dictionary<string, double>(StringComparer.OrdinalIgnoreCase);
        string? line;
        int lineNo = 0;
        while ((line = reader.ReadLine()) != null)
        {
            lineNo++;
            line = line.TrimEnd('\r');
            if (line.Length == 0 || (lineNo == 1 && line.StartsWith("Genus", StringComparison.OrdinalIgnoreCase)))
                continue;

            string[] cols = line.Split('\t');
            if (cols.Length < 5 || string.IsNullOrWhiteSpace(cols[4]))
                continue; // genus with < 10 genomes: no published threshold

            string genus = NormalizeSkewIGenus(cols[0]);
            if (!double.TryParse(cols[4].Trim(), System.Globalization.NumberStyles.Float,
                    System.Globalization.CultureInfo.InvariantCulture, out double threshold)
                || !double.IsFinite(threshold))
                throw new FormatException($"Line {lineNo}: malformed SkewI threshold '{cols[4]}'.");
            if (genus.Length == 0 || !map.TryAdd(genus, threshold))
                throw new FormatException($"Line {lineNo}: empty or duplicate genus '{cols[0]}'.");
        }

        return new System.Collections.ObjectModel.ReadOnlyDictionary<string, double>(map);
    }

    /// <summary>
    /// Looks up a genus in a SkewIT threshold table (see <see cref="ParseSkewIGenusThresholds"/>):
    /// exact genus name, case-insensitive (ordinal); an optional <c>g__</c> prefix and surrounding
    /// white space are ignored.
    /// </summary>
    /// <returns>True when the table holds a threshold for the genus.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="thresholds"/> is null.</exception>
    public static bool TryGetSkewIThreshold(
        IReadOnlyDictionary<string, double> thresholds, string genus, out double threshold)
    {
        ArgumentNullException.ThrowIfNull(thresholds);
        threshold = 0;
        if (string.IsNullOrWhiteSpace(genus))
            return false;

        string key = NormalizeSkewIGenus(genus);
        if (thresholds.TryGetValue(key, out threshold))
            return true;

        // Caller-built maps may use a case-sensitive comparer.
        foreach (var kv in thresholds)
        {
            if (string.Equals(NormalizeSkewIGenus(kv.Key), key, StringComparison.OrdinalIgnoreCase))
            {
                threshold = kv.Value;
                return true;
            }
        }

        threshold = 0;
        return false;
    }

    /// <summary>
    /// SkewIT significance rule (Lu &amp; Salzberg 2020): a genome whose SkewI is <b>below</b> the
    /// threshold (strict &lt;) is flagged as atypical / potentially mis-assembled.
    /// </summary>
    /// <param name="sequence">Complete chromosome sequence.</param>
    /// <param name="threshold">SkewI threshold, e.g. a genus value from SkewIT's table (Escherichia 0.7110).</param>
    /// <param name="windowSize">SkewIT window k. SkewIT's tables use its default k = 20 000; other
    /// values give a SkewI not comparable to the published thresholds.</param>
    /// <returns>True/false, or null when <see cref="CalculateSkewIndex(string,int)"/> is null.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="threshold"/> is not finite, or
    /// <paramref name="windowSize"/> is less than 1.</exception>
    public static bool? IsSkewIBelowThreshold(string sequence, double threshold, int windowSize = DefaultSkewIndexWindow)
    {
        if (!double.IsFinite(threshold))
            throw new ArgumentOutOfRangeException(nameof(threshold), threshold, "Threshold must be finite.");

        double? skewI = CalculateSkewIndex(sequence, windowSize);
        return skewI.HasValue ? skewI.Value < threshold : null;
    }

    /// <summary>
    /// <see cref="IsSkewIBelowThreshold"/> with the threshold of <paramref name="genus"/> taken from a
    /// SkewIT table (<see cref="ParseSkewIGenusThresholds"/>, <see cref="TryGetSkewIThreshold"/>).
    /// </summary>
    /// <returns>Null when the genus has no threshold in the table or SkewI is null.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="thresholds"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is less than 1.</exception>
    public static bool? IsSkewIBelowGenusThreshold(
        string sequence, string genus, IReadOnlyDictionary<string, double> thresholds,
        int windowSize = DefaultSkewIndexWindow)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        if (!TryGetSkewIThreshold(thresholds, genus, out double threshold))
            return null;

        return IsSkewIBelowThreshold(sequence, threshold, windowSize);
    }

    private static string NormalizeSkewIGenus(string genus)
    {
        string g = genus.Trim();
        return g.StartsWith("g__", StringComparison.OrdinalIgnoreCase) ? g[3..] : g;
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

    /// <summary>
    /// Comprehensive GC analysis with Biopython <c>gc_fraction(seq, ambiguous=…)</c> IUPAC-ambiguity
    /// handling for the GC content (overall and windowed): both are scored by the canonical
    /// <see cref="SequenceExtensions.CalculateGcFraction(ReadOnlySpan{char},SequenceExtensions.GcAmbiguityMode)"/>.
    /// <c>Remove</c> (Biopython default): S counts as GC, A/C/G/T/U/S/W form the denominator, other
    /// codes are excluded (<c>"GGSW"</c> → 0.75, where the mode-less overloads give 2/2 = 1.0);
    /// <c>Ignore</c>: denominator = length; <c>Weighted</c>: each ambiguity code adds its mean GC
    /// (N = 0.5). The GC/AT skews and their window profile are unchanged — they count only G/C and A/T.
    /// </summary>
    /// <remarks>
    /// Same windows, guards, variance definition and null/empty result as
    /// <see cref="AnalyzeGcContent(string,int,int,bool)"/> (windowed GC = the windows of
    /// <see cref="CalculateWindowedGcContent(string,int,int,bool,SequenceExtensions.GcAmbiguityMode)"/>).
    /// Case is folded ASCII-only for the GC count (as Biopython's literal <c>"CGScgs"</c> counting), so a
    /// non-ASCII letter such as U+017F is never counted as S. No <see cref="DnaSequence"/> overload:
    /// a <see cref="DnaSequence"/> holds only A/C/G/T, for which all three modes equal the default.
    /// </remarks>
    /// <param name="sequence">Nucleotide sequence, IUPAC codes allowed.</param>
    /// <param name="windowSize">Sliding-window length (≥ 1).</param>
    /// <param name="stepSize">Step between window starts (≥ 1).</param>
    /// <param name="fraction">True → GC content in [0,1]; false → percentage [0,100].</param>
    /// <param name="ambiguityMode">Biopython <c>ambiguous</c> mode.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (validated eagerly, before the null/empty check).</exception>
    public static GcAnalysisResult AnalyzeGcContent(
        string sequence,
        int windowSize,
        int stepSize,
        bool fraction,
        SequenceExtensions.GcAmbiguityMode ambiguityMode)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        if (string.IsNullOrEmpty(sequence))
            return new GcAnalysisResult(0, 0, 0, 0, 0, Array.Empty<GcSkewPoint>(), Array.Empty<GcContentPoint>(), 0);

        return AnalyzeGcContentCore(sequence.ToUpperInvariant(), windowSize, stepSize, fraction, ambiguityMode, sequence);
    }

    // seq: upper-cased input for the G/C/A/T skew counts. mode == null → default A/C/G/T/U GC counting
    // on seq; otherwise Biopython gc_fraction ambiguity handling on gcSeq (the raw input; the canonical
    // mode counter folds case ASCII-only).
    private static GcAnalysisResult AnalyzeGcContentCore(
        string seq, int windowSize, int stepSize, bool fraction,
        SequenceExtensions.GcAmbiguityMode? mode = null, string? gcSeq = null)
    {
        gcSeq ??= seq;
        var windowedSkew = CalculateWindowedGcSkewCore(seq, windowSize, stepSize).ToList();
        var windowedContent = CalculateWindowedGcContentCore(gcSeq, windowSize, stepSize, fraction, mode).ToList();

        // Opt-in Biopython convention: fraction == true reports GC content in [0,1]
        // (Bio.SeqUtils.gc_fraction); the default (false) keeps the existing percentage [0,100].
        double overallGcContent = CalculateGcContent(gcSeq, fraction, mode);
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

    /// <summary>
    /// Sliding-window GC-content profile: the GC content of every complete window of length
    /// <paramref name="windowSize"/> starting at 0, <paramref name="stepSize"/>, 2·<paramref name="stepSize"/>, …
    /// (while start + windowSize ≤ length; a trailing partial window is not reported). This is the
    /// single sliding-GC driver behind <see cref="AnalyzeGcContent(string,int,int,bool)"/>'s
    /// <see cref="GcAnalysisResult.WindowedGcContent"/>.
    /// </summary>
    /// <remarks>
    /// Per-window GC is the canonical <see cref="SequenceExtensions.CalculateGcFraction(ReadOnlySpan{char})"/>
    /// (ASCII case-insensitive; G+C over A+C+G+T+U; every other symbol excluded from both counts; 0 for a
    /// window with no valid base) — for an A/C/G/T/U window this equals Biopython 1.88
    /// <c>Bio.SeqUtils.gc_fraction(seq[i:i+windowSize])</c> (default <c>ambiguous="remove"</c>).
    /// <c>Position</c> is the 0-based window centre <c>WindowStart + windowSize / 2</c> (integer division),
    /// <c>WindowEnd = WindowStart + windowSize − 1</c> (inclusive), matching
    /// <see cref="CalculateWindowedGcSkew(string,int,int)"/>. Returns an empty sequence for null/empty
    /// input or when <paramref name="windowSize"/> exceeds the length.
    /// </remarks>
    /// <param name="sequence">Nucleotide sequence (DNA or RNA).</param>
    /// <param name="windowSize">Window length (≥ 1; default 1000).</param>
    /// <param name="stepSize">Step between window starts (≥ 1; default 100).</param>
    /// <param name="fraction">When true, GC content is reported in [0,1] (Biopython <c>gc_fraction</c>);
    /// default false reports a percentage in [0,100].</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (validated eagerly; a zero step would never terminate).</exception>
    public static IEnumerable<GcContentPoint> CalculateWindowedGcContent(
        string sequence,
        int windowSize = 1000,
        int stepSize = 100,
        bool fraction = false)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        if (string.IsNullOrEmpty(sequence))
            return Enumerable.Empty<GcContentPoint>();

        return CalculateWindowedGcContentCore(sequence, windowSize, stepSize, fraction);
    }

    /// <summary>
    /// Sliding-window GC-content profile with Biopython <c>gc_fraction(window, ambiguous=…)</c>
    /// IUPAC-ambiguity handling: each complete window is scored by the canonical
    /// <see cref="SequenceExtensions.CalculateGcFraction(ReadOnlySpan{char},SequenceExtensions.GcAmbiguityMode)"/>
    /// (<c>Remove</c>: S counts as GC, S/W in the denominator, other codes excluded; <c>Ignore</c>:
    /// denominator = window length; <c>Weighted</c>: ambiguity codes add their mean GC, e.g. N = 0.5).
    /// Window/step/position semantics and guards are those of
    /// <see cref="CalculateWindowedGcContent(string,int,int,bool)"/>.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1.</exception>
    public static IEnumerable<GcContentPoint> CalculateWindowedGcContent(
        string sequence,
        int windowSize,
        int stepSize,
        bool fraction,
        SequenceExtensions.GcAmbiguityMode ambiguityMode)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        if (string.IsNullOrEmpty(sequence))
            return Enumerable.Empty<GcContentPoint>();

        return CalculateWindowedGcContentCore(sequence, windowSize, stepSize, fraction, ambiguityMode);
    }

    /// <summary>
    /// <see cref="DnaSequence"/> overload of <see cref="CalculateWindowedGcContent(string,int,int,bool)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1.</exception>
    public static IEnumerable<GcContentPoint> CalculateWindowedGcContent(
        DnaSequence sequence,
        int windowSize = 1000,
        int stepSize = 100,
        bool fraction = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        return CalculateWindowedGcContentCore(sequence.Sequence, windowSize, stepSize, fraction);
    }

    // Single sliding-GC driver (complete windows only). mode == null → default A/C/G/T/U counting;
    // otherwise Biopython gc_fraction ambiguity handling.
    private static IEnumerable<GcContentPoint> CalculateWindowedGcContentCore(
        string seq,
        int windowSize,
        int stepSize,
        bool fraction = false,
        SequenceExtensions.GcAmbiguityMode? mode = null)
    {
        double scale = fraction ? 1.0 : PercentScale;
        for (int i = 0; windowSize <= seq.Length - i; i += stepSize)
        {
            var window = seq.AsSpan(i, windowSize);
            double gcFraction = mode is { } m ? window.CalculateGcFraction(m) : window.CalculateGcFraction();

            yield return new GcContentPoint(
                Position: i + windowSize / 2,
                GcContent: gcFraction * scale,
                WindowStart: i,
                WindowEnd: i + windowSize - 1);

            // Guard against int overflow of i + stepSize on huge steps.
            if (stepSize > seq.Length - i)
                yield break;
        }
    }

    // GC content as a percentage of all bases: GC% = (G+C)/(A+T+G+C)·100
    // per Madigan & Martinko, Brock Biology of Microorganisms (via Wikipedia "GC-content").
    private const double PercentScale = 100.0;

    // Delegates to the canonical SequenceExtensions.CalculateGcFraction (case-insensitive;
    // G/C over A/C/G/T/U — U is the RNA counterpart of T, as in Biopython gc_fraction
    // "remove", whose denominator counts ATWU; every other symbol is excluded from both counts).
    // mode != null → Biopython gc_fraction ambiguity handling (canonical CalculateGcFraction(mode)).
    private static double CalculateGcContent(
        ReadOnlySpan<char> seq, bool fraction = false, SequenceExtensions.GcAmbiguityMode? mode = null)
    {
        // Opt-in Biopython convention: fraction == true reports [0,1] (Bio.SeqUtils.gc_fraction);
        // default (false) keeps the percentage GC% = fraction·100.
        double gcFraction = mode is { } m ? seq.CalculateGcFraction(m) : seq.CalculateGcFraction();
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
