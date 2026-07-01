using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for GC skew calculations.
/// Verifies skew range, complement negation, and windowed consistency.
///
/// Test Unit: SEQ-GCSKEW-001 (Property Extension), SEQ-ATSKEW-001, SEQ-REPLICATION-001, SEQ-GC-ANALYSIS-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class GcSkewProperties
{
    private static Arbitrary<string> DnaArbitrary(int minLen = 10) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// GC skew is in [-1, 1].
    /// Evidence: (G-C)/(G+C), bounded by definition.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcSkew_InRange()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            double skew = GcSkewCalculator.CalculateGcSkew(seq);
            return (skew >= -1.0 - 0.0001 && skew <= 1.0 + 0.0001)
                .Label($"GcSkew={skew:F4} must be in [-1, 1]");
        });
    }

    /// <summary>
    /// AT skew is in [-1, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AtSkew_InRange()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            double skew = GcSkewCalculator.CalculateAtSkew(seq);
            return (skew >= -1.0 - 0.0001 && skew <= 1.0 + 0.0001)
                .Label($"AtSkew={skew:F4} must be in [-1, 1]");
        });
    }

    /// <summary>
    /// Complement negates GC skew: skew(complement) == -skew(original).
    /// Replacing G↔C inverts the (G-C)/(G+C) ratio.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_NegatesGcSkew()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            // Swap G↔C only to keep (G+C) constant
            string swapped = new(seq.Select(c => c switch
            {
                'G' => 'C',
                'C' => 'G',
                _ => c
            }).ToArray());

            double skew1 = GcSkewCalculator.CalculateGcSkew(seq);
            double skew2 = GcSkewCalculator.CalculateGcSkew(swapped);
            return (Math.Abs(skew1 + skew2) < 0.0001)
                .Label($"skew(orig)={skew1:F4}, skew(swapped)={skew2:F4}, sum={skew1 + skew2:F4}");
        });
    }

    /// <summary>
    /// SEQ-GCSKEW-001 (P): the cumulative GC-skew diagram has exactly one point per nucleotide
    /// when accumulated position-by-position (windowSize = 1) — its length equals the sequence length.
    /// Evidence: the cumulative skew S_i = Σ_{j≤i} sign(G/C) is defined per position
    /// (Grigoriev 1998; Rosalind BA1F minimum-skew diagram).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CumulativeGcSkew_LengthEqualsSequenceLength()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            int count = GcSkewCalculator.CalculateCumulativeGcSkew(seq, windowSize: 1).Count();
            return (count == seq.Length)
                .Label($"cumulative points={count}, expected sequence length={seq.Length}");
        });
    }

    #region SEQ-ATSKEW-001: R: AT skew ∈ [-1,1]; S: complement reverses sign; D: deterministic

    // AT skew = (A−T)/(A+T) ∈ [-1,1]. The DNA complement swaps A↔T, so it negates the AT skew
    // (Lobry 1996). (Range is also covered by AtSkew_InRange above.)

    /// <summary>
    /// INV-1 (S): complementing the strand (A↔T, C↔G) negates the AT skew.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Complement_NegatesAtSkew()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            string complement = new(seq.Select(c => c switch
            {
                'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c
            }).ToArray());
            double s1 = GcSkewCalculator.CalculateAtSkew(seq);
            double s2 = GcSkewCalculator.CalculateAtSkew(complement);
            return (Math.Abs(s1 + s2) < 1e-9).Label($"AtSkew(orig)={s1:F4}, AtSkew(complement)={s2:F4}");
        });
    }

    /// <summary>
    /// INV-2 (D): AT skew is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AtSkew_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
            (GcSkewCalculator.CalculateAtSkew(seq) == GcSkewCalculator.CalculateAtSkew(seq))
                .Label("CalculateAtSkew must be deterministic"));
    }

    #endregion

    #region SEQ-REPLICATION-001: R: origin index ∈ [0,len]; P: at cumulative-skew extremum; D: deterministic

    // PredictReplicationOrigin returns the prefix index of the minimum cumulative GC skew (origin) and
    // maximum (terminus) — Grigoriev (1998); Rosalind BA1F minimum-skew problem.

    /// <summary>
    /// INV-1 (R + P): the predicted origin/terminus are prefix indices in [0,len] sitting at the first
    /// global minimum/maximum of the cumulative GC-skew diagram (independently recomputed).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReplicationOrigin_AtCumulativeSkewExtremum()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            var pred = GcSkewCalculator.PredictReplicationOrigin(seq);

            // Independent cumulative skew diagram (Skew_0 = 0; +1 G, −1 C, 0 otherwise).
            int cum = 0, min = 0, max = 0, minPos = 0, maxPos = 0;
            for (int i = 0; i < seq.Length; i++)
            {
                if (seq[i] == 'G') cum++; else if (seq[i] == 'C') cum--;
                if (cum < min) { min = cum; minPos = i + 1; }
                if (cum > max) { max = cum; maxPos = i + 1; }
            }
            bool ok = pred.PredictedOrigin >= 0 && pred.PredictedOrigin <= seq.Length
                      && pred.PredictedTerminus >= 0 && pred.PredictedTerminus <= seq.Length
                      && pred.PredictedOrigin == minPos && pred.OriginSkew == min
                      && pred.PredictedTerminus == maxPos && pred.TerminusSkew == max;
            return ok.Label($"origin/terminus not at the skew extrema (origin={pred.PredictedOrigin}, min@{minPos})");
        });
    }

    /// <summary>
    /// INV-2 (D): Replication-origin prediction is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReplicationOrigin_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
            (GcSkewCalculator.PredictReplicationOrigin(seq) == GcSkewCalculator.PredictReplicationOrigin(seq))
                .Label("PredictReplicationOrigin must be deterministic"));
    }

    #endregion

    #region SEQ-GC-ANALYSIS-001: R: GC% ∈ [0,100]; P: windows tile sequence; D: deterministic

    // AnalyzeGcContent reports overall GC% and a windowed GC-content/skew profile. With step =
    // windowSize the windows are non-overlapping consecutive tiles.

    /// <summary>
    /// INV-1 (R): overall and per-window GC% are in [0,100], skews in [-1,1], variances ≥ 0, and the
    /// reported sequence length matches the input.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcAnalysis_ValuesInRange()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var r = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 5, stepSize: 5);
            bool ok = r.OverallGcContent is >= 0.0 and <= 100.0
                      && r.OverallGcSkew is >= -1.0 and <= 1.0 && r.OverallAtSkew is >= -1.0 and <= 1.0
                      && r.GcContentVariance >= 0 && r.GcSkewVariance >= 0
                      && r.SequenceLength == seq.Length
                      && r.WindowedGcContent.All(w => w.GcContent is >= 0.0 and <= 100.0);
            return ok.Label($"GC analysis out of range (overall={r.OverallGcContent})");
        });
    }

    /// <summary>
    /// INV-2 (P): with step = windowSize the GC-content windows tile the sequence — consecutive,
    /// non-overlapping, in-bounds windows whose count is ⌊len/windowSize⌋.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcAnalysis_Windows_TileSequence()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            const int w = 5;
            var windows = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: w, stepSize: w).WindowedGcContent;
            bool tiles = windows.Count == seq.Length / w
                         && windows.Select((win, i) => win.WindowStart == i * w && win.WindowEnd == i * w + w - 1).All(x => x)
                         && windows.All(win => win.WindowEnd < seq.Length);
            return tiles.Label($"windows do not tile: count={windows.Count}, expected {seq.Length / w}");
        });
    }

    /// <summary>
    /// INV-3 (D): GC-content analysis is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GcAnalysis_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var a = GcSkewCalculator.AnalyzeGcContent(seq, 5, 5);
            var b = GcSkewCalculator.AnalyzeGcContent(seq, 5, 5);
            return (a.OverallGcContent == b.OverallGcContent
                    && a.WindowedGcContent.Count == b.WindowedGcContent.Count)
                .Label("AnalyzeGcContent must be deterministic");
        });
    }

    #endregion

}
