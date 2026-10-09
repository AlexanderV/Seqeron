// SEQ-GCSKEW-001 — windowed / cumulative GC skew, trailing partial window (includePartialWindow)
// TestSpec: tests/TestSpecs/SEQ-GCSKEW-001.md
// Reference: Biopython 1.88 Bio.SeqUtils.GC_skew(seq, window):
//   for i in range(0, len(seq), window): s = seq[i:i+window]; skew = (g-c)/(g+c) or 0.0 on ZeroDivision.
// Expected values below were produced by running Biopython 1.88 (GC_skew + itertools.accumulate),
// and for stepSize != windowSize by direct Python slicing seq[i:i+window] for i in range(0, len, step).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GcSkewCalculator_PartialWindow_Tests
{
    private const double Tol = 1e-12;

    // Biopython: GC_skew("GGGCACGTGGCCCCATG", 4) = [0.5, 0.0, 0.0, -1.0, 1.0] (last window "G", length 1).
    [Test]
    public void Windowed_IncludePartial_StepEqualsWindow_MatchesBiopythonGcSkew()
    {
        var pts = GcSkewCalculator.CalculateWindowedGcSkew("GGGCACGTGGCCCCATG", 4, 4, includePartialWindow: true).ToList();

        Assert.That(pts.Select(p => p.GcSkew), Is.EqualTo(new[] { 0.5, 0.0, 0.0, -1.0, 1.0 }).Within(Tol));
        Assert.Multiple(() =>
        {
            Assert.That(pts[^1].WindowStart, Is.EqualTo(16));
            Assert.That(pts[^1].WindowEnd, Is.EqualTo(16));
            Assert.That(pts[^1].Position, Is.EqualTo(16), "partial window: start + actualLength/2 = 16 + 1/2");
            Assert.That(pts[0].Position, Is.EqualTo(2), "complete window: start + windowSize/2");
        });
    }

    [Test]
    public void Windowed_Default_DropsTrailingPartialWindow_Unchanged()
    {
        var pts = GcSkewCalculator.CalculateWindowedGcSkew("GGGCACGTGGCCCCATG", 4, 4).ToList();
        Assert.That(pts.Select(p => p.GcSkew), Is.EqualTo(new[] { 0.5, 0.0, 0.0, -1.0 }).Within(Tol));

        var explicitFalse = GcSkewCalculator.CalculateWindowedGcSkew("GGGCACGTGGCCCCATG", 4, 4, includePartialWindow: false).ToList();
        Assert.That(explicitFalse, Is.EqualTo(pts));
    }

    // Biopython: GC_skew("ggGCAcgtNNatGCCCa", 5) = [0.5, 0.0, -0.3333333333333333, -1.0] — lowercase and N ignored.
    [Test]
    public void Windowed_IncludePartial_LowercaseAndAmbiguous_MatchesBiopython()
    {
        var skews = GcSkewCalculator.CalculateWindowedGcSkew("ggGCAcgtNNatGCCCa", 5, 5, true).Select(p => p.GcSkew);
        Assert.That(skews, Is.EqualTo(new[] { 0.5, 0.0, -1.0 / 3.0, -1.0 }).Within(Tol));
    }

    // Random sequences (Python random.seed(2026)); Biopython 1.88 GC_skew values.
    [TestCase("AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTA", 8,
        new[] { -0.3333333333333333, 1.0, 0.6, 0.6, 1.0 },
        new[] { -0.3333333333333333, 0.6666666666666667, 1.2666666666666666, 1.8666666666666667, 2.8666666666666667 })]
    [TestCase("TTATTTGTTACCAATTCTCATTGTGTTTCGGAACTTGCGTTTTAGGTATGTCT", 10,
        new[] { 1.0, -1.0, 0.5, 0.2, 1.0, -1.0 },
        new[] { 1.0, 0.0, 0.5, 0.7, 1.7, 0.7 })]
    public void IncludePartial_RandomSequences_MatchBiopython(string seq, int window, double[] skew, double[] cumulative)
    {
        var w = GcSkewCalculator.CalculateWindowedGcSkew(seq, window, window, true).Select(p => p.GcSkew);
        var c = GcSkewCalculator.CalculateCumulativeGcSkew(seq, window, true).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(w, Is.EqualTo(skew).Within(Tol));
            Assert.That(c.Select(p => p.GcSkew), Is.EqualTo(skew).Within(Tol));
            Assert.That(c.Select(p => p.CumulativeGcSkew), Is.EqualTo(cumulative).Within(1e-9));
        });
    }

    // All-A/T trailing window: Biopython GC_skew("GGGGAT", 4) -> [1.0, 0.0] (ZeroDivision -> 0.0); ours is also 0.
    [Test]
    public void IncludePartial_AllAtWindow_IsZero_LikeBiopython()
    {
        var skews = GcSkewCalculator.CalculateWindowedGcSkew("GGGGAT", 4, 4, true).Select(p => p.GcSkew);
        Assert.That(skews, Is.EqualTo(new[] { 1.0, 0.0 }).Within(Tol));
    }

    // step != window: every start i = 0,3,6,... < 17 is emitted; the last ([15,16]) is truncated.
    // Direct Python: [(0,3,2,0.5),(3,6,5,-1/3),(6,9,8,1.0),(9,12,11,-0.5),(12,15,14,-1.0),(15,16,16,1.0)].
    [Test]
    public void Windowed_IncludePartial_StepLessThanWindow_EmitsTruncatedTail()
    {
        var pts = GcSkewCalculator.CalculateWindowedGcSkew("GGGCACGTGGCCCCATG", 4, 3, true).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(pts.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 3, 6, 9, 12, 15 }));
            Assert.That(pts.Select(p => p.WindowEnd), Is.EqualTo(new[] { 3, 6, 9, 12, 15, 16 }));
            Assert.That(pts.Select(p => p.Position), Is.EqualTo(new[] { 2, 5, 8, 11, 14, 16 }));
            Assert.That(pts.Select(p => p.GcSkew), Is.EqualTo(new[] { 0.5, -1.0 / 3.0, 1.0, -0.5, -1.0, 1.0 }).Within(Tol));
        });
        // Default keeps only complete windows (starts 0..12).
        Assert.That(GcSkewCalculator.CalculateWindowedGcSkew("GGGCACGTGGCCCCATG", 4, 3).Count(), Is.EqualTo(5));
    }

    // Biopython GC_skew("GC", 4) = [0.0]: a sequence shorter than the window yields one partial window.
    [Test]
    public void IncludePartial_SequenceShorterThanWindow_OnePoint()
    {
        Assert.Multiple(() =>
        {
            Assert.That(GcSkewCalculator.CalculateWindowedGcSkew("GC", 4, 4).Count(), Is.EqualTo(0));
            var p = GcSkewCalculator.CalculateWindowedGcSkew("GC", 4, 4, true).Single();
            Assert.That(p, Is.EqualTo(new GcSkewPoint(1, 0.0, 0, 1)));
            var c = GcSkewCalculator.CalculateCumulativeGcSkew("GC", 4, true).Single();
            Assert.That(c, Is.EqualTo(new CumulativeGcSkewPoint(1, 0.0, 0.0)));
        });
    }

    // Biopython: GC_skew(..., 4) = [0.5,0,0,-1,1] -> accumulate = [0.5,0.5,0.5,-0.5,0.5].
    [Test]
    public void Cumulative_IncludePartial_DnaSequence_MatchesBiopythonAccumulate()
    {
        var pts = GcSkewCalculator.CalculateCumulativeGcSkew(new DnaSequence("GGGCACGTGGCCCCATG"), 4, true).ToList();
        Assert.That(pts.Select(p => p.CumulativeGcSkew), Is.EqualTo(new[] { 0.5, 0.5, 0.5, -0.5, 0.5 }).Within(Tol));
        Assert.That(GcSkewCalculator.CalculateCumulativeGcSkew(new DnaSequence("GGGCACGTGGCCCCATG"), 4).Count(), Is.EqualTo(4));
    }

    [Test]
    public void Windowed_IncludePartial_DnaSequence_EqualsStringOverload()
    {
        var a = GcSkewCalculator.CalculateWindowedGcSkew(new DnaSequence("GGGCACGTGGCCCCATG"), 4, 3, true).ToList();
        var b = GcSkewCalculator.CalculateWindowedGcSkew("gggcacgtggccccatg", 4, 3, true).ToList();
        Assert.That(a, Is.EqualTo(b));
    }

    [Test]
    public void IncludePartial_Guards()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.CalculateWindowedGcSkew((DnaSequence)null!, 4, 4, true));
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.CalculateCumulativeGcSkew((DnaSequence)null!, 4, true));
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateWindowedGcSkew("GC", 0, 4, true));
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateWindowedGcSkew("GC", 4, 0, true));
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateCumulativeGcSkew("GC", 0, true));
            Assert.That(GcSkewCalculator.CalculateWindowedGcSkew("", 4, 4, true), Is.Empty);
            Assert.That(GcSkewCalculator.CalculateCumulativeGcSkew((string)null!, 4, true), Is.Empty);
            // Huge step must not overflow the window start.
            Assert.That(GcSkewCalculator.CalculateWindowedGcSkew("GGC", 2, int.MaxValue, true).Count(), Is.EqualTo(1));
            Assert.That(GcSkewCalculator.CalculateWindowedGcSkew("GGC", int.MaxValue, int.MaxValue, true).Single().GcSkew,
                Is.EqualTo(1.0 / 3.0).Within(Tol));
        });
    }
}
