// SEQ-ATSKEW-001 — windowed and cumulative AT skew
// Evidence: docs/Evidence/SEQ-ATSKEW-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-ATSKEW-001.md
// Source: AT skew = (A−T)/(A+T), Charneski CA et al. (2011) PLoS Genet 7(9):e1002283; Lobry JR (1996)
// Mol Biol Evol 13(5):660-665. Expected values from a direct numpy computation over
// seq[i:i+window] for i in range(0, len, step) (A/T counted case-insensitively; A+T = 0 ⇒ 0.0).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GcSkewCalculator_WindowedAtSkew_Tests
{
    private const double Tol = 1e-12;

    // numpy: "AAATTTAAGCAT", w=5 step=5 -> (0,4,2,0.2),(5,9,7,1/3),(10,11,11,0.0)
    [Test]
    public void WindowedAtSkew_CompleteWindowsByDefault_AndPartialOnRequest()
    {
        var full = GcSkewCalculator.CalculateWindowedAtSkew("AAATTTAAGCAT", 5, 5).ToList();
        var part = GcSkewCalculator.CalculateWindowedAtSkew("AAATTTAAGCAT", 5, 5, includePartialWindow: true).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(full, Is.EqualTo(new[] { new AtSkewPoint(2, 0.2, 0, 4), new AtSkewPoint(7, 1.0 / 3.0, 5, 9) }));
            Assert.That(part.Select(p => p.AtSkew), Is.EqualTo(new[] { 0.2, 1.0 / 3.0, 0.0 }).Within(Tol));
            Assert.That(part[^1], Is.EqualTo(new AtSkewPoint(11, 0.0, 10, 11)));
        });
    }

    // numpy: "aaTtGCgcAtT", w=4 step=3 -> (0,3,2,0.0),(3,6,5,-1.0),(6,9,8,0.0),(9,10,10,-1.0); GC-only window -> 0.
    [Test]
    public void WindowedAtSkew_LowercaseOverlapping_MatchesDirectComputation()
    {
        var pts = GcSkewCalculator.CalculateWindowedAtSkew("aaTtGCgcAtT", 4, 3, true).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(pts.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 3, 6, 9 }));
            Assert.That(pts.Select(p => p.WindowEnd), Is.EqualTo(new[] { 3, 6, 9, 10 }));
            Assert.That(pts.Select(p => p.Position), Is.EqualTo(new[] { 2, 5, 8, 10 }));
            Assert.That(pts.Select(p => p.AtSkew), Is.EqualTo(new[] { 0.0, -1.0, 0.0, -1.0 }).Within(Tol));
        });
        Assert.That(GcSkewCalculator.CalculateWindowedAtSkew("aaTtGCgcAtT", 4, 3).Count(), Is.EqualTo(3));
    }

    // numpy: "AAATTTAAGCAT", w=4 -> skews [0.5, 0.0, 0.0], cumulative [0.5, 0.5, 0.5];
    //        w=5 (+partial) -> [0.2, 1/3, 0.0], cumulative [0.2, 0.5333333333333333, 0.5333333333333333].
    [Test]
    public void CumulativeAtSkew_MatchesDirectComputation()
    {
        var c4 = GcSkewCalculator.CalculateCumulativeAtSkew(new DnaSequence("AAATTTAAGCAT"), 4).ToList();
        var c5 = GcSkewCalculator.CalculateCumulativeAtSkew("AAATTTAAGCAT", 5, includePartialWindow: true).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(c4.Select(p => p.AtSkew), Is.EqualTo(new[] { 0.5, 0.0, 0.0 }).Within(Tol));
            Assert.That(c4.Select(p => p.CumulativeAtSkew), Is.EqualTo(new[] { 0.5, 0.5, 0.5 }).Within(Tol));
            Assert.That(c4.Select(p => p.Position), Is.EqualTo(new[] { 2, 6, 10 }));
            Assert.That(c5.Select(p => p.CumulativeAtSkew), Is.EqualTo(new[] { 0.2, 0.5333333333333333, 0.5333333333333333 }).Within(1e-9));
            Assert.That(GcSkewCalculator.CalculateCumulativeAtSkew("AAATTTAAGCAT", 5).Count(), Is.EqualTo(2));
        });
    }

    // Each windowed point must equal the canonical whole-sequence CalculateAtSkew of its window (no re-implementation drift).
    [Test]
    public void WindowedAtSkew_EachPointEqualsScalarAtSkewOfWindow()
    {
        const string seq = "ATTAGCAAATCGTTTAAGCTA";
        foreach (var p in GcSkewCalculator.CalculateWindowedAtSkew(new DnaSequence(seq), 6, 4, true))
        {
            string w = seq.Substring(p.WindowStart, p.WindowEnd - p.WindowStart + 1);
            Assert.That(p.AtSkew, Is.EqualTo(GcSkewCalculator.CalculateAtSkew(w)).Within(Tol), w);
        }
    }

    [Test]
    public void AtSkewWindows_Guards()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.CalculateWindowedAtSkew((DnaSequence)null!));
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.CalculateCumulativeAtSkew((DnaSequence)null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateWindowedAtSkew("AT", 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateWindowedAtSkew("AT", 1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateCumulativeAtSkew("AT", 0));
            // Eager validation even for empty input.
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateWindowedAtSkew("", 0, 1));
            Assert.That(GcSkewCalculator.CalculateWindowedAtSkew("", 4, 4, true), Is.Empty);
            Assert.That(GcSkewCalculator.CalculateCumulativeAtSkew((string)null!, 4, true), Is.Empty);
        });
    }
}
