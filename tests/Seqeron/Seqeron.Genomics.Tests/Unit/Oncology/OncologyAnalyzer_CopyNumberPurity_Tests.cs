// ONCO-CNA-001 / ONCO-CNA-003 — CNVkit purity path (B24 F29).
// Evidence: docs/Evidence/ONCO-CNA-001-Evidence.md, docs/Evidence/ONCO-CNA-003-Evidence.md
// Source: CNVkit cnvlib/call.py do_call(method="threshold", purity=p): absolute_clonal → _log2_ratio_to_absolute
//         (n = max(0,(r·2^v − x(1−p))/p)) → log2_ratios (log2(max(n/ploidy,1e-3))) → absolute_threshold.
// Reference values: CNVkit 0.9.14 do_call on a chr1 CopyNumArray (ploidy 2 unless stated, default thresholds).

using CopyNumberArmSegment = Seqeron.Genomics.Oncology.OncologyAnalyzer.CopyNumberArmSegment;
using CopyNumberState = Seqeron.Genomics.Oncology.OncologyAnalyzer.CopyNumberState;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_CopyNumberPurity_Tests
{
    // (log2, purity, ploidy, CNVkit cn)
    [TestCase(-1.0, 0.7, 2.0, 0)]   // pure → 1; abs 0.571 → log2 −1.807 ⇒ DEL
    [TestCase(-0.2, 0.7, 2.0, 1)]   // pure → 2
    [TestCase(0.15, 0.7, 2.0, 3)]   // pure → 2
    [TestCase(0.8, 0.7, 2.0, 5)]
    [TestCase(1.5, 0.7, 2.0, 8)]
    [TestCase(-0.6, 0.3, 2.0, 0)]   // clamp: abs 0 → log2(1e-3)
    [TestCase(-0.25, 0.3, 2.0, 1)]
    [TestCase(0.0, 0.3, 2.0, 2)]    // abs 2.0000000000000004 → log2 3.2e-16
    [TestCase(0.3, 0.3, 2.0, 4)]
    [TestCase(1.0, 0.3, 2.0, 9)]
    [TestCase(-0.4, 0.5, 2.0, 1)]
    [TestCase(1.0, 0.5, 2.0, 6)]    // abs 6 → log2 1.58496… → ceil(2·2^…) = 6 (round-trip)
    [TestCase(-1.1, 1.0, 2.0, 0)]   // purity 1 = pure path
    [TestCase(-1.0, 1.0, 2.0, 1)]
    [TestCase(1.5, 1.0, 2.0, 6)]
    [TestCase(-1.0, 0.6, 3.0, 0)]
    [TestCase(-0.3, 0.6, 3.0, 1)]
    [TestCase(0.4, 0.6, 3.0, 3)]
    [TestCase(1.0, 0.6, 3.0, 8)]
    public void CallCopyNumber_Purity_MatchesCnvkitDoCall(double log2, double purity, double ploidy, int expected)
    {
        Assert.That(OncologyAnalyzer.CallCopyNumber(log2, null, ploidy, purity), Is.EqualTo(expected));
    }

    [Test]
    public void CallCopyNumber_Purity_CustomThresholds_MatchesCnvkit()
    {
        double[] thr = { -1.5, -0.3, 0.3, 1.0 };
        // CNVkit do_call(purity=0.7, thresholds=(-1.5,-0.3,0.3,1.0)) on [-0.5, 0.1, 0.35] → [1, 2, 3].
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CallCopyNumber(-0.5, thr, 2.0, 0.7), Is.EqualTo(1));
            Assert.That(OncologyAnalyzer.CallCopyNumber(0.1, thr, 2.0, 0.7), Is.EqualTo(2));
            Assert.That(OncologyAnalyzer.CallCopyNumber(0.35, thr, 2.0, 0.7), Is.EqualTo(3));
        });
    }

    [Test]
    public void CallCopyNumber_Purity_DeletionBoundary_MatchesCnvkit()
    {
        // v* = log2((0.7·2·2^−1.1 + 2·0.3)/2) = −0.6744718626824432 maps to rescaled log2 −1.1 (CNVkit).
        // CNVkit: v*−1e-9 → rescaled −1.1000000019 → CN 0; v*+1e-9 → −1.0999999981 → CN 1.
        const double vStar = -0.6744718626824432;
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CallCopyNumber(vStar - 1e-9, null, 2.0, 0.7), Is.EqualTo(0));
            Assert.That(OncologyAnalyzer.CallCopyNumber(vStar + 1e-9, null, 2.0, 0.7), Is.EqualTo(1));
            Assert.That(OncologyAnalyzer.CallCopyNumber(vStar + 1e-9, null, 2.0, 1.0), Is.EqualTo(1));
        });
    }

    [Test]
    public void CallCopyNumber_Purity_NaNAndNegativeInfinity_MatchCnvkit()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CallCopyNumber(double.NaN, null, 2.0, 0.7), Is.EqualTo(2));
            Assert.That(OncologyAnalyzer.CallCopyNumber(double.NegativeInfinity, null, 2.0, 0.7), Is.EqualTo(0));
            Assert.That(OncologyAnalyzer.CallCopyNumber(double.PositiveInfinity, null, 2.0, 0.7), Is.EqualTo(int.MaxValue));
        });
    }

    [Test]
    public void PurityOne_IsIdenticalToPurityLessOverloads()
    {
        foreach (double v in new[] { -3.0, -1.1, -1.0, -0.25, -0.2, 0.0, 0.2, 0.7, 0.8, 1.5, double.NaN })
        {
            Assert.That(OncologyAnalyzer.CallCopyNumber(v, null, 2.0, 1.0), Is.EqualTo(OncologyAnalyzer.CallCopyNumber(v)));
            Assert.That(OncologyAnalyzer.ClassifyCopyNumber(v, null, 2.0, 1.0), Is.EqualTo(OncologyAnalyzer.ClassifyCopyNumber(v)));
            if (!double.IsNaN(v))
            {
                Assert.That(OncologyAnalyzer.Log2RatioToCopyNumber(v, 2.0, 1.0), Is.EqualTo(OncologyAnalyzer.Log2RatioToCopyNumber(v)));
            }
        }
    }

    [Test]
    public void Log2RatioToCopyNumber_Purity_MatchesCnvkitAbsolute()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.Log2RatioToCopyNumber(-1.0, 2.0, 0.7), Is.EqualTo(0.5714285714285713).Within(1e-15));
            Assert.That(OncologyAnalyzer.Log2RatioToCopyNumber(0.8, 2.0, 0.7), Is.EqualTo(4.117431790263566).Within(1e-14));
            Assert.That(OncologyAnalyzer.Log2RatioToCopyNumber(-1.5, 2.0, 0.5), Is.EqualTo(0.0));
            Assert.That(OncologyAnalyzer.Log2RatioToCopyNumber(1.0, 3.0, 0.6), Is.EqualTo(8.0).Within(1e-14));
        });
    }

    [Test]
    public void ClassifyCopyNumber_Purity_ReportsAbsoluteAndState()
    {
        var call = OncologyAnalyzer.ClassifyCopyNumber(0.15, null, 2.0, 0.7);
        Assert.Multiple(() =>
        {
            Assert.That(call.Log2Ratio, Is.EqualTo(0.15));
            Assert.That(call.AbsoluteCopyNumber, Is.EqualTo(2.313055634479557).Within(1e-14));
            Assert.That(call.IntegerCopyNumber, Is.EqualTo(3));
            Assert.That(call.State, Is.EqualTo(CopyNumberState.Gain));
        });

        var calls = OncologyAnalyzer.ClassifyCopyNumbers(new[] { -1.0, 0.0, 1.5 }, null, 2.0, 0.7);
        Assert.That(calls.Select(c => c.IntegerCopyNumber), Is.EqualTo(new[] { 0, 2, 8 }));
    }

    [TestCase(0.0)]
    [TestCase(1.5)]
    [TestCase(double.NaN)]
    public void PurityOverloads_InvalidPurity_Throw(double purity)
    {
        var seg = new CopyNumberArmSegment("9p", 0, 100, 1000, -1.0);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.Log2RatioToCopyNumber(0.0, 2.0, purity));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CallCopyNumber(0.0, null, 2.0, purity));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyCopyNumbers(Array.Empty<double>(), null, 2.0, purity));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.IsHomozygousDeletion(seg, null, 2.0, purity));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.DetectHomozygousDeletions(Array.Empty<CopyNumberArmSegment>(), null, 2.0, purity));
        });
    }

    [Test]
    public void IsHomozygousDeletion_Purity_MatchesCnvkit()
    {
        var segMinus1 = new CopyNumberArmSegment("9p", 0, 100, 1000, -1.0);
        var segMinus04 = new CopyNumberArmSegment("10q", 0, 100, 1000, -0.4);
        Assert.Multiple(() =>
        {
            // Pure: −1.0 → CN 1 (not HD). Purity 0.7 / 0.5: CN 0 (CNVkit).
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(segMinus1), Is.False);
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(segMinus1, null, 2.0, 1.0), Is.False);
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(segMinus1, null, 2.0, 0.7), Is.True);
            // −0.4: CN 1 at purity 0.5, CN 0 at purity 0.3 (CNVkit).
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(segMinus04, null, 2.0, 0.5), Is.False);
            Assert.That(OncologyAnalyzer.IsHomozygousDeletion(segMinus04, null, 2.0, 0.3), Is.True);
        });

        var hd = OncologyAnalyzer.DetectHomozygousDeletions(new[] { segMinus1, segMinus04 }, null, 2.0, 0.3);
        Assert.That(hd, Is.EqualTo(new[] { segMinus1, segMinus04 }));
        Assert.That(OncologyAnalyzer.DetectHomozygousDeletions(new[] { segMinus1, segMinus04 }, null, 2.0, 0.5),
            Is.EqualTo(new[] { segMinus1 }));
    }
}
