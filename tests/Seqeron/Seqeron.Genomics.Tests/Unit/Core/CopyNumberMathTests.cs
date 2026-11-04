// CopyNumberMath — canonical log2 ↔ absolute copy-number conversions (B24 F28, ONCO-CNA-001).
// Source: CNVkit cnvlib/call.py (_log2_ratio_to_absolute_pure, _log2_ratio_to_absolute, log2_ratios).
// Reference values: CNVkit 0.9.14 (pip) — call._log2_ratio_to_absolute(_pure)(...) printed with repr().

using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class CopyNumberMathTests
{
    [TestCase(0.3, 2.0, 2.4622888266898326)]
    [TestCase(-0.7, 1.5, 0.9233583100086873)]
    [TestCase(1.0, 2.0, 4.0)]
    [TestCase(-1.0, 2.0, 1.0)]
    public void Log2RatioToAbsolute_Pure_MatchesCnvkit(double v, double reference, double expected)
    {
        Assert.That(CopyNumberMath.Log2RatioToAbsolute(v, reference), Is.EqualTo(expected));
    }

    [TestCase(0.3, 2.0, 2.0, 0.45, 3.0273085037551835)]
    [TestCase(-1.0, 2.0, 2.0, 0.7, 0.5714285714285713)]
    [TestCase(1.5, 2.0, 2.0, 0.7, 7.224077499274829)]
    [TestCase(1.0, 2.0, 2.0, 0.3, 8.666666666666668)]
    [TestCase(-0.3, 3.0, 3.0, 0.6, 2.0612619817811773)]
    [TestCase(0.4, 3.0, 3.0, 0.6, 4.597539553864471)]
    [TestCase(0.5, 1.0, 2.0, 0.6, 1.0236892706218252)] // chrX: male reference (r=1), female sample (x=2)
    public void Log2RatioToAbsolute_Purity_MatchesCnvkit(double v, double r, double x, double p, double expected)
    {
        Assert.That(CopyNumberMath.Log2RatioToAbsolute(v, r, x, p), Is.EqualTo(expected).Within(1e-15 * expected));
    }

    [TestCase(-3.0, 0.7)]
    [TestCase(-1.0, 0.5)]
    [TestCase(-0.6, 0.3)]
    public void Log2RatioToAbsolute_BelowContaminationFloor_ClampsToZero(double v, double p)
    {
        // CNVkit #503 clamp: (2·2^v − 2(1−p))/p < 0 ⇒ 0 (CNVkit 0.9.14 returns 0.0 for all three).
        Assert.That(CopyNumberMath.Log2RatioToAbsolute(v, 2.0, 2.0, p), Is.EqualTo(0.0));
        Assert.That(CopyNumberMath.Log2RatioToAbsolute(-2.0, 1.0, 2.0, 0.6), Is.EqualTo(0.0));
    }

    [Test]
    public void Log2RatioToAbsolute_PurityOne_IsPurePath()
    {
        Assert.Multiple(() =>
        {
            // CNVkit: _log2_ratio_to_absolute(0.5, 1, 2, 1.0) == _pure(0.5, 1) == 1.4142135623730951 (expect ignored).
            Assert.That(CopyNumberMath.Log2RatioToAbsolute(0.5, 1.0, 2.0, 1.0), Is.EqualTo(1.4142135623730951));
            foreach (double v in new[] { -3.0, -1.1, 0.0, 0.37, 2.5 })
            {
                Assert.That(CopyNumberMath.Log2RatioToAbsolute(v, 2.0, 2.0, 1.0),
                    Is.EqualTo(CopyNumberMath.Log2RatioToAbsolute(v, 2.0)));
            }
        });
    }

    [Test]
    public void Log2RatioToAbsolute_NaN_Propagates()
    {
        Assert.That(CopyNumberMath.Log2RatioToAbsolute(double.NaN, 2.0, 2.0, 0.45), Is.NaN);
        Assert.That(CopyNumberMath.Log2RatioToAbsolute(double.NaN, 2.0), Is.NaN);
    }

    [TestCase(0.0)]
    [TestCase(-0.1)]
    [TestCase(1.0000001)]
    [TestCase(double.NaN)]
    public void Log2RatioToAbsolute_PurityOutsideUnitInterval_Throws(double p)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => CopyNumberMath.Log2RatioToAbsolute(0.0, 2.0, 2.0, p));
    }

    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void Log2RatioToAbsolute_InvalidCopies_Throws(double copies)
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => CopyNumberMath.Log2RatioToAbsolute(0.0, copies));
            Assert.Throws<ArgumentOutOfRangeException>(() => CopyNumberMath.Log2RatioToAbsolute(0.0, 2.0, copies, 0.5));
        });
    }

    [Test]
    public void AbsoluteToLog2Ratio_MatchesCnvkitLog2Ratios()
    {
        Assert.Multiple(() =>
        {
            // CNVkit do_call(purity=0.7) rescaled log2 column: absolute 0.5714285714285713 → −1.8073549220576044.
            Assert.That(CopyNumberMath.AbsoluteToLog2Ratio(0.5714285714285713, 2.0), Is.EqualTo(-1.8073549220576044).Within(1e-15));
            // Zero copies floored at 1e-3: log2(1e-3) = −9.965784284662087 (CNVkit).
            Assert.That(CopyNumberMath.AbsoluteToLog2Ratio(0.0, 2.0), Is.EqualTo(-9.965784284662087).Within(1e-15));
            Assert.That(CopyNumberMath.AbsoluteToLog2Ratio(6.0, 2.0), Is.EqualTo(1.584962500721156).Within(1e-15));
            Assert.That(CopyNumberMath.AbsoluteToLog2Ratio(double.NaN, 2.0), Is.NaN);
        });
    }

    [TestCase(0.0)]
    [TestCase(double.NaN)]
    public void AbsoluteToLog2Ratio_InvalidPloidy_Throws(double ploidy)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => CopyNumberMath.AbsoluteToLog2Ratio(1.0, ploidy));
    }

    // Shared guards (B24 F66 / E7): the oncology CNA methods call these instead of private copies; same exception type,
    // parameter name and message as before. CNVkit purity_value accepts (0, 1]; ploidy must be finite and > 0.
    [TestCase(0.0)]
    [TestCase(-0.1)]
    [TestCase(1.0000000000000002)]
    [TestCase(double.NaN)]
    public void ValidatePurity_OutsideUnitInterval_Throws(double purity)
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => CopyNumberMath.ValidatePurity(purity));
        Assert.That(ex!.ParamName, Is.EqualTo("purity"));
    }

    [TestCase(0.0)]
    [TestCase(-2.0)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NaN)]
    public void ValidatePloidy_NotFinitePositive_Throws(double ploidy)
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => CopyNumberMath.ValidatePloidy(ploidy));
        Assert.That(ex!.ParamName, Is.EqualTo("ploidy"));
    }

    [Test]
    public void ValidatePurityPloidy_Boundaries_Accepted()
    {
        Assert.DoesNotThrow(() => CopyNumberMath.ValidatePurity(1.0));
        Assert.DoesNotThrow(() => CopyNumberMath.ValidatePurity(double.Epsilon));
        Assert.DoesNotThrow(() => CopyNumberMath.ValidatePloidy(double.Epsilon));
    }
}
