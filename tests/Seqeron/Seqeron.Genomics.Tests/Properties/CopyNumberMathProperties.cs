using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// FIN-B24 heavy tier — property tests for the canonical <see cref="CopyNumberMath"/> conversions (B24 F28) and the
/// CNVkit purity path of the oncology copy-number calls (B24 F29). Source: CNVkit <c>cnvlib/call.py</c>
/// (<c>_log2_ratio_to_absolute_pure</c>: n = r·2^v; <c>_log2_ratio_to_absolute</c>: n = max(0, (r·2^v − x(1−p))/p);
/// <c>log2_ratios</c>: log2(max(n/ploidy, 1e-3))). Oracles are the closed forms restated here.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Oncology")]
public class CopyNumberMathProperties
{
    private static Arbitrary<(double Log2, double Ploidy, double Purity, double Expected)> CnArbitrary() =>
        (from vi in Gen.Choose(-6000, 6000)
         from ploidyTenths in Gen.Choose(10, 80)
         from purityPct in Gen.Choose(1, 100)
         from expectedTenths in Gen.Choose(0, 60)
         select (vi / 1000.0, ploidyTenths / 10.0, purityPct / 100.0, expectedTenths / 10.0)).ToArbitrary();

    /// <summary>
    /// F28 round trip: AbsoluteToLog2Ratio(Log2RatioToAbsolute(v, r), r) = v whenever 2^v ≥ the 1e-3 floor (log2 of
    /// r·2^v / r), and v ↦ n is non-negative and non-decreasing.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PurePath_RoundTrips_AndIsMonotone()
    {
        return Prop.ForAll(CnArbitrary(), t =>
        {
            double n = CopyNumberMath.Log2RatioToAbsolute(t.Log2, t.Ploidy);
            double back = CopyNumberMath.AbsoluteToLog2Ratio(n, t.Ploidy);
            double nNext = CopyNumberMath.Log2RatioToAbsolute(t.Log2 + 0.001, t.Ploidy);
            return (n >= 0 && nNext >= n && Math.Abs(back - t.Log2) <= 1e-12)
                .Label($"n={n}, back={back}, n(v+0.001)={nNext}");
        });
    }

    /// <summary>
    /// F28/F29 purity path: p = 1 is bit-identical to the pure path; for p &lt; 1 the result is ≥ 0, non-decreasing in
    /// v, and inverts the mixture 2^v = (p·n + (1 − p)·x)/r — feeding v = log2((p·n + (1 − p)x)/r) recovers n.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PurityPath_PureLimit_Monotone_InvertsMixture()
    {
        return Prop.ForAll(CnArbitrary(), t =>
        {
            double pure = CopyNumberMath.Log2RatioToAbsolute(t.Log2, t.Ploidy);
            double atOne = CopyNumberMath.Log2RatioToAbsolute(t.Log2, t.Ploidy, t.Expected, 1.0);
            double n = CopyNumberMath.Log2RatioToAbsolute(t.Log2, t.Ploidy, t.Expected, t.Purity);
            double nNext = CopyNumberMath.Log2RatioToAbsolute(t.Log2 + 0.001, t.Ploidy, t.Expected, t.Purity);

            double trueCn = Math.Abs(t.Log2) * 2; // any non-negative tumour copy number in [0, 12]
            double v = Math.Log2((t.Purity * trueCn + (1 - t.Purity) * t.Expected) / t.Ploidy);
            double recovered = double.IsNegativeInfinity(v)
                ? 0.0
                : CopyNumberMath.Log2RatioToAbsolute(v, t.Ploidy, t.Expected, t.Purity);
            return (atOne == pure && n >= 0 && nNext >= n && Math.Abs(recovered - trueCn) <= 1e-9 * Math.Max(1, trueCn) / t.Purity)
                .Label($"pure={pure}, atOne={atOne}, n={n}, nNext={nNext}, trueCn={trueCn}, recovered={recovered}");
        });
    }

    /// <summary>
    /// F29 (CNVkit <c>do_call</c> with <c>--purity</c>): the purity-aware call/classify overloads at p = 1 are identical
    /// to the purity-less ones; the purity-corrected absolute CN and the integer call are non-decreasing in log2; and at
    /// fixed log2 lowering the purity moves the corrected CN away from the ploidy (normal contamination dilutes the
    /// deviation, so the correction amplifies it: n − ψ = (ψ·2^v − ψ)/p).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property OncologyPurityOverloads_PureLimit_Monotone_PurityAmplifiesDeviation()
    {
        return Prop.ForAll(CnArbitrary(), t =>
        {
            var pureCall = OncologyAnalyzer.ClassifyCopyNumber(t.Log2, null, t.Ploidy);
            var atOne = OncologyAnalyzer.ClassifyCopyNumber(t.Log2, null, t.Ploidy, 1.0);
            int call = OncologyAnalyzer.CallCopyNumber(t.Log2, null, t.Ploidy, t.Purity);
            int callNext = OncologyAnalyzer.CallCopyNumber(t.Log2 + 0.01, null, t.Ploidy, t.Purity);
            double abs = OncologyAnalyzer.Log2RatioToCopyNumber(t.Log2, t.Ploidy, t.Purity);
            double absNext = OncologyAnalyzer.Log2RatioToCopyNumber(t.Log2 + 0.01, t.Ploidy, t.Purity);
            double pureAbs = OncologyAnalyzer.Log2RatioToCopyNumber(t.Log2, t.Ploidy);
            bool deviation = abs == 0.0 || Math.Abs(abs - t.Ploidy) >= Math.Abs(pureAbs - t.Ploidy) - 1e-12;
            return (atOne == pureCall && callNext >= call && absNext >= abs && deviation)
                .Label($"call={call}, next={callNext}, abs={abs}, absNext={absNext}, pureAbs={pureAbs}");
        });
    }
}
