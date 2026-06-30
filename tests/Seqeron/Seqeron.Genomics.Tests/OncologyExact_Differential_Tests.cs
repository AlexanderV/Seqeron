// 08_DIFFERENTIAL_TESTING rows 88, 92, 115, 116 (Oncology exact metrics). Independent closed-form oracles
// for VAF = alt/total, TMB = mutations/Mb, the CCF formula, and the MATH score (100·1.4826·MAD/median).

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyExact_Differential_Tests
{
    private const double Tol = 1e-9;

    // ---- Row 88: ONCO-VAF-001 — VAF = alt / total ----

    [Test]
    [Category("ONCO-VAF-001")]
    [TestCase(3, 10)]
    [TestCase(7, 10)]
    [TestCase(1, 4)]
    [TestCase(0, 5)]
    [TestCase(0, 0)]   // no coverage -> 0
    public void Vaf_MatchesAltOverTotal(int alt, int total)
    {
        double expected = total == 0 ? 0.0 : (double)alt / total;
        Assert.That(OncologyAnalyzer.CalculateVAF(alt, total), Is.EqualTo(expected).Within(Tol));
    }

    // ---- Row 92: ONCO-TMB-001 — TMB = mutations / Mb ----

    [Test]
    [Category("ONCO-TMB-001")]
    [TestCase(100, 2.0)]
    [TestCase(10, 1.0)]
    [TestCase(3, 1.5)]
    [TestCase(0, 30.0)]
    public void Tmb_MatchesMutationsPerMb(int mutations, double mb)
    {
        Assert.That(OncologyAnalyzer.CalculateTMB(mutations, mb), Is.EqualTo(mutations / mb).Within(Tol));
    }

    // ---- Row 115: ONCO-CCF-001 — CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m), capped at 1 ----

    [Test]
    [Category("ONCO-CCF-001")]
    [TestCase(0.25, 0.5, 2, 1)]
    [TestCase(0.1, 0.8, 3, 1)]
    [TestCase(0.6, 0.5, 2, 1)]   // raw CCF 2.4 -> capped to 1.0
    public void Ccf_MatchesFormula(double vaf, double purity, int cn, int mult)
    {
        double totalDna = purity * cn + 2.0 * (1.0 - purity);
        double rawCcf = vaf * totalDna / (purity * mult);
        double capped = Math.Min(1.0, rawCcf);

        var est = OncologyAnalyzer.EstimateCcf(vaf, purity, cn, mult);
        Assert.That(est.RawCcf, Is.EqualTo(rawCcf).Within(Tol), "raw CCF");
        Assert.That(est.Ccf, Is.EqualTo(capped).Within(Tol), "capped CCF");
    }

    // ---- Row 116: ONCO-HETERO-001 — MATH = 100·1.4826·MAD(VAF)/median(VAF) ----

    private static double Median(double[] v)
    {
        var s = (double[])v.Clone();
        Array.Sort(s);
        int n = s.Length, mid = n / 2;
        return n % 2 == 1 ? s[mid] : (s[mid - 1] + s[mid]) / 2.0;
    }

    [Test]
    [Category("ONCO-HETERO-001")]
    public void Math_MatchesMadOverMedianFormula()
    {
        foreach (var vafs in new[]
        {
            new[] { 0.1, 0.2, 0.3, 0.4, 0.5 },
            new[] { 0.2, 0.4 },
            new[] { 0.3, 0.3, 0.31, 0.6, 0.15 },
        })
        {
            double median = Median(vafs);
            double mad = Median(vafs.Select(x => Math.Abs(x - median)).ToArray());
            double expected = 100.0 * 1.4826 * mad / median;
            Assert.That(OncologyAnalyzer.CalculateITH(vafs), Is.EqualTo(expected).Within(Tol));
        }
    }
}
