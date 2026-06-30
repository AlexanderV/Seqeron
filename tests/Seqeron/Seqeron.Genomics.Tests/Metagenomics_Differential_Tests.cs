// 08_DIFFERENTIAL_TESTING rows 55, 56 (Metagenomics diversity). Independent oracles: the Shannon/Simpson/
// Pielou alpha-diversity formulas (with the log2-vs-ln proportionality relation) and the Bray-Curtis /
// Jaccard beta-diversity formulas, each cross-checked against identity/disjoint anchors. (Rows 53/54/57 —
// classification / profiling / binning — are deferred: loose "correlated" comparisons.)

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class Metagenomics_Differential_Tests
{
    private const double Tol = 1e-12;

    // ---- Row 55: META-ALPHA-001 — Shannon (ln) / Simpson / Pielou + log2-proportionality ----

    [Test]
    [Category("META-ALPHA-001")]
    public void AlphaDiversity_MatchesClosedFormAndLog2Proportionality()
    {
        var abund = new Dictionary<string, double> { ["a"] = 10, ["b"] = 20, ["c"] = 30, ["d"] = 40 };
        double sum = abund.Values.Sum();
        var p = abund.Values.Select(v => v / sum).ToList();

        double shannonLn = -p.Sum(pi => pi * Math.Log(pi));
        double simpson = p.Sum(pi => pi * pi);
        double pielou = shannonLn / Math.Log(abund.Count);
        double shannonLog2 = -p.Sum(pi => pi * Math.Log2(pi));

        var ad = MetagenomicsAnalyzer.CalculateAlphaDiversity(abund);
        Assert.That(ad.ShannonIndex, Is.EqualTo(shannonLn).Within(Tol), "Shannon (natural log)");
        Assert.That(ad.SimpsonIndex, Is.EqualTo(simpson).Within(Tol), "Simpson");
        Assert.That(ad.PielouEvenness, Is.EqualTo(pielou).Within(Tol), "Pielou evenness");
        Assert.That(ad.ObservedSpecies, Is.EqualTo(4));

        // The ln- and log2-based Shannon indices are proportional by ln(2).
        Assert.That(ad.ShannonIndex, Is.EqualTo(shannonLog2 * Math.Log(2)).Within(1e-9), "ln = log2 * ln(2)");
    }

    // ---- Row 56: META-BETA-001 — Bray-Curtis / Jaccard vs independent formulas + anchors ----

    private static (double bc, double jac) BetaOracle(
        Dictionary<string, double> s1, Dictionary<string, double> s2)
    {
        var all = s1.Keys.Union(s2.Keys).ToList();
        double sumMin = 0, sumTotal = 0;
        int shared = 0, u1 = 0, u2 = 0;
        foreach (var sp in all)
        {
            double a1 = s1.GetValueOrDefault(sp, 0), a2 = s2.GetValueOrDefault(sp, 0);
            sumMin += Math.Min(a1, a2);
            sumTotal += a1 + a2;
            bool in1 = a1 > 0, in2 = a2 > 0;
            if (in1 && in2) shared++; else if (in1) u1++; else if (in2) u2++;
        }
        double bc = sumTotal > 0 ? 1 - 2 * sumMin / sumTotal : 0;
        int total = shared + u1 + u2;
        double jac = total > 0 ? 1.0 - (double)shared / total : 0;
        return (bc, jac);
    }

    [Test]
    [Category("META-BETA-001")]
    public void BetaDiversity_MatchesIndependentFormulas()
    {
        var s1 = new Dictionary<string, double> { ["a"] = 1, ["b"] = 2, ["c"] = 3 };
        var s2 = new Dictionary<string, double> { ["b"] = 2, ["c"] = 1, ["d"] = 4 };

        var bd = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", s1, "s2", s2);
        var (bc, jac) = BetaOracle(s1, s2);
        Assert.That(bd.BrayCurtis, Is.EqualTo(bc).Within(Tol));
        Assert.That(bd.JaccardDistance, Is.EqualTo(jac).Within(Tol));
        Assert.That(bd.BrayCurtis, Is.InRange(0.0, 1.0));
        Assert.That(bd.JaccardDistance, Is.InRange(0.0, 1.0));
        Assert.That(bd.SharedSpecies, Is.EqualTo(2)); // b, c
    }

    [Test]
    [Category("META-BETA-001")]
    public void BetaDiversity_IdentityAndDisjointAnchors()
    {
        var s = new Dictionary<string, double> { ["a"] = 5, ["b"] = 3 };
        var same = MetagenomicsAnalyzer.CalculateBetaDiversity("s", s, "s", s);
        Assert.That(same.BrayCurtis, Is.EqualTo(0.0).Within(Tol), "identical samples -> BC 0");
        Assert.That(same.JaccardDistance, Is.EqualTo(0.0).Within(Tol), "identical samples -> Jaccard 0");

        var d1 = new Dictionary<string, double> { ["a"] = 1 };
        var d2 = new Dictionary<string, double> { ["z"] = 1 };
        var disjoint = MetagenomicsAnalyzer.CalculateBetaDiversity("d1", d1, "d2", d2);
        Assert.That(disjoint.BrayCurtis, Is.EqualTo(1.0).Within(Tol), "disjoint -> BC 1");
        Assert.That(disjoint.JaccardDistance, Is.EqualTo(1.0).Within(Tol), "disjoint -> Jaccard 1");
    }
}
