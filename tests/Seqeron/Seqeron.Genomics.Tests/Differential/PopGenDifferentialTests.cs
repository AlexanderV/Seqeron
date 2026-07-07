// 08_DIFFERENTIAL_TESTING rows 43-47 (Population genetics). Each statistic is recomputed by an
// INDEPENDENT oracle from its published definition (manual allele counting; brute pairwise diversity +
// closed-form Watterson/Tajima; manual HWE expected genotypes + chi-square; the Fst variance ratio; and
// an independent Pearson-correlation LD) and cross-checked against hand-known anchor values.

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class PopGenDifferentialTests
{
    private const double Tol = 1e-12;

    // ---- Row 43: POP-FREQ-001 — allele frequencies vs manual count/total ----

    [Test]
    [Category("POP-FREQ-001")]
    [TestCase(10, 5, 5)]
    [TestCase(0, 0, 7)]
    [TestCase(3, 4, 1)]
    public void AlleleFrequencies_MatchManualCount(int homMajor, int het, int homMinor)
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(homMajor, het, homMinor);
        int total = 2 * (homMajor + het + homMinor);
        double majOracle = total == 0 ? 0 : (2.0 * homMajor + het) / total;
        double minOracle = total == 0 ? 0 : (2.0 * homMinor + het) / total;
        Assert.That(major, Is.EqualTo(majOracle).Within(Tol));
        Assert.That(minor, Is.EqualTo(minOracle).Within(Tol));
    }

    [Test]
    [Category("POP-FREQ-001")]
    public void Maf_MatchesManualCount()
    {
        var genos = new[] { 0, 1, 2, 1, 0, 2, 1 };          // alt alleles = sum = 7, 2N = 14
        double altFreq = (double)genos.Sum() / (2 * genos.Length);
        Assert.That(PopulationGeneticsAnalyzer.CalculateMAF(genos),
            Is.EqualTo(Math.Min(altFreq, 1 - altFreq)).Within(Tol));
    }

    // ---- Row 44: POP-DIV-001 — π (brute), Watterson θ, Tajima's D vs closed form ----

    private static readonly string[] Sample = { "AAAACCCCGG", "AAAACCTCGG", "AAGACCCCGT", "TAAACCCCGG" };

    [Test]
    [Category("POP-DIV-001")]
    public void NucleotideDiversity_Watterson_Tajima_MatchOracles()
    {
        var seqs = Sample.Select(s => (IReadOnlyList<char>)s.ToCharArray()).ToList();
        int n = Sample.Length, L = Sample[0].Length;

        // Brute pairwise π.
        double totalDiff = 0; int pairs = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                int d = 0;
                for (int k = 0; k < L; k++) if (Sample[i][k] != Sample[j][k]) d++;
                totalDiff += d; pairs++;
            }
        double pi = totalDiff / (pairs * L);
        Assert.That(PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(seqs), Is.EqualTo(pi).Within(Tol));
        Assert.That(pi, Is.GreaterThanOrEqualTo(0));

        // Segregating sites + Watterson θ = S / (a1 * L).
        int S = 0;
        for (int k = 0; k < L; k++)
            if (Sample.Select(s => s[k]).Distinct().Count() > 1) S++;
        double a1 = Enumerable.Range(1, n - 1).Sum(i => 1.0 / i);
        double theta = S / (a1 * L);
        Assert.That(PopulationGeneticsAnalyzer.CalculateWattersonTheta(S, n, L), Is.EqualTo(theta).Within(Tol));
        Assert.That(theta, Is.GreaterThanOrEqualTo(0));

        // Tajima's D closed form (Tajima 1989).
        double kHat = totalDiff / pairs;   // average pairwise differences = π * L
        double a2 = Enumerable.Range(1, n - 1).Sum(i => 1.0 / (i * i));
        double b1 = (n + 1.0) / (3 * (n - 1));
        double b2 = 2.0 * (n * n + n + 3) / (9.0 * n * (n - 1));
        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);
        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);
        double variance = e1 * S + e2 * S * (S - 1);
        double expectedD = variance <= 0 ? 0 : (kHat - S / a1) / Math.Sqrt(variance);
        Assert.That(PopulationGeneticsAnalyzer.CalculateTajimasD(kHat, S, n), Is.EqualTo(expectedD).Within(1e-9));
    }

    // ---- Row 45: POP-HW-001 — chi-square vs manual expected genotypes ----

    [Test]
    [Category("POP-HW-001")]
    [TestCase(25, 50, 25, 0.0)]    // perfect HWE -> chi2 = 0
    [TestCase(40, 20, 40, 36.0)]   // p=0.5, expected 25/50/25 -> chi2 = 9+18+9
    public void HardyWeinberg_ChiSquare_MatchesManualExpected(int aa, int ab, int bb, double expectedChi2)
    {
        int total = aa + ab + bb;
        double p = (2.0 * aa + ab) / (2.0 * total), q = 1 - p;
        double eAA = p * p * total, eAa = 2 * p * q * total, eaa = q * q * total;
        double chi2 = 0;
        if (eAA > 0) chi2 += Math.Pow(aa - eAA, 2) / eAA;
        if (eAa > 0) chi2 += Math.Pow(ab - eAa, 2) / eAa;
        if (eaa > 0) chi2 += Math.Pow(bb - eaa, 2) / eaa;

        var r = PopulationGeneticsAnalyzer.TestHardyWeinberg("v", aa, ab, bb);
        Assert.That(r.ExpectedAA, Is.EqualTo(eAA).Within(Tol));
        Assert.That(r.ExpectedAa, Is.EqualTo(eAa).Within(Tol));
        Assert.That(r.Expectedaa, Is.EqualTo(eaa).Within(Tol));
        Assert.That(r.ChiSquare, Is.EqualTo(chi2).Within(1e-9));
        Assert.That(chi2, Is.EqualTo(expectedChi2).Within(1e-9));   // hand-known anchor
        Assert.That(r.PValue, Is.InRange(0.0, 1.0));
        if (expectedChi2 == 0.0) Assert.That(r.InEquilibrium, Is.True, "chi2=0 -> p=1 -> in equilibrium");
    }

    // ---- Row 46: POP-FST-001 — Fst variance ratio vs independent formula + anchors ----

    private static double FstOracle((double p, int n)[] pop1, (double p, int n)[] pop2)
    {
        double num = 0, den = 0;
        for (int i = 0; i < pop1.Length; i++)
        {
            double p1 = pop1[i].p, p2 = pop2[i].p;
            int n1 = pop1[i].n, n2 = pop2[i].n;
            double pBar = (n1 * p1 + n2 * p2) / (n1 + n2);
            double var = ((p1 - pBar) * (p1 - pBar) * n1 + (p2 - pBar) * (p2 - pBar) * n2) / (n1 + n2);
            num += var;
            den += pBar * (1 - pBar);
        }
        return den > 0 ? num / den : 0;
    }

    [Test]
    [Category("POP-FST-001")]
    public void Fst_MatchesVarianceRatioAndAnchors()
    {
        var a = new[] { (0.2, 50), (0.6, 50) };
        var b = new[] { (0.5, 50), (0.4, 50) };
        double fst = PopulationGeneticsAnalyzer.CalculateFst(
            a.Select(x => (x.Item1, x.Item2)), b.Select(x => (x.Item1, x.Item2)));
        Assert.That(fst, Is.EqualTo(FstOracle(a, b)).Within(Tol));
        Assert.That(fst, Is.InRange(0.0, 1.0));

        // Anchors: identical populations -> 0; fully divergent (0 vs 1, equal n) -> 1.
        var same = new[] { (0.3, 40) };
        Assert.That(PopulationGeneticsAnalyzer.CalculateFst(same.Select(x => (x.Item1, x.Item2)), same.Select(x => (x.Item1, x.Item2))),
            Is.EqualTo(0.0).Within(Tol));
        var zero = new[] { (0.0, 30) }; var one = new[] { (1.0, 30) };
        Assert.That(PopulationGeneticsAnalyzer.CalculateFst(zero.Select(x => (x.Item1, x.Item2)), one.Select(x => (x.Item1, x.Item2))),
            Is.EqualTo(1.0).Within(Tol));
    }

    // ---- Row 47: POP-LD-001 — D' / r² vs independent Pearson-correlation LD ----

    [Test]
    [Category("POP-LD-001")]
    public void Ld_MatchesIndependentPearsonAndAnchors()
    {
        var g = new[] { (0, 0), (1, 1), (2, 2), (1, 0), (0, 1), (2, 1) };
        var ld = PopulationGeneticsAnalyzer.CalculateLD("v1", "v2", g, distance: 100);

        int n = g.Length;
        double p1 = g.Sum(x => x.Item1) / (2.0 * n), p2 = g.Sum(x => x.Item2) / (2.0 * n);
        double m1 = g.Average(x => (double)x.Item1), m2 = g.Average(x => (double)x.Item2);
        double cov = g.Sum(x => (x.Item1 - m1) * (x.Item2 - m2)) / n;
        double v1 = g.Sum(x => Math.Pow(x.Item1 - m1, 2)) / n;
        double v2 = g.Sum(x => Math.Pow(x.Item2 - m2, 2)) / n;
        double r2 = (v1 > 0 && v2 > 0) ? (cov * cov) / (v1 * v2) : 0;
        double d = cov / 2;
        double dMax = d >= 0 ? Math.Min(p1 * (1 - p2), (1 - p1) * p2) : Math.Min(p1 * p2, (1 - p1) * (1 - p2));
        double dPrime = dMax > 1e-10 ? Math.Min(Math.Abs(d) / dMax, 1.0) : 0;

        Assert.That(ld.RSquared, Is.EqualTo(r2).Within(1e-12));
        Assert.That(ld.DPrime, Is.EqualTo(dPrime).Within(1e-12));

        // Anchor: perfectly correlated genotypes -> r² = 1, D' = 1.
        var perfect = new[] { (0, 0), (1, 1), (2, 2), (0, 0), (2, 2) };
        var ldp = PopulationGeneticsAnalyzer.CalculateLD("v1", "v2", perfect, 1);
        Assert.That(ldp.RSquared, Is.EqualTo(1.0).Within(1e-9));
        Assert.That(ldp.DPrime, Is.EqualTo(1.0).Within(1e-9));
    }
}
