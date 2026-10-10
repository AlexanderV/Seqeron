// ONCO-PURITY-001 — CNAqc BMix mixture peaks and n_bootstrap (FIN-B24 F63):
//   OncologyAnalyzer.FitBinomialMixture (BMix bmixfit, K.Binomials = 1:4: R kmeans Hartigan-Wong nstart = 100, sample,
//   runif jitter, one EM step, ICL), PurityPeakOptions.FitMixturePeaks / BootstrapCount / Seed in AnalyzePurityPeaks,
//   BootstrapCount in AnalyzeComplexKaryotypePeaks / AnalyzeSubclonalPurityPeaks.
// Evidence: docs/Evidence/ONCO-PURITY-001-Evidence.md (§ "BMix and n_bootstrap")
// TestSpec: tests/TestSpecs/ONCO-PURITY-001.md (§ 5.10)
// Reference: caravagnalab/BMix R/bmixfit.R, bmixfit_EM.R, lib.R; CNAqc 1.1.5 R/peak_algorithms.R (combined /
//            mixture / simple peak detectors) run in R 4.3.3 after set.seed(seed); outputs in TestData/CNAqc/
//            cnaqc_bmix_R.txt (matprod "default" = OpenBLAS and "internal" = long-double dot products),
//            cnaqc_common_boot_R.txt, cnaqc_general_boot_R.txt, cnaqc_subclonal_boot_R.txt.

using System.Globalization;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Oncology;
using Seqeron.Genomics.Tests.Helpers;
using static Seqeron.Genomics.Tests.Unit.Oncology.OncologyAnalyzer_CnaqcGeneralSubclonalPeaks_Tests;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_CnaqcBMixBootstrap_Tests
{
    private static List<OncologyAnalyzer.PurityPeakMutation> Mutations(string dataset) =>
        CnaqcTestData.Load(dataset)
            .Select(r => new OncologyAnalyzer.PurityPeakMutation(r.Vaf, r.Major, r.Minor) { AlternateReads = r.Nv, Depth = r.Dp })
            .ToList();

    private static List<string> Reference(string file, string tag) =>
        ReferenceBlock(file, tag).Where(l => !l.StartsWith("  RNG-after", StringComparison.Ordinal)).ToList();

    // Line-by-line equality with numeric tokens compared to 1e−12 relative (1e−14 absolute): KDE heights are exact only
    // up to the direct-convolution vs FFT rounding of StatisticsHelper.GaussianKernelDensity (F33), and the KDE grid x
    // of a snapped BMix peak by the last bit of R's seq().
    private static void AssertLinesMatch(IReadOnlyList<string> actual, IReadOnlyList<string> expected)
    {
        Assert.That(actual.Count, Is.EqualTo(expected.Count), "line count");
        var split = new Regex("([ =,|])");
        for (int i = 0; i < expected.Count; i++)
        {
            string[] a = split.Split(actual[i]), e = split.Split(expected[i]);
            Assert.That(a.Length, Is.EqualTo(e.Length), $"line {i}: {actual[i]} vs {expected[i]}");
            for (int t = 0; t < e.Length; t++)
            {
                if (a[t] == e[t]) continue;
                bool numbers = double.TryParse(a[t], NumberStyles.Float, CultureInfo.InvariantCulture, out double x)
                               & double.TryParse(e[t], NumberStyles.Float, CultureInfo.InvariantCulture, out double y);
                Assert.That(numbers, Is.True, $"line {i}: '{a[t]}' vs '{e[t]}'");
                Assert.That(Math.Abs(x - y) <= Math.Max(1e-14, 1e-12 * Math.Max(Math.Abs(x), Math.Abs(y))), Is.True,
                    $"line {i}: {x:R} vs {y:R}");
            }
        }
    }

    #region BMix bmixfit

    private static IEnumerable<TestCaseData> BMixCases()
    {
        foreach ((string ds, int major, int minor, int seed) in new[]
                 { ("D1", 2, 1, 7), ("D1", 1, 1, 7), ("D1", 2, 2, 11), ("D3", 2, 0, 7), ("D3", 1, 1, 3), ("G1", 3, 1, 42), ("G1", 1, 1, 5) })
        {
            yield return new TestCaseData(ds, major, minor, seed).SetName($"FitBinomialMixture_BMixR_{ds}_{major}_{minor}_seed{seed}");
        }
    }

    // Every grid fit (NLL, BIC, ICL) and the ICL-best means / proportions equal R with matprod = "internal" bit for bit
    // (long-double dot products, emulated exactly); R's default OpenBLAS dot product moves the means by ≤ 1e−15 relative.
    [TestCaseSource(nameof(BMixCases))]
    public void FitBinomialMixture_MatchesBMixR(string dataset, int major, int minor, int seed)
    {
        var rows = CnaqcTestData.Load(dataset).Where(r => r.Major == major && r.Minor == minor).ToList();
        var (best, grid) = OncologyAnalyzer.FitBinomialMixture(rows.Select(r => r.Nv).ToList(), rows.Select(r => r.Dp).ToList(), seed);
        var actual = new List<string> { $"== {dataset} {major}:{minor} seed={seed} matprod=internal n={rows.Count}" };
        foreach (var g in grid)
            actual.Add($"  G K={g!.Components} NLL={F(g.NegativeLogLikelihood)} BIC={F(g.Bic)} ICL={F(g.Icl)}");
        actual.Add($"  BEST K={best.Components} B={string.Join(",", best.Means.Select(F))} pi={string.Join(",", best.MixingProportions.Select(F))}");
        Assert.That(actual, Is.EqualTo(Reference("cnaqc_bmix_R.txt", $"{dataset} {major}:{minor} seed={seed} matprod=internal")));

        var defaults = Reference("cnaqc_bmix_R.txt", $"{dataset} {major}:{minor} seed={seed} matprod=default");
        AssertLinesMatch(actual.Skip(1).ToList(), defaults.Skip(1).ToList());
    }

    // BMix's EM stops after one iteration (its convergence test subtracts the initial NLL .Machine$integer.max), so the
    // 2:1 fit of D1 is k-means + jitter + one E/M step; the ICL picks K = 2 with means 0.517 (m = 2) and 0.263 (m = 1)
    // around the expected 0.519 / 0.259 at π 0.7.
    [Test]
    public void FitBinomialMixture_D1_2_1_SelectsTwoComponentsNearExpectedPeaks()
    {
        var rows = CnaqcTestData.Load("D1").Where(r => r.Major == 2 && r.Minor == 1).ToList();
        var (best, grid) = OncologyAnalyzer.FitBinomialMixture(rows.Select(r => r.Nv).ToList(), rows.Select(r => r.Dp).ToList(), 7);
        Assert.That(grid.Select(g => g!.Components), Is.EqualTo(new[] { 1, 1, 2, 2, 3, 3, 4, 4 }));
        Assert.That(best.Components, Is.EqualTo(2));
        Assert.That(best.Icl, Is.EqualTo(1901.1530982123516));
        Assert.That(best.Means[0], Is.EqualTo(0.51667932020615859).Within(1e-15));
        Assert.That(best.Means[1], Is.EqualTo(0.26289651855238316).Within(1e-15));
    }

    // R kmeans() errors for K > number of distinct frequencies ("more cluster centers than distinct data points");
    // BMix's runner retries, easypar drops the task: K = 3, 4 are missing from the grid and K ≤ 2 still fit.
    [Test]
    public void FitBinomialMixture_TooFewDistinctFrequencies_DropsLargeK()
    {
        int[] nv = { 10, 10, 10, 30, 30, 30 }, dp = { 100, 100, 100, 100, 100, 100 };
        var (best, grid) = OncologyAnalyzer.FitBinomialMixture(nv, dp, 1);
        Assert.That(grid.Select(g => g is null), Is.EqualTo(new[] { false, false, false, false, true, true, true, true }));
        Assert.That(best.Components, Is.EqualTo(2));
    }

    [Test]
    public void FitBinomialMixture_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.FitBinomialMixture(null!, new[] { 10 }, 1));
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.FitBinomialMixture(new[] { 1, 2 }, new[] { 10 }, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.FitBinomialMixture(new[] { 11 }, new[] { 10 }, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.FitBinomialMixture(new[] { 1 }, new[] { 10 }, 1, new[] { 0 }));
    }

    #endregion

    #region analyze_peaks with BMix and n_bootstrap

    private static IEnumerable<TestCaseData> CommonCases()
    {
        static TestCaseData C(string tag, string ds, double p, OncologyAnalyzer.PurityPeakOptions o) =>
            new TestCaseData(tag, ds, p, o).SetName($"AnalyzePurityPeaks_BMixBootstrap_CNAqcR_{tag.Replace(' ', '_')}");
        yield return C("D1 p=0.7 bmix seed=7 nb=1", "D1", 0.7, new() { FitMixturePeaks = true, Seed = 7 });
        yield return C("D3 p=0.45 bmix seed=7 nb=1", "D3", 0.45, new() { FitMixturePeaks = true, Seed = 7 });
        yield return C("D4 p=0.3 bmix seed=11 nb=1", "D4", 0.3, new() { FitMixturePeaks = true, Seed = 11 });
        yield return C("D1 p=0.7 bmix seed=7 nb=3", "D1", 0.7, new() { FitMixturePeaks = true, Seed = 7, BootstrapCount = 3 });
        yield return C("D3 p=0.45 bmix seed=21 nb=5", "D3", 0.45, new() { FitMixturePeaks = true, Seed = 21, BootstrapCount = 5 });
        yield return C("D4 p=0.3 bmix seed=11 nb=2", "D4", 0.3, new() { FitMixturePeaks = true, Seed = 11, BootstrapCount = 2 });
        yield return C("D1 p=0.6 bmix seed=3 nb=4 adj=0.5", "D1", 0.6, new() { FitMixturePeaks = true, Seed = 3, BootstrapCount = 4, KernelAdjust = 0.5 });
        yield return C("D3 p=0.45 kde seed=5 nb=5", "D3", 0.45, new() { Seed = 5, BootstrapCount = 5 });
        yield return C("D1 p=0.7 kde seed=9 nb=10 adj=0.5", "D1", 0.7, new() { Seed = 9, BootstrapCount = 10, KernelAdjust = 0.5 });
    }

    [TestCaseSource(nameof(CommonCases))]
    public void AnalyzePurityPeaks_BMixAndBootstrap_MatchCnaqcR(string tag, string dataset, double purity, OncologyAnalyzer.PurityPeakOptions options)
    {
        var r = OncologyAnalyzer.AnalyzePurityPeaks(Mutations(dataset), purity, options);
        static string Src(OncologyAnalyzer.PurityDataPeak p) => p.Source == OncologyAnalyzer.PurityPeakSource.Kde ? "KDE" : "BMix";
        var actual = new List<string> { $"== {tag} score={F(r.Score)} QC={(r.Pass == true ? "PASS" : "FAIL")}" };
        foreach (var k in r.Karyotypes)
        {
            foreach (var m in k.Matches)
                actual.Add($"  M {m.MajorCopyNumber}:{m.MinorCopyNumber} m={m.Multiplicity} x={F(m.MatchedPeak.X)} y={F(m.MatchedPeak.Y)} " +
                           $"cpb={m.MatchedPeak.CountsPerBin?.ToString(CultureInfo.InvariantCulture) ?? "NA"} from={Src(m.MatchedPeak)} off={F(m.Offset)} " +
                           $"matched={B(m.Matched)} QC={(k.Pass ? "PASS" : "FAIL")}");
        }

        foreach (var k in r.Karyotypes)
        {
            foreach (var d in k.Peaks)
                actual.Add($"  P {k.MajorCopyNumber}:{k.MinorCopyNumber} x={F(d.X)} y={F(d.Y)} " +
                           $"cpb={d.CountsPerBin?.ToString(CultureInfo.InvariantCulture) ?? "NA"} disc={B(d.Discarded)} from={Src(d)}");
        }

        AssertLinesMatch(actual, Reference("cnaqc_common_boot_R.txt", tag));
    }

    // Fitted BMix peaks reproduce the F34 run that supplied R's B.params as MixturePeaks (D1, set.seed 7): same λ.
    [Test]
    public void AnalyzePurityPeaks_FitMixturePeaks_EqualsSuppliedRMeans()
    {
        var fitted = OncologyAnalyzer.AnalyzePurityPeaks(Mutations("D1"), 0.7, new() { FitMixturePeaks = true, Seed = 7 });
        Assert.That(fitted.Score, Is.EqualTo(0.0045358591466179345).Within(1e-15));
        Assert.That(fitted.Pass, Is.True);
        Assert.That(fitted.Karyotypes.SelectMany(k => k.Peaks).Count(p => p.Source == OncologyAnalyzer.PurityPeakSource.Mixture), Is.GreaterThan(0));
    }

    // Defaults are unchanged: no read counts needed, no randomness (F34 results bit-identical).
    [Test]
    public void AnalyzePurityPeaks_DefaultOptions_IgnoreSeedAndReads()
    {
        var plain = OncologyAnalyzer.AnalyzePurityPeaks(
            CnaqcTestData.Load("D1").Select(r => new OncologyAnalyzer.PurityPeakMutation(r.Vaf, r.Major, r.Minor)), 0.7);
        var seeded = OncologyAnalyzer.AnalyzePurityPeaks(Mutations("D1"), 0.7, new() { Seed = 123 });
        Assert.That(seeded.Score, Is.EqualTo(plain.Score));
        Assert.That(plain.Score, Is.EqualTo(0.0023756289876209163));
    }

    [Test]
    public void AnalyzePurityPeaks_BMixBootstrapOptions_InvalidArguments_Throw()
    {
        var withReads = Mutations("D4");
        var withoutReads = CnaqcTestData.Load("D4").Select(r => new OncologyAnalyzer.PurityPeakMutation(r.Vaf, r.Major, r.Minor)).ToList();
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzePurityPeaks(withReads, 0.3, new() { BootstrapCount = 0 }));
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.AnalyzePurityPeaks(withoutReads, 0.3, new() { FitMixturePeaks = true }));
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.AnalyzePurityPeaks(withReads, 0.3, new()
        {
            FitMixturePeaks = true,
            MixturePeaks = new Dictionary<(int Major, int Minor), IReadOnlyList<double>> { [(1, 1)] = new[] { 0.15 } },
        }));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(withReads, 0.3, new() { BootstrapCount = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(
            CnaqcTestData.Segments("S1"), 0.7, new() { BootstrapCount = 0 }));
    }

    private static IEnumerable<TestCaseData> GeneralBootCases()
    {
        yield return new TestCaseData("G1 p=0.6 seed=4 nb=3", "G1", 0.6, new OncologyAnalyzer.PurityPeakOptions { Seed = 4, BootstrapCount = 3 })
            .SetName("AnalyzeComplexKaryotypePeaks_Bootstrap_CNAqcR_G1_seed4_nb3");
        yield return new TestCaseData("G1 p=0.6 seed=8 nb=10 adj=0.5", "G1", 0.6, new OncologyAnalyzer.PurityPeakOptions { Seed = 8, BootstrapCount = 10, KernelAdjust = 0.5 })
            .SetName("AnalyzeComplexKaryotypePeaks_Bootstrap_CNAqcR_G1_seed8_nb10_adj0.5");
        yield return new TestCaseData("G2 p=0.4 seed=1 nb=5", "G2", 0.4, new OncologyAnalyzer.PurityPeakOptions { Seed = 1, BootstrapCount = 5 })
            .SetName("AnalyzeComplexKaryotypePeaks_Bootstrap_CNAqcR_G2_seed1_nb5");
    }

    [TestCaseSource(nameof(GeneralBootCases))]
    public void AnalyzeComplexKaryotypePeaks_Bootstrap_MatchesCnaqcR(string tag, string dataset, double purity, OncologyAnalyzer.PurityPeakOptions options)
    {
        var r = OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(Mutations(dataset), purity, options);
        AssertLinesMatch(FormatGeneral(tag, r), Reference("cnaqc_general_boot_R.txt", tag));
    }

    private static IEnumerable<TestCaseData> SubclonalBootCases()
    {
        yield return new TestCaseData("S1 p=0.7 seed=13 nb=3", "S1", 0.7, new OncologyAnalyzer.SubclonalPeakOptions { Seed = 13, BootstrapCount = 3 })
            .SetName("AnalyzeSubclonalPurityPeaks_Bootstrap_CNAqcR_S1_seed13_nb3");
        yield return new TestCaseData("S2 p=0.5 seed=2 nb=6 adj=0.5", "S2", 0.5, new OncologyAnalyzer.SubclonalPeakOptions { Seed = 2, BootstrapCount = 6, KernelAdjust = 0.5 })
            .SetName("AnalyzeSubclonalPurityPeaks_Bootstrap_CNAqcR_S2_seed2_nb6_adj0.5");
    }

    // Bootstrap resamples are drawn after the model identifiers on the same R stream (all expectations first).
    [TestCaseSource(nameof(SubclonalBootCases))]
    public void AnalyzeSubclonalPurityPeaks_Bootstrap_MatchesCnaqcR(string tag, string dataset, double purity, OncologyAnalyzer.SubclonalPeakOptions options)
    {
        var r = OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(CnaqcTestData.Segments(dataset), purity, options);
        AssertLinesMatch(FormatSubclonal(tag, r), Reference("cnaqc_subclonal_boot_R.txt", tag));
    }

    // Bootstrapping only adds peaks: the full-data peaks come first, unchanged.
    [Test]
    public void AnalyzePurityPeaks_Bootstrap_KeepsFullDataPeaksFirst()
    {
        var single = OncologyAnalyzer.AnalyzePurityPeaks(Mutations("D3"), 0.45);
        var boot = OncologyAnalyzer.AnalyzePurityPeaks(Mutations("D3"), 0.45, new() { Seed = 5, BootstrapCount = 5 });
        for (int k = 0; k < single.Karyotypes.Count; k++)
        {
            var p1 = single.Karyotypes[k].Peaks;
            Assert.That(boot.Karyotypes[k].Peaks.Take(p1.Count), Is.EqualTo(p1));
            Assert.That(boot.Karyotypes[k].Peaks.Count, Is.GreaterThanOrEqualTo(p1.Count));
        }
    }

    #endregion
}
