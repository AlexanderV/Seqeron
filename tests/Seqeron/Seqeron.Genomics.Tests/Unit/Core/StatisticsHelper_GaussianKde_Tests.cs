// StatisticsHelper.BandwidthNrd0 / GaussianKernelDensity / PeakPick — generic ports added for CNAqc's peak-based purity
// QC (ONCO-PURITY-001, FIN-B24 F33). Reference values: R 4.3.3 stats::bw.nrd0 and stats::density.default (old
// coordinates), R-trunk src/library/stats/R/density.R density.default (old.coords = FALSE, the R ≥ 4.4 default) and
// peakPick 0.11 peakpick (cran/peakPick), printed with sprintf("%.17g"). Inputs: the seeded CNAqc test datasets
// (TestData/CNAqc). Evidence: docs/Evidence/ONCO-PURITY-001-Evidence.md.

using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Tests.Helpers;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_GaussianKde_Tests
{
    private static readonly int[] GridIndices = { 1, 2, 100, 200, 256, 300, 400, 511, 512 }; // R 1-based

    // ── bw.nrd0 ──

    [TestCase("D1", 1, 1, 0.012856625638346407)]
    [TestCase("D3", 1, 1, 0.018444382013876173)]
    [TestCase("D1", 2, 1, 0.038758432750663083)]
    public void BandwidthNrd0_MatchesR(string dataset, int major, int minor, double expected)
    {
        double bw = StatisticsHelper.BandwidthNrd0(CnaqcTestData.Vafs(dataset, major, minor));
        Assert.That(bw, Is.EqualTo(expected).Within(1e-14 * expected));
    }

    [Test]
    public void BandwidthNrd0_SmallSamples_MatchR()
    {
        // R 4.3.3 sprintf("%.17g", bw.nrd0(...)):
        // c(1,2,3,4,10): min(sd = 3.54, IQR/1.34 = 2/1.34) → 0.9·(2/1.34)·5^(−1/5)
        Assert.That(StatisticsHelper.BandwidthNrd0(new[] { 1.0, 2, 3, 4, 10 }), Is.EqualTo(0.97358462285063574).Within(1e-15));
        // c(5,5,5): zero spread → |x₁| fallback, 0.9·5·3^(−1/5)
        Assert.That(StatisticsHelper.BandwidthNrd0(new[] { 5.0, 5, 5 }), Is.EqualTo(3.6123370279210381).Within(1e-15));
        // c(0,0): → 1 fallback, 0.9·2^(−1/5)
        Assert.That(StatisticsHelper.BandwidthNrd0(new[] { 0.0, 0 }), Is.EqualTo(0.78349550696651171).Within(1e-15));
        // c(1,1,1,1,9): IQR = 0 but sd = √12.8 > 0 → sd
        Assert.That(StatisticsHelper.BandwidthNrd0(new[] { 1.0, 1, 1, 1, 9 }), Is.EqualTo(2.3337454992375779).Within(1e-15));
    }

    [Test]
    public void BandwidthNrd0_FewerThanTwoPoints_Throws()
    {
        Assert.Throws<ArgumentException>(() => StatisticsHelper.BandwidthNrd0(new[] { 0.3 }));
        Assert.Throws<ArgumentNullException>(() => StatisticsHelper.BandwidthNrd0(null!));
    }

    // ── density.default (gaussian, n = 512, cut = 3) ──

    private static IEnumerable<TestCaseData> DensityCases()
    {
        // D1 1:1 (n = 400), new coordinates (R ≥ 4.4 default), adjust 1.
        yield return new TestCaseData("D1", 1, 1, 1.0, false, 0.012856625638346407,
            new[] { 0.18572918850552153, 0.18640389725753814, 0.25252535495516548, 0.31999623015682599, 0.35777992026975591, 0.38746710535848655, 0.45493798056014711, 0.52983065203399027, 0.53050536078600696 },
            new[] { 0.00095861337313726954, 0.0011257419735323385, 0.86931357010833754, 7.0056211434408713, 7.8057017230652654, 5.925178018820767, 1.1278329130976437, 0.002005470733908801, 0.0017138350430965449 })
            .SetName("Density_D1_1:1_NewCoords_Adjust1");
        // Same data, R 4.3.3 native (old.coords): identical grid, densities ≈ 0.999× apart.
        yield return new TestCaseData("D1", 1, 1, 1.0, true, 0.012856625638346407,
            new[] { 0.18572918850552153, 0.18640389725753814, 0.25252535495516548, 0.31999623015682599, 0.35777992026975591, 0.38746710535848655, 0.45493798056014711, 0.52983065203399027, 0.53050536078600696 },
            new[] { 0.00096753848631085074, 0.0011358885100311421, 0.87046261983246254, 7.0120985461430658, 7.8128770632778197, 5.930523020174487, 1.1290139079606947, 0.0020226310572655859, 0.0017290294062646317 })
            .SetName("Density_D1_1:1_OldCoords_Adjust1");
        // adjust = 0.5 (CNAqc kernel_adjust).
        yield return new TestCaseData("D1", 1, 1, 0.5, false, 0.0064283128191732035,
            new[] { 0.20501412696304114, 0.20561335650387372, 0.26433785150546618, 0.32426080558872372, 0.35781765987534797, 0.38418375967198137, 0.44410671375523891, 0.51062119278765483, 0.51122042232848741 },
            new[] { 0.0017633510738935855, 0.0023262423381932296, 1.2547462115904831, 7.197370044310027, 8.0782179161373886, 6.7205916201413238, 1.6317535165310431, 0.0044983725870440501, 0.003409292319568916 })
            .SetName("Density_D1_1:1_NewCoords_Adjust0.5");
        // D3 1:1 (n = 700, bimodal: clonal + subclonal tail; grid starts below 0).
        yield return new TestCaseData("D3", 1, 1, 1.0, false, 0.018444382013876173,
            new[] { -0.041047431755914236, -0.040152144693470684, 0.047585987425997185, 0.13711469367035214, 0.18725076916719094, 0.22664339991470711, 0.31617210615906211, 0.41554897009029612, 0.41644425715273964 },
            new[] { 0.00059395677465531177, 0.00069474298609031144, 2.1812747684738212, 2.2120686808528052, 4.691219100473865, 5.7357833537818435, 0.96824420037119274, 0.00046227934780752791, 0.00039856389401291491 })
            .SetName("Density_D3_1:1_NewCoords_Adjust1");
        yield return new TestCaseData("D3", 1, 1, 0.5, true, 0.0092221910069380866,
            new[] { -0.013380858735099975, -0.012593855715790728, 0.064532440176515493, 0.14323274210744022, 0.18730491118875806, 0.22193304403836495, 0.30063334596928964, 0.38799068111261609, 0.38877768413192537 },
            new[] { 0.00074140515086296957, 0.00095847331180325666, 3.3600873996651703, 2.1217991448603124, 4.6698992525609491, 5.8651870018292005, 1.5650865181902316, 0.00090957793343726147, 0.0007071262822027838 })
            .SetName("Density_D3_1:1_OldCoords_Adjust0.5");
        // D1 2:1 (n = 250, two multiplicity peaks).
        yield return new TestCaseData("D1", 2, 1, 1.0, false, 0.038758432750663083,
            new[] { 0.0097751219160779751, 0.011161169054708769, 0.14699378864052665, 0.28559850250360613, 0.36321714226693064, 0.42420321636668556, 0.56280793022976505, 0.71665916261778329, 0.71804520975641406 },
            new[] { 0.00048465958546188543, 0.00054092241038430867, 0.51480828288663949, 4.0050278964229742, 1.0689448698667177, 0.82248497867643522, 1.839986286399125, 0.0035858324206598335, 0.0032057670642611858 })
            .SetName("Density_D1_2:1_NewCoords_Adjust1");
    }

    [TestCaseSource(nameof(DensityCases))]
    public void GaussianKernelDensity_MatchesRDensityDefault(
        string dataset, int major, int minor, double adjust, bool legacy, double bw, double[] xs, double[] ys)
    {
        KernelDensityEstimate kde = StatisticsHelper.GaussianKernelDensity(
            CnaqcTestData.Vafs(dataset, major, minor), adjust, legacyCoordinates: legacy);

        Assert.That(kde.X, Has.Count.EqualTo(512));
        Assert.That(kde.Y, Has.Count.EqualTo(512));
        Assert.That(kde.Bandwidth, Is.EqualTo(bw).Within(1e-14 * bw));
        for (int k = 0; k < GridIndices.Length; k++)
        {
            int i = GridIndices[k] - 1;
            Assert.That(kde.X[i], Is.EqualTo(xs[k]).Within(1e-14), $"x[{GridIndices[k]}]");
            // R convolves by FFT; the direct sum agrees to ≈ 1e−15 absolute (1e−12 relative on the peak).
            Assert.That(kde.Y[i], Is.EqualTo(ys[k]).Within(1e-12 * Math.Max(1.0, ys[k])), $"y[{GridIndices[k]}]");
        }
    }

    [Test]
    public void GaussianKernelDensity_IntegratesToOne()
    {
        KernelDensityEstimate kde = StatisticsHelper.GaussianKernelDensity(CnaqcTestData.Vafs("D1", 2, 1));
        double dx = kde.X[1] - kde.X[0], area = 0;
        for (int i = 1; i < kde.Y.Count; i++) area += 0.5 * (kde.Y[i] + kde.Y[i - 1]) * dx;
        Assert.That(area, Is.EqualTo(1.0).Within(0.01), "±3 bw grid holds ≈ 99.7 % of a Gaussian KDE's mass");
    }

    [Test]
    public void GaussianKernelDensity_InvalidArguments_Throw()
    {
        double[] v = { 0.1, 0.2, 0.3 };
        Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.GaussianKernelDensity(v, adjust: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.GaussianKernelDensity(v, points: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.GaussianKernelDensity(v, cut: -1));
        Assert.Throws<ArgumentException>(() => StatisticsHelper.GaussianKernelDensity(new[] { 0.1 }));
    }

    // ── peakPick::peakpick ──

    [TestCase("D1", 1, 1, 1.0, false, new[] { 243 })]
    [TestCase("D1", 1, 1, 0.5, false, new[] { 193, 234, 254, 291, 362, 479 })]
    [TestCase("D1", 1, 1, 0.5, true, new[] { 193, 234, 253, 291, 362, 479 })] // old coords move one peak by a grid step
    [TestCase("D3", 1, 1, 1.0, false, new[] { 140, 292 })]
    [TestCase("D3", 1, 1, 0.5, true, new[] { 102, 140, 285, 327, 441 })]
    [TestCase("D1", 2, 1, 1.0, false, new[] { 181, 368 })]
    public void PeakPick_OnRDensity_MatchesPeakPick(string dataset, int major, int minor, double adjust, bool legacy, int[] expected1Based)
    {
        KernelDensityEstimate kde = StatisticsHelper.GaussianKernelDensity(
            CnaqcTestData.Vafs(dataset, major, minor), adjust, legacyCoordinates: legacy);
        for (int nl = 1; nl <= 5; nl++)
        {
            bool[] mask = StatisticsHelper.PeakPick(kde.Y, nl);
            int[] peaks = Enumerable.Range(0, mask.Length).Where(i => mask[i]).Select(i => i + 1).ToArray();
            Assert.That(peaks, Is.EqualTo(expected1Based), $"neighlim = {nl}");
        }
    }

    private static double[] Synthetic() =>
        Enumerable.Range(0, 200).Select(i => (0.05 * Math.Sin(i / 4.0)) + (0.03 * Math.Sin(i / 1.7)) + (0.0004 * i)).ToArray();

    // R: i <- 0:199; v <- 0.05*sin(i/4) + 0.03*sin(i/1.7) + 0.0004*i; which(peakpick(v, neighlim = nl)[,1]) - 1
    [TestCase(0, new[] { 12, 13, 25, 26, 33, 34, 56, 57, 78, 79, 86, 87, 100, 101, 108, 109, 131, 132, 153, 154, 161, 162, 175, 183, 184 })]
    [TestCase(1, new[] { 12, 25, 34, 56, 78, 87, 100, 109, 131, 153, 162, 175, 184 })]
    [TestCase(8, new[] { 12, 25, 34, 56, 78, 87, 100, 109, 131, 153, 162, 175, 184 })]
    [TestCase(12, new[] { 12, 34, 56, 78, 109, 131, 153, 184 })]
    [TestCase(20, new[] { 12, 34, 56, 78, 109, 131, 153, 184 })]
    public void PeakPick_NeighborLimit_MatchesPeakPick(int neighlim, int[] expected0Based)
    {
        bool[] mask = StatisticsHelper.PeakPick(Synthetic(), neighlim);
        Assert.That(Enumerable.Range(0, mask.Length).Where(i => mask[i]).ToArray(), Is.EqualTo(expected0Based));
    }

    [Test]
    public void PeakPick_NonDefaultParameters_MatchPeakPick()
    {
        int[] Peaks(bool[] m) => Enumerable.Range(0, m.Length).Where(i => m[i]).ToArray();
        int[] wide = { 4, 12, 25, 34, 46, 56, 66, 78, 87, 100, 109, 121, 131, 140, 153, 162, 175, 184, 195 };
        // peak.npos = 3, peak.min.sd = 0
        Assert.That(Peaks(StatisticsHelper.PeakPick(Synthetic(), 1, peakMinSd: 0, peakPositions: 3)), Is.EqualTo(wide));
        // peak.npos = 0: sd() of one value is NA in R, so no candidate is deleted by the small-peak filter.
        Assert.That(Peaks(StatisticsHelper.PeakPick(Synthetic(), 1, peakPositions: 0)), Is.EqualTo(wide));
        // deriv.lim = 0.01
        Assert.That(Peaks(StatisticsHelper.PeakPick(Synthetic(), 1, derivativeLimit: 0.01)),
            Is.EqualTo(new[] { 12, 25, 34, 56, 78, 87, 100, 109, 131, 153, 162, 175, 184 }));
    }

    [Test]
    public void PeakPick_ShortOrInvalid()
    {
        Assert.That(StatisticsHelper.PeakPick(new[] { 1.0, 2.0 }, 1), Is.EqualTo(new[] { false, false }));
        Assert.Throws<ArgumentOutOfRangeException>(() => StatisticsHelper.PeakPick(new[] { 1.0, 2, 1 }, -1));
        Assert.Throws<ArgumentNullException>(() => StatisticsHelper.PeakPick(null!, 1));
    }
}
