// ONCO-PURITY-001 — CNAqc peak-based purity QC (FIN-B24 F34): OncologyAnalyzer.AnalyzePurityPeaks.
// Evidence: docs/Evidence/ONCO-PURITY-001-Evidence.md (§ "CNAqc analyze_peaks — peak-based purity QC")
// TestSpec: tests/TestSpecs/ONCO-PURITY-001.md (§ AnalyzePurityPeaks)
// Reference: caravagnalab/CNAqc 1.1.5 (commit 4b7cea4a) R/peak_algorithms.R analyze_peaks_common, R/vaf_functions.R,
//            R/equations.R (delta_vaf_karyo, compute_delta_purity), R/peak_detector.R (overlap_bands; legacy rightmost
//            matching), peakPick 0.11, BMix (caravagnalab, bmixfit) — run in R 4.3.3 with the R ≥ 4.4 density lattice
//            (R-trunk density.R) unless "old"; values printed with sprintf("%.17g"). Data: TestData/CNAqc (seeded).

using Seqeron.Genomics.Tests.Helpers;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_AnalyzePurityPeaks_Tests
{
    // Expected match row (CNAqc peaks_analysis$matches) and data peak (fits[[k]]$xy_peaks).
    public sealed record E(int Major, int Minor, int M, double Peak, double DeltaVaf, double X, double Y, int Counts,
        bool Mixture, double OffsetVaf, double Offset, double Weight, bool Matched, bool KaryotypePass);

    public sealed record P(int Major, int Minor, double X, double Y, int Counts, bool Discarded, bool Mixture);

    private static List<OncologyAnalyzer.PurityPeakMutation> Mutations(string dataset) =>
        CnaqcTestData.Load(dataset).Select(r => new OncologyAnalyzer.PurityPeakMutation(r.Vaf, r.Major, r.Minor)).ToList();

    private static IEnumerable<TestCaseData> CnaqcReferenceCases()
    {
        yield return new TestCaseData("D1", 0.7, new OncologyAnalyzer.PurityPeakOptions {  }, 0.0023756289876209163, true,
            new[]
            {
                new E(1, 0, 1, 0.53846153846153844, 0.059171597633136092, 0.56000000000000005, 8.6300000000000008, 8, false, -0.021538461538461617, -0.017700905274869835, 0.13043478260869565, true, true),
                new E(1, 1, 1, 0.34999999999999998, 0.025000000000000001, 0.34999999999999998, 7.9400000000000004, 36, false, 0, 0, 0.43478260869565216, true, true),
                new E(2, 1, 1, 0.25925925925925924, 0.013717421124828532, 0.26000000000000001, 4.4500000000000002, 16, false, -0.00074074074074076401, -0.0027054081108136012, 0.27173913043478259, true, true),
                new E(2, 1, 2, 0.51851851851851849, 0.027434842249657063, 0.52000000000000002, 2.3799999999999999, 6, false, -0.001481481481481528, -0.0027054081108136012, 0.27173913043478259, true, true),
                new E(2, 2, 1, 0.20588235294117646, 0.0086505190311418692, 0.20000000000000001, 3.6499999999999999, 11, false, 0.0058823529411764497, 0.032679738562091387, 0.16304347826086957, true, true),
                new E(2, 2, 2, 0.41176470588235292, 0.017301038062283738, 0.40999999999999998, 3.4399999999999999, 7, false, 0.001764705882352946, 0.0050695371512581031, 0.16304347826086957, true, true),
            },
            new[]
            {
                new P(1, 0, 0.56, 8.63, 8, false, false),
                new P(1, 1, 0.35, 7.94, 36, false, false),
                new P(2, 1, 0.26, 4.45, 16, false, false),
                new P(2, 1, 0.52, 2.38, 6, false, false),
                new P(2, 2, 0.2, 3.65, 11, false, false),
                new P(2, 2, 0.41, 3.44, 7, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D1_p0.7_new");
        yield return new TestCaseData("D1", 0.7, new OncologyAnalyzer.PurityPeakOptions { LegacyDensityCoordinates = true }, 0.0023756289876209163, true,
            new[]
            {
                new E(1, 0, 1, 0.53846153846153844, 0.059171597633136092, 0.56000000000000005, 8.6400000000000006, 8, false, -0.021538461538461617, -0.017700905274869835, 0.13043478260869565, true, true),
                new E(1, 1, 1, 0.34999999999999998, 0.025000000000000001, 0.34999999999999998, 7.9500000000000002, 36, false, 0, 0, 0.43478260869565216, true, true),
                new E(2, 1, 1, 0.25925925925925924, 0.013717421124828532, 0.26000000000000001, 4.4500000000000002, 16, false, -0.00074074074074076401, -0.0027054081108136012, 0.27173913043478259, true, true),
                new E(2, 1, 2, 0.51851851851851849, 0.027434842249657063, 0.52000000000000002, 2.3799999999999999, 6, false, -0.001481481481481528, -0.0027054081108136012, 0.27173913043478259, true, true),
                new E(2, 2, 1, 0.20588235294117646, 0.0086505190311418692, 0.20000000000000001, 3.6499999999999999, 11, false, 0.0058823529411764497, 0.032679738562091387, 0.16304347826086957, true, true),
                new E(2, 2, 2, 0.41176470588235292, 0.017301038062283738, 0.40999999999999998, 3.4399999999999999, 7, false, 0.001764705882352946, 0.0050695371512581031, 0.16304347826086957, true, true),
            },
            new[]
            {
                new P(1, 0, 0.56, 8.64, 8, false, false),
                new P(1, 1, 0.35, 7.95, 36, false, false),
                new P(2, 1, 0.26, 4.45, 16, false, false),
                new P(2, 1, 0.52, 2.38, 6, false, false),
                new P(2, 2, 0.2, 3.65, 11, false, false),
                new P(2, 2, 0.41, 3.44, 7, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D1_p0.7_old");
        yield return new TestCaseData("D1", 0.5, new OncologyAnalyzer.PurityPeakOptions {  }, -0.29645326597409133, false,
            new[]
            {
                new E(1, 0, 1, 0.33333333333333331, 0.044444444444444446, 0.56000000000000005, 8.6300000000000008, 8, false, -0.22666666666666674, -0.1862809555117248, 0.13043478260869565, false, false),
                new E(1, 1, 1, 0.25, 0.025000000000000001, 0.34999999999999998, 7.9400000000000004, 36, false, -0.099999999999999978, -0.19999999999999996, 0.43478260869565216, false, false),
                new E(2, 1, 1, 0.20000000000000001, 0.016, 0.26000000000000001, 4.4500000000000002, 16, false, -0.059999999999999998, -0.21913805697589481, 0.27173913043478259, false, false),
                new E(2, 1, 2, 0.40000000000000002, 0.032000000000000001, 0.52000000000000002, 2.3799999999999999, 6, false, -0.12, -0.21913805697589481, 0.27173913043478259, false, false),
                new E(2, 2, 1, 0.16666666666666666, 0.011111111111111112, 0.20000000000000001, 3.6499999999999999, 11, false, -0.033333333333333354, -0.18518518518518531, 0.16304347826086957, false, false),
                new E(2, 2, 2, 0.33333333333333331, 0.022222222222222223, 0.40999999999999998, 3.4399999999999999, 7, false, -0.076666666666666661, -0.22024322512687919, 0.16304347826086957, false, false),
            },
            new[]
            {
                new P(1, 0, 0.56, 8.63, 8, false, false),
                new P(1, 1, 0.35, 7.94, 36, false, false),
                new P(2, 1, 0.26, 4.45, 16, false, false),
                new P(2, 1, 0.52, 2.38, 6, false, false),
                new P(2, 2, 0.2, 3.65, 11, false, false),
                new P(2, 2, 0.41, 3.44, 7, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D1_p0.5_new");
        yield return new TestCaseData("D1", 0.7, new OncologyAnalyzer.PurityPeakOptions { KernelAdjust = 0.5 }, 0.023857243903529196, true,
            new[]
            {
                new E(1, 0, 1, 0.53846153846153844, 0.059171597633136092, 0.56000000000000005, 10.779999999999999, 8, false, -0.021538461538461617, -0.017700905274869835, 0.13043478260869565, true, true),
                new E(1, 1, 1, 0.34999999999999998, 0.025000000000000001, 0.34000000000000002, 8.1899999999999995, 36, false, 0.0099999999999999534, 0.019999999999999907, 0.43478260869565216, true, true),
                new E(2, 1, 1, 0.25925925925925924, 0.013717421124828532, 0.25, 5.5599999999999996, 16, false, 0.0092592592592592449, 0.032921810699588425, 0.27173913043478259, true, true),
                new E(2, 1, 2, 0.51851851851851849, 0.027434842249657063, 0.53000000000000003, 2.8199999999999998, 11, false, -0.011481481481481537, -0.021253147265457056, 0.27173913043478259, true, true),
                new E(2, 2, 1, 0.20588235294117646, 0.0086505190311418692, 0.19, 4.6200000000000001, 10, false, 0.015882352941176459, 0.082634510620064813, 0.16304347826086957, true, true),
                new E(2, 2, 2, 0.41176470588235292, 0.017301038062283738, 0.40999999999999998, 4.2300000000000004, 7, false, 0.001764705882352946, 0.0050695371512581031, 0.16304347826086957, true, true),
            },
            new[]
            {
                new P(1, 0, 0.45, 3.48, 4, false, false),
                new P(1, 0, 0.5, 6.09, 9, false, false),
                new P(1, 0, 0.56, 10.78, 8, false, false),
                new P(1, 0, 0.65, 1.45, 2, false, false),
                new P(1, 1, 0.32, 7.33, 30, false, false),
                new P(1, 1, 0.34, 8.19, 36, false, false),
                new P(1, 1, 0.36, 8.08, 30, false, false),
                new P(1, 1, 0.38, 6.99, 27, false, false),
                new P(1, 1, 0.42, 3.15, 12, false, false),
                new P(1, 1, 0.49, 0.31, 0, true, false),
                new P(2, 1, 0.25, 5.56, 16, false, false),
                new P(2, 1, 0.53, 2.82, 11, false, false),
                new P(2, 2, 0.19, 4.62, 10, false, false),
                new P(2, 2, 0.41, 4.23, 7, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D1_p0.7_new_adjust0.5");
        yield return new TestCaseData("D1", 0.7, new OncologyAnalyzer.PurityPeakOptions { PurityError = 0.01, VafTolerance = 0.005 }, 0.0023756289876209163, true,
            new[]
            {
                new E(1, 0, 1, 0.53846153846153844, 0.011834319526627219, 0.56000000000000005, 8.6300000000000008, 8, false, -0.021538461538461617, -0.017700905274869835, 0.13043478260869565, false, false),
                new E(1, 1, 1, 0.34999999999999998, 0.0050000000000000001, 0.34999999999999998, 7.9400000000000004, 36, false, 0, 0, 0.43478260869565216, true, true),
                new E(2, 1, 1, 0.25925925925925924, 0.0027434842249657062, 0.26000000000000001, 4.4500000000000002, 16, false, -0.00074074074074076401, -0.0027054081108136012, 0.27173913043478259, true, true),
                new E(2, 1, 2, 0.51851851851851849, 0.0054869684499314125, 0.52000000000000002, 2.3799999999999999, 6, false, -0.001481481481481528, -0.0027054081108136012, 0.27173913043478259, true, true),
                new E(2, 2, 1, 0.20588235294117646, 0.001730103806228374, 0.20000000000000001, 3.6499999999999999, 11, false, 0.0058823529411764497, 0.032679738562091387, 0.16304347826086957, true, true),
                new E(2, 2, 2, 0.41176470588235292, 0.0034602076124567479, 0.40999999999999998, 3.4399999999999999, 7, false, 0.001764705882352946, 0.0050695371512581031, 0.16304347826086957, true, true),
            },
            new[]
            {
                new P(1, 0, 0.56, 8.63, 8, false, false),
                new P(1, 1, 0.35, 7.94, 36, false, false),
                new P(2, 1, 0.26, 4.45, 16, false, false),
                new P(2, 1, 0.52, 2.38, 6, false, false),
                new P(2, 2, 0.2, 3.65, 11, false, false),
                new P(2, 2, 0.41, 3.44, 7, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D1_p0.7_new_purity_error0.01_tol0.005");
        yield return new TestCaseData("D3", 0.45, new OncologyAnalyzer.PurityPeakOptions {  }, 0.014615384615384629, true,
            new[]
            {
                new E(1, 1, 1, 0.22500000000000001, 0.025000000000000001, 0.22, 5.7800000000000002, 45, false, 0.0050000000000000044, 0.010000000000000009, 0.76923076923076927, true, true),
                new E(2, 0, 1, 0.22500000000000001, 0.025000000000000001, 0.20999999999999999, 1.21, 2, false, 0.015000000000000013, 0.030000000000000027, 0.23076923076923078, true, true),
                new E(2, 0, 2, 0.45000000000000001, 0.050000000000000003, 0.45000000000000001, 2.8500000000000001, 9, false, 0, 0, 0.23076923076923078, true, true),
            },
            new[]
            {
                new P(1, 1, 0.08, 3.21, 18, false, false),
                new P(1, 1, 0.22, 5.78, 45, false, false),
                new P(2, 0, 0.21, 1.21, 2, false, false),
                new P(2, 0, 0.45, 2.85, 9, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D3_p0.45_new");
        yield return new TestCaseData("D4", 0.3, new OncologyAnalyzer.PurityPeakOptions {  }, 0.0039996443291489122, true,
            new[]
            {
                new E(1, 0, 1, 0.17647058823529413, 0.034602076124567477, 0.17000000000000001, 8.2899999999999991, 11, false, 0.0064705882352941169, 0.0094537047779883372, 0.42307692307692307, true, true),
                new E(1, 1, 1, 0.14999999999999999, 0.025000000000000001, 0.14999999999999999, 8.1300000000000008, 11, false, 0, 0, 0.57692307692307687, true, true),
            },
            new[]
            {
                new P(1, 0, 0.17, 8.29, 11, false, false),
                new P(1, 1, 0.15, 8.13, 11, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D4_p0.3_new");
        yield return new TestCaseData("D4", 0.3, new OncologyAnalyzer.PurityPeakOptions { MinKaryotypeSize = 0.35 }, 0, true,
            new[]
            {
                new E(1, 1, 1, 0.14999999999999999, 0.025000000000000001, 0.14999999999999999, 8.1300000000000008, 11, false, 0, 0, 1, true, true),
            },
            new[]
            {
                new P(1, 1, 0.15, 8.13, 11, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D4_p0.3_new_size0.35");
        yield return new TestCaseData("D4", 0.3, new OncologyAnalyzer.PurityPeakOptions { MinAbsoluteKaryotypeMutations = 80 }, -0.010660215188956853, true,
            new[]
            {
                new E(1, 0, 1, 0.17647058823529413, 0.034602076124567477, 0.17000000000000001, 8.2899999999999991, 11, false, 0.0064705882352941169, 0.0094537047779883372, 0.31428571428571428, true, true),
                new E(1, 1, 1, 0.14999999999999999, 0.025000000000000001, 0.14999999999999999, 8.1300000000000008, 11, false, 0, 0, 0.42857142857142855, true, true),
                new E(2, 1, 1, 0.13043478260869565, 0.018903591682419663, 0.14999999999999999, 4.21, 3, false, -0.019565217391304346, -0.054159771325409958, 0.25714285714285712, true, true),
                new E(2, 1, 2, 0.2608695652173913, 0.037807183364839327, 0.26000000000000001, 4.1399999999999997, 4, false, 0.00086956521739128823, 0.0011488508619253378, 0.25714285714285712, true, true),
            },
            new[]
            {
                new P(1, 0, 0.17, 8.29, 11, false, false),
                new P(1, 1, 0.15, 8.13, 11, false, false),
                new P(2, 1, 0.15, 4.21, 3, false, false),
                new P(2, 1, 0.26, 4.14, 4, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D4_p0.3_new_minabs80");
        yield return new TestCaseData("D3", 0.45, new OncologyAnalyzer.PurityPeakOptions { MinVaf = 0.12 }, 0.030000000000000027, true,
            new[]
            {
                new E(1, 1, 1, 0.22500000000000001, 0.025000000000000001, 0.20999999999999999, 8.1699999999999999, 40, false, 0.015000000000000013, 0.030000000000000027, 0.71034482758620687, true, true),
                new E(2, 0, 1, 0.22500000000000001, 0.025000000000000001, 0.20999999999999999, 1.21, 2, false, 0.015000000000000013, 0.030000000000000027, 0.28965517241379313, true, true),
                new E(2, 0, 2, 0.45000000000000001, 0.050000000000000003, 0.45000000000000001, 2.8500000000000001, 9, false, 0, 0, 0.28965517241379313, true, true),
            },
            new[]
            {
                new P(1, 1, 0.21, 8.17, 40, false, false),
                new P(2, 0, 0.21, 1.21, 2, false, false),
                new P(2, 0, 0.45, 2.85, 9, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D3_p0.45_new_minVAF0.12");
        yield return new TestCaseData("D3", 0.45, new OncologyAnalyzer.PurityPeakOptions { MatchingStrategy = OncologyAnalyzer.PurityPeakMatchingStrategy.Rightmost }, 0.014615384615384629, true,
            new[]
            {
                new E(1, 1, 1, 0.22500000000000001, 0.025000000000000001, 0.22, 5.7800000000000002, 45, false, 0.0050000000000000044, 0.010000000000000009, 0.76923076923076927, true, true),
                new E(2, 0, 2, 0.45000000000000001, 0.050000000000000003, 0.45000000000000001, 2.8500000000000001, 9, false, 0, 0, 0.23076923076923078, true, true),
                new E(2, 0, 1, 0.22500000000000001, 0.025000000000000001, 0.20999999999999999, 1.21, 2, false, 0.015000000000000013, 0.030000000000000027, 0.23076923076923078, true, true),
            },
            new[]
            {
                new P(1, 1, 0.08, 3.21, 18, false, false),
                new P(1, 1, 0.22, 5.78, 45, false, false),
                new P(2, 0, 0.21, 1.21, 2, false, false),
                new P(2, 0, 0.45, 2.85, 9, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D3_p0.45_new_rightmost");
        yield return new TestCaseData("D1", 0.7, new OncologyAnalyzer.PurityPeakOptions { MatchingStrategy = OncologyAnalyzer.PurityPeakMatchingStrategy.Rightmost }, 0.0023756289876209163, true,
            new[]
            {
                new E(1, 0, 1, 0.53846153846153844, 0.059171597633136092, 0.56000000000000005, 8.6300000000000008, 8, false, -0.021538461538461617, -0.017700905274869835, 0.13043478260869565, true, true),
                new E(1, 1, 1, 0.34999999999999998, 0.025000000000000001, 0.34999999999999998, 7.9400000000000004, 36, false, 0, 0, 0.43478260869565216, true, true),
                new E(2, 1, 2, 0.51851851851851849, 0.027434842249657063, 0.52000000000000002, 2.3799999999999999, 6, false, -0.001481481481481528, -0.0027054081108136012, 0.27173913043478259, true, true),
                new E(2, 1, 1, 0.25925925925925924, 0.013717421124828532, 0.26000000000000001, 4.4500000000000002, 16, false, -0.00074074074074076401, -0.0027054081108136012, 0.27173913043478259, true, true),
                new E(2, 2, 2, 0.41176470588235292, 0.017301038062283738, 0.40999999999999998, 3.4399999999999999, 7, false, 0.001764705882352946, 0.0050695371512581031, 0.16304347826086957, true, true),
                new E(2, 2, 1, 0.20588235294117646, 0.0086505190311418692, 0.20000000000000001, 3.6499999999999999, 11, false, 0.0058823529411764497, 0.032679738562091387, 0.16304347826086957, true, true),
            },
            new[]
            {
                new P(1, 0, 0.56, 8.63, 8, false, false),
                new P(1, 1, 0.35, 7.94, 36, false, false),
                new P(2, 1, 0.26, 4.45, 16, false, false),
                new P(2, 1, 0.52, 2.38, 6, false, false),
                new P(2, 2, 0.2, 3.65, 11, false, false),
                new P(2, 2, 0.41, 3.44, 7, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D1_p0.7_new_rightmost");
        yield return new TestCaseData("D1", 0.62, new OncologyAnalyzer.PurityPeakOptions {  }, -0.11362481021020236, true,
            new[]
            {
                new E(1, 0, 1, 0.44927536231884063, 0.052509976895610176, 0.56000000000000005, 8.6300000000000008, 8, false, -0.11072463768115942, -0.090996579290893664, 0.13043478260869565, false, false),
                new E(1, 1, 1, 0.31, 0.025000000000000001, 0.34999999999999998, 7.9400000000000004, 36, false, -0.03999999999999998, -0.07999999999999996, 0.43478260869565216, true, true),
                new E(2, 1, 1, 0.23664122137404578, 0.014567915622632713, 0.26000000000000001, 4.4500000000000002, 16, false, -0.02335877862595423, -0.085313289357027866, 0.27173913043478259, true, true),
                new E(2, 1, 2, 0.47328244274809156, 0.029135831245265427, 0.52000000000000002, 2.3799999999999999, 6, false, -0.046717557251908459, -0.085313289357027866, 0.27173913043478259, false, true),
                new E(2, 2, 1, 0.19135802469135801, 0.0095259868922420356, 0.20000000000000001, 3.6499999999999999, 11, false, -0.008641975308641997, -0.048010973936899987, 0.16304347826086957, true, true),
                new E(2, 2, 2, 0.38271604938271603, 0.019051973784484071, 0.40999999999999998, 3.4399999999999999, 7, false, -0.027283950617283947, -0.078379634062866826, 0.16304347826086957, true, true),
            },
            new[]
            {
                new P(1, 0, 0.56, 8.63, 8, false, false),
                new P(1, 1, 0.35, 7.94, 36, false, false),
                new P(2, 1, 0.26, 4.45, 16, false, false),
                new P(2, 1, 0.52, 2.38, 6, false, false),
                new P(2, 2, 0.2, 3.65, 11, false, false),
                new P(2, 2, 0.41, 3.44, 7, false, false),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D1_p0.62_new");
        yield return new TestCaseData("D1", 0.7, new OncologyAnalyzer.PurityPeakOptions { MixturePeaks = new Dictionary<(int Major, int Minor), IReadOnlyList<double>> { [(1, 0)] = new[] { 0.53894927536231885 }, [(1, 1)] = new[] { 0.35226677286160424 }, [(2, 1)] = new[] { 0.26285730659499063, 0.5165650460721285 }, [(2, 2)] = new[] { 0.41368298792735059, 0.19704338051414069 } } }, 0.0045358591466179345, true,
            new[]
            {
                new E(1, 0, 1, 0.53846153846153844, 0.059171597633136092, 0.53921175884249484, 7.3772290573174262, 9, true, -0.0007502203809564012, -0.00063331840442059511, 0.13043478260869565, true, true),
                new E(1, 1, 1, 0.34999999999999998, 0.025000000000000001, 0.34999999999999998, 7.9400000000000004, 36, false, 0, 0, 0.43478260869565216, true, true),
                new E(2, 1, 1, 0.25925925925925924, 0.013717421124828532, 0.26000000000000001, 4.4500000000000002, 16, false, -0.00074074074074076401, -0.0027054081108136012, 0.27173913043478259, true, true),
                new E(2, 1, 2, 0.51851851851851849, 0.027434842249657063, 0.51706837465494881, 2.3799330581481675, 6, true, 0.0014501438635696751, 0.0026377208132664556, 0.27173913043478259, true, true),
                new E(2, 2, 1, 0.20588235294117646, 0.0086505190311418692, 0.20000000000000001, 3.6499999999999999, 11, false, 0.0058823529411764497, 0.032679738562091387, 0.16304347826086957, true, true),
                new E(2, 2, 2, 0.41176470588235292, 0.017301038062283738, 0.41322467581493705, 3.4361817783104756, 7, true, -0.0014599699325841264, -0.004240335576719669, 0.16304347826086957, true, true),
            },
            new[]
            {
                new P(1, 0, 0.56, 8.63, 8, false, false),
                new P(1, 0, 0.539211758842495, 7.37722905731743, 9, false, true),
                new P(1, 1, 0.35, 7.94, 36, false, false),
                new P(1, 1, 0.352382250253623, 7.91997299371671, 36, false, true),
                new P(2, 1, 0.26, 4.45, 16, false, false),
                new P(2, 1, 0.52, 2.38, 6, false, false),
                new P(2, 1, 0.263421748285513, 4.43099152204969, 16, false, true),
                new P(2, 1, 0.517068374654949, 2.37993305814817, 6, false, true),
                new P(2, 2, 0.2, 3.65, 11, false, false),
                new P(2, 2, 0.41, 3.44, 7, false, false),
                new P(2, 2, 0.413224675814937, 3.43618177831048, 7, false, true),
                new P(2, 2, 0.196509650749174, 3.65123689868539, 11, false, true),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D1_p0.7_new_bmix");
        yield return new TestCaseData("D3", 0.45, new OncologyAnalyzer.PurityPeakOptions { MixturePeaks = new Dictionary<(int Major, int Minor), IReadOnlyList<double>> { [(1, 1)] = new[] { 0.08283459223334233, 0.22238674826433397 }, [(2, 0)] = new[] { 0.61731564204021494, 0.22663755735414512, 0.44882962304120372 } } }, 0.0037889790909820253, true,
            new[]
            {
                new E(1, 1, 1, 0.22500000000000001, 0.025000000000000001, 0.22216696460248936, 5.7754176545467608, 45, true, 0.002833035397510647, 0.0056660707950212941, 0.76923076923076927, true, true),
                new E(2, 0, 1, 0.22500000000000001, 0.025000000000000001, 0.22623399662790777, 1.1916068509032687, 3, true, -0.0012339966279077685, -0.0024679932558155371, 0.23076923076923078, true, true),
                new E(2, 0, 2, 0.45000000000000001, 0.050000000000000003, 0.45000000000000001, 2.8500000000000001, 9, false, 0, 0, 0.23076923076923078, true, true),
            },
            new[]
            {
                new P(1, 1, 0.08, 3.21, 18, false, false),
                new P(1, 1, 0.22, 5.78, 45, false, false),
                new P(1, 1, 0.0825021828612956, 3.2141006931568, 18, false, true),
                new P(1, 1, 0.222166964602489, 5.77541765454676, 45, false, true),
                new P(2, 0, 0.21, 1.21, 2, false, false),
                new P(2, 0, 0.45, 2.85, 9, false, false),
                new P(2, 0, 0.617959913926028, 1.96379706312039, 8, false, true),
                new P(2, 0, 0.226233996627908, 1.19160685090327, 3, false, true),
                new P(2, 0, 0.449614891616092, 2.8461656236106, 9, false, true),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D3_p0.45_new_bmix");
        yield return new TestCaseData("D3", 0.45, new OncologyAnalyzer.PurityPeakOptions { MixturePeaks = new Dictionary<(int Major, int Minor), IReadOnlyList<double>> { [(1, 1)] = new[] { 0.08283459223334233, 0.22238674826433397 }, [(2, 0)] = new[] { 0.61731564204021494, 0.22663755735414512, 0.44882962304120372 } }, MatchingStrategy = OncologyAnalyzer.PurityPeakMatchingStrategy.Rightmost }, -0.13824761798675933, true,
            new[]
            {
                new E(1, 1, 1, 0.22500000000000001, 0.025000000000000001, 0.22216696460248936, 5.7754176545467608, 45, true, 0.002833035397510647, 0.0056660707950212941, 0.76923076923076927, true, true),
                new E(2, 0, 2, 0.45000000000000001, 0.050000000000000003, 0.61795991392602811, 1.9637970631203869, 8, true, -0.1679599139260281, -0.1679599139260281, 0.23076923076923078, false, false),
                new E(2, 0, 1, 0.22500000000000001, 0.025000000000000001, 0.45000000000000001, 2.8500000000000001, 9, false, -0.22500000000000001, -0.45000000000000001, 0.23076923076923078, false, false),
            },
            new[]
            {
                new P(1, 1, 0.08, 3.21, 18, false, false),
                new P(1, 1, 0.22, 5.78, 45, false, false),
                new P(1, 1, 0.0825021828612956, 3.2141006931568, 18, false, true),
                new P(1, 1, 0.222166964602489, 5.77541765454676, 45, false, true),
                new P(2, 0, 0.21, 1.21, 2, false, false),
                new P(2, 0, 0.45, 2.85, 9, false, false),
                new P(2, 0, 0.617959913926028, 1.96379706312039, 8, false, true),
                new P(2, 0, 0.226233996627908, 1.19160685090327, 3, false, true),
                new P(2, 0, 0.449614891616092, 2.8461656236106, 9, false, true),
            }).SetName("AnalyzePurityPeaks_CNAqcR_D3_p0.45_new_bmix_rightmost");

    }

    [TestCaseSource(nameof(CnaqcReferenceCases))]
    public void AnalyzePurityPeaks_MatchesCnaqcR(
        string dataset, double purity, OncologyAnalyzer.PurityPeakOptions options, double score, bool pass, E[] matches, P[] peaks)
    {
        OncologyAnalyzer.PurityPeakAnalysis result = OncologyAnalyzer.AnalyzePurityPeaks(Mutations(dataset), purity, options);

        Assert.That(result.Pass, Is.EqualTo(pass), "sample QC (weight-majority of karyotype PASS/FAIL)");
        Assert.That(result.Score, Is.EqualTo(score).Within(1e-12), "λ = Σ weight·offset");
        Assert.That(result.Matches, Has.Count.EqualTo(matches.Length));
        for (int i = 0; i < matches.Length; i++)
        {
            E e = matches[i];
            OncologyAnalyzer.PurityPeakMatch a = result.Matches[i];
            string at = $"match {i} ({e.Major}:{e.Minor}, m = {e.M})";
            Assert.That((a.MajorCopyNumber, a.MinorCopyNumber, a.Multiplicity), Is.EqualTo((e.Major, e.Minor, e.M)), at);
            Assert.That(a.ExpectedPeak, Is.EqualTo(e.Peak).Within(1e-15), at + " expected peak");
            Assert.That(a.DeltaVaf, Is.EqualTo(e.DeltaVaf).Within(1e-15), at + " δ");
            Assert.That(a.MatchedPeak.X, Is.EqualTo(e.X).Within(1e-12), at + " x");
            Assert.That(a.MatchedPeak.Y, Is.EqualTo(e.Y).Within(1e-12 * Math.Max(1, e.Y)), at + " y");
            Assert.That(a.MatchedPeak.CountsPerBin, Is.EqualTo(e.Counts), at + " counts_per_bin");
            Assert.That(a.MatchedPeak.Source == OncologyAnalyzer.PurityPeakSource.Mixture, Is.EqualTo(e.Mixture), at + " from");
            Assert.That(a.OffsetVaf, Is.EqualTo(e.OffsetVaf).Within(1e-12), at + " offset_VAF");
            Assert.That(a.Offset, Is.EqualTo(e.Offset).Within(1e-11), at + " offset");
            Assert.That(a.Weight, Is.EqualTo(e.Weight).Within(1e-15), at + " weight");
            Assert.That(a.Matched, Is.EqualTo(e.Matched), at + " matched");
            OncologyAnalyzer.PurityPeakKaryotype k = result.Karyotypes.Single(
                q => q.MajorCopyNumber == e.Major && q.MinorCopyNumber == e.Minor);
            Assert.That(k.Pass, Is.EqualTo(e.KaryotypePass), at + " karyotype QC");
        }

        var actualPeaks = result.Karyotypes.SelectMany(k => k.Peaks.Select(p => (k.MajorCopyNumber, k.MinorCopyNumber, p))).ToList();
        Assert.That(actualPeaks, Has.Count.EqualTo(peaks.Length), "peak table size");
        for (int i = 0; i < peaks.Length; i++)
        {
            P e = peaks[i];
            (int major, int minor, OncologyAnalyzer.PurityDataPeak p) = actualPeaks[i];
            string at = $"peak {i} ({e.Major}:{e.Minor})";
            Assert.That((major, minor), Is.EqualTo((e.Major, e.Minor)), at);
            // R prints mixture-peak coordinates with 15 significant digits.
            Assert.That(p.X, Is.EqualTo(e.X).Within(1e-12), at + " x");
            Assert.That(p.Y, Is.EqualTo(e.Y).Within(1e-12 * Math.Max(1, e.Y)), at + " y");
            Assert.That(p.CountsPerBin, Is.EqualTo(e.Counts), at + " counts");
            Assert.That(p.Discarded, Is.EqualTo(e.Discarded), at + " discarded");
            Assert.That(p.Source == OncologyAnalyzer.PurityPeakSource.Mixture, Is.EqualTo(e.Mixture), at + " from");
        }
    }

    // CNAqc: with min_absolute_karyotype_mutations = 200 no D4 karyotype qualifies → "No karyotypes satisfy input data
    // filters", the object is returned without peaks_analysis.
    [Test]
    public void AnalyzePurityPeaks_NoKaryotypePassesFilters_ReturnsNotAnalysed()
    {
        var result = OncologyAnalyzer.AnalyzePurityPeaks(Mutations("D4"), 0.3,
            new OncologyAnalyzer.PurityPeakOptions { MinAbsoluteKaryotypeMutations = 200 });
        Assert.That(result.Pass, Is.Null);
        Assert.That(result.Score, Is.NaN);
        Assert.That(result.Karyotypes, Is.Empty);
        Assert.That(result.Matches, Is.Empty);
    }

    // CNAqc's QC is the weight-majority of karyotype verdicts, not |λ| ≤ ε: at π = 0.62 the 1:0 karyotype FAILs and the
    // 2:1 m = 2 band misses, yet the sample PASSes with λ = −0.1136 (|λ| > ε = 0.05). (R-locked in the case source.)
    [Test]
    public void AnalyzePurityPeaks_ScoreIsNotThresholded()
    {
        var result = OncologyAnalyzer.AnalyzePurityPeaks(Mutations("D1"), 0.62);
        Assert.That(Math.Abs(result.Score), Is.GreaterThan(OncologyAnalyzer.PurityPeakOptions.Default.PurityError));
        Assert.That(result.Pass, Is.True);
        Assert.That(result.Karyotypes.Single(k => k.MajorCopyNumber == 1 && k.MinorCopyNumber == 0).Pass, Is.False);
    }

    [Test]
    public void AnalyzePurityPeaks_KaryotypeScoresSumToSampleScore()
    {
        var result = OncologyAnalyzer.AnalyzePurityPeaks(Mutations("D1"), 0.5);
        Assert.That(result.Karyotypes.Sum(k => k.Score), Is.EqualTo(result.Score).Within(1e-15));
        Assert.That(result.Karyotypes.Sum(k => k.Weight), Is.EqualTo(1.0).Within(1e-15));
        Assert.That(result.Karyotypes.Select(k => k.MutationCount), Is.EqualTo(new[] { 120, 400, 250, 150 }));
    }

    // Non-simple karyotypes (3:1, 60 mutations in D1) are not analysed but count towards N for MinKaryotypeSize:
    // 1:0 has 120/980 = 0.1224 of all mutations.
    [Test]
    public void AnalyzePurityPeaks_MinKaryotypeSize_UsesAllKaryotypesAsDenominator()
    {
        var keep = OncologyAnalyzer.AnalyzePurityPeaks(Mutations("D1"), 0.7,
            new OncologyAnalyzer.PurityPeakOptions { MinKaryotypeSize = 120.0 / 980 });
        var drop = OncologyAnalyzer.AnalyzePurityPeaks(Mutations("D1"), 0.7,
            new OncologyAnalyzer.PurityPeakOptions { MinKaryotypeSize = 120.0 / 980 + 1e-9 });
        Assert.That(keep.Karyotypes.Any(k => k.MajorCopyNumber == 1 && k.MinorCopyNumber == 0), Is.True);
        Assert.That(drop.Karyotypes.Any(k => k.MajorCopyNumber == 1 && k.MinorCopyNumber == 0), Is.False);
    }

    [Test]
    public void AnalyzePurityPeaks_InvalidArguments_Throw()
    {
        var m = Mutations("D4");
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.AnalyzePurityPeaks(null!, 0.3));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzePurityPeaks(m, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzePurityPeaks(m, 1.01));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzePurityPeaks(m, 0.3,
            new OncologyAnalyzer.PurityPeakOptions { PurityError = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzePurityPeaks(m, 0.3,
            new OncologyAnalyzer.PurityPeakOptions { MinKaryotypeSize = 1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzePurityPeaks(m, 0.3,
            new OncologyAnalyzer.PurityPeakOptions { KernelAdjust = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzePurityPeaks(m, 0.3,
            new OncologyAnalyzer.PurityPeakOptions { MinAbsoluteKaryotypeMutations = -1 }));
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.AnalyzePurityPeaks(m, 0.3,
            new OncologyAnalyzer.PurityPeakOptions { Karyotypes = new[] { (3, 1) } }));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzePurityPeaks(
            new[] { new OncologyAnalyzer.PurityPeakMutation(1.2, 1, 1) }, 0.3));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzePurityPeaks(
            new[] { new OncologyAnalyzer.PurityPeakMutation(0.2, -1, 1) }, 0.3));
    }
}
