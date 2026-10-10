// ONCO-PURITY-001 — CNAqc complex-karyotype and subclonal peak QC (FIN-B24 F62):
//   OncologyAnalyzer.AnalyzeComplexKaryotypePeaks (CNAqc analyze_peaks_general) and
//   OncologyAnalyzer.AnalyzeSubclonalPurityPeaks (CNAqc analyze_peaks_subclonal + expectations_subclonal).
// Evidence: docs/Evidence/ONCO-PURITY-001-Evidence.md (§ "CNAqc analyze_peaks_general / analyze_peaks_subclonal")
// TestSpec: tests/TestSpecs/ONCO-PURITY-001.md (§ 5.9)
// Reference: caravagnalab/CNAqc 1.1.5 (commit 4b7cea4a) R/peak_algorithms.R, R/equations.R, R/analyze_peaks.R run in
//            R 4.3.3 (R >= 4.4 density lattice); outputs printed with sprintf("%.17g") into
//            TestData/CNAqc/cnaqc_general_R.txt and cnaqc_subclonal_R.txt (segment ids reduced to the chromosome). The C#
//            results are formatted the same way and must equal the R text line by line (bit-identical doubles, the
//            same random 8-letter mutation identifiers from set.seed).

using System.Globalization;
using Seqeron.Genomics.Oncology;
using Seqeron.Genomics.Tests.Helpers;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_CnaqcGeneralSubclonalPeaks_Tests
{
    private static string F(double v) => v.ToString("G17", CultureInfo.InvariantCulture);

    private static string B(bool b) => b ? "TRUE" : "FALSE";

    private static List<OncologyAnalyzer.PurityPeakMutation> Mutations(string dataset) =>
        CnaqcTestData.Load(dataset).Select(r => new OncologyAnalyzer.PurityPeakMutation(r.Vaf, r.Major, r.Minor)).ToList();

    // The block of the R reference file that starts with "== <tag>".
    private static List<string> ReferenceBlock(string file, string tag)
    {
        var lines = CnaqcTestData.Text(file).Split('\n').Select(l => l.TrimEnd('\r')).ToList();
        int start = lines.FindIndex(l => l == $"== {tag}" || l.StartsWith($"== {tag} ", StringComparison.Ordinal));
        Assert.That(start, Is.GreaterThanOrEqualTo(0), $"tag {tag} not in {file}");
        int end = lines.FindIndex(start + 1, l => l.StartsWith("== ", StringComparison.Ordinal));
        return lines.GetRange(start, (end < 0 ? lines.Count : end) - start).Where(l => l.Length > 0).ToList();
    }

    #region analyze_peaks_general

    private static IEnumerable<TestCaseData> GeneralCases()
    {
        static TestCaseData C(string tag, string ds, double p, OncologyAnalyzer.PurityPeakOptions o) =>
            new TestCaseData(tag, ds, p, o).SetName($"AnalyzeComplexKaryotypePeaks_CNAqcR_{tag.Replace(' ', '_')}");
        yield return C("G1 p=0.6", "G1", 0.6, new());
        yield return C("G1 p=0.6 eps=0.01", "G1", 0.6, new() { PurityError = 0.01 });
        yield return C("G1 p=0.6 adj=0.5", "G1", 0.6, new() { KernelAdjust = 0.5 });
        yield return C("G1 p=0.75", "G1", 0.75, new());
        yield return C("G2 p=0.4", "G2", 0.4, new());
        yield return C("G2 p=0.55 minVAF=0.05", "G2", 0.55, new() { MinVaf = 0.05 });
        yield return C("G3 p=0.85", "G3", 0.85, new());
        yield return C("G3 p=0.85 minabs=101", "G3", 0.85, new() { MinAbsoluteKaryotypeMutations = 101 });
        yield return C("G3 p=0.85 minabs=150", "G3", 0.85, new() { MinAbsoluteKaryotypeMutations = 150 });
        yield return C("G3 p=0.85 minabs=60 k=1:1", "G3", 0.85, new() { MinAbsoluteKaryotypeMutations = 60, Karyotypes = new[] { (1, 1) } });
    }

    [TestCaseSource(nameof(GeneralCases))]
    public void AnalyzeComplexKaryotypePeaks_MatchesCnaqcR(string tag, string dataset, double purity, OncologyAnalyzer.PurityPeakOptions options)
    {
        var r = OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(Mutations(dataset), purity, options);
        var actual = new List<string>();
        if (!r.Ran)
        {
            actual.Add($"== {tag} NOT-RUN");
        }
        else
        {
            var analysis = r.ExpectedPeaks.Select(e => $"{e.MajorCopyNumber}:{e.MinorCopyNumber}").Distinct().ToList();
            actual.Add($"== {tag} analysis={string.Join(",", analysis)}");
            foreach (var e in r.ExpectedPeaks)
            {
                int n = r.Karyotypes.Single(k => k.MajorCopyNumber == e.MajorCopyNumber && k.MinorCopyNumber == e.MinorCopyNumber).MutationCount;
                actual.Add($"  E {e.MajorCopyNumber}:{e.MinorCopyNumber} m={e.Multiplicity} peak={F(e.ExpectedPeak)} matched={B(e.Matched)} n={n}");
            }

            // R data_peaks follow group_split(karyotype) order (C collation).
            foreach (var k in r.Karyotypes.OrderBy(x => $"{x.MajorCopyNumber}:{x.MinorCopyNumber}", StringComparer.Ordinal))
            {
                foreach (var d in k.Peaks)
                    actual.Add($"  D {k.MajorCopyNumber}:{k.MinorCopyNumber} x={F(d.X)} y={F(d.Y)} cpb={d.CountsPerBin?.ToString(CultureInfo.InvariantCulture) ?? "NA"} disc={B(d.Discarded)}");
            }

            foreach (var s in r.Karyotypes)
                actual.Add($"  S {s.MajorCopyNumber}:{s.MinorCopyNumber} n={s.MutationCount} matched={s.MatchedPeaks} mismatched={s.MismatchedPeaks} prop={F(s.MatchedProportion)}");
        }

        Assert.That(actual, Is.EqualTo(ReferenceBlock("cnaqc_general_R.txt", tag)));
    }

    // G1, π = 0.6 (R): 3:0 expects 3 peaks (m = 1..3, minor 0 → Major multiplicities), 3:1 → 1..3, 4:1 → 1..4; the
    // 3:2 karyotype (80 mutations < 100) and the simple 1:1 are not analysed. prop = 1, 2/3, 1/2 → QC PASS for all three.
    [Test]
    public void AnalyzeComplexKaryotypePeaks_G1_SummaryAndExpectations()
    {
        var r = OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(Mutations("G1"), 0.6);
        Assert.That(r.Ran, Is.True);
        Assert.That(r.Karyotypes.Select(k => (k.MajorCopyNumber, k.MinorCopyNumber)), Is.EqualTo(new[] { (3, 0), (3, 1), (4, 1) }));
        Assert.That(r.Karyotypes.Select(k => k.MatchedProportion), Is.EqualTo(new[] { 1.0, 0.66666666666666663, 0.5 }));
        Assert.That(r.Karyotypes.Select(k => k.Pass), Is.EqualTo(new[] { true, true, true }));
        Assert.That(r.Karyotypes.Select(k => k.MutationCount), Is.EqualTo(new[] { 200, 250, 150 }));
        var k31 = r.ExpectedPeaks.Where(e => e.MajorCopyNumber == 3 && e.MinorCopyNumber == 1).ToList();
        Assert.That(k31.Select(e => e.ExpectedPeak), Is.EqualTo(new[] { 0.18749999999999997, 0.37499999999999994, 0.56249999999999989 }));
        Assert.That(k31.Select(e => e.Matched), Is.EqualTo(new[] { true, false, true }));
    }

    // ε only widens the strict |x − peak| < ε window: at ε = 0.01 the 4:1 m = 4 peak (0.6316 vs data 0.62) is lost → 0.25 → FAIL.
    [Test]
    public void AnalyzeComplexKaryotypePeaks_SmallEpsilon_FailsKaryotype()
    {
        var r = OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(Mutations("G1"), 0.6, new() { PurityError = 0.01 });
        var k41 = r.Karyotypes.Single(k => k.MajorCopyNumber == 4 && k.MinorCopyNumber == 1);
        Assert.That(k41.MatchedProportion, Is.EqualTo(0.25));
        Assert.That(k41.Pass, Is.False);
    }

    // Gate (analyze_peaks): strict n > 100 for a karyotype outside `karyotypes`; analysis then uses n ≥ 100 (G3 4:0 with
    // exactly 100 mutations is analysed). Raising the minimum to 150 closes the gate.
    [Test]
    public void AnalyzeComplexKaryotypePeaks_GateIsStrict_AnalysisInclusive()
    {
        var r = OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(Mutations("G3"), 0.85);
        Assert.That(r.Karyotypes.Any(k => k.MajorCopyNumber == 4 && k.MinorCopyNumber == 0 && k.MutationCount == 100), Is.True);
        var closed = OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(Mutations("G3"), 0.85, new() { MinAbsoluteKaryotypeMutations = 150 });
        Assert.That(closed.Ran, Is.False);
        Assert.That(closed.Karyotypes, Is.Empty);
    }

    // Only simple karyotypes: the gate opens when `karyotypes` omits one with n > 100 (R analyze_peaks), but the general
    // analysis has nothing to analyse — CNAqc stops with an R error (`1:nrow(NULL)`).
    [Test]
    public void AnalyzeComplexKaryotypePeaks_GateOpenWithoutComplexKaryotype_Throws()
    {
        var simpleOnly = Mutations("D1").Where(m => m.MajorCopyNumber + m.MinorCopyNumber <= 4 && m.MajorCopyNumber <= 2).ToList();
        Assert.Throws<InvalidOperationException>(() => OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(simpleOnly, 0.7,
            new() { Karyotypes = new[] { (1, 1) } }));
        Assert.That(OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(simpleOnly, 0.7).Ran, Is.False);
    }

    [Test]
    public void AnalyzeComplexKaryotypePeaks_InvalidArguments_Throw()
    {
        var m = Mutations("G1");
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(null!, 0.6));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(m, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(m, 0.6, new() { PurityError = 1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(
            new[] { new OncologyAnalyzer.PurityPeakMutation(-0.1, 3, 1) }, 0.6));
    }

    #endregion

    #region analyze_peaks_subclonal

    private static IEnumerable<TestCaseData> SubclonalCases()
    {
        static TestCaseData C(string tag, string ds, double p, OncologyAnalyzer.SubclonalPeakOptions o, int take = int.MaxValue) =>
            new TestCaseData(tag, ds, p, o, take).SetName($"AnalyzeSubclonalPurityPeaks_CNAqcR_{tag.Replace(' ', '_')}");
        yield return C("S1 p=0.7 seed=1", "S1", 0.7, new() { Seed = 1 });
        yield return C("S1 p=0.7 seed=99 eps=0.02", "S1", 0.7, new() { Seed = 99, Epsilon = 0.02 });
        yield return C("S1 p=0.7 seed=5 nmin=150 adj=0.5", "S1", 0.7, new() { Seed = 5, MinMutations = 150, KernelAdjust = 0.5 });
        yield return C("S2 p=0.5 seed=2", "S2", 0.5, new() { Seed = 2 });
        yield return C("S2 p=0.5 seed=3 start=2:2", "S2", 0.5, new() { Seed = 3, StartingState = (2, 2) }, 2);
        yield return C("S3 p=0.8 seed=4", "S3", 0.8, new() { Seed = 4 });
        yield return C("S3 p=0.8 seed=6 start=1:0", "S3", 0.8, new() { Seed = 6, StartingState = (1, 0) });
        yield return C("S3 p=0.8 seed=4 nmin=300", "S3", 0.8, new() { Seed = 4, MinMutations = 300 });
    }

    [TestCaseSource(nameof(SubclonalCases))]
    public void AnalyzeSubclonalPurityPeaks_MatchesCnaqcR(string tag, string dataset, double purity, OncologyAnalyzer.SubclonalPeakOptions options, int take)
    {
        var r = OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(CnaqcTestData.Segments(dataset).Take(take), purity, options);
        var actual = new List<string>();
        if (r.Count == 0)
        {
            actual.Add($"== {tag} NO-SUBCLONAL");
        }
        else
        {
            static string K(string? g)
            {
                if (g is null) return "NA";
                int a = g.Count(c => c == 'A'), b = g.Count(c => c == 'B');
                return $"{Math.Max(a, b)}:{Math.Min(a, b)}";
            }

            static string M(OncologyAnalyzer.SubclonalEvolutionModel m) => m == OncologyAnalyzer.SubclonalEvolutionModel.Linear ? "linear" : "branching";
            actual.Add($"== {tag}");
            foreach (var s in r)
            {
                foreach (var e in s.ExpectedPeaks)
                    actual.Add($"  E {s.Segment.Chromosome}|{M(e.Model)}|{e.ModelId}|{e.MutationId} n1={e.FirstCloneCopies} n2={e.SecondCloneCopies} " +
                               $"k1={K(e.FirstCloneGenotype)} g1={e.FirstCloneGenotype ?? "NA"} k2={K(e.SecondCloneGenotype)} g2={e.SecondCloneGenotype ?? "NA"} " +
                               $"role={(e.Shared ? "shared" : "private")} peak={F(e.ExpectedPeak)} matched={B(e.Matched)}");
            }

            foreach (var s in r)
            {
                foreach (var d in s.Peaks)
                    actual.Add($"  D {s.Segment.Chromosome} x={F(d.X)} y={F(d.Y)} cpb={d.CountsPerBin?.ToString(CultureInfo.InvariantCulture) ?? "NA"} disc={B(d.Discarded)}");
            }

            foreach (var s in r)
            {
                foreach (var b in s.BestModels)
                    actual.Add($"  S {s.Segment.Chromosome}|{b.ModelId}|{M(b.Model)} prop={F(b.MatchedProportion)}");
            }
        }

        Assert.That(actual, Is.EqualTo(ReferenceBlock("cnaqc_subclonal_R.txt", tag)));
    }

    // S1 chr1 (2:1 at CCF 0.6 / 1:1 at 0.4, π 0.7): three models, the branching A1B1 → A1A2B1 | A1B1 and the linear
    // A1B1 → A1A2B1 → A2B1 tie at 3/4 matched peaks; peaks 0.1157 (private to the 1:1 clone), 0.1736, 0.2893, 0.4628.
    [Test]
    public void AnalyzeSubclonalPurityPeaks_S1Chr1_ModelsAndPeaks()
    {
        var r = OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(CnaqcTestData.Segments("S1"), 0.7, new() { Seed = 1 });
        Assert.That(r.Select(s => s.Segment.Chromosome), Is.EqualTo(new[] { "chr1", "chr3", "chr5" })); // chr7: n = 100, not > 100
        var chr1 = r[0];
        Assert.That(chr1.Rankings.Select(m => (m.ModelId, m.MatchedProportion)), Is.EqualTo(new[]
        {
            ("A1B1 -> A1A2B1 -> A2B1", 0.75),
            ("A1B1 -> A1A2B1 | A1B1", 0.75),
            ("A1B1 -> A1B1 -> A1A2B1", 2.0 / 3),
        }));
        Assert.That(chr1.BestModels.Select(m => m.Model), Is.EqualTo(new[]
        {
            OncologyAnalyzer.SubclonalEvolutionModel.Linear, OncologyAnalyzer.SubclonalEvolutionModel.Branching,
        }));
        var branching = chr1.ExpectedPeaks.Where(e => e.Model == OncologyAnalyzer.SubclonalEvolutionModel.Branching).ToList();
        Assert.That(branching.Select(e => e.ExpectedPeak), Is.EqualTo(new[]
            { 0.11570247933884296, 0.17355371900826447, 0.28925619834710742, 0.46280991735537186 }));
        Assert.That(branching.Select(e => e.Shared), Is.EqualTo(new[] { false, false, true, true }));
    }

    // Peaks, matches and decisions do not depend on the seed (only the mutation identifiers do).
    [Test]
    public void AnalyzeSubclonalPurityPeaks_Seed_OnlyChangesIdentifiers()
    {
        var a = OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(CnaqcTestData.Segments("S2"), 0.5, new() { Seed = 2 });
        var b = OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(CnaqcTestData.Segments("S2"), 0.5, new() { Seed = 12345 });
        for (int i = 0; i < a.Count; i++)
        {
            Assert.That(b[i].ExpectedPeaks.Select(e => (e.ModelId, e.ExpectedPeak, e.Matched)),
                Is.EqualTo(a[i].ExpectedPeaks.Select(e => (e.ModelId, e.ExpectedPeak, e.Matched))));
            Assert.That(b[i].BestModels, Is.EqualTo(a[i].BestModels));
        }

        Assert.That(b[0].ExpectedPeaks[0].MutationId, Is.Not.EqualTo(a[0].ExpectedPeaks[0].MutationId));
    }

    // A LOH starting state (1:0) cannot reach a segment without LOH: CNAqc aborts that segment's expectations, so it
    // keeps its data peaks but has no model (R: S3 chr12 2:2/1:1 and chr17 1:1/2:0 are missing from the summary).
    [Test]
    public void AnalyzeSubclonalPurityPeaks_LohStart_UnreachableSegmentHasNoModels()
    {
        var r = OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(CnaqcTestData.Segments("S3"), 0.8, new() { Seed = 6, StartingState = (1, 0) });
        Assert.That(r.Select(s => s.ExpectedPeaks.Count > 0), Is.EqualTo(new[] { true, false, false }));
        Assert.That(r[1].Peaks, Is.Not.Empty);
        Assert.That(r[1].BestModels, Is.Empty);
    }

    // CNAqc's evolve() never terminates for 2:2 → 1:0 (every deletion exceeds the ploidy cap 2·1): reported, not hung.
    [Test]
    public void AnalyzeSubclonalPurityPeaks_UnreachableUnderCap_ThrowsInsteadOfLooping()
    {
        var s2 = CnaqcTestData.Segments("S2"); // chr8 is 1:0 / 2:0
        Assert.Throws<InvalidOperationException>(() => OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(s2, 0.5, new() { StartingState = (2, 2) }));
    }

    [Test]
    public void AnalyzeSubclonalPurityPeaks_InvalidArguments_Throw()
    {
        var s = CnaqcTestData.Segments("S1");
        var seg = s[0];
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(null!, 0.7));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(s, 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(s, 0.7, new() { Epsilon = 0 }));
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(s, 0.7, new() { StartingState = (1, 2) }));
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(new[] { seg with { MajorCopyNumber = 3 } }, 0.7));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(new[] { seg with { Ccf = 1.0 } }, 0.7));
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(new[] { seg with { Vafs = new[] { 1.5 } } }, 0.7));
    }

    #endregion
}
