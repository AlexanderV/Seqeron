// ONCO-CNA-002 — GISTIC2 ziggurat deconstruction (B24 F50–F52)
// Evidence: docs/Evidence/ONCO-CNA-002-Evidence.md (§ "GISTIC2 ziggurat deconstruction — Octave runs")
// TestSpec: tests/TestSpecs/ONCO-CNA-002.md (C16+)
// Source: broadinstitute/gistic2 master 26c590bd (+ snputil cf3172b8): make_sample_B.m, normalize_by_arm_length.m,
//         deconstruct_sample.m, deconstruct_chr.m, prepare_B.m, merge_adj_segs.m, atomic_zigg_deconstruction.m,
//         add_broad_levels_to_zigg.m. Expected values = GNU Octave 8.4 running the ORIGINAL .m files (printf %.17g).

using System.Globalization;
using ZChr = Seqeron.Genomics.Oncology.OncologyAnalyzer.ZigguratChromosome;
using Row = Seqeron.Genomics.Oncology.OncologyAnalyzer.GisticZiggRow;
using ZSeg = Seqeron.Genomics.Oncology.OncologyAnalyzer.ZigguratSegment;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_ZigguratDeconstruction_Tests
{
    // Octave layout [10 10; 8 12]: chr "1" = 10 p + 10 q markers, chr "2" = 8 p + 12 q markers (global 21..40).
    private static readonly ZChr[] Layout = { new("1", 10, 10), new("2", 8, 12) };

    private static OncologyAnalyzer.GisticMarkerLayout MarkerLayout() => new(Layout);

    private static ZSeg S(string chr, int start, int end, double value) => new(chr, start, end, value);

    private static double[][] Parse(string block) =>
        block.Split('\n', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
            .Select(line => line.Split(' ', StringSplitOptions.RemoveEmptyEntries)
                .Select(t => double.Parse(t, CultureInfo.InvariantCulture)).ToArray())
            .ToArray();

    private static double[] BColumns(Row r) =>
        new double[] { r.Chromosome, r.Start, r.End, r.Amplitude, r.Sample, r.Fraction };

    private static double[] ZColumns(Row r) =>
        new double[] { r.Chromosome, r.Start, r.End, r.Amplitude, r.Sample, r.StartLevel, r.EndLevel, r.Fraction };

    private static void AssertRows(IReadOnlyList<Row> actual, Func<Row, double[]> columns, string expected, string what)
    {
        double[][] want = expected.Length == 0 ? Array.Empty<double[]>() : Parse(expected);
        Assert.That(actual, Has.Count.EqualTo(want.Length), $"{what}: row count (Octave)");
        for (int i = 0; i < want.Length; i++)
        {
            Assert.That(columns(actual[i]), Is.EqualTo(want[i]), $"{what}: row {i + 1} must equal Octave bit-for-bit");
        }
    }

    #region F50 — make_sample_B / deconstruct_sample (Octave-locked)

    // Case A — interrupted arm gain 0.5 | 1.5 | 0.5 on 1p, 1q flat 0, chr 2 flat 0.2.
    [Test]
    public void MakeSampleB_And_DeconstructSample_InterruptedArmGain_MatchOctave()
    {
        var layout = MarkerLayout();
        var segments = new[] { S("1", 1, 3, 0.5), S("1", 4, 6, 1.5), S("1", 7, 10, 0.5), S("1", 11, 20, 0), S("2", 1, 20, 0.2) };

        List<Row> b = OncologyAnalyzer.GisticMakeSampleB(layout, segments, 1);
        AssertRows(b, BColumns, """
            1 1 3 0.5 1 0.29999999999999999
            1 4 6 1.5 1 0.29999999999999999
            1 7 10 0.5 1 0.40000000000000002
            1 11 20 0 1 1
            2 21 40 0.20000000000000001 1 2
            """, "A.B");

        var (za0, zd0) = OncologyAnalyzer.GisticDeconstructSample(b, new double[4], new[] { 20, 40 });
        AssertRows(za0, ZColumns, """
            1 4 6 1 1 0.5 1.5 0.29999999999999999
            1 1 10 0.5 1 0 0.5 1
            2 21 40 0.20000000000000001 1 0 0.20000000000000001 2
            """, "A0.ZA");
        AssertRows(zd0, ZColumns, "", "A0.ZD");

        // Against broad levels 1p 0.5 / 1q 0 (breakpoint = marker 10) and chr 2 0.2: only the +1.0 focal step remains.
        var (za1, zd1) = OncologyAnalyzer.GisticDeconstructSample(b, new[] { 0.5, 0, 0.2, 0.2 }, new[] { 10, 40 });
        AssertRows(za1, ZColumns, "1 4 6 1 1 0.5 1.5 0.29999999999999999", "A1.ZA");
        AssertRows(zd1, ZColumns, "", "A1.ZD");
    }

    // Case B — deletion inside a gain on 1q (0.6 | −0.8 | 0.6), chr 2 flat −0.3; includes the reference quirk that a
    // breakpoint on the chromosome's first row deconstructs every row against the q level.
    [Test]
    public void MakeSampleB_And_DeconstructSample_DeletionInsideGain_MatchOctave()
    {
        var layout = MarkerLayout();
        var segments = new[] { S("1", 1, 10, 0), S("1", 11, 14, 0.6), S("1", 15, 16, -0.8), S("1", 17, 20, 0.6), S("2", 1, 20, -0.3) };

        List<Row> b = OncologyAnalyzer.GisticMakeSampleB(layout, segments, 1);
        AssertRows(b, BColumns, """
            1 1 10 0 1 1
            1 11 14 0.59999999999999998 1 0.40000000000000002
            1 15 16 -0.80000000000000004 1 0.20000000000000001
            1 17 20 0.59999999999999998 1 0.40000000000000002
            2 21 40 -0.29999999999999999 1 2
            """, "B.B");

        var (za0, zd0) = OncologyAnalyzer.GisticDeconstructSample(b, new double[4], new[] { 20, 40 });
        AssertRows(za0, ZColumns, """
            1 11 14 0.59999999999999998 1 0 0.59999999999999998 0.40000000000000002
            1 17 20 0.59999999999999998 1 0 0.59999999999999998 0.40000000000000002
            """, "B0.ZA");
        AssertRows(zd0, ZColumns, """
            1 15 16 0.80000000000000004 1 0 0.80000000000000004 0.20000000000000001
            2 21 40 0.29999999999999999 1 0 0.29999999999999999 2
            """, "B0.ZD");

        var (za1, zd1) = OncologyAnalyzer.GisticDeconstructSample(b, new[] { 0, 0.6, -0.3, -0.3 }, new[] { 10, 40 });
        AssertRows(za1, ZColumns, "", "B1.ZA");
        AssertRows(zd1, ZColumns, """
            1 15 16 1.3999999999999999 1 0.59999999999999998 2 0.20000000000000001
            1 1 10 0.59999999999999998 1 0.59999999999999998 1.2 1
            """, "B1.ZD");
    }

    // Case C — centromere-spanning gain on chr 2 (local markers 6..12: 3 of 8 p + 4 of 12 q = 0.7083…), sample 2;
    // whole-chromosome segments have fraction 2 (p + q, ref_length 2).
    [Test]
    public void MakeSampleB_CentromereSpanningSegment_SumsArmFractions_MatchOctave()
    {
        var layout = MarkerLayout();
        var segments = new[] { S("1", 1, 20, 0.3), S("2", 1, 5, 0), S("2", 6, 12, 1.0), S("2", 13, 20, 0) };

        List<Row> b = OncologyAnalyzer.GisticMakeSampleB(layout, segments, 2);
        AssertRows(b, BColumns, """
            1 1 20 0.29999999999999999 2 2
            2 21 25 0 2 0.625
            2 26 32 1 2 0.70833333333333326
            2 33 40 0 2 0.66666666666666663
            """, "C.B");

        var (za0, zd0) = OncologyAnalyzer.GisticDeconstructSample(b, new double[4], new[] { 20, 40 });
        AssertRows(za0, ZColumns, """
            1 1 20 0.29999999999999999 2 0 0.29999999999999999 2
            2 26 32 1 2 0 1 0.70833333333333326
            """, "C0.ZA");
        AssertRows(zd0, ZColumns, "", "C0.ZD");
    }

    // Case D — breakpoint on the first row of chr 2 (0.4 | 1.2 | −0.5) with p 0.7 / q 0.4: all rows go to q.
    [Test]
    public void DeconstructSample_BreakpointOnFirstRow_TreatsChromosomeAsQ_MatchOctave()
    {
        var layout = MarkerLayout();
        var segments = new[] { S("1", 1, 20, 0.1), S("2", 1, 4, 0.4), S("2", 5, 9, 1.2), S("2", 10, 20, -0.5) };

        List<Row> b = OncologyAnalyzer.GisticMakeSampleB(layout, segments, 1);
        AssertRows(b, BColumns, """
            1 1 20 0.10000000000000001 1 2
            2 21 24 0.40000000000000002 1 0.5
            2 25 29 1.2 1 0.58333333333333337
            2 30 40 -0.5 1 0.91666666666666663
            """, "D.B");

        var (za, zd) = OncologyAnalyzer.GisticDeconstructSample(b, new[] { 0.1, 0.1, 0.7, 0.4 }, new[] { 20, 24 });
        AssertRows(za, ZColumns, "2 25 29 0.79999999999999993 1 0.40000000000000002 1.2 0.58333333333333337", "D1.ZA");
        AssertRows(zd, ZColumns, "2 30 40 0.90000000000000002 1 0.40000000000000002 1.3 0.91666666666666663", "D1.ZD");
    }

    // make_sample_B merges adjacent input segments with equal values (diff(D.dat) == 0 ⇒ no breakpoint).
    [Test]
    public void MakeSampleB_EqualAdjacentValues_Merged()
    {
        List<Row> b = OncologyAnalyzer.GisticMakeSampleB(
            MarkerLayout(), new[] { S("1", 1, 4, 0.5), S("1", 5, 12, 0.5), S("1", 13, 20, 0), S("2", 1, 20, 0) }, 1);

        Assert.That(b.Select(r => (r.Start, r.End, r.Fraction)), Is.EqualTo(new[] { (1, 12, 1.2), (13, 20, 0.8), (21, 40, 2.0) }),
            "Equal adjacent values form one B row (1..12 spans the centromere: 10/10 + 2/10).");
    }

    [Test]
    public void MakeSampleB_InvalidTiling_Throws()
    {
        var layout = MarkerLayout();
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.GisticMakeSampleB(layout, new[] { S("1", 1, 20, 0) }, 1),
                NUnit.Framework.Throws.ArgumentException, "chromosome 2 not covered");
            Assert.That(() => OncologyAnalyzer.GisticMakeSampleB(layout, new[] { S("1", 1, 10, 0), S("1", 12, 20, 0), S("2", 1, 20, 0) }, 1),
                NUnit.Framework.Throws.ArgumentException, "gap");
            Assert.That(() => OncologyAnalyzer.GisticMakeSampleB(layout, new[] { S("1", 1, 21, 0), S("2", 1, 20, 0) }, 1),
                NUnit.Framework.Throws.ArgumentException, "beyond chromosome end");
            Assert.That(() => OncologyAnalyzer.GisticMakeSampleB(layout, new[] { S("1", 1, 20, double.NaN), S("2", 1, 20, 0) }, 1),
                NUnit.Framework.Throws.ArgumentException, "NaN value");
            Assert.That(() => OncologyAnalyzer.GisticMakeSampleB(layout, new[] { S("1", 1, 20, 0), S("2", 1, 20, 0), S("3", 1, 1, 0) }, 1),
                NUnit.Framework.Throws.ArgumentException, "unknown chromosome");
            Assert.That(() => new OncologyAnalyzer.GisticMarkerLayout(new[] { new ZChr("1", 0, 0) }),
                NUnit.Framework.Throws.ArgumentException, "chromosome without markers");
            Assert.That(() => new OncologyAnalyzer.GisticMarkerLayout(new[] { new ZChr("1", 1, 1), new ZChr("1", 2, 2) }),
                NUnit.Framework.Throws.ArgumentException, "duplicate chromosome");
        });
    }

    #endregion

    #region F51 — generate_2d_hists / find_max_broad_level_by_table (Octave-locked)

    // 3-sample cohort (copy number − 2 units) on the [10 10; 8 12] layout: S1 interrupted 1p gain + chr2 0.2;
    // S2 deletion inside a 1q gain + centromere-spanning chr2 gain 0.9 over −0.3; S3 chr1 −0.4 + focal 2.4 on 2p.
    private static readonly ZSeg[][] Cohort =
    {
        new[] { S("1", 1, 3, 0.5), S("1", 4, 6, 1.5), S("1", 7, 10, 0.5), S("1", 11, 20, 0), S("2", 1, 20, 0.2) },
        new[] { S("1", 1, 10, 0), S("1", 11, 14, 0.6), S("1", 15, 16, -0.8), S("1", 17, 20, 0.6),
                S("2", 1, 5, -0.3), S("2", 6, 12, 0.9), S("2", 13, 20, -0.3) },
        new[] { S("1", 1, 20, -0.4), S("2", 1, 8, 0.25), S("2", 9, 10, 2.4), S("2", 11, 20, 0.25) },
    };

    // Octave: initial deconstruction (levels 0, perform_deconstruction step 1) → generate_2d_hists(QA,QD,[],[],.01,1),
    // then find_max_broad_level_by_table(Bt(1:i,:)) / (Bt(i+1:end,:)) for every sample, chromosome and breakpoint i.
    private const string F51Octave = """
        bg -9.4415314548696934
        cell 51 5 -2.714898452106028
        cell 16 6 -2.714898452106028
        cell 38 8 -2.714898452106028
        cell 33 11 -2.0223507320495973
        cell 22 16 -2.714898452106028
        cell 22 17 -2.714898452106028
        cell 37 18 -2.714898452106028
        cell 32 26 -2.714898452106028
        cell 20 51 -2.714898452106028
        cell 28 51 -2.714898452106028
        cell 29 51 -2.714898452106028
        fm 1 1 1 1 level 0.5 score -9.4415314548696934 n 1 rows 1
          1  1  3  0.5  1  0  0.5  0.29999999999999999  0
        fm 1 1 1 2 level 0 score -21.597961361845414 n 3 rows 2
          1  4  6  1  1  0.5  1.5  0.29999999999999999  -2.714898452106028
          1  4  10  0.5  1  0  0.5  0.69999999999999996  -9.4415314548696934
        fm 1 1 2 1 level 0.5 score -12.15642990697572 n 2 rows 2
          1  4  6  1  1  0.5  1.5  0.29999999999999999  -2.714898452106028
          1  1  6  0.5  1  0  0.5  0.59999999999999998  -9.4415314548696934
        fm 1 1 2 2 level 0 score -18.883062909739387 n 2 rows 1
          1  7  10  0.5  1  0  0.5  0.40000000000000002  -9.4415314548696934
        fm 1 1 3 1 level 0.5 score -12.15642990697572 n 2 rows 2
          1  4  6  1  1  0.5  1.5  0.29999999999999999  -2.714898452106028
          1  1  10  0.5  1  0  0.5  1  -9.4415314548696934
        fm 1 1 3 2 level 0 score -9.4415314548696934 n 1 rows 1
          1  11  20  0  1  0  0  1  0
        fm 1 1 4 1 level 0 score -21.597961361845414 n 3 rows 2
          1  4  6  1  1  0.5  1.5  0.29999999999999999  -2.714898452106028
          1  1  10  0.5  1  0  0.5  1  -9.4415314548696934
        fm 1 1 4 2 level 0 score 0 n 0 rows 0
        fm 1 2 1 1 level 0.20000000000000001 score -9.4415314548696934 n 1 rows 1
          2  21  40  0.20000000000000001  1  0  0.20000000000000001  2  0
        fm 1 2 1 2 level 0 score 0 n 0 rows 0
        fm 2 1 1 1 level 0 score -9.4415314548696934 n 1 rows 1
          1  1  10  0  2  0  0  1  0
        fm 2 1 1 2 level 0.59999999999999998 score -18.883062909739387 n 2 rows 2
          1  15  16  -1.3999999999999999  2  0.59999999999999998  -0.79999999999999993  0.20000000000000001  -9.4415314548696934
          1  11  20  0.59999999999999998  2  0  0.59999999999999998  1  -9.4415314548696934
        fm 2 1 2 1 level 0 score -18.883062909739387 n 2 rows 1
          1  11  14  0.59999999999999998  2  0  0.59999999999999998  0.40000000000000002  -9.4415314548696934
        fm 2 1 2 2 level -0.80000000000000004 score -18.883062909739387 n 2 rows 2
          1  17  20  1.3999999999999999  2  -0.80000000000000004  0.59999999999999987  0.40000000000000002  -9.4415314548696934
          1  15  20  -0.80000000000000004  2  0  -0.80000000000000004  0.60000000000000009  -9.4415314548696934
        fm 2 1 3 1 level -0.80000000000000004 score -28.32459436460908 n 3 rows 3
          1  11  14  0.59999999999999987  2  0  0.59999999999999987  0.40000000000000002  -9.4415314548696934
          1  1  14  0.80000000000000004  2  -0.80000000000000004  0  1.3999999999999999  -9.4415314548696934
          1  1  16  -0.80000000000000004  2  0  -0.80000000000000004  1.5999999999999999  -9.4415314548696934
        fm 2 1 3 2 level 0.59999999999999998 score -9.4415314548696934 n 1 rows 1
          1  17  20  0.59999999999999998  2  0  0.59999999999999998  0.40000000000000002  0
        fm 2 1 4 1 level 0.59999999999999998 score -28.32459436460908 n 3 rows 3
          1  15  16  -1.3999999999999999  2  0.59999999999999998  -0.79999999999999993  0.20000000000000001  -9.4415314548696934
          1  1  10  -0.59999999999999998  2  0.59999999999999998  -0  1  -9.4415314548696934
          1  1  20  0.59999999999999998  2  0  0.59999999999999998  2  -9.4415314548696934
        fm 2 1 4 2 level 0 score 0 n 0 rows 0
        fm 2 2 1 1 level -0.29999999999999999 score -2.714898452106028 n 1 rows 1
          2  21  25  -0.29999999999999999  2  0  -0.29999999999999999  0.625  0
        fm 2 2 1 2 level -0.29999999999999999 score -18.883062909739387 n 2 rows 2
          2  26  32  1.2  2  -0.29999999999999999  0.89999999999999991  0.70833333333333326  -9.4415314548696934
          2  26  40  -0.29999999999999999  2  0  -0.29999999999999999  1.375  -9.4415314548696934
        fm 2 2 2 1 level -0.29999999999999999 score -18.883062909739387 n 2 rows 2
          2  26  32  1.2  2  -0.29999999999999999  0.89999999999999991  0.70833333333333326  -9.4415314548696934
          2  21  32  -0.29999999999999999  2  0  -0.29999999999999999  1.3333333333333333  -9.4415314548696934
        fm 2 2 2 2 level -0.29999999999999999 score -2.714898452106028 n 1 rows 1
          2  33  40  -0.29999999999999999  2  0  -0.29999999999999999  0.66666666666666663  0
        fm 2 2 3 1 level -0.29999999999999999 score -18.883062909739387 n 2 rows 2
          2  26  32  1.2  2  -0.29999999999999999  0.89999999999999991  0.70833333333333326  -9.4415314548696934
          2  21  40  -0.29999999999999999  2  0  -0.29999999999999999  2  -9.4415314548696934
        fm 2 2 3 2 level 0 score 0 n 0 rows 0
        fm 3 1 1 1 level -0.40000000000000002 score -9.4415314548696934 n 1 rows 1
          1  1  20  -0.40000000000000002  3  0  -0.40000000000000002  2  0
        fm 3 1 1 2 level 0 score 0 n 0 rows 0
        fm 3 2 1 1 level 0.25 score -9.4415314548696934 n 1 rows 1
          2  21  28  0.25  3  0  0.25  1  0
        fm 3 2 1 2 level 0.25 score -12.15642990697572 n 2 rows 2
          2  29  30  2.1499999999999999  3  0.25  2.3999999999999999  0.16666666666666666  -2.714898452106028
          2  29  40  0.25  3  0  0.25  1  -9.4415314548696934
        fm 3 2 2 1 level 0.25 score -12.15642990697572 n 2 rows 2
          2  29  30  2.1499999999999999  3  0.25  2.3999999999999999  0.16666666666666666  -2.714898452106028
          2  21  30  0.25  3  0  0.25  1.1666666666666667  -9.4415314548696934
        fm 3 2 2 2 level 0.25 score -9.4415314548696934 n 1 rows 1
          2  31  40  0.25  3  0  0.25  0.83333333333333337  0
        fm 3 2 3 1 level 0.25 score -12.15642990697572 n 2 rows 2
          2  29  30  2.1499999999999999  3  0.25  2.3999999999999999  0.16666666666666666  -2.714898452106028
          2  21  40  0.25  3  0  0.25  2  -9.4415314548696934
        fm 3 2 3 2 level 0 score 0 n 0 rows 0
        """;

    private static (OncologyAnalyzer.GisticLengthAmplitudeTable Table, List<Row>[] B) CohortInitialTable()
    {
        var layout = MarkerLayout();
        var qa = new List<Row>();
        var qd = new List<Row>();
        var bs = new List<Row>[Cohort.Length];
        for (int j = 0; j < Cohort.Length; j++)
        {
            bs[j] = OncologyAnalyzer.GisticMakeSampleB(layout, Cohort[j], j + 1);
            var (za, zd) = OncologyAnalyzer.GisticDeconstructSample(bs[j], new double[4], new[] { 20, 40 });
            qa.AddRange(za);
            qd.AddRange(zd.Select(r => { r.Amplitude = -r.Amplitude; return r; }));
        }

        return (OncologyAnalyzer.GisticGenerate2dHistogram(qa, qd), bs);
    }

    [Test]
    public void Generate2dHistogram_CohortInitialEvents_MatchOctaveLogTable()
    {
        var (table, _) = CohortInitialTable();
        string[] lines = F51Octave.Split('\n', StringSplitOptions.TrimEntries | StringSplitOptions.RemoveEmptyEntries);
        double background = double.Parse(lines[0].Split(' ')[1], CultureInfo.InvariantCulture);
        var cells = lines.Where(l => l.StartsWith("cell", StringComparison.Ordinal))
            .Select(l => l.Split(' '))
            .ToDictionary(p => (int.Parse(p[1], CultureInfo.InvariantCulture) - 1, int.Parse(p[2], CultureInfo.InvariantCulture) - 1),
                          p => double.Parse(p[3], CultureInfo.InvariantCulture));

        Assert.That(cells, Has.Count.EqualTo(11), "Octave: 11 occupied (amplitude, length) bins");
        Assert.That(table.LogDensity.GetLength(0), Is.EqualTo(51));
        Assert.That(table.LogDensity.GetLength(1), Is.EqualTo(51));
        for (int i = 0; i < 51; i++)
        {
            for (int j = 0; j < 51; j++)
            {
                double expected = cells.TryGetValue((i, j), out double v) ? v : background;
                Assert.That(table.LogDensity[i, j], Is.EqualTo(expected), $"log_hd({i + 1},{j + 1}) must equal Octave bit-for-bit");
            }
        }

        Assert.That(table.AmplitudeEdges[30], Is.EqualTo(0.39999999999999991), "xamp(31) = Octave -2:.08:2 element");
        Assert.That(table.AmplitudeEdges[^1], Is.EqualTo(2.0));
        Assert.That(table.LengthEdges[^1], Is.EqualTo(2.0));
    }

    [Test]
    public void FindMaxBroadLevelByTable_EveryArmSplitOfCohort_MatchesOctave()
    {
        var (table, bs) = CohortInitialTable();
        string[] lines = F51Octave.Split('\n', StringSplitOptions.TrimEntries | StringSplitOptions.RemoveEmptyEntries);
        int checkedCalls = 0;
        for (int li = 0; li < lines.Length; li++)
        {
            if (!lines[li].StartsWith("fm ", StringComparison.Ordinal))
            {
                continue;
            }

            string[] h = lines[li].Split(' ', StringSplitOptions.RemoveEmptyEntries);
            int sample = int.Parse(h[1], CultureInfo.InvariantCulture);
            int chr = int.Parse(h[2], CultureInfo.InvariantCulture);
            int i = int.Parse(h[3], CultureInfo.InvariantCulture);
            bool pPart = h[4] == "1";
            double level = double.Parse(h[6], CultureInfo.InvariantCulture);
            double score = double.Parse(h[8], CultureInfo.InvariantCulture);
            int n = int.Parse(h[10], CultureInfo.InvariantCulture);
            int rows = int.Parse(h[12], CultureInfo.InvariantCulture);

            List<Row> bt = bs[sample - 1].FindAll(r => r.Chromosome == chr);
            List<Row> part = pPart ? bt.GetRange(0, i) : bt.GetRange(i, bt.Count - i);
            double fract = 0;
            foreach (Row r in part)
            {
                fract += r.Fraction;
            }

            var choice = OncologyAnalyzer.GisticFindMaxBroadLevel(part, table, fract);
            string where = $"sample {sample} chr {chr} bpt {i} {(pPart ? "p" : "q")}";
            Assert.That(choice.BroadLevel, Is.EqualTo(level), $"{where}: broad level");
            Assert.That(choice.Score, Is.EqualTo(score), $"{where}: max score");
            Assert.That(choice.LevelCount, Is.EqualTo(n), $"{where}: num_levels");
            Assert.That(choice.Events, Has.Count.EqualTo(rows), $"{where}: rows of max_Q");
            for (int r = 0; r < rows; r++)
            {
                double[] want = Parse(lines[li + 1 + r])[0];
                Row e = choice.Events[r];
                double[] got = { e.Chromosome, e.Start, e.End, e.Amplitude, e.Sample, e.StartLevel, e.EndLevel, e.Fraction, e.Score };
                Assert.That(got, Is.EqualTo(want), $"{where}: max_Q row {r + 1}");
            }

            checkedCalls++;
        }

        Assert.That(checkedCalls, Is.EqualTo(32), "every (sample, chromosome, breakpoint, arm) call of the Octave run");
    }

    #endregion

    #region F52 — perform_deconstruction / make_final_Qs / DeconstructZiggurat (Octave-locked)

    // Expected rows: "type sample chr startMarker endMarker amplitude cn_st cn_en fract score armLevel" with
    // type 0 = amp, 1 = del, 2 = aod, 3 = doa (Qs field order), sample 0-based, markers local, amplitude = Qs column 12.
    private static double[] EventColumns(OncologyAnalyzer.ZigguratEvent e) =>
        new[] { (double)(int)e.Type, e.Sample, double.Parse(e.Chromosome, CultureInfo.InvariantCulture), e.StartMarker, e.EndMarker,
                e.Amplitude, e.StartLevel, e.EndLevel, e.ArmFraction, e.Score, e.ArmLevel };

    private static void AssertEvents(OncologyAnalyzer.ZigguratDeconstruction result, string expected, string what)
    {
        double[][] want = Parse(expected);
        Assert.That(result.Events, Has.Count.EqualTo(want.Length), $"{what}: event count (Octave)");
        for (int i = 0; i < want.Length; i++)
        {
            double[] got = EventColumns(result.Events[i]);
            Assert.That(got, Is.EqualTo(want[i]), $"{what}: event {i + 1} must equal Octave bit-for-bit");
            for (int c = 5; c < got.Length; c++)
            {
                Assert.That(double.IsNegative(got[c]), Is.EqualTo(double.IsNegative(want[i][c])), $"{what}: event {i + 1} column {c + 1} sign (incl. −0)");
            }
        }
    }

    private static readonly OncologyAnalyzer.ZigguratOptions CopyNumberInput = new() { InputIsLog2 = false, Cap = null };

    // R1 — one sample, 1p 0.5 | 1.5 | 0.5 (copy number − 2 units): one broad 0.5 event over 1p (fraction 1) and one
    // focal +1.0 event (markers 4–6, fraction 0.3) — the B24 LIMITATIONS repro, resolved.
    [Test]
    public void DeconstructZiggurat_InterruptedArmGain_OneBroadOneFocal_MatchesOctave()
    {
        var sample = new[] { S("1", 1, 3, 0.5), S("1", 4, 6, 1.5), S("1", 7, 10, 0.5), S("1", 11, 20, 0), S("2", 1, 20, 0) };
        string expected = """
            0 0 1 4 6 1 0.5 1.5 0.29999999999999999 -0.9241382834507913 0
            0 0 1 1 10 0.5 0 0.5 1 -9.4415314548696951 0
            """;

        var result = OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { sample }, CopyNumberInput);
        AssertEvents(result, expected, "R1 niters 1");
        AssertEvents(OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { sample }, CopyNumberInput with { Iterations = 2 }),
            """
            0 0 1 4 6 1 0.5 1.5 0.29999999999999999 -0.9241382834507913 0
            0 0 1 1 10 0.5 0 0.5 1 -9.4415314548696951 0
            """, "R1 niters 2");

        var focal = OncologyAnalyzer.DetectFocalAmplificationEvents(result);
        Assert.That(focal, Has.Count.EqualTo(1), "GISTIC2 focal filter on the events: only the +1.0 step is focal");
        Assert.That((focal[0].StartMarker, focal[0].EndMarker, focal[0].Amplitude, focal[0].StartLevel, focal[0].EndLevel),
            Is.EqualTo((4, 6, 1.0, 0.5, 1.5)), "focal +1.0 over the 0.5 level");
        var broad = result.Events.Where(e => !e.IsFocal()).ToList();
        Assert.That(broad.Select(e => (e.StartMarker, e.EndMarker, e.Amplitude)), Is.EqualTo(new[] { (1, 10, 0.5) }),
            "one broad 0.5 event over the whole p arm (fraction 1 ≥ 0.98)");
        Assert.That(result.ForSample(0), Has.Count.EqualTo(2));

        // Contrast: the per-segment filter (unchanged) reports both flanks and the peak as focal.
        var perSegment = OncologyAnalyzer.DetectFocalAmplifications(new[]
        {
            new OncologyAnalyzer.CopyNumberArmSegment("1p", 0, 3, 10, 0.5) { MarkerCount = 3, ArmMarkerCount = 10 },
            new OncologyAnalyzer.CopyNumberArmSegment("1p", 3, 6, 10, 1.5) { MarkerCount = 3, ArmMarkerCount = 10 },
            new OncologyAnalyzer.CopyNumberArmSegment("1p", 6, 10, 10, 0.5) { MarkerCount = 4, ArmMarkerCount = 10 },
        });
        Assert.That(perSegment, Has.Count.EqualTo(3), "per-segment overload unchanged");
    }

    // R2 — deletion inside a 1q gain (0.6 | −0.8 | 0.6) and chr 2 −0.3: the −1.4 step crosses 0 and is split
    // (make_final_Qs) into a deletion 0 → −0.8 and a deletion-over-amplification 0.6 → 0; signed zeros as in Octave.
    [Test]
    public void DeconstructZiggurat_DeletionInsideGain_SplitAtZero_MatchesOctave()
    {
        var sample = new[] { S("1", 1, 10, 0), S("1", 11, 14, 0.6), S("1", 15, 16, -0.8), S("1", 17, 20, 0.6), S("2", 1, 20, -0.3) };
        AssertEvents(OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { sample }, CopyNumberInput), """
            0 0 1 1 20 0.59999999999999998 0 0.59999999999999998 2 -9.4415314548696916 0.59999999999999998
            1 0 2 1 20 0.29999999999999999 0 -0.29999999999999999 2 0 -0.29999999999999999
            1 0 1 15 16 0.79999999999999993 0 -0.79999999999999993 0.20000000000000001 -9.4415314548696916 0.59999999999999998
            3 0 1 1 10 0.59999999999999998 0.59999999999999998 -0 1 -9.4415314548696916 0.59999999999999998
            3 0 1 15 16 0.59999999999999998 0.59999999999999998 0 0.20000000000000001 -9.4415314548696916 0.59999999999999998
            """, "R2 niters 1");
        AssertEvents(OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { sample }, CopyNumberInput with { Iterations = 2 }),
            """
            0 0 1 1 20 0.59999999999999998 0 0.59999999999999998 2 -9.4415314548696951 0.59999999999999998
            1 0 2 1 20 0.29999999999999999 0 -0.29999999999999999 2 0 -0.29999999999999999
            1 0 1 15 16 0.79999999999999993 0 -0.79999999999999993 0.20000000000000001 -9.4415314548696951 0.59999999999999998
            3 0 1 1 10 0.59999999999999998 0.59999999999999998 -0 1 -9.4415314548696951 0.59999999999999998
            3 0 1 15 16 0.59999999999999998 0.59999999999999998 0 0.20000000000000001 -9.4415314548696951 0.59999999999999998
            """, "R2 niters 2");
    }

    // R4 — centromere-spanning gain on chr 2 (local markers 6..12): one event, fraction 3/8 + 4/12.
    [Test]
    public void DeconstructZiggurat_CentromereSpanningGain_MatchesOctave()
    {
        var sample = new[] { S("1", 1, 20, 0), S("2", 1, 5, 0), S("2", 6, 12, 1.0), S("2", 13, 20, 0) };
        var result = OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { sample }, CopyNumberInput);
        AssertEvents(result, """
            0 0 2 6 12 1 0 1 0.70833333333333326 -0.23109108789317978 0
            """, "R4");
        Assert.That(OncologyAnalyzer.DetectFocalAmplificationEvents(result), Has.Count.EqualTo(1),
            "0.7083 (p + q) < 0.98 ⇒ focal");
    }

    // R3 — 3-sample cohort (the F51 cohort), niters 1 and 2 (the table is rebuilt from the final events).
    [Test]
    public void DeconstructZiggurat_ThreeSampleCohort_MatchesOctave()
    {
        AssertEvents(OncologyAnalyzer.DeconstructZiggurat(Layout, Cohort, CopyNumberInput), """
            0 0 1 4 6 1 0.5 1.5 0.29999999999999999 -2.714898452106028 0
            0 0 1 1 10 0.5 0 0.5 1 -9.4415314548696934 0
            0 1 1 1 20 0.59999999999999998 0 0.59999999999999998 2 -9.4415314548696934 0.59999999999999998
            0 0 2 1 20 0.20000000000000001 0 0.20000000000000001 2 0 0.20000000000000001
            0 2 2 9 10 2.1499999999999999 0.25 2.3999999999999999 0.16666666666666666 -2.714898452106028 0.25
            0 2 2 1 20 0.25 0 0.25 2 -9.4415314548696934 0.25
            0 1 2 6 12 0.89999999999999991 0 0.89999999999999991 0.70833333333333326 -9.4415314548696934 -0.29999999999999999
            1 2 1 1 20 0.40000000000000002 0 -0.40000000000000002 2 0 -0.40000000000000002
            1 1 2 1 20 0.29999999999999999 0 -0.29999999999999999 2 -9.4415314548696934 -0.29999999999999999
            1 1 1 15 16 0.79999999999999993 0 -0.79999999999999993 0.20000000000000001 -9.4415314548696934 0.59999999999999998
            2 1 2 6 12 0.29999999999999999 -0.29999999999999999 0 0.70833333333333326 -9.4415314548696934 -0.29999999999999999
            3 1 1 1 10 0.59999999999999998 0.59999999999999998 -0 1 -9.4415314548696934 0.59999999999999998
            3 1 1 15 16 0.59999999999999998 0.59999999999999998 0 0.20000000000000001 -9.4415314548696934 0.59999999999999998
            """, "R3 niters 1");
        AssertEvents(OncologyAnalyzer.DeconstructZiggurat(Layout, Cohort, CopyNumberInput with { Iterations = 2 }),
            """
            0 0 1 4 6 1 0.5 1.5 0.29999999999999999 -2.5327766755544729 0
            0 0 1 1 10 0.5 0 0.5 1 -9.4415314548696934 0
            0 1 1 1 20 0.59999999999999998 0 0.59999999999999998 2 -9.4415314548696934 0.59999999999999998
            0 0 2 1 20 0.20000000000000001 0 0.20000000000000001 2 0 0.20000000000000001
            0 2 2 9 10 2.1499999999999999 0.25 2.3999999999999999 0.16666666666666666 -2.5327766755544729 0.25
            0 2 2 1 20 0.25 0 0.25 2 -9.4415314548696934 0.25
            0 1 2 6 12 0.89999999999999991 0 0.89999999999999991 0.70833333333333326 -9.4415314548696934 -0.29999999999999999
            1 2 1 1 20 0.40000000000000002 0 -0.40000000000000002 2 0 -0.40000000000000002
            1 1 2 1 20 0.29999999999999999 0 -0.29999999999999999 2 -9.4415314548696934 -0.29999999999999999
            1 1 1 15 16 0.79999999999999993 0 -0.79999999999999993 0.20000000000000001 -9.4415314548696934 0.59999999999999998
            2 1 2 6 12 0.29999999999999999 -0.29999999999999999 0 0.70833333333333326 -9.4415314548696934 -0.29999999999999999
            3 1 1 1 10 0.59999999999999998 0.59999999999999998 -0 1 -9.4415314548696934 0.59999999999999998
            3 1 1 15 16 0.59999999999999998 0.59999999999999998 0 0.20000000000000001 -9.4415314548696934 0.59999999999999998
            """, "R3 niters 2");
    }

    // R5 — full perform_ziggurat_deconstruction: log2 input, default cap 1.5 (1.9 and 1.6 capped, −∞ capped to −1.5),
    // 2^(x+1) − 2 conversion; Qs.amp/del/aod/doa with column 12 = amplitude.
    [Test]
    public void DeconstructZiggurat_Log2InputWithDefaultCap_MatchesOctavePerformZigguratDeconstruction()
    {
        var samples = new[]
        {
            new[] { S("1", 1, 4, 0.3), S("1", 5, 6, 1.9), S("1", 7, 20, 0.3), S("2", 1, 20, 0) },
            new[] { S("1", 1, 20, -0.2), S("2", 1, 3, 0), S("2", 4, 15, double.NegativeInfinity), S("2", 16, 20, 0) },
            new[] { S("1", 1, 12, 0.584962500721156), S("1", 13, 20, -1), S("2", 1, 9, 0.1), S("2", 10, 11, 1.6), S("2", 12, 20, 0.1) },
        };
        AssertEvents(OncologyAnalyzer.DeconstructZiggurat(Layout, samples), """
            0 0 1 5 6 3.194565422802548 0.46228882668983262 3.6568542494923806 0.20000000000000001 -2.309832944402781 0.46228882668983262
            0 0 1 1 20 0.46228882668983262 0 0.46228882668983262 2 -9.4415314548696916 0.46228882668983262
            0 2 2 10 11 3.5133073244197943 0.14354692507258626 3.6568542494923806 0.16666666666666666 -2.309832944402781 0.14354692507258626
            0 2 2 1 20 0.14354692507258626 0 0.14354692507258626 2 -9.4415314548696916 0.14354692507258626
            0 2 1 1 12 1 0 1 1.2 -9.4415314548696916 -1
            1 1 1 1 20 0.25889887340775175 0 -0.25889887340775175 2 0 -0.25889887340775175
            1 2 1 1 20 1 0 -1 2 -9.4415314548696916 -1
            1 1 2 4 15 1.2928932188134525 0 -1.2928932188134525 1.2083333333333335 -2.309832944402781 0
            2 2 1 1 12 1 -1 0 1.2 -9.4415314548696916 -1
            """, "R5");
    }

    // The GISTIC2 focal filter on events (reconstruct_genomes.m): amp + aod, fraction < cutoff, amplitude ≥ t_amp.
    [Test]
    public void DetectFocalAmplificationEvents_AppliesReconstructGenomesFilter()
    {
        var result = OncologyAnalyzer.DeconstructZiggurat(Layout, Cohort, CopyNumberInput);
        var focal = OncologyAnalyzer.DetectFocalAmplificationEvents(result);
        Assert.That(focal.Select(e => (e.Type, e.Sample, e.Chromosome, e.StartMarker, e.EndMarker)), Is.EqualTo(new[]
        {
            (OncologyAnalyzer.ZigguratEventType.Amplification, 0, "1", 4, 6),
            (OncologyAnalyzer.ZigguratEventType.Amplification, 2, "2", 9, 10),
            (OncologyAnalyzer.ZigguratEventType.Amplification, 1, "2", 6, 12),
            (OncologyAnalyzer.ZigguratEventType.AmplificationOverDeletion, 1, "2", 6, 12),
        }), "focal amp + aod events of the cohort (fraction < 0.98, amplitude ≥ 0.1)");

        // ≥ t_amp (reconstruct_genomes.m): an event of amplitude exactly t_amp passes.
        var exact = OncologyAnalyzer.DetectFocalAmplificationEvents(result, new OncologyAnalyzer.FocalAmplificationThresholds(1.0, 0.98));
        Assert.That(exact.Select(e => e.Amplitude), Is.EqualTo(new[] { 1.0, 2.1499999999999999 }), "amplitude 1.0 ≥ t_amp 1.0 is kept");
        Assert.That(() => OncologyAnalyzer.DetectFocalAmplificationEvents(result, new OncologyAnalyzer.FocalAmplificationThresholds(double.NaN, 0.98)),
            NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
        Assert.That(() => OncologyAnalyzer.DetectFocalAmplificationEvents(null!), NUnit.Framework.Throws.ArgumentNullException);
    }

    // All-neutral cohort: every chromosome is one zero segment ⇒ no events (Octave perform_deconstruction: empty Q arrays).
    [Test]
    public void DeconstructZiggurat_AllNeutral_NoEvents()
    {
        var neutral = new[] { S("1", 1, 20, 0), S("2", 1, 4, 0), S("2", 5, 20, 0) };
        Assert.That(OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { neutral, neutral }, CopyNumberInput).Events, Is.Empty);
    }

    [Test]
    public void DeconstructZiggurat_InvalidArguments_Throw()
    {
        var ok = new[] { S("1", 1, 20, 0.5), S("2", 1, 20, 0) };
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(null!, new[] { ok }), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(Layout, null!), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { ok, null! }), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(Layout, Array.Empty<ZSeg[]>()), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { ok }, new() { Cap = 0 }),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { ok }, new() { Cap = double.NaN }),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { ok }, new() { Iterations = 0 }),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(Layout, new[] { new[] { S("1", 1, 20, double.PositiveInfinity), S("2", 1, 20, 0) } },
                new() { Cap = null }), NUnit.Framework.Throws.ArgumentException, "+∞ log2 without a cap is not finite after conversion");
        });
    }

    #endregion

    #region F54 — remove_noisy_samples (Octave-locked)

    // FIN-B24 WP26 / F54. Log2 input, default cap 1.5. Octave: SegArray branch of remove_noisy_samples.m
    // (cap_vals → getbpt_counts, keepers = count <= max) gives counts [2 5 2 6]: S0 0.4 runs across the chr1/chr2
    // boundary (one segment), S2 2.0 | 2.2 both capped to 1.5 (annealed), S3 two equal 0.5 segments joined. The
    // uncompressed branch would give [1 4 1 5] (diff ~= 0) — not followed (GISTIC2 default use_segarray = 1).
    // Events: perform_ziggurat_deconstruction on the kept columns, sample = original index (Qs.sdesc).
    private static readonly ZSeg[][] NoisyCohort =
    {
        new[] { S("1", 1, 20, 0.4), S("2", 1, 8, 0.4), S("2", 9, 20, -0.5) },
        new[] { S("1", 1, 5, 0.3), S("1", 6, 8, 0.9), S("1", 9, 20, 0.3), S("2", 1, 10, -0.4), S("2", 11, 20, 0.2) },
        new[] { S("1", 1, 4, 2.0), S("1", 5, 9, 2.2), S("1", 10, 20, 0.5), S("2", 1, 20, 0.5) },
        new[] { S("1", 1, 6, 0.5), S("1", 7, 12, 0.5), S("1", 13, 20, 1.0), S("2", 1, 5, 0.2), S("2", 6, 9, 0.6), S("2", 10, 15, 0.2), S("2", 16, 20, -0.6) },
    };

    private const string NoisyCohortAllKept = """
            0 0 1 1 20 0.6390158215457884 0 0.6390158215457884 2 0 0.6390158215457884
            0 1 1 6 8 1.269843139457397 0.46228882668983262 1.7321319661472296 0.29999999999999999 -2.9377424078719847 0.46228882668983262
            0 1 1 1 20 0.46228882668983262 0 0.46228882668983262 2 -9.4415314548696934 0.46228882668983262
            0 2 1 1 9 2.8284271247461903 0.82842712474619029 3.6568542494923806 0.90000000000000002 -2.9377424078719847 0.82842712474619029
            0 2 1 1 20 0.82842712474619029 0 0.82842712474619029 2 -9.4415314548696934 0.82842712474619029
            0 3 1 13 20 1.1715728752538097 0.82842712474619029 2 0.80000000000000004 -9.4415314548696934 0.82842712474619029
            0 3 1 1 20 0.82842712474619029 0 0.82842712474619029 2 -9.4415314548696934 0.82842712474619029
            0 1 2 11 20 0.29739670999406975 0 0.29739670999406975 0.83333333333333337 0 0.29739670999406975
            0 2 2 1 20 0.82842712474619029 0 0.82842712474619029 2 0 0.82842712474619029
            0 3 2 6 9 0.73403642302672667 0.29739670999406975 1.0314331330207964 0.45833333333333331 -2.9377424078719847 0.29739670999406975
            0 3 2 1 15 0.29739670999406975 0 0.29739670999406975 1.5833333333333333 -2.9377424078719847 0.29739670999406975
            0 0 2 1 8 0.6390158215457884 0 0.6390158215457884 1 -9.4415314548696934 -0.58578643762690485
            1 0 2 1 20 0.58578643762690485 0 -0.58578643762690485 2 -9.4415314548696934 -0.58578643762690485
            1 1 2 1 10 0.48428343348960201 0 -0.48428343348960201 1.1666666666666667 0 -0.48428343348960201
            1 3 2 16 20 0.6804920892271058 0 -0.6804920892271058 0.41666666666666669 0 -0.6804920892271058
            2 0 2 1 8 0.58578643762690485 -0.58578643762690485 0 1 -9.4415314548696934 -0.58578643762690485
            """;

    [Test]
    public void DeconstructZiggurat_NoiseFilter_SegmentCountsAndDefault2500_MatchOctave()
    {
        Assert.That(OncologyAnalyzer.ZigguratOptions.Default.MaxSegmentsPerSample, Is.EqualTo(2500), "gistic2_param_defaults.m l.136");
        var result = OncologyAnalyzer.DeconstructZiggurat(Layout, NoisyCohort);
        Assert.That(result.SegmentCounts, Is.EqualTo(new[] { 2, 5, 2, 6 }), "Octave getbpt_counts(cap_vals(SegArray(D.dat)))");
        Assert.That(result.RemovedSamples, Is.Empty);
        AssertEvents(result, NoisyCohortAllKept, "max 2500");
        AssertEvents(OncologyAnalyzer.DeconstructZiggurat(Layout, NoisyCohort, new() { MaxSegmentsPerSample = 6 }), NoisyCohortAllKept, "max 6 (count 6 kept: <=)");

        var unfiltered = OncologyAnalyzer.DeconstructZiggurat(Layout, NoisyCohort, new() { MaxSegmentsPerSample = null });
        Assert.That(unfiltered.SegmentCounts, Is.Empty, "filter disabled");
        AssertEvents(unfiltered, NoisyCohortAllKept, "no filter");
    }

    [Test]
    public void DeconstructZiggurat_NoiseFilter_Max5_RemovesSample3_MatchesOctave()
    {
        var result = OncologyAnalyzer.DeconstructZiggurat(Layout, NoisyCohort, new() { MaxSegmentsPerSample = 5 });
        Assert.That(result.RemovedSamples, Is.EqualTo(new[] { 3 }), "count 6 > 5 removed, count 5 kept");
        AssertEvents(result, """
            0 0 1 1 20 0.6390158215457884 0 0.6390158215457884 2 0 0.6390158215457884
            0 1 1 6 8 1.269843139457397 0.46228882668983262 1.7321319661472296 0.29999999999999999 -2.5327766755544725 0.46228882668983262
            0 1 1 1 20 0.46228882668983262 0 0.46228882668983262 2 -9.4415314548696934 0.46228882668983262
            0 2 1 1 9 2.8284271247461903 0.82842712474619029 3.6568542494923806 0.90000000000000002 -2.5327766755544725 0.82842712474619029
            0 2 1 1 20 0.82842712474619029 0 0.82842712474619029 2 -9.4415314548696934 0.82842712474619029
            0 1 2 11 20 0.29739670999406975 0 0.29739670999406975 0.83333333333333337 0 0.29739670999406975
            0 2 2 1 20 0.82842712474619029 0 0.82842712474619029 2 0 0.82842712474619029
            0 0 2 1 8 0.6390158215457884 0 0.6390158215457884 1 -9.4415314548696934 -0.58578643762690485
            1 0 2 1 20 0.58578643762690485 0 -0.58578643762690485 2 -9.4415314548696934 -0.58578643762690485
            1 1 2 1 10 0.48428343348960201 0 -0.48428343348960201 1.1666666666666667 0 -0.48428343348960201
            2 0 2 1 8 0.58578643762690485 -0.58578643762690485 0 1 -9.4415314548696934 -0.58578643762690485
            """, "max 5");
    }

    [Test]
    public void DeconstructZiggurat_NoiseFilter_Max2_KeepsSamples0And2_OriginalIndices_MatchesOctave()
    {
        var result = OncologyAnalyzer.DeconstructZiggurat(Layout, NoisyCohort, new() { MaxSegmentsPerSample = 2 });
        Assert.That(result.RemovedSamples, Is.EqualTo(new[] { 1, 3 }));
        Assert.That(result.SegmentCounts, Is.EqualTo(new[] { 2, 5, 2, 6 }));
        AssertEvents(result, """
            0 0 1 1 20 0.6390158215457884 0 0.6390158215457884 2 0 0.6390158215457884
            0 2 1 1 9 2.8284271247461903 0.82842712474619029 3.6568542494923806 0.90000000000000002 -2.0223507320495968 0.82842712474619029
            0 2 1 1 20 0.82842712474619029 0 0.82842712474619029 2 -9.4415314548696916 0.82842712474619029
            0 2 2 1 20 0.82842712474619029 0 0.82842712474619029 2 0 0.82842712474619029
            0 0 2 1 8 0.6390158215457884 0 0.6390158215457884 1 -9.4415314548696916 -0.58578643762690485
            1 0 2 1 20 0.58578643762690485 0 -0.58578643762690485 2 -9.4415314548696916 -0.58578643762690485
            2 0 2 1 8 0.58578643762690485 -0.58578643762690485 0 1 -9.4415314548696916 -0.58578643762690485
            """, "max 2");
    }

    [Test]
    public void DeconstructZiggurat_NoiseFilter_AllRemoved_OrNegative_Throws()
    {
        Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(Layout, NoisyCohort, new() { MaxSegmentsPerSample = 1 }),
            NUnit.Framework.Throws.ArgumentException, "GISTIC2 snp:perform_ziggurat_deconstruction:all_data_removed");
        Assert.That(() => OncologyAnalyzer.DeconstructZiggurat(Layout, NoisyCohort, new() { MaxSegmentsPerSample = -1 }),
            NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
    }

    #endregion
}
