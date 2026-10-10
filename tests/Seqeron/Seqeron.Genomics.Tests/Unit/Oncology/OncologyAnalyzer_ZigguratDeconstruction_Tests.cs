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
}
