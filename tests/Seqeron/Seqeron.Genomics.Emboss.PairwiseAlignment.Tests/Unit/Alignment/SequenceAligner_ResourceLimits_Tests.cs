// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-MULTI-001 (audit A3-2, A3-3; B09 F131, F132): the resource limits of the program ports (one <see cref="ArgumentException"/> style, 64-bit arithmetic, raised
/// BEFORE anything large is allocated: every limit test asserts the exception AND that the call allocated less than 64 MB) and the NaN / Infinity / overflowing penalty contract of the EMBOSS ports.
/// The programs' own limit is "none": installed needle / water die with "Sequences too big. Try 'stretcher'" / "Try 'matcher' or 'supermatcher'" when the allocation fails (measured under a 4 GB
/// address-space limit with 30 000 x 30 000 residues), installed stretcher wraps its 32-bit score and warns "stretcher CheckScore failed" (1e9 penalties: score 1 949 672 837).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ResourceLimits_Tests
{

    // ------------------------------------------------------------------ A3-2 limits

    [Test]
    public void Needle_Above2Pow27Cells_IsRefusedWithoutAllocating()
    {
        string a = Dna(12_000, 1), b = Dna(12_000, 2);              // 144 000 000 cells > 134 217 728
        var ex = Refused<ArgumentException>(() => EmbossPairwisePrograms.NeedleAlign(a, b));
        Assert.That(ex.Message, Does.StartWith("EMBOSS needle keeps").And.Contain("144,000,000 cells").And.Contain("SequenceAligner.NeedleMaxCells").And.Contain("StretcherAlign"));
        Assert.That(SequenceAligner.NeedleMaxCells, Is.EqualTo(134_217_728L));
    }

    /// <summary>A4-15: the guard messages format numbers with the invariant culture (thousands separator "," whatever the current culture: here "." as in de-DE).</summary>
    [Test]
    public void LimitMessages_AreCultureInvariant()
    {
        var odd = (System.Globalization.CultureInfo)System.Globalization.CultureInfo.InvariantCulture.Clone();
        odd.NumberFormat.NumberGroupSeparator = ".";
        odd.NumberFormat.NumberDecimalSeparator = ",";
        var saved = System.Globalization.CultureInfo.CurrentCulture;
        string needle, merger;
        try
        {
            System.Globalization.CultureInfo.CurrentCulture = odd;
            Assert.That(144000000L.ToString("N0"), Is.EqualTo("144.000.000"), "the test culture really differs");
            string a = Dna(12_000, 1), b = Dna(12_000, 2);
            needle = Refused<ArgumentException>(() => EmbossPairwisePrograms.NeedleAlign(a, b)).Message;
            merger = Refused<ArgumentException>(() => EmbossPairwisePrograms.Merger(a, b)).Message;
        }
        finally { System.Globalization.CultureInfo.CurrentCulture = saved; }
        Assert.That(needle, Does.Contain("144,000,000 cells").And.Not.Contain("144.000.000"));
        Assert.That(merger, Does.Contain("12,000 x 12,000 residues").And.Not.Contain("12.000"));
    }

    [Test]
    public void Water_Above2Pow28Cells_IsRefusedWithoutAllocating()
    {
        string a = Dna(17_000, 3), b = Dna(17_000, 4);              // 289 000 000 cells > 268 435 456
        var ex = Refused<ArgumentException>(() => EmbossPairwisePrograms.WaterAlign(a, b));
        Assert.That(ex.Message, Does.StartWith("EMBOSS water keeps").And.Contain("SequenceAligner.WaterMaxCells").And.Contain("Supermatcher"));
    }

    [Test]
    public void Supermatcher_HugeWidth_IsRefusedWithoutAllocating()
    {
        var seq = new[] { new EmbossSequence("a", Dna(200, 5)) };
        var ex = Refused<ArgumentException>(() => EmbossPairwisePrograms.Supermatcher(seq, seq, width: 1_000_000_000));
        Assert.That(ex.Message, Does.StartWith("EMBOSS supermatcher keeps").And.Contain("SequenceAligner.SupermatcherMaxPathCells"));
    }

    // ------------------------------------------------------------------ A3-3 penalties

    [Test]
    public void Needle_Water_NaNAndInfinitePenalties_AreRejected()
    {
        double[] bad = [double.NaN, double.PositiveInfinity, double.NegativeInfinity];
        Assert.Multiple(() =>
        {
            foreach (double v in bad)
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAlign("ACGT", "ACT", null, v, 0.5), $"needle gapOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAlign("ACGT", "ACT", null, 10, v), $"needle gapExtend {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAlign("ACGT", "ACT", null, 10, 0.5, true, v, 0.5), $"needle endOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAlign("ACGT", "ACT", null, 10, 0.5, true, 10, v), $"needle endExtend {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.WaterAlign("ACGT", "ACT", null, v, 0.5), $"water gapOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.WaterAlign("ACGT", "ACT", null, 10, v), $"water gapExtend {v}");
                var one = new[] { new EmbossSequence("a", "ACGTACGT") };
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAll(one, one, null, v), $"needleall gapOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAll(one, one, null, 10, 0.5, false, 10, 0.5, v), $"needleall minscore {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.Supermatcher(one, one, null, v), $"supermatcher gapOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.Supermatcher(one, one, null, 10, v), $"supermatcher gapExtend {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.Supermatcher(one, one, null, 10, 0.5, 6, 16, v), $"supermatcher minscore {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, v, 0.5), $"MatrixScoring.Emboss open {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, v), $"MatrixScoring.Emboss extend {v}");
            }

            // a NaN used to pass `x < min || x > max`: needle("ACGT", "ACT", gapOpen NaN) returned score 0, water(gapExtend NaN) returned 10
            Assert.That(EmbossPairwisePrograms.NeedleAlign("ACGT", "ACT").Score, Is.EqualTo(EmbossPairwisePrograms.NeedleAlign("ACGT", "ACT", null, 10, 0.5).Score));
        });
    }

    [Test]
    public void Stretcher_OverflowingPenalties_AreRefusedLikeTheScalarContract()
    {
        string a = "ACGT" + new string('A', 50), b = "ACT" + new string('C', 50);
        // installed stretcher -gapopen 1000000000 -gapextend 1000000000: "stretcher CheckScore failed", Score: 1949672837 (a 32-bit wrap; int.MaxValue gives 1687091105 or 2147483643)
        var ex = Assert.Throws<ArgumentException>(() => EmbossPairwisePrograms.StretcherAlign(a, b, null, 1_000_000_000, 1_000_000_000));
        Assert.That(ex!.Message, Does.Contain("Score range too large").And.Contain("CheckScore"));
        Assert.Throws<ArgumentException>(() => EmbossPairwisePrograms.StretcherAlign(a, b, null, int.MaxValue, int.MaxValue));
        Assert.Throws<ArgumentException>(() => EmbossPairwisePrograms.StretcherAlign(a, b, null, 100_000, 1_000_000_000));
        // large but exact penalties are fine: the score is the exact integer (no wrap)
        Assert.That(EmbossPairwisePrograms.StretcherAlign(a, b, null, 1_000_000, 1_000_000).Score, Is.LessThan(0));
        Assert.That(EmbossPairwisePrograms.StretcherAlign("ACGT", "ACGT").Score, Is.EqualTo(20), "installed stretcher 6.6.0 -asequence ACGT -bsequence ACGT: # Score: 20");
    }

    [Test]
    public void TopNSim_MaxAlignmentsIntMaxValue_MeansEveryAlignment_LikeTheWatermanEggertOrder()
    {
        // the vertex list used to be allocated up front (new Vertex[int.MaxValue]: OutOfMemoryException) while the Waterman-Eggert order accepted int.MaxValue (23 alignments)
        string a = Dna(60, 21), b = Dna(60, 22);
        var scoring = new ScoringMatrix(5, -4, -10, -4);
        var sim = EmbossPairwisePrograms.LocalAlignAffineTopN(a, b, scoring, int.MaxValue, 0, TopNTieOrder.Sim);
        var simBig = EmbossPairwisePrograms.LocalAlignAffineTopN(a, b, scoring, 100_000, 0, TopNTieOrder.Sim);
        Assert.That(sim, Has.Count.GreaterThan(1));
        Assert.That(sim.Select(x => (x.Score, x.AlignedSequence1, x.AlignedSequence2)), Is.EqualTo(simBig.Select(x => (x.Score, x.AlignedSequence1, x.AlignedSequence2))));
        // A5-4: int.MaxValue is exactly "every alignment": the same rows as 100 000 in every order and overload
        static IEnumerable<(int, string, string)> Rows(IEnumerable<AlignmentResult> r) => r.Select(x => (x.Score, x.AlignedSequence1, x.AlignedSequence2));
        Assert.That(Rows(EmbossPairwisePrograms.LocalAlignAffineTopN(a, b, scoring, int.MaxValue, 0, TopNTieOrder.WatermanEggert)), Is.EqualTo(Rows(EmbossPairwisePrograms.LocalAlignAffineTopN(a, b, scoring, 100_000, 0, TopNTieOrder.WatermanEggert))).And.Not.Empty);
        Assert.That(Rows(EmbossPairwisePrograms.LocalAlignTopN(a, b, scoring, int.MaxValue, 0, TopNTieOrder.Sim)), Is.EqualTo(Rows(EmbossPairwisePrograms.LocalAlignTopN(a, b, scoring, 100_000, 0, TopNTieOrder.Sim))).And.Not.Empty);
        Assert.That(Rows(EmbossPairwisePrograms.LocalAlignTopN(a, b, scoring, int.MaxValue, 0, TopNTieOrder.WatermanEggert)), Is.EqualTo(Rows(EmbossPairwisePrograms.LocalAlignTopN(a, b, scoring, 100_000, 0, TopNTieOrder.WatermanEggert))).And.Not.Empty);
    }

    [Test]
    public void TopNSim_VertexListIsBounded_ByTheDocumentedLimit()
    {
        string a = Dna(200, 23), b = Dna(200, 24);
        var ex = Refused<ArgumentException>(() => EmbossPairwisePrograms.LocalTopNSimCore(a, b, new ScoringMatrix(5, -4, -10, -4), affine: true, int.MaxValue, 0, maxVertices: 10));
        Assert.That(ex.Message, Does.StartWith("EMBOSS matcher (Huang-Miller SIM) keeps a list of candidate start points").And.Contain("SequenceAligner.TopNSimMaxVertices").And.Contain("TopNTieOrder.WatermanEggert"));
        Assert.That(ex.ParamName, Is.EqualTo("maxAlignments"));
        Assert.That(SequenceAligner.TopNSimMaxVertices, Is.EqualTo(4_194_304L));
    }

    [Test]
    public void DotMatcher_SegmentList_IsBounded_WhenTheThresholdLetsEveryMatchBeASegment()
    {
        // window 1, threshold 0: every matching cell of a 6 000 x 6 000 plot (EDNAFULL +5 / -4) starts a segment: 0.1875 x 36 million = 6.75 million segments (> 2^22)
        string a = Dna(6_000, 41), b = Dna(6_000, 42);
        long before = GC.GetAllocatedBytesForCurrentThread();
        var ex = Assert.Throws<ArgumentException>(() => DotMatcher.Compute(a, b, windowSize: 1, threshold: 0))!;
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(512L << 20), "the list stops growing at the limit (64 MiB of segments)");
        Assert.That(ex.Message, Does.StartWith("EMBOSS dotmatcher keeps the segments of the plot").And.Contain("SequenceAligner.DotMatcherMaxSegments"));
        Assert.That(SequenceAligner.DotMatcherMaxSegments, Is.EqualTo(4_194_304L));
        // A5-4: the plot of an ordinary pair: 397 segments, the first one (316,318)-(326,328) as the first "Line" of installed dotmatcher 6.6.0 -windowsize 10 -threshold 23 (the full plot comparison is DotMatcher_Tests, 600 / 600)
        var ordinary = DotMatcher.Compute(a[..500], b[..500], windowSize: 10, threshold: 23);
        Assert.That(ordinary, Has.Count.EqualTo(397));
        Assert.That((ordinary[0].Start1, ordinary[0].Start2, ordinary[0].End1, ordinary[0].End2), Is.EqualTo((318, 316, 328, 326)));
    }

    private static string Dna(int n, int seed)
    {
        var r = new Random(seed);
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = "ACGT"[r.Next(4)];
        return new string(c);
    }

    /// <summary>The exception of <paramref name="call"/>, asserting it allocated (almost) nothing.</summary>
    private static T Refused<T>(Action call) where T : Exception
    {
        long before = GC.GetAllocatedBytesForCurrentThread();
        var ex = Assert.Throws<T>(() => call());
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(allocated, Is.LessThan(MaxAllocatedBytes), "the limit is checked before the matrices are allocated");
        return ex!;
    }

    private const long MaxAllocatedBytes = 64L << 20;
}
