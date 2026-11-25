// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using Seqeron.Genomics.TCoffee;
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
    private const long MaxAllocatedBytes = 64L << 20;

    /// <summary>The exception of <paramref name="call"/>, asserting it allocated (almost) nothing.</summary>
    private static T Refused<T>(Action call) where T : Exception
    {
        long before = GC.GetAllocatedBytesForCurrentThread();
        var ex = Assert.Throws<T>(() => call());
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(allocated, Is.LessThan(MaxAllocatedBytes), "the limit is checked before the matrices are allocated");
        return ex!;
    }

    private static string Dna(int n, int seed)
    {
        var r = new Random(seed);
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = "ACGT"[r.Next(4)];
        return new string(c);
    }

    [Test]
    public void ProbCons_LongPair_IsRefusedWithProbConsWording_NotTCoffees()
    {
        string a = string.Concat(Enumerable.Repeat("ACDEFGHIKL", 740)), b = string.Concat(Enumerable.Repeat("LKIHGFEDCA", 740));      // 7 400 x 7 400: 5 x 7 401^2 floats > 2^28
        var ex = Refused<ArgumentException>(() => Seqeron.Genomics.TCoffee.ProbConsAligner.MultipleAlignProbCons([a, b]));
        Assert.That(ex.Message, Does.StartWith("The ProbCons pair-HMM keeps").And.Contain("ProbCons 1.12 itself").And.Contain("TCoffeeAligner.TCoffeeProbaMaxMatrixFloats"));
        Assert.That(ex.Message, Does.Not.Contain("T-Coffee").And.Not.Contain("t_coffee").And.Not.Contain("proba_pair"));
        var tc = Assert.Throws<ArgumentException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignTCoffeeProgram([a, b]));
        Assert.That(tc!.Message, Does.StartWith("The T-Coffee proba_pair pair-HMM keeps").And.Contain("slow_pair"), "the T-Coffee caller keeps its own wording");
    }

    [Test]
    public void Mafft_Above2Pow28Cells_IsRefusedWithoutAllocating()
    {
        string a = Dna(17_000, 6), b = Dna(17_000, 7);
        var ex = Refused<ArgumentException>(() => SequenceAligner.MultipleAlignMafft([a, b]));
        Assert.That(ex.Message, Does.StartWith("MAFFT keeps").And.Contain("SequenceAligner.MafftMaxPairCells"));
    }

    [Test]
    public void Ksw2_LargeTraceback_IsRefusedWithoutAllocating()
    {
        var q = new byte[40_000];
        var t = new byte[40_000];
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, 4, 1);
        var scalar = Refused<ArgumentException>(() => Ksw2Aligner.ExtendZ(q, t, 5, mat, 4, 2, -1, -1));
        Assert.That(scalar.Message, Does.StartWith("ksw2 ExtendZ keeps a traceback matrix").And.Contain("Ksw2Aligner.MaxTracebackBytes").And.Contain("ScoreOnly"));
        Refused<ArgumentException>(() => Ksw2Aligner.ExtendD(q, t, 5, mat, 4, 2, 24, 1, -1, -1));
        Refused<ArgumentException>(() => Ksw2Aligner.ExtendZ2(q, t, 5, mat, 4, 2, -1, -1, 0));
        Refused<ArgumentException>(() => Ksw2Aligner.ExtendD2(q, t, 5, mat, 4, 2, 24, 1, -1, -1, 0));
        // score-only needs no traceback: the same sizes are not refused (a narrow band keeps the run short)
        Assert.DoesNotThrow(() => Ksw2Aligner.ExtendZ(q, t, 5, mat, 4, 2, 10, -1, Ksw2Option.ScoreOnly));
    }

    [Test]
    public void Ksw2_BandWiderThanTheMatrix_IsNoBand_AndIntMaxValueDoesNotOverflow()
    {
        var q = Ksw2Aligner.EncodeNucleotides("ACGTACGTTAGCTAGCATCG");
        var t = Ksw2Aligner.EncodeNucleotides("ACGTACGATAGCTAGCTCG");
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, 4, 1);
        var full = Ksw2Aligner.ExtendZ(q, t, 5, mat, 4, 2, -1, -1);
        Assert.Multiple(() =>
        {
            var wide = Ksw2Aligner.ExtendZ(q, t, 5, mat, 4, 2, int.MaxValue, -1);
            Assert.That(wide.Score, Is.EqualTo(full.Score));
            Assert.That(wide.Cigar, Is.EqualTo(full.Cigar));
            Assert.That(Ksw2Aligner.ExtendD(q, t, 5, mat, 4, 2, 24, 1, int.MaxValue, -1).Score, Is.EqualTo(Ksw2Aligner.ExtendD(q, t, 5, mat, 4, 2, 24, 1, -1, -1).Score));
            Assert.That(Ksw2Aligner.ExtendZ2(q, t, 5, mat, 4, 2, int.MaxValue, -1, 0).Score, Is.EqualTo(Ksw2Aligner.ExtendZ2(q, t, 5, mat, 4, 2, -1, -1, 0).Score));
        });
    }

    [Test]
    public void AlignmentFormats_HeaderCounts_AreNotTrusted()
    {
        long before = GC.GetAllocatedBytesForCurrentThread();
        var phylip = Assert.Throws<FormatException>(() => AlignmentFormats.Parse("999999999 999999999\nA ACGT\n", AlignmentFileFormat.Phylip));
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(allocated, Is.LessThan(MaxAllocatedBytes));
        Assert.That(phylip!.Message, Does.Contain("999999999 sequences").And.Contain("only 1 lines"));
        Assert.Throws<FormatException>(() => AlignmentFormats.Parse("-3 4\nA ACGT\n", AlignmentFileFormat.Phylip));
        Assert.Throws<FormatException>(() => AlignmentFormats.Parse("#NEXUS\nbegin data;\ndimensions ntax=99999999999 nchar=4;\nformat datatype=dna;\nmatrix\nA ACGT\n;\nend;\n", AlignmentFileFormat.Nexus));
        Assert.Throws<FormatException>(() => AlignmentFormats.Parse("# STOCKHOLM 1.0\n#=GF RN [99999999999999999999999]\nA ACGT\n//\n", AlignmentFileFormat.Stockholm));
        Assert.Throws<FormatException>(() => AlignmentFormats.Parse("# STOCKHOLM 1.0\n#=GF SQ 99999999999999999999999\nA ACGT\n//\n", AlignmentFileFormat.Stockholm));
        Assert.Throws<FormatException>(() => AlignmentFormats.Parse("# STOCKHOLM 1.0\n#=GF SQ 99999999999\nA ACGT\n//\n", AlignmentFileFormat.Stockholm), "Biopython: Inconsistent number of sequences in alignment");
        // a 64-bit reference number is read (Python ints are unbounded; Biopython reads 99999999999)
        Assert.That(AlignmentFormats.Parse("# STOCKHOLM 1.0\n#=GF RN [99999999999]\nA ACGT\n//\n", AlignmentFileFormat.Stockholm)[0].Records.Count, Is.EqualTo(1));
        // a consistent header is read as before
        Assert.That(AlignmentFormats.Parse("2 4\nA         ACGT\nB         ACGA\n", AlignmentFileFormat.Phylip)[0].Records.Count, Is.EqualTo(2));
    }

    [Test]
    public void DualAffine_IntMinValue_IsAnArgumentException_NotOverflow()
    {
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignDualAffine("ACGT", "ACGT", new DualAffineScoring(2, -4, int.MinValue, -2, -24, -1)));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignDualAffine("ACGT", "ACGT", new DualAffineScoring(int.MinValue, -4, -4, -2, -24, -1)));
    }

    [Test]
    public void Ksw2_GapPenalties_FollowTheHeaderOfKsw2()
    {
        var q = Ksw2Aligner.EncodeNucleotides("ACGTACGT");
        var t = Ksw2Aligner.EncodeNucleotides("ACGACGT");
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, 4, 1);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.ExtendZ(q, t, 5, mat, -1, 2, -1, -1), "q < 0");
            Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.ExtendZ(q, t, 5, mat, 4, 0, -1, -1), "e = 0");
            Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.ExtendZ(q, t, 5, mat, 200, 2, -1, -1), "q beyond int8_t");
            Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.ExtendD(q, t, 5, mat, 4, 2, -24, 1, -1, -1), "q2 < 0");
            Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.ExtendZ2(q, t, 5, mat, 4, 0, -1, -1, 0), "e = 0 (SIMD)");
            Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.ExtendD2(q, t, 5, mat, 4, 2, 24, -1, -1, -1, 0), "e2 < 0 (SIMD)");
            Assert.DoesNotThrow(() => Ksw2Aligner.ExtendD2(q, t, 5, mat, 4, 2, 24, 0, -1, -1, 0), "e2 = 0: minimap2 accepts -E 2,0 (a long gap at a fixed cost); the second piece may have a zero extension penalty");
            Assert.DoesNotThrow(() => Ksw2Aligner.ExtendZ(q, t, 5, mat, 0, 1, -1, -1), "q = 0, e = 1 is the smallest valid pair");
        });
    }

    // ------------------------------------------------------------------ A4-2: the allocations added with the EMBOSS / minimap2 / nucmer ports (B09 F158)

    private sealed class FakeAnchors(int count) : IReadOnlyList<MinimapAnchor>
    {
        public int Count => count;

        public MinimapAnchor this[int index] => throw new InvalidOperationException("the limit is checked before an anchor is read");

        public IEnumerator<MinimapAnchor> GetEnumerator() => throw new InvalidOperationException("the limit is checked before the anchors are enumerated");

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
    }

    [Test]
    public void Minimap2CollectAnchors_PolyAAndTandemQueries_AreRefusedBeforeTheAnchorsAreStored()
    {
        // a 100 kb poly-A query against itself: every position is a minimizer and mid_occ is 100 001 - 10^10 anchors (OutOfMemoryException after 10 s before A4-2)
        string polyA = new('A', 100_000);
        var ex = Refused<ArgumentException>(() => Minimap2Chaining.CollectAnchors([polyA], polyA));      // A5-4: the deterministic proxy is the allocation bound of Refused (10^10 anchors were an OutOfMemoryException after 10 s)
        Assert.That(ex.Message, Does.StartWith("minimap2 seeding stores one anchor").And.Contain("SequenceAligner.MinimapMaxAnchors").And.Contain("mm_seed_select").And.Contain("8,388,608"));
        Assert.That(ex.ParamName, Is.EqualTo("query"));
        Assert.That(SequenceAligner.MinimapMaxAnchors, Is.EqualTo(8_388_608L));
        // a tandem repeat (period 2, 60 kb) against itself: hundreds of millions of anchors (the audit's 57 kb tandem repeat: 45 million, 5.2 GB, 10.8 s)
        string tandem = string.Concat(Enumerable.Repeat("AC", 30_000));
        Refused<ArgumentException>(() => Minimap2Chaining.CollectAnchors([tandem], tandem));
        // the anchors of an ordinary pair are untouched
        string target = Dna(20_000, 11);
        Assert.That(Minimap2Chaining.CollectAnchors([target], target[5_000..15_000]), Has.Count.EqualTo(1866), "A5-4: the number of 'CN' lines of `minimap2 --print-seeds t.fa q.fa` (2.26) on these inputs");
    }

    [Test]
    public void Minimap2Chain_AboveTheAnchorLimit_IsRefusedBeforeAnyAnchorIsRead()
    {
        var ex = Refused<ArgumentException>(() => Minimap2Chaining.Chain(new FakeAnchors(8_388_609)));
        Assert.That(ex.Message, Does.StartWith("minimap2 chaining keeps four arrays").And.Contain("SequenceAligner.MinimapMaxAnchors"));
        Assert.That(ex.ParamName, Is.EqualTo("anchors"));
    }

    [Test]
    public void NucmerCluster_OfMoreThan100000Matches_IsRefusedBeforeTheQuadraticProgram()
    {
        // 100 001 collinear matches, 30 apart: ONE cluster whose chain dynamic program is O(c^2) (the program's own loop)
        var matches = Enumerable.Range(0, 100_001).Select(i => ((long)(1 + 30 * i), (long)(1 + 30 * i), 20L)).ToList();
        long before = GC.GetAllocatedBytesForCurrentThread();
        var ex = Assert.Throws<ArgumentException>(() => NucmerClusters.Cluster(matches))!;      // A5-4: refused before the O(c^2) program runs; the allocation bound below replaces the wall-clock one
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(512L << 20));
        Assert.That(ex.Message, Does.StartWith("mgaps extracts the chains of a cluster").And.Contain("SequenceAligner.NucmerMaxClusterMatches"));
        Assert.That(SequenceAligner.NucmerMaxClusterMatches, Is.EqualTo(100_000L));
        var small = NucmerClusters.Cluster(matches.Take(2_000).ToList());
        Assert.That(small.Select(c => c.Matches.Count), Is.EqualTo(new[] { 2_000 }), "a cluster of 2 000 collinear matches is untouched: one cluster holding all of them");
    }
}
