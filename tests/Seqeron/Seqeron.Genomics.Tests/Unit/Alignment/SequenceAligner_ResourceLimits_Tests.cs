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

    // ------------------------------------------------------------------ A3-2 limits

    [Test]
    public void Needle_Above2Pow27Cells_IsRefusedWithoutAllocating()
    {
        string a = Dna(12_000, 1), b = Dna(12_000, 2);              // 144 000 000 cells > 134 217 728
        var ex = Refused<ArgumentException>(() => SequenceAligner.NeedleAlign(a, b));
        Assert.That(ex.Message, Does.StartWith("EMBOSS needle keeps").And.Contain("144,000,000 cells").And.Contain("SequenceAligner.NeedleMaxCells").And.Contain("StretcherAlign"));
        Assert.That(SequenceAligner.NeedleMaxCells, Is.EqualTo(134_217_728L));
    }

    [Test]
    public void Water_Above2Pow28Cells_IsRefusedWithoutAllocating()
    {
        string a = Dna(17_000, 3), b = Dna(17_000, 4);              // 289 000 000 cells > 268 435 456
        var ex = Refused<ArgumentException>(() => SequenceAligner.WaterAlign(a, b));
        Assert.That(ex.Message, Does.StartWith("EMBOSS water keeps").And.Contain("SequenceAligner.WaterMaxCells").And.Contain("Supermatcher"));
    }

    [Test]
    public void Supermatcher_HugeWidth_IsRefusedWithoutAllocating()
    {
        var seq = new[] { new EmbossSequence("a", Dna(200, 5)) };
        var ex = Refused<ArgumentException>(() => SequenceAligner.Supermatcher(seq, seq, width: 1_000_000_000));
        Assert.That(ex.Message, Does.StartWith("EMBOSS supermatcher keeps").And.Contain("SequenceAligner.SupermatcherMaxPathCells"));
    }

    [Test]
    public void Poa_Above2Pow27Cells_IsRefusedWithoutAllocating()
    {
        string a = new('A', 12_000), b = new('C', 12_000);
        var ex = Refused<ArgumentException>(() => SequenceAligner.MultipleAlignPoa([a, b]));
        Assert.That(ex.Message, Does.StartWith("POA keeps").And.Contain("SequenceAligner.PoaMaxMatrixCells").And.Contain("2^31"));
    }

    [Test]
    public void ProbCons_LongPair_IsRefusedWithProbConsWording_NotTCoffees()
    {
        string a = string.Concat(Enumerable.Repeat("ACDEFGHIKL", 740)), b = string.Concat(Enumerable.Repeat("LKIHGFEDCA", 740));      // 7 400 x 7 400: 5 x 7 401^2 floats > 2^28
        var ex = Refused<ArgumentException>(() => SequenceAligner.MultipleAlignProbCons([a, b]));
        Assert.That(ex.Message, Does.StartWith("The ProbCons pair-HMM keeps").And.Contain("ProbCons 1.12 itself").And.Contain("SequenceAligner.TCoffeeProbaMaxMatrixFloats"));
        Assert.That(ex.Message, Does.Not.Contain("T-Coffee").And.Not.Contain("t_coffee").And.Not.Contain("proba_pair"));
        var tc = Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignTCoffeeProgram([a, b]));
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

    // ------------------------------------------------------------------ A3-3 penalties

    [Test]
    public void Needle_Water_NaNAndInfinitePenalties_AreRejected()
    {
        double[] bad = [double.NaN, double.PositiveInfinity, double.NegativeInfinity];
        Assert.Multiple(() =>
        {
            foreach (double v in bad)
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAlign("ACGT", "ACT", null, v, 0.5), $"needle gapOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAlign("ACGT", "ACT", null, 10, v), $"needle gapExtend {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAlign("ACGT", "ACT", null, 10, 0.5, true, v, 0.5), $"needle endOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAlign("ACGT", "ACT", null, 10, 0.5, true, 10, v), $"needle endExtend {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.WaterAlign("ACGT", "ACT", null, v, 0.5), $"water gapOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.WaterAlign("ACGT", "ACT", null, 10, v), $"water gapExtend {v}");
                var one = new[] { new EmbossSequence("a", "ACGTACGT") };
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAll(one, one, null, v), $"needleall gapOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAll(one, one, null, 10, 0.5, false, 10, 0.5, v), $"needleall minscore {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.Supermatcher(one, one, null, v), $"supermatcher gapOpen {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.Supermatcher(one, one, null, 10, v), $"supermatcher gapExtend {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.Supermatcher(one, one, null, 10, 0.5, 6, 16, v), $"supermatcher minscore {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, v, 0.5), $"MatrixScoring.Emboss open {v}");
                Assert.Throws<ArgumentOutOfRangeException>(() => MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, v), $"MatrixScoring.Emboss extend {v}");
            }

            // a NaN used to pass `x < min || x > max`: needle("ACGT", "ACT", gapOpen NaN) returned score 0, water(gapExtend NaN) returned 10
            Assert.That(SequenceAligner.NeedleAlign("ACGT", "ACT").Score, Is.EqualTo(SequenceAligner.NeedleAlign("ACGT", "ACT", null, 10, 0.5).Score));
        });
    }

    [Test]
    public void Stretcher_OverflowingPenalties_AreRefusedLikeTheScalarContract()
    {
        string a = "ACGT" + new string('A', 50), b = "ACT" + new string('C', 50);
        // installed stretcher -gapopen 1000000000 -gapextend 1000000000: "stretcher CheckScore failed", Score: 1949672837 (a 32-bit wrap; int.MaxValue gives 1687091105 or 2147483643)
        var ex = Assert.Throws<ArgumentException>(() => SequenceAligner.StretcherAlign(a, b, null, 1_000_000_000, 1_000_000_000));
        Assert.That(ex!.Message, Does.Contain("Score range too large").And.Contain("CheckScore"));
        Assert.Throws<ArgumentException>(() => SequenceAligner.StretcherAlign(a, b, null, int.MaxValue, int.MaxValue));
        Assert.Throws<ArgumentException>(() => SequenceAligner.StretcherAlign(a, b, null, 100_000, 1_000_000_000));
        // large but exact penalties are fine: the score is the exact integer (no wrap)
        Assert.That(SequenceAligner.StretcherAlign(a, b, null, 1_000_000, 1_000_000).Score, Is.LessThan(0));
        Assert.That(SequenceAligner.StretcherAlign("ACGT", "ACGT").Score, Is.GreaterThan(0));
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
            Assert.Throws<ArgumentOutOfRangeException>(() => Ksw2Aligner.ExtendD2(q, t, 5, mat, 4, 2, 24, 0, -1, -1, 0), "e2 = 0 (SIMD)");
            Assert.DoesNotThrow(() => Ksw2Aligner.ExtendZ(q, t, 5, mat, 0, 1, -1, -1), "q = 0, e = 1 is the smallest valid pair");
        });
    }
}
