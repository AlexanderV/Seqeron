using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A5-2, B09 F168): the quadratic pair dynamic programs of the T-Coffee port (<c>biphasic_pair_wise</c>, <c>glocal_pair_wise</c>, <c>gotoh_pair_wise_lgp</c>, <c>gotoh_pair_wise</c> /
/// <c>_sw</c>, <c>viterbi_pair_wise</c>, the FASTA band of <c>fast_pair</c>) refuse an input whose matrices exceed <see cref="SequenceAligner.TCoffeePairDpMaxBytes"/> with ONE
/// <see cref="ArgumentException"/>, in 64-bit arithmetic BEFORE anything is allocated (before: a bare <see cref="OutOfMemoryException"/> after seconds, or a 32-bit overflow); every one polls the cancellation token
/// once per row; <c>-lalign_n_top</c> = <see cref="int.MaxValue"/> means "every alignment" (the vertex list is lazily grown). The program has no limit of its own: the installed t_coffee 13.41 stops with
/// "FAILED TO ALLOCATE REQUIRED MEMORY (vmalloc)" (exit status 1; a segmentation fault for gotoh_pair_wise at 40 000 x 40 000) under a 4 GB limit.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_TCoffeeLimits_Tests
{
    private const long MaxAllocatedBytes = 64L << 20;

    private static string Dna(int n, int seed)
    {
        var r = new Random(seed);
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = "ACGT"[r.Next(4)];
        return new string(c);
    }

    private static readonly Func<int, int, int> Cost = (i, j) => (i * 7 + j * 13) % 11 - 4;

    /// <summary>The exception of <paramref name="call"/>, asserting it allocated (almost) nothing.</summary>
    private static ArgumentException Refused(Action call)
    {
        long before = GC.GetAllocatedBytesForCurrentThread();
        var ex = Assert.Throws<ArgumentException>(() => call())!;
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(MaxAllocatedBytes), "the limit is checked before the matrices are allocated");
        Assert.That(ex.Message, Does.Contain("SequenceAligner.TCoffeePairDpMaxBytes").And.Contain("2,147,483,648").And.Contain("FAILED TO ALLOCATE REQUIRED MEMORY (vmalloc)"));
        return ex;
    }

    // ------------------------------------------------------------------ limits

    [Test]
    public void BiphasicPair_Above2GiB_IsRefused_ThroughEveryEntryPoint()
    {
        Assert.That(SequenceAligner.TCoffeePairDpMaxBytes, Is.EqualTo(2_147_483_648L));
        string a = Dna(7_000, 1), b = Dna(7_000, 2);                // 7 001 x 7 001 x 48 bytes = 2.35 GB
        var ex = Refused(() => SequenceAligner.TCoffeeBiphasicPairAlign(a, b));
        Assert.That(ex.Message, Does.StartWith("T-Coffee biphasic_pair_wise keeps twelve int matrices").And.Contain("7,000 x 7,000 positions").And.Contain("slow_pair"));
        // the program: -method biphasic_pair (20 000 x 20 000 residues was an OutOfMemoryException after 1.4 s)
        Refused(() => SequenceAligner.MultipleAlignTCoffeeProgram([a, b], TCoffeeProgramOptions.FromCommandLine("-method biphasic_pair")));
        // 64-bit arithmetic: no overflow at the largest int lengths, nothing is allocated
        Refused(() => SequenceAligner.TcBiphasicMasks(int.MaxValue - 1, int.MaxValue - 1, Cost, -10, -1));
        // an ordinary pair is untouched
        var (ra, rb) = SequenceAligner.TCoffeeBiphasicPairAlign("ACGTACGTTAGC", "ACGTCGTTAGC");
        Assert.That(ra.Replace("-", ""), Is.EqualTo("ACGTACGTTAGC"));
        Assert.That(rb.Replace("-", ""), Is.EqualTo("ACGTCGTTAGC"));
    }

    [Test]
    public void FastPair_FlaggedDiagonalBand_IsRefused()
    {
        // two random 20 000-base sequences flag nearly every diagonal: five int matrices of 40 001 rows x 40 000 columns (OutOfMemoryException after 3 s)
        string a = Dna(20_000, 3), b = Dna(20_000, 4);
        var ex = Refused(() => SequenceAligner.MultipleAlignTCoffeeProgram([a, b], TCoffeeProgramOptions.FromCommandLine("-method fast_pair")));
        Assert.That(ex.Message, Does.StartWith("T-Coffee fasta_pair_wise (the banded Gotoh of fast_pair) keeps five int matrices").And.Contain("flagged diagonals").And.Contain("-ndiag"));
        // the band itself, direct: 60 000 flagged diagonals over 30 000 x 30 000 positions
        var list = new int[60_001];
        list[0] = 60_000;
        Refused(() => SequenceAligner.TcMakeFastaGotoh(30_000, 30_000, list, Cost, -100, -10, 0));
    }

    [Test]
    public void GotohPairWise_TracebackMatrix_IsRefused_AndItsSwAndLgpAndGlocalForms()
    {
        var gotoh = Refused(() => SequenceAligner.TcGotohProfile(25_000, 25_000, Cost, -10, -1, 0, false));      // 25 001^2 x 4 bytes = 2.5 GB (the 20 000 x 20 000 case allocated 1.6 GB unchecked)
        Assert.That(gotoh.Message, Does.StartWith("T-Coffee gotoh_pair_wise keeps an int trace-back matrix"));
        var sw = Refused(() => SequenceAligner.TcGotohSwMasks(25_000, 25_000, Cost, -10, -1, 0));
        Assert.That(sw.Message, Does.StartWith("T-Coffee gotoh_pair_wise_sw keeps an int trace-back matrix"));
        var glocal = Refused(() => SequenceAligner.TcGlocalMasks(7_500, 7_500, Cost, -10, -1));                  // 7 501^2 x 40 bytes = 2.25 GB
        Assert.That(glocal.Message, Does.StartWith("T-Coffee glocal_pair_wise keeps ten int matrices"));
        var lgp = Refused(() => SequenceAligner.TcRequireLgpFits(9_000, 9_000));                                  // 9 001^2 x 32 bytes = 2.59 GB
        Assert.That(lgp.Message, Does.StartWith("T-Coffee gotoh_pair_wise_lgp keeps eight int matrices"));
        // a pair just below the limit of the traceback matrix is accepted by the guard: the largest int lengths are refused without an overflow
        Refused(() => SequenceAligner.TcGotohProfile(int.MaxValue - 1, int.MaxValue - 1, Cost, -10, -1, 0, false));
    }

    [Test]
    public void ViterbiPairWise_IsRefused_WithoutTheIntOverflowOfTheProgram()
    {
        // the program computes nmt x (I + 1) x (J + 1) in a 32-bit int (it wraps beyond 20 700 x 20 700) and 400 multiplications per cell of the match probability; a 25 000 x 25 000 call ran for minutes
        var ex = Refused(() => SequenceAligner.TcRequireViterbiFits(7_500, 7_500));                               // 7 501^2 x 44 bytes = 2.47 GB
        Assert.That(ex.Message, Does.StartWith("T-Coffee viterbi_pair_wise keeps the match probability").And.Contain("7,500 x 7,500 columns"));
        Refused(() => SequenceAligner.TcRequireViterbiFits(25_000, 25_000));
        Refused(() => SequenceAligner.TcRequireViterbiFits(int.MaxValue - 1, int.MaxValue - 1));
        Assert.DoesNotThrow(() => SequenceAligner.TcRequireViterbiFits(6_900, 6_900));                            // 6 901^2 x 44 = 2.095 GB < 2^31
    }

    [Test]
    public void LalignTopAlignments_IntMaxValue_MeansEveryAlignment()
    {
        // new Vertex[int.MaxValue] was an OutOfMemoryException ("Array dimensions exceeded"); the program calloc()s K vertices and stops with "FAILED TO ALLOCATE REQUIRED MEMORY (vmalloc)" for 2147483647 and 100000000
        string a = Dna(120, 5), b = Dna(120, 6);
        long before = GC.GetAllocatedBytesForCurrentThread();
        var every = SequenceAligner.TCoffeeLalignTopAlignments(a, b, int.MaxValue);
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(MaxAllocatedBytes));
        var big = SequenceAligner.TCoffeeLalignTopAlignments(a, b, 100_000);
        Assert.That(every, Is.Not.Empty);
        Assert.That(every.Select(x => (x.RowA, x.RowB, x.Score)), Is.EqualTo(big.Select(x => (x.RowA, x.RowB, x.Score))));
        Assert.That(SequenceAligner.TCoffeeLalignTopAlignments(a, b, int.MaxValue, protein: false).Count, Is.EqualTo(every.Count));
    }

    [Test]
    public void FracMinHashSketch_KAboveTheSequenceLength_IsAnEmptySketch_NotAnOutOfMemory()
    {
        // KmerAnalyzer.cs (B06's file, B09 for B06): new byte[k] for k = int.MaxValue
        const string s = "ACGTACGTTAGCACGTAGCTAGCTAGGATC";
        long before = GC.GetAllocatedBytesForCurrentThread();
        var huge = KmerAnalyzer.CreateFracMinHashSketch(s, int.MaxValue, 1);
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(MaxAllocatedBytes));
        Assert.That(huge.Hashes, Is.Empty);
        Assert.That(KmerAnalyzer.CreateFracMinHashSketch(s, s.Length + 1, 1).Hashes, Is.Empty);
        Assert.That(KmerAnalyzer.CreateFracMinHashSketch(s, s.Length, 1).Hashes, Has.Count.EqualTo(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.CreateFracMinHashSketch(s, 0, 1));
    }

    // ------------------------------------------------------------------ cancellation: a pre-cancelled token, and a cancellation at the middle poll of each method (PollCancellation: no sleep)

    private static readonly CancellationToken Cancelled = new(true);

    [Test]
    public void BiphasicPair_Cancels_PreCancelled_AndAtAMidRunPoll()
    {
        string a = Dna(200, 11), b = Dna(200, 12);
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.TCoffeeBiphasicPairAlign(a, b)));
        PollCancellation.CancelInside(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.TCoffeeBiphasicPairAlign(a, b)), "TcBiphasicMasks");
    }

    [Test]
    public void GlocalGotohSwAndGotohProfile_Cancel_PreCancelled_AndAtAMidRunPoll()
    {
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.TcGlocalMasks(100, 100, Cost, -10, -1)));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.TcGotohSwMasks(100, 100, Cost, -10, -1, 0)));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.TcGotohProfile(100, 100, Cost, -10, -1, 0, false)));
        PollCancellation.CancelInside(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.TcGlocalMasks(150, 150, Cost, -10, -1)), "TcGlocalMasks");
        PollCancellation.CancelInside(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.TcGotohSwMasks(150, 150, Cost, -10, -1, 0)), "TcGotohSwMasks");
        PollCancellation.CancelInside(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.TcGotohProfile(150, 150, Cost, -10, -1, 0, false)), "TcGotohProfile");
    }

    [Test]
    public void FastPair_Cancels_PreCancelled_AndAtAMidRunPoll_InsideTheBandAndTheDiagonalScan()
    {
        string a = Dna(400, 13), b = Dna(400, 14);
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignTCoffeeProgram([a, b], TCoffeeProgramOptions.FromCommandLine("-method fast_pair") with { CancellationToken = Cancelled }));
        PollCancellation.CancelInside(ct => SequenceAligner.MultipleAlignTCoffeeProgram([a, b], TCoffeeProgramOptions.FromCommandLine("-method fast_pair") with { CancellationToken = ct }), "TcMakeFastaGotoh");
        // a poly-A pair gives offsets1 x offsets2 pairs per k-tuple: the scan of the diagonals polls per occurrence
        string polyA = new('A', 300);
        PollCancellation.CancelInside(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.TcDiagonalsByKtup(polyA, polyA, 5, Cost)), "TcDiagonalsByKtup");
    }

    [Test]
    public void LgpAndViterbiDpModes_Cancel_PreCancelled_AndAtAMidRunPoll()
    {
        var seqs = new[] { Dna(150, 15), Dna(150, 16), Dna(150, 17) };
        foreach (var (mode, method) in new[] { ("gotoh_pair_wise_lgp", "TcGotohLgpMasks"), ("viterbi_pair_wise", "TcViterbiMasks") })
        {
            var options = TCoffeeProgramOptions.FromCommandLine($"-method slow_pair -dp_mode {mode}");
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(seqs, options with { CancellationToken = Cancelled }), mode);
            PollCancellation.CancelInside(ct => SequenceAligner.MultipleAlignTCoffeeProgram(seqs, options with { CancellationToken = ct }), method);
        }
    }

    [Test]
    public void LalignAndSlowPair_Cancel_AtAMidRunPoll()
    {
        string a = Dna(250, 18), b = Dna(250, 19);
        PollCancellation.CancelInside(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.TCoffeeLalignTopAlignments(a, b, 5)), "TcSim");
        PollCancellation.CancelInside(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.TCoffeeSlowPairAlign(a, b)), "MyersMillerState");
    }
}
