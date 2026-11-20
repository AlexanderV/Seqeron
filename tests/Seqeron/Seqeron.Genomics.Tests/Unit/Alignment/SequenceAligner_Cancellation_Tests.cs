namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 / ALIGN-MULTI-001 / GENOMIC-* (audit A2-40, B09 F86): cooperative cancellation of the long-running APIs.
///
/// Contract: a call whose token is cancelled (before it starts, or while it runs) throws <see cref="OperationCanceledException"/> promptly - the dynamic
/// programs poll the token at their row / iteration boundaries - and returns no partial result; the token never leaks into a later call. The token is
/// supplied through <c>MsaOptions / ClustalWOptions / ClustalWProteinOptions / MuscleOptions / TCoffeeProgramOptions.CancellationToken</c> (the convention of
/// <see cref="System.Threading.Tasks.ParallelOptions.CancellationToken"/>), through <see cref="SequenceAligner.WithCancellation{T}(CancellationToken, Func{T})"/> for
/// every other aligner, and through the trailing <c>cancellationToken</c> parameter of the <see cref="GenomicAnalyzer"/> scans.
/// </summary>
[TestFixture]
// Audit A3-15: runs alone. The fixture builds 3 M-base suffix trees (~0.5 GB each) and measures how long cancelled calls take to stop; next to seven other workers (LevelOfParallelism 8 on 4 cores) the
// uninterruptible tree builds took longer than the limit once (flake "the cancelled call did not stop within the limit") and the sum of the fixtures' working sets once ended in an OOM kill.
[NonParallelizable]
public class SequenceAligner_Cancellation_Tests
{
    private static readonly CancellationToken Cancelled = new(true);

    private static string Random(int n, int seed, string alphabet = "ACGT")
    {
        var rng = new System.Random(seed);
        return new string(Enumerable.Range(0, n).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
    }

    private static T[] Pair<T>(Func<string, T> make, int n, string alphabet = "ACGT") =>
        new[] { make(Random(n, 1, alphabet)), make(Random(n, 2, alphabet)) };

    /// <summary>
    /// A5-4: the call runs ONCE with a token that is cancelled at the third poll (of <paramref name="pollInside"/> when named) - <see cref="PollCancellation"/>: deterministic, no sleep, no thread, no time bound (the
    /// earlier version cancelled after a wall-clock delay, which under CPU load landed in another stage). The call must stop with <see cref="OperationCanceledException"/> raised by that poll and run no later poll.
    /// </summary>
    private static void AssertCancelsPromptly(Func<CancellationToken, object?> work, string pollInside, int nth = 3)
    {
        string stack = PollCancellation.CancelInside(work, pollInside, nth);
        Assert.That(stack, Does.Contain(pollInside), "the cancelling poll is inside " + pollInside);
    }

    // ------------------------------------------------------------------ pre-cancelled

    [Test]
    public void PreCancelled_MultipleAlignFamily_Throws()
    {
        var dna = Pair(s => new DnaSequence(s), 60).Append(new DnaSequence(Random(60, 3))).ToArray();
        var o = new MsaOptions { CancellationToken = Cancelled };
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlign(dna, null, o));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignAnchored(dna, null, o));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignProgressive(dna, null, o));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignIterative(dna, null, o));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignConsistency(dna, null, o));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignExact(dna, null, o));
        var protein = new[] { "MKVLAAGIVG", "MKVLAGIVG", "MKILAAGVG" };
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlign(protein, new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -1), o));
    }

    [Test]
    public void PreCancelled_ProgramPorts_Throw()
    {
        var dna = new[] { "ACGTACGTAC", "ACGTACGAAC", "ACGGACGTAC" };
        var protein = new[] { "MKVLAAGIVG", "MKVLAGIVG", "MKILAAGVG" };
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignClustalW(dna, new ClustalWOptions { CancellationToken = Cancelled }));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignClustalW(protein, new ClustalWProteinOptions { CancellationToken = Cancelled }));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignMuscle(dna, new MuscleOptions { CancellationToken = Cancelled }));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlignTCoffeeProgram(dna, new TCoffeeProgramOptions { CancellationToken = Cancelled }));
    }

    [Test]
    public void PreCancelled_WithCancellation_CoversEveryOtherAligner()
    {
        const string a = "ACGTACGTTAGCAC", b = "ACGTCGTTAGCAC";
        var matrix = new MatrixScoring(SubstitutionMatrix.EdnaFull, -5, -2);
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.LocalAlignTopN(a, b, null, 3, 1)));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.LocalAlignAffineCoOptimal(a, b, matrix, 5)));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.GlobalAlignAffine(a, b, matrix)));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.LocalAlign(a, b)));
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.GlobalAlignGeneralGap(a, b, k => -2.0 - Math.Sqrt(k))));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.WithCancellation<int>(Cancelled, null!));
    }

    [Test]
    public void PreCancelled_GenomicAnalyzerScans_Throw()
    {
        var s = new DnaSequence(Random(5000, 7));
        var t = new DnaSequence(Random(5000, 8));
        Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindMaximalRepeats(s, 4, Cancelled).ToList());
        Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindRepeats(s, 4, Cancelled).ToList());
        Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindCommonRegions(s, t, 5, Cancelled).ToList());
        Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindKnownMotifHits(s, new[] { "ACGT" }, true, false, 1, Cancelled));
        Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindKnownMotifHits(s, new[] { "ACGT" }, false, false, 0, Cancelled));
        Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindOpenReadingFrames(s, 30, OrfSearchOptions.Default, Cancelled).ToList());
        Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindOpenReadingFrames(s.Sequence, 30, OrfSearchOptions.Default, Cancelled).ToList());
    }

    [Test]
    public void UncancelledToken_GivesTheSameResultAsNoToken_AndTheScopeDoesNotLeak()
    {
        var dna = new[] { "ACGTACGTAC", "ACGTACGAAC", "ACGGACGTAC", "ACGTTCGTAC" }.Select(x => new DnaSequence(x)).ToArray();
        var plain = SequenceAligner.MultipleAlign(dna);
        using var cts = new CancellationTokenSource();
        var withToken = SequenceAligner.MultipleAlign(dna, null, new MsaOptions { CancellationToken = cts.Token });
        Assert.That(withToken.AlignedSequences, Is.EqualTo(plain.AlignedSequences));
        Assert.That(withToken.TotalScore, Is.EqualTo(plain.TotalScore));

        // a cancelled call must not poison the next call of the same thread
        Assert.Catch<OperationCanceledException>(() => SequenceAligner.MultipleAlign(dna, null, new MsaOptions { CancellationToken = Cancelled }));
        Assert.That(SequenceAligner.MultipleAlign(dna).AlignedSequences, Is.EqualTo(plain.AlignedSequences));
        Assert.That(AlignmentCancellation.Token.CanBeCanceled, Is.False);
        Assert.That(GenomicAnalyzer.FindMaximalRepeats(dna[0], 2, cts.Token).Select(r => r.Sequence),
            Is.EqualTo(GenomicAnalyzer.FindMaximalRepeats(dna[0], 2).Select(r => r.Sequence)));
    }

    // ------------------------------------------------------------------ cancelled while running (A5-4: at the third poll, see PollCancellation; the inputs only have to be long enough to reach it)

    [Test]
    public void CancelledWhileRunning_StarMultipleAlign_StopsPromptly()
    {
        var dna = Enumerable.Range(0, 5).Select(i => new DnaSequence(Random(600, 10 + i))).ToArray();
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlign(dna, null, new MsaOptions { GapModel = MsaGapModel.Affine, CancellationToken = ct }), pollInside: "StarAlignCore");
    }

    [Test]
    public void CancelledWhileRunning_ProgressiveAndIterative_StopPromptly()
    {
        var dna = Enumerable.Range(0, 4).Select(i => new DnaSequence(Random(500, 20 + i))).ToArray();
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignProgressive(dna, null, new MsaOptions { CancellationToken = ct }), pollInside: "PairwiseIdentityDistance");
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignIterative(dna, null, new MsaOptions { CancellationToken = ct }), pollInside: "PairwiseIdentityDistance");
    }

    [Test]
    public void CancelledWhileRunning_ClustalW_StopsPromptly()
    {
        var seqs = Pair(s => s, 3000);
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWOptions { CancellationToken = ct }), pollInside: "ForwardPass");
    }

    [Test]
    public void CancelledWhileRunning_Muscle_StopsPromptly()
    {
        var seqs = Pair(s => s, 2000);
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignMuscle(seqs, new MuscleOptions { CancellationToken = ct }), pollInside: "NwSmall");
    }

    [Test]
    public void CancelledWhileRunning_TCoffeeProbaPair_StopsPromptly()
    {
        var seqs = Pair(s => s, 400, "ACDEFGHIKLMNPQRSTVWY");
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignTCoffeeProgram(seqs, new TCoffeeProgramOptions { CancellationToken = ct }), pollInside: "TcProbaForward");
    }

    [Test]
    public void CancelledWhileRunning_TopNCoOptimalAndGeneralGap_StopPromptly()
    {
        string a = Random(600, 31), b = Random(600, 32);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.LocalAlignTopN(a, b, null, 5, 1)), pollInside: "LocalDeclumper");
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.LocalAlignAffine(a, b)), pollInside: "AffineFill");
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.GlobalAlignAffine(Random(800, 33), Random(800, 34))), pollInside: "AffineFill");
        string g1 = Random(1200, 35), g2 = Random(1200, 36);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.GlobalAlignGeneralGap(g1, g2, k => -2.0 - 3.0 * Math.Sqrt(k))), pollInside: "GeneralGapCore");
    }

    [Test]
    public void CancelledWhileRunning_GenomicAnalyzerScans_StopPromptly()
    {
        var big = new DnaSequence(Random(300_000, 41));
        var other = new DnaSequence(Random(30_000, 42));
        // the ORF scan is the lazy OrfScanner (Core: it polls the token itself, every 65 536 codons): cancel the token after the first ORF has been returned, the enumeration must then stop at the next poll
        using (var cts = new CancellationTokenSource())
        {
            int returned = 0;
            Assert.Catch<OperationCanceledException>(() =>
            {
                foreach (var orf in GenomicAnalyzer.FindOpenReadingFrames(big, 30, OrfSearchOptions.Default, cts.Token))
                {
                    if (returned++ == 0)
                        cts.Cancel();
                }
            });
            Assert.That(returned, Is.GreaterThanOrEqualTo(1));
        }

        // the suffix tree of a sequence is built once, uninterruptibly (DnaSequence.SuffixTree), on the first call that needs it: the cancelled call stops at the poll after the build
        AssertCancelsPromptly(ct => GenomicAnalyzer.FindKnownMotifHits(big, Enumerable.Range(0, 200).Select(i => Random(8, 100 + i)), true, false, 3, ct), "FindKnownMotifHits", nth: 1);
        AssertCancelsPromptly(ct => GenomicAnalyzer.FindCommonRegions(big, other, 3, ct).Count(), "FindCommonRegions", nth: 1);
    }

    // ------------------------------------------------------------------ A3-5: the pairwise programs, ksw2, T-Coffee DPs, MUSCLE refinement, anchors, distmat (B09 F133)

    [Test]
    public void PreCancelled_A3_5_EveryPortedProgram_Throws()
    {
        const string a = "ACGTACGTTAGCACGTAGCTAGCTAGGATC", b = "ACGTCGTTAGCACGTAGCTTAGCTAGGATC";
        var one = new[] { new EmbossSequence("a", a) };
        var q = Ksw2Aligner.EncodeNucleotides(a);
        var t = Ksw2Aligner.EncodeNucleotides(b);
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, 4, 1);
        Assert.Multiple(() =>
        {
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.NeedleAlign(a, b)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.WaterAlign(a, b)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.StretcherAlign(a, b)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.Supermatcher(one, one)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.NeedleAll(one, one)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => Ksw2Aligner.ExtendZ(q, t, 5, mat, 4, 2, -1, -1)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => Ksw2Aligner.ExtendD(q, t, 5, mat, 4, 2, 24, 1, -1, -1)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => Ksw2Aligner.ExtendZ2(q, t, 5, mat, 4, 2, -1, -1, 0)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => Ksw2Aligner.ExtendD2(q, t, 5, mat, 4, 2, 24, 1, -1, -1, 0)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.GlobalAlignDualAffine(a, b, DualAffineScoring.Minimap2)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.TCoffeeSlowPairAlign(a, b)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => MultipleAlignmentStatistics.CalculateDistanceMatrix(new[] { "ACGT", "ACGA", "ACCA" }, new MsaDistanceOptions())));
        });
    }

    [Test]
    public void CancelledWhileRunning_Needle_Water_StopPromptly()
    {
        string a = Random(600, 51), b = Random(600, 52);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.NeedleAlign(a, b)), pollInside: "NeedleCore");
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.WaterAlign(a, b)), pollInside: "WaterCore");
    }

    [Test]
    public void CancelledWhileRunning_Stretcher_StopsPromptly()
    {
        string a = Random(3000, 53), b = Random(3000, 54);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.StretcherAlign(a, b)), pollInside: "StretcherEngine");
    }

    [Test]
    public void CancelledWhileRunning_Supermatcher_StopsPromptly()
    {
        var seq = new[] { new EmbossSequence("a", Random(3000, 55)) };
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.Supermatcher(seq, seq, width: 600)), pollInside: "SwPathCalcFast");
    }

    [Test]
    public void CancelledWhileRunning_NeedleAll_StopsPromptly()
    {
        var set = Enumerable.Range(0, 60).Select(i => new EmbossSequence("s" + i, Random(600, 200 + i))).ToArray();
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.NeedleAll(set, set)), pollInside: "NeedleAll");
    }

    [Test]
    public void CancelledWhileRunning_Ksw2AndDualAffine_StopPromptly()
    {
        var q = Ksw2Aligner.EncodeNucleotides(Random(3000, 56));
        var t = Ksw2Aligner.EncodeNucleotides(Random(3000, 57));
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, 4, 1);
        const Ksw2Option so = Ksw2Option.ScoreOnly;
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => Ksw2Aligner.ExtendZ(q, t, 5, mat, 4, 2, -1, -1, so)), pollInside: "ExtendZ");
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => Ksw2Aligner.ExtendD(q, t, 5, mat, 4, 2, 24, 1, -1, -1, so)), pollInside: "ExtendD");
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => Ksw2Aligner.ExtendZ2(q, t, 5, mat, 4, 2, -1, -1, 0, so)), pollInside: "ExtendZ2");
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => Ksw2Aligner.ExtendD2(q, t, 5, mat, 4, 2, 24, 1, -1, -1, 0, so)), pollInside: "ExtendD2");
        string a = Random(600, 58), b = Random(600, 59);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.GlobalAlignDualAffine(a, b, DualAffineScoring.Minimap2)), pollInside: "DualAffineCore");
    }

    [Test]
    public void CancelledWhileRunning_TCoffeePrimaryLibrary_StopsPromptly()
    {
        string a = Random(2000, 60), b = Random(2000, 61);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.TCoffeeSlowPairAlign(a, b)), pollInside: "MyersMillerState");
    }

    [Test]
    public void CancelledWhileRunning_TCoffeeProgressiveDp_StopsPromptly()
    {
        // TcGotohProfile (TCoffeeDp.cs) is the profile-profile Gotoh program of the progressive stage; the library stage would dominate any end-to-end run, so the DP is called directly
        Func<int, int, int> cost = (i, j) => (i * 7 + j * 13) % 11 - 4;
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.TcGotohProfile(600, 600, cost, -10, -1, 0, false)), pollInside: "TcGotohProfile");
    }

    [Test]
    public void CancelledWhileRunning_MuscleRefinement_StopsPromptly()
    {
        var seqs = Enumerable.Range(0, 24).Select(i => Random(300, 300 + i, "ACDEFGHIKLMNPQRSTVWY")).ToArray();
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignMuscle(seqs, new MuscleOptions { MaxIterations = 16, SequenceType = MuscleSequenceType.Protein, CancellationToken = ct }), pollInside: "Refine");
    }

    [Test]
    public void CancelledWhileRunning_AnchoredAlignment_StopsPromptly()
    {
        // 3 x 150 kb with a substitution every 40 bases: thousands of anchors
        var baseSeq = Random(150_000, 81).ToCharArray();
        var other = (char[])baseSeq.Clone();
        for (int i = 20; i < other.Length; i += 40) other[i] = other[i] == 'A' ? 'C' : 'A';
        var dna = new[] { new DnaSequence(new string(baseSeq)), new DnaSequence(new string(other)), new DnaSequence(new string(baseSeq)) };
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignAnchored(dna, null, new MsaOptions { CancellationToken = ct }), pollInside: "ChainAnchors");
    }

    [Test]
    public void CancelledWhileRunning_DistanceMatrix_StopsPromptly()
    {
        var rows = Enumerable.Range(0, 40).Select(i => Random(800, 400 + i)).ToArray();
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => MultipleAlignmentStatistics.CalculateDistanceMatrix(rows, new MsaDistanceOptions())), pollInside: "CalcMatch");
    }

    // ------------------------------------------------------------------ A4-1: the ports added in WP36 - WP38 poll the ambient token (one test per class: cancelled while running, the poll named)

    [Test]
    public void PreCancelled_A4_1_EveryPortAddedInWp36ToWp38_Throws()
    {
        string a = Random(300, 71), b = Random(300, 72);
        var anchors = new[] { new MinimapAnchor(0, false, 100, 100, 15), new MinimapAnchor(0, false, 120, 120, 15) };
        var matches = new (long, long, long)[] { (1, 1, 20), (31, 31, 20) };
        Assert.Multiple(() =>
        {
            Assert.Catch<OperationCanceledException>(() => DotMatcher.Compute(a, b, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => Minimap2Chaining.Chain(anchors, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => Minimap2Chaining.CollectAnchors([a], b, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => Minimap2Chaining.MapChains([a], b, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => NucmerClusters.Cluster(matches, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.Merger(a, b, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.LocalAlignAffineTopN(a, b, new ScoringMatrix(5, -4, -10, -4), 3, 0, TopNTieOrder.Sim)));
            Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindTranslatedPatternHits(a, "L-x-S", "F", 0, Cancelled));
            Assert.Catch<OperationCanceledException>(() => MultipleAlignmentDisplay.ShowAlign(["a", "b"], [a, b], new ShowAlignOptions { CancellationToken = Cancelled }));
            Assert.Catch<OperationCanceledException>(() => MotifFinder.GenerateEmbossAmbiguousConsensus([a, b], cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindNucmerClusters(a, b, 8, cancellationToken: Cancelled));
        });
    }

    [Test]
    public void CancelledWhileRunning_DotMatcher_StopsPromptly()
    {
        string a = Random(3000, 73), b = Random(3000, 74);
        AssertCancelsPromptly(ct => DotMatcher.Compute(a, b, cancellationToken: ct), pollInside: "Diagonal");
    }

    [Test]
    public void CancelledWhileRunning_Minimap2Chain_StopsPromptly()
    {
        // 100 000 anchors on one diagonal, 20 apart: about 250 predecessors inside the 5 000-base gap limit for every anchor
        var anchors = new MinimapAnchor[100_000];
        for (int i = 0; i < anchors.Length; i++)
            anchors[i] = new MinimapAnchor(0, false, 100 + 20 * i, 100 + 20 * i, 15);
        AssertCancelsPromptly(ct => Minimap2Chaining.Chain(anchors, cancellationToken: ct), pollInside: "Chain");
    }

    [Test]
    public void CancelledWhileRunning_Minimap2CollectAnchors_StopsPromptly()
    {
        string target = Random(300_000, 75), query = Random(300_000, 76);
        AssertCancelsPromptly(ct => Minimap2Chaining.CollectAnchors([target], query, cancellationToken: ct), pollInside: "CollectAnchors");
    }

    [Test]
    public void CancelledWhileRunning_NucmerClustersCluster_StopsPromptly()
    {
        // 3 000 collinear matches in ONE cluster: the program's own O(c^2) chain dynamic program (300 000 of them did not return in two minutes)
        var matches = Enumerable.Range(0, 3_000).Select(i => ((long)(1 + 30 * i), (long)(1 + 30 * i), 20L)).ToList();
        AssertCancelsPromptly(ct => NucmerClusters.Cluster(matches, cancellationToken: ct), pollInside: "ProcessCluster");
    }

    [Test]
    public void CancelledWhileRunning_Merger_StopsPromptly()
    {
        string a = Random(600, 77), b = Random(600, 78);
        AssertCancelsPromptly(ct => SequenceAligner.Merger(a, b, cancellationToken: ct), pollInside: "MergerPathCalc");
    }

    [Test]
    public void CancelledWhileRunning_TopNSim_StopsPromptly()
    {
        string a = Random(600, 79), b = Random(600, 80);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.LocalAlignAffineTopN(a, b, new ScoringMatrix(5, -4, -10, -4), 40, 0, TopNTieOrder.Sim)), pollInside: "SimMatcher");
    }

    [Test]
    public void CancelledWhileRunning_MegaMergeMinimalMatchSet_StopsPromptly()
    {
        // 3 000 matches that all survive (the same diagonal, 100 apart): every kept match is compared with every later one
        var list = Enumerable.Range(0, 3_000).Select(i => new GenomicAnalyzer.MegaMatch(100 * i, 100 * i, 30)).ToList();
        AssertCancelsPromptly(ct => { GenomicAnalyzer.MinimalMatchSet(list, 20, ct); return list; }, pollInside: "MinimalMatchSet");
    }

    [Test]
    public void CancelledWhileRunning_FuzztranTranslatedPatternHits_StopsPromptly()
    {
        string nt = Random(600_000, 81);
        AssertCancelsPromptly(ct => GenomicAnalyzer.FindTranslatedPatternHits(nt, "L-x-S", "F", 0, ct), pollInside: "FindTranslatedPatternHits");
    }

    [Test]
    public void CancelledWhileRunning_FuzzproProteinPatternHits_StopsPromptly_InsideThePositionLoop()
    {
        string protein = Random(300_000, 86, "ACDEFGHIKLMNPQRSTVWY");
        AssertCancelsPromptly(ct => GenomicAnalyzer.FindProteinPatternHits(protein, "L-x-S", ct), pollInside: "ProteinPatternHitsCore");
    }

    [Test]
    public void CancelledWhileRunning_ShowAlign_StopsPromptly()
    {
        // the consensus of 150 000 columns is computed first (its own poll), then 150 000 one-column blocks are written (the poll of Output)
        string row1 = Random(150_000, 82), row2 = Random(150_000, 83);
        AssertCancelsPromptly(ct => MultipleAlignmentDisplay.ShowAlign(["a", "b"], [row1, row2], new ShowAlignOptions { Width = 1, CancellationToken = ct }), pollInside: "Output");
        AssertCancelsPromptly(ct => MultipleAlignmentDisplay.ShowAlign(["a", "b"], [row1, row2], new ShowAlignOptions { Width = 1, CancellationToken = ct }), pollInside: "EmbossConsensusCore");
    }

    [Test]
    public void CancelledWhileRunning_Consambig_StopsPromptly()
    {
        var rows = Enumerable.Range(0, 20).Select(i => Random(40_000, 90 + i)).ToList();
        AssertCancelsPromptly(ct => MotifFinder.GenerateEmbossAmbiguousConsensus(rows, cancellationToken: ct), pollInside: "GenerateEmbossAmbiguousConsensus");
    }
}
