using System.Diagnostics;

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
    /// Runs <paramref name="work"/> on a dedicated thread and cancels it from the test thread after <paramref name="delayMs"/> (a sleep, not a timer callback: the
    /// pool may be saturated by parallel fixtures); the call must stop with OperationCanceledException well before it would finish.
    /// </summary>
    private static void AssertCancelsPromptly(Func<CancellationToken, object> work, int delayMs = 100, int limitMs = 20000)
    {
        using var cts = new CancellationTokenSource();
        var token = cts.Token;
        var task = Task.Factory.StartNew(() => work(token), CancellationToken.None, TaskCreationOptions.LongRunning, TaskScheduler.Default);
        Thread.Sleep(delayMs);
        var sw = Stopwatch.StartNew();
        cts.Cancel();
        bool finished = task.ContinueWith(static _ => { }, TaskScheduler.Default).Wait(limitMs);   // waits without rethrowing the task's exception
        Assert.That(finished, Is.True, "the cancelled call did not stop within the limit");
        Assert.That(task.IsFaulted || task.IsCanceled, Is.True, "the call finished before it was cancelled: enlarge the input");
        var inner = task.Exception?.Flatten().InnerException;
        Assert.That(inner, Is.InstanceOf<OperationCanceledException>(), "the call must stop with OperationCanceledException");
        Assert.That(sw.ElapsedMilliseconds, Is.LessThan(limitMs));
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

    // ------------------------------------------------------------------ cancelled while running (the same calls would run for many seconds)

    [Test]
    public void CancelledWhileRunning_StarMultipleAlign_StopsPromptly()
    {
        var dna = Enumerable.Range(0, 5).Select(i => new DnaSequence(Random(6000, 10 + i))).ToArray();
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlign(dna, null, new MsaOptions { GapModel = MsaGapModel.Affine, CancellationToken = ct }));
    }

    [Test]
    public void CancelledWhileRunning_ProgressiveAndIterative_StopPromptly()
    {
        var dna = Enumerable.Range(0, 4).Select(i => new DnaSequence(Random(4000, 20 + i))).ToArray();
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignProgressive(dna, null, new MsaOptions { CancellationToken = ct }));
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignIterative(dna, null, new MsaOptions { CancellationToken = ct }));
    }

    [Test]
    public void CancelledWhileRunning_ClustalW_StopsPromptly()
    {
        // two 30 kb sequences take minutes in Debug
        var seqs = Pair(s => s, 30000);
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWOptions { CancellationToken = ct }));
    }

    [Test]
    public void CancelledWhileRunning_Muscle_StopsPromptly()
    {
        var seqs = Pair(s => s, 20000);
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignMuscle(seqs, new MuscleOptions { CancellationToken = ct }));
    }

    [Test]
    public void CancelledWhileRunning_TCoffeeProbaPair_StopsPromptly()
    {
        var seqs = Pair(s => s, 2500, "ACDEFGHIKLMNPQRSTVWY");
        AssertCancelsPromptly(ct => SequenceAligner.MultipleAlignTCoffeeProgram(seqs, new TCoffeeProgramOptions { CancellationToken = ct }));
    }

    [Test]
    public void CancelledWhileRunning_TopNCoOptimalAndGeneralGap_StopPromptly()
    {
        string a = Random(4000, 31), b = Random(4000, 32);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.LocalAlignTopN(a, b, null, 5, 1)));
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.LocalAlignAffine(a, b)));
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.GlobalAlignAffine(Random(6000, 33), Random(6000, 34))));
        string g1 = Random(1200, 35), g2 = Random(1200, 36);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => SequenceAligner.GlobalAlignGeneralGap(g1, g2, k => -2.0 - 3.0 * Math.Sqrt(k))));
    }

    [Test]
    public void CancelledWhileRunning_GenomicAnalyzerScans_StopPromptly()
    {
        var big = new DnaSequence(Random(3_000_000, 41));
        var other = new DnaSequence(Random(300_000, 42));
        AssertCancelsPromptly(ct => GenomicAnalyzer.FindOpenReadingFrames(big, 30, OrfSearchOptions.Default, ct).Count(), 1);
        AssertCancelsPromptly(ct => GenomicAnalyzer.FindKnownMotifHits(big, Enumerable.Range(0, 200).Select(i => Random(8, 100 + i)), true, false, 3, ct), 20);
        AssertCancelsPromptly(ct => GenomicAnalyzer.FindCommonRegions(big, other, 3, ct).Count(), 20);
    }
}
