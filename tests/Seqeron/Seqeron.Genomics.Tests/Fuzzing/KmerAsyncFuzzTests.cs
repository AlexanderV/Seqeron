namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the K-mer area — ASYNCHRONOUS k-mer counting (KMER-ASYNC-001):
/// the cooperatively cancelable, progress-reporting, thread-pool variant of
/// synchronous k-mer counting.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang / deadlock on await, no
/// state corruption, no nonsense output, no swallowed cancellation, and no *unhandled*
/// runtime exception (IndexOutOfRangeException from k &gt; len, OutOfMemoryException on a
/// large input). Every input must resolve to EITHER a well-defined, theory-correct
/// result, OR a *documented, intentional* exception — for the async surface that means
/// the awaited Task either completes with the documented result, surfaces the documented
/// ArgumentOutOfRangeException (k ≤ 0), or surfaces the documented
/// OperationCanceledException on a signaled token. A raw runtime exception, a hang on
/// await, a swallowed cancellation, or an async result differing from the synchronous
/// reference is a bug, not a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8
/// "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-ASYNC-001 — asynchronous k-mer counting
/// Checklist: docs/checklists/03_FUZZING.md, row 156.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty
///          (docs/checklists/03_FUZZING.md §Description). The checklist row names the
///          BE targets precisely: empty, k &gt; len, cancellation, huge input.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The async-k-mer-counting contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Asynchronous k-mer counting computes, off the calling thread, the count of every
/// length-k substring of a sequence. The "asynchronous" qualifier concerns EXECUTION
/// (thread-pool offload + cancellation), not the mathematical definition of the count:
/// the numeric result is byte-for-byte identical to the synchronous reference
/// (Asynchronous_K-mer_Counting.md §1, §2.2, §5.2). The API entry under test is
///   KmerAnalyzer.CountKmersAsync(string sequence, int k, CancellationToken, IProgress&lt;double&gt;?)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines 93–100),
/// a thread-pool wrapper that queues the synchronous cancellation/progress overload
///   KmerAnalyzer.CountKmers(string, int, CancellationToken, IProgress&lt;double&gt;?)
///   (KmerAnalyzer.cs lines 52–88) via `Task.Run(..., token)` (KmerAnalyzer.cs line 99).
///
/// THE FOUR DOCUMENTED INVARIANTS (Asynchronous_K-mer_Counting.md §2.4):
///   • INV-01 — `await CountKmersAsync(S, k)` equals `CountKmers(S, k)` (same keys AND
///     counts): the async method runs the identical synchronous algorithm via Task.Run.
///   • INV-02 — Σ counts = L − k + 1 for 1 ≤ k ≤ L: one k-mer per sliding-window start.
///   • INV-03 — a signaled cancellation token ⇒ awaiting the task throws
///     OperationCanceledException (TaskCanceledException is the wrapped/concrete form;
///     §2.2): cooperative cancellation via `ThrowIfCancellationRequested` /
///     pre-start cancellation through the Task.Run token.
///   • INV-04 — empty/null S or k &gt; L ⇒ empty result; k ≤ 0 ⇒
///     ArgumentOutOfRangeException (surfaced through the awaited task).
///
/// Documented parameter contract (Asynchronous_K-mer_Counting.md §3.1, §3.3, §6.1):
///   • Null / empty sequence → empty dictionary (L = 0 ⇒ no k-mers). The underlying
///     synchronous overload's `string.IsNullOrEmpty` guard runs BEFORE k is validated
///     (KmerAnalyzer.cs lines 58–59), so empty/null wins even with a degenerate k.
///   • k &gt; sequence.Length → empty dictionary (L − k + 1 ≤ 0; no window fits). The loop
///     bound `i ≤ L − k` is negative so the loop never runs and no Substring/Slice is
///     taken past the end — never an IndexOutOfRange on the async path (KmerAnalyzer.cs
///     lines 64–65, 73).
///   • k ≤ 0 with non-empty input → ArgumentOutOfRangeException(nameof(k)), surfaced
///     unchanged through the awaited Task (KmerAnalyzer.cs lines 61–62; §6.1).
///   • Token signaled BEFORE the call → Task.Run cancels work not yet started; awaiting
///     throws OperationCanceledException (§6.1, [3]).
///   • Token signaled DURING the run → observed at the next ThrowIfCancellationRequested
///     checkpoint (every CheckInterval = 1000 windows); awaiting throws
///     OperationCanceledException (§6.1, [2]).
///   • Lowercase / mixed case → same counts as uppercase (input is ToUpperInvariant'd).
/// The alphabet is NOT restricted — non-ACGT symbols are counted literally (§3.3, §6.2);
/// these tests exercise only the BE targets of THIS fuzz row (empty, k &gt; len,
/// cancellation, huge input).
///
/// The four checklist targets map to these documented behaviours:
///   • empty        → awaited Task completes with the empty dictionary, no crash, no hang
///                    (and empty + degenerate k still empty, NOT a throw).
///   • k &gt; len      → awaited Task completes with the empty dictionary — no window fits,
///                    no IndexOutOfRange / negative-length Slice on the async path.
///   • cancellation → a pre-signaled OR mid-run signaled token ⇒ the awaited Task throws
///                    OperationCanceledException; it does not hang and does not return a
///                    partial-garbage dictionary (INV-03).
///   • huge input   → a large sequence completes, terminates, does NOT overflow, and the
///                    awaited result matches the synchronous reference EXACTLY (INV-01);
///                    [CancelAfter] guards every awaited path against a deadlock/hang.
/// A positive-sanity test pins the cross-surface equality (INV-01) on a known sequence —
/// `await CountKmersAsync("ATGG", 3)` returns {ATG:1, TGG:1}, the §7.1 worked example,
/// equal to the synchronous CountKmers — and a pre-cancelled token throws the documented
/// OperationCanceledException.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class KmerAsyncFuzzTests
{
    #region Helpers

    /// <summary>
    /// Synchronous <see cref="IProgress{T}"/> — invokes the handler inline on the reporting
    /// thread instead of marshaling it asynchronously the way <see cref="Progress{T}"/> does.
    /// Used where a progress callback must run deterministically during the operation (e.g. to
    /// signal cancellation mid-run) and complete before the caller disposes its resources.
    /// </summary>
    private sealed class SynchronousProgress<T> : IProgress<T>
    {
        private readonly Action<T> _handler;
        public SynchronousProgress(Action<T> handler) => _handler = handler;
        public void Report(T value) => _handler(value);
    }

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Well-formed-result helper: a non-empty async count result is a true sliding-window
    /// k-mer tally — every key is a length-k window, every count is strictly positive, and
    /// the counts sum to the window count L − k + 1 (INV-02, Asynchronous_K-mer_Counting.md
    /// §2.4). An empty result (L = 0 or k &gt; L) is trivially well-formed.
    /// </summary>
    private static void AssertWellFormed(IReadOnlyDictionary<string, int> counts, int length, int k)
    {
        if (counts.Count == 0)
            return; // empty result for L = 0 or k > L is well-formed by definition.

        counts.Keys.Should().OnlyContain(key => key.Length == k,
            $"every emitted key must be a length-{k} window — no shorter/longer fragments");
        counts.Values.Should().OnlyContain(c => c > 0, "every observed k-mer has a strictly positive count");
        counts.Values.Sum().Should().Be(length - k + 1,
            "INV-02: the sum of all counts equals the sliding-window count L − k + 1");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  KMER-ASYNC-001 — asynchronous k-mer counting : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-ASYNC-001 — asynchronous k-mer counting

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The underlying synchronous
    /// overload's `string.IsNullOrEmpty` guard short-circuits to the empty dictionary
    /// BEFORE k is validated (KmerAnalyzer.cs lines 58–59), so empty/null input awaited
    /// through the Task NEVER throws and NEVER hangs — even when k is itself degenerate
    /// (k = 0). We pin that empty, null, and empty-with-degenerate-k all complete with the
    /// empty dictionary (INV-04; Asynchronous_K-mer_Counting.md §6.1: "Empty / null
    /// sequence → Empty dictionary"). [CancelAfter] guards the awaits against a deadlock.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public async Task CountKmersAsync_EmptyOrNullSequence_AwaitsEmptyAndDoesNotThrow()
    {
        var empty = await KmerAnalyzer.CountKmersAsync(string.Empty, 3);
        var fromNull = await KmerAnalyzer.CountKmersAsync(null!, 3);
        var emptyDegenerateK = await KmerAnalyzer.CountKmersAsync(string.Empty, 0);

        empty.Should().BeEmpty("an empty sequence has no windows; the guard returns an empty dictionary, not a crash");
        fromNull.Should().BeEmpty("null input is treated as empty, not as an error, by the underlying guard");
        emptyDegenerateK.Should().BeEmpty(
            "the empty/null guard runs BEFORE k is validated, so empty input wins even with a degenerate k = 0");
    }

    #endregion

    #region BE — Boundary: k > sequence length

    /// <summary>
    /// BE: k far larger than the sequence length must complete with the EMPTY dictionary,
    /// never an IndexOutOfRange on the async path. No length-k window fits, so the window
    /// count L − k + 1 is ≤ 0 and the loop bound `i ≤ L − k` is negative — the loop never
    /// runs, no Slice is taken past the end (KmerAnalyzer.cs lines 64–65, 73;
    /// Asynchronous_K-mer_Counting.md §6.1: "k &gt; L → Empty dictionary"). We pin the
    /// awaited result empty at a huge k AND at the exact off-by-one k = L + 1 boundary.
    /// [CancelAfter] guards the awaits against a deadlock.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public async Task CountKmersAsync_KGreaterThanSequenceLength_AwaitsEmpty()
    {
        var hugeK = await KmerAnalyzer.CountKmersAsync("ACGT", 1000);
        hugeK.Should().BeEmpty(
            "no length-1000 window fits a 4-base sequence; the awaited result is empty, never an out-of-range Slice");

        // k = L + 1 is the exact off-by-one boundary one past the last fitting window.
        var offByOne = await KmerAnalyzer.CountKmersAsync("ACGT", 5);
        offByOne.Should().BeEmpty(
            "k = L + 1 is one past the last fitting window; still empty, never an IndexOutOfRange on the async path");

        // int.MaxValue is the BE numeric ceiling — must not overflow the loop bound.
        var maxK = await KmerAnalyzer.CountKmersAsync("ACGT", int.MaxValue);
        maxK.Should().BeEmpty("k = int.MaxValue is the numeric ceiling; no window fits, so the awaited result is empty");
    }

    #endregion

    #region BE — Cancellation: pre-signaled token

    /// <summary>
    /// BE / KEY (INV-03): a token already signaled BEFORE the call cancels the work that
    /// has not yet started — Task.Run observes the canceled token and the task transitions
    /// to Canceled, so awaiting it throws OperationCanceledException (its concrete wrapped
    /// form is TaskCanceledException, a subclass; Asynchronous_K-mer_Counting.md §2.2,
    /// §6.1 "Token signaled before call", [3]). We pin that the await THROWS the documented
    /// exception — cancellation is surfaced, NOT swallowed into a partial-garbage result —
    /// and does not hang. [CancelAfter] guards against a deadlock if cancellation were
    /// silently swallowed and the task instead ran to completion or stalled.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void CountKmersAsync_TokenSignaledBeforeCall_AwaitThrowsOperationCanceled()
    {
        using var cts = new CancellationTokenSource();
        cts.Cancel(); // signaled before the operation starts.

        // A modestly large input so, were cancellation swallowed, the work would actually run.
        string seq = RandomDna(50_000, seed: 156_001);

        var act = async () => await KmerAnalyzer.CountKmersAsync(seq, 6, cts.Token);

        act.Should().ThrowAsync<OperationCanceledException>(
                "INV-03: a token signaled before the call cancels the not-yet-started work; awaiting the canceled task throws")
            .GetAwaiter().GetResult();
    }

    #endregion

    #region BE — Cancellation: token signaled during the run

    /// <summary>
    /// BE / KEY (INV-03): a token signaled WHILE the scan is running must be observed at
    /// the next ThrowIfCancellationRequested checkpoint (every CheckInterval = 1000
    /// windows; KmerAnalyzer.cs lines 71–78), so the awaited task throws
    /// OperationCanceledException — it does NOT hang and does NOT return a
    /// partial-garbage dictionary (Asynchronous_K-mer_Counting.md §6.1 "Token signaled
    /// during run", [2]). We drive cancellation deterministically via a progress reporter:
    /// the FIRST progress callback (reported at the i = 0 checkpoint, before the bulk of
    /// the scan) signals the token, guaranteeing a still-running cancellation without a
    /// fragile sleep/race. A huge input ensures many checkpoints remain after the signal.
    /// [CancelAfter] guards against a deadlock if the mid-run cancellation were swallowed.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void CountKmersAsync_TokenSignaledDuringRun_AwaitThrowsOperationCanceled()
    {
        using var cts = new CancellationTokenSource();

        // Cancel from the FIRST progress report — that fires at the i = 0 checkpoint, so the
        // very next checkpoint (i = CheckInterval) observes the signal mid-scan.
        // NB: use a SYNCHRONOUS reporter, not Progress<double>. Progress<T> marshals the
        // callback asynchronously (SynchronizationContext.Post / thread pool), so under load
        // cts.Cancel() could run AFTER this test returns and `using` disposes cts — throwing
        // ObjectDisposedException on a background thread and crashing the test host. A
        // synchronous reporter cancels inline on the scanning thread, during the run (exactly
        // the invariant under test), and always completes before cts is disposed.
        var cancelOnFirstProgress = new SynchronousProgress<double>(_ => cts.Cancel());

        // Large enough that thousands of checkpoints remain after the i = 0 report.
        string seq = RandomDna(2_000_000, seed: 156_002);

        var act = async () => await KmerAnalyzer.CountKmersAsync(seq, 6, cts.Token, cancelOnFirstProgress);

        act.Should().ThrowAsync<OperationCanceledException>(
                "INV-03: a token signaled mid-run is observed at the next ThrowIfCancellationRequested checkpoint; awaiting throws, it does not hang or return a partial result")
            .GetAwaiter().GetResult();
    }

    #endregion

    #region BE — Boundary: k ≤ 0 surfaces through the awaited task

    /// <summary>
    /// BE: k ≤ 0 is the degenerate floor and a meaningless k-mer length — there is no
    /// "length-0 substring" to tally. The underlying synchronous overload rejects k ≤ 0 on
    /// non-empty input with ArgumentOutOfRangeException(nameof(k)) (KmerAnalyzer.cs lines
    /// 61–62), and that rejection is surfaced UNCHANGED through the awaited Task (INV-04;
    /// Asynchronous_K-mer_Counting.md §6.1 "k ≤ 0 → ArgumentOutOfRangeException (via
    /// awaited task)"). We pin that both k = 0 and a negative k throw on await and carry
    /// the documented "k" parameter name — the rejection boundary is exactly k ≤ 0, so a
    /// non-positive length can never slip into the window loop via the async surface.
    /// [CancelAfter] guards the awaits against a deadlock.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void CountKmersAsync_NonPositiveK_AwaitThrowsArgumentOutOfRange()
    {
        var zero = async () => await KmerAnalyzer.CountKmersAsync("ACGTACGT", 0);
        var negative = async () => await KmerAnalyzer.CountKmersAsync("ACGTACGT", -3);

        zero.Should().ThrowAsync<ArgumentOutOfRangeException>(
                "a 0-length k-mer is meaningless; the awaited task surfaces the underlying k <= 0 rejection on non-empty input")
            .Result.Which.ParamName.Should().Be("k");

        negative.Should().ThrowAsync<ArgumentOutOfRangeException>(
                "a negative k-mer length is nonsensical; the contract rejects all k <= 0 on non-empty input")
            .Result.Which.ParamName.Should().Be("k");
    }

    #endregion

    #region BE — Boundary: huge input (terminates + matches the synchronous reference)

    /// <summary>
    /// BE / KEY (INV-01 + INV-02): a large sequence must COMPLETE (terminate, no overflow,
    /// no hang) and the awaited result must equal the synchronous reference EXACTLY — the
    /// async surface runs the identical synchronous algorithm via Task.Run, so the count is
    /// byte-for-byte identical (Asynchronous_K-mer_Counting.md §1, §2.4 INV-01, §5.2). We
    /// scan several k values on a 100k-base sequence, asserting per k that (a) the result
    /// is well-formed (every key length k, counts sum to L − k + 1, INV-02), and (b) the
    /// awaited async map equals the synchronous CountKmers map key-for-key and count-for-
    /// count. [CancelAfter] is the load-bearing deadlock guard for the huge-input awaits.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public async Task CountKmersAsync_HugeInput_TerminatesAndMatchesSynchronousReference()
    {
        const int length = 100_000;
        string seq = RandomDna(length, seed: 156_003);

        foreach (int k in new[] { 1, 2, 3, 5, 8, 13, 21 })
        {
            var asyncResult = await KmerAnalyzer.CountKmersAsync(seq, k);
            var sync = KmerAnalyzer.CountKmers(seq, k); // synchronous reference oracle

            AssertWellFormed(asyncResult, length, k);
            asyncResult.Should().BeEquivalentTo(sync,
                $"INV-01: at k = {k} the awaited async result equals the synchronous reference key-for-key and count-for-count");
        }
    }

    #endregion

    #region Positive sanity — cross-surface equality and the documented cancellation throw

    /// <summary>
    /// Positive sanity (INV-01): the §7.1 worked example — `await CountKmersAsync("ATGG", 3)`
    /// returns {ATG:1, TGG:1} (Wikipedia: "ATGG has two 3-mers: ATG and TGG") — equal to
    /// the synchronous CountKmers. This pins the load-bearing cross-surface equality on a
    /// known sequence so the boundary/cancellation hardening never comes at the cost of the
    /// core async count silently breaking. [CancelAfter] guards the await against a hang.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public async Task CountKmersAsync_KnownSequence_EqualsSynchronousReference()
    {
        var asyncResult = await KmerAnalyzer.CountKmersAsync("ATGG", 3);

        asyncResult.Should().HaveCount(2, "'ATGG' has exactly two 3-mers: ATG and TGG");
        asyncResult.Should().ContainKey("ATG").WhoseValue.Should().Be(1, "'ATG' starts at position 0");
        asyncResult.Should().ContainKey("TGG").WhoseValue.Should().Be(1, "'TGG' starts at position 1");
        asyncResult.Values.Sum().Should().Be("ATGG".Length - 3 + 1, "INV-02: Σ counts = L − k + 1 = 2");

        asyncResult.Should().BeEquivalentTo(KmerAnalyzer.CountKmers("ATGG", 3),
            "INV-01: the awaited async result equals the synchronous reference");
    }

    /// <summary>
    /// Positive sanity (INV-01 + INV-03): a NON-cancelled run on a known sequence completes
    /// and matches the synchronous reference (the happy path), while a pre-cancelled token
    /// on the SAME inputs throws the documented OperationCanceledException — pinning both
    /// halves of the cancellation contract side by side so a swallowed cancellation (which
    /// would let the cancelled call also "succeed") cannot pass. [CancelAfter] guards both
    /// awaits against a hang.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public async Task CountKmersAsync_CompletesWithoutToken_ButThrowsOnPreCancelledToken()
    {
        const string seq = "ACGTACGTACGT";
        const int k = 4;

        // Happy path: no cancellation → completes and matches the synchronous reference.
        var ok = await KmerAnalyzer.CountKmersAsync(seq, k, CancellationToken.None);
        ok.Should().BeEquivalentTo(KmerAnalyzer.CountKmers(seq, k),
            "INV-01: an uncancelled async run equals the synchronous reference");

        // Same inputs, pre-cancelled token → the documented OperationCanceledException.
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        var act = async () => await KmerAnalyzer.CountKmersAsync(seq, k, cts.Token);
        act.Should().ThrowAsync<OperationCanceledException>(
                "INV-03: a pre-cancelled token throws on await — cancellation is surfaced, not swallowed into a result")
            .GetAwaiter().GetResult();
    }

    #endregion

    #endregion
}
