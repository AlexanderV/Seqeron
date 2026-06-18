// KMER-ASYNC-001 — Asynchronous K-mer Counting
// Evidence: docs/Evidence/KMER-ASYNC-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-ASYNC-001.md
// Source: Wikipedia — K-mer (https://en.wikipedia.org/wiki/K-mer);
//         Microsoft Learn — Task Cancellation (.NET); Task.Run Method.

using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for KMER-ASYNC-001: Asynchronous K-mer Counting.
///
/// Canonical method: KmerAnalyzer.CountKmersAsync (thread-pool wrapper over the
/// synchronous CountKmers; result MUST be identical to the synchronous reference).
/// Evidence: Wikipedia K-mer (L − k + 1, GTAGAGCTGT table, ATGG example);
/// Microsoft Learn Task Cancellation / Task.Run (cooperative cancellation contract).
/// All async tests await completion — no race-dependent assertions.
/// </summary>
[TestFixture]
public class KmerAnalyzer_CountKmersAsync_Tests
{
    // Wikipedia K-mer worked sequence (L = 10).
    private const string Gtagagctgt = "GTAGAGCTGT";

    #region Determinism — async result equals synchronous reference

    /// <summary>
    /// M1: Async result is identical to the synchronous reference (GTAGAGCTGT, k=3).
    /// Evidence: Wikipedia table — k=3 yields 8 total k-mers; the 8 windows
    /// GTA,TAG,AGA,GAG,AGC,GCT,CTG,TGT are all distinct (8 distinct, hand-derived); INV-1.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_Gtagagctgt_K3_EqualsSynchronousReference()
    {
        var expected = KmerAnalyzer.CountKmers(Gtagagctgt, 3);

        var actual = await KmerAnalyzer.CountKmersAsync(Gtagagctgt, 3);

        Assert.Multiple(() =>
        {
            Assert.That(actual, Is.EqualTo(expected),
                "Async k-mer counts must equal the synchronous reference (INV-1).");
            Assert.That(actual.Values.Sum(), Is.EqualTo(8),
                "GTAGAGCTGT (L=10), k=3 has L−k+1 = 8 total k-mers (Wikipedia formula).");
            // k=3 distinct k-mers (derived from GTAGAGCTGT): GTA,TAG,AGA,GAG,AGC,GCT,CTG,TGT — all 8 distinct.
            Assert.That(actual.Count, Is.EqualTo(8),
                "GTAGAGCTGT, k=3 has 8 distinct 3-mers (none repeat; derived from the sequence).");
        });
    }

    /// <summary>
    /// M2: Determinism across multiple k (GTAGAGCTGT, k=2 and k=4).
    /// Evidence: Wikipedia L−k+1 formula — k=2 → 9 total, k=4 → 7 total; INV-1, INV-2.
    /// Distinct counts derived directly from the sequence: k=2 has GT and AG repeated
    /// (GT,TA,AG,GA,GC,CT,TG = 7 distinct); k=4 windows are all distinct (7).
    /// </summary>
    [Test]
    public async Task CountKmersAsync_Gtagagctgt_K2AndK4_EqualSynchronousReference()
    {
        var sync2 = KmerAnalyzer.CountKmers(Gtagagctgt, 2);
        var sync4 = KmerAnalyzer.CountKmers(Gtagagctgt, 4);

        var async2 = await KmerAnalyzer.CountKmersAsync(Gtagagctgt, 2);
        var async4 = await KmerAnalyzer.CountKmersAsync(Gtagagctgt, 4);

        Assert.Multiple(() =>
        {
            Assert.That(async2, Is.EqualTo(sync2), "k=2 async must equal sync (INV-1).");
            Assert.That(async2.Values.Sum(), Is.EqualTo(9), "k=2 total = L−k+1 = 9 (Wikipedia formula).");
            Assert.That(async2.Count, Is.EqualTo(7), "k=2 distinct = 7 (GT and AG each occur twice).");
            Assert.That(async2["GT"], Is.EqualTo(2), "GT occurs at positions 0 and 8.");
            Assert.That(async2["AG"], Is.EqualTo(2), "AG occurs at positions 2 and 4.");
            Assert.That(async4, Is.EqualTo(sync4), "k=4 async must equal sync (INV-1).");
            Assert.That(async4.Values.Sum(), Is.EqualTo(7), "k=4 total = L−k+1 = 7 (Wikipedia formula).");
            Assert.That(async4.Count, Is.EqualTo(7), "k=4 distinct = 7 (all 4-mer windows distinct).");
        });
    }

    /// <summary>
    /// M3: Wikipedia worked example — "ATGG has two 3-mers: ATG and TGG".
    /// Evidence: Wikipedia K-mer lead example.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_Atgg_K3_ReturnsAtgAndTgg()
    {
        var counts = await KmerAnalyzer.CountKmersAsync("ATGG", 3);

        Assert.Multiple(() =>
        {
            Assert.That(counts.Count, Is.EqualTo(2), "ATGG has two distinct 3-mers (Wikipedia).");
            Assert.That(counts["ATG"], Is.EqualTo(1), "ATG occurs once (Wikipedia).");
            Assert.That(counts["TGG"], Is.EqualTo(1), "TGG occurs once (Wikipedia).");
        });
    }

    /// <summary>
    /// M4: Total-count invariant Σcounts = L − k + 1 for k = 1..10 on GTAGAGCTGT.
    /// Evidence: Wikipedia "a sequence of length L will have L − k + 1 k-mers"; INV-2.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_Gtagagctgt_TotalCountInvariant_ForAllK()
    {
        int length = Gtagagctgt.Length; // 10

        for (int k = 1; k <= length; k++)
        {
            var counts = await KmerAnalyzer.CountKmersAsync(Gtagagctgt, k);
            int expectedTotal = length - k + 1; // 10,9,...,1

            Assert.That(counts.Values.Sum(), Is.EqualTo(expectedTotal),
                $"k={k}: total k-mers must equal L−k+1 = {expectedTotal} (Wikipedia, INV-2).");
        }
    }

    /// <summary>
    /// M4b: k=1 distinct count = 4 on GTAGAGCTGT (the four bases G,T,A,C).
    /// Evidence: Wikipedia GTAGAGCTGT table states k=1 → "4 distinct k-mers" — the only
    /// distinct figure the source gives for this sequence; carried through the async path.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_Gtagagctgt_K1_HasFourDistinctBases()
    {
        var counts = await KmerAnalyzer.CountKmersAsync(Gtagagctgt, 1);

        Assert.Multiple(() =>
        {
            Assert.That(counts.Count, Is.EqualTo(4),
                "k=1 distinct = 4 (bases G,T,A,C) — Wikipedia GTAGAGCTGT table.");
            Assert.That(counts.Values.Sum(), Is.EqualTo(10),
                "k=1 total = L−k+1 = 10 (Wikipedia formula, INV-2).");
        });
    }

    #endregion

    #region Edge Cases — empty, null, boundary, invalid k

    /// <summary>
    /// M5: Empty sequence returns an empty dictionary.
    /// Evidence: Wikipedia (L=0 ⇒ no k-mers); INV-4.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_EmptySequence_ReturnsEmptyDictionary()
    {
        var counts = await KmerAnalyzer.CountKmersAsync("", 4);
        Assert.That(counts, Is.Empty, "Empty sequence has no k-mers (Wikipedia, INV-4).");
    }

    /// <summary>
    /// M6: Null sequence returns an empty dictionary.
    /// Evidence: synchronous contract preserved through the wrapper; INV-4.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_NullSequence_ReturnsEmptyDictionary()
    {
        string? nullSequence = null;
        var counts = await KmerAnalyzer.CountKmersAsync(nullSequence!, 4);
        Assert.That(counts, Is.Empty, "Null sequence yields an empty result (INV-4).");
    }

    /// <summary>
    /// M7: k greater than sequence length returns an empty dictionary.
    /// Evidence: Wikipedia (L−k+1 ≤ 0); INV-4.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_KLargerThanSequence_ReturnsEmptyDictionary()
    {
        var counts = await KmerAnalyzer.CountKmersAsync("ACG", 4);
        Assert.That(counts, Is.Empty, "k > L means no k-mers exist (Wikipedia, INV-4).");
    }

    /// <summary>
    /// M8: k ≤ 0 surfaces ArgumentOutOfRangeException via the awaited task.
    /// Evidence: synchronous contract (KMER-COUNT-001) preserved through Task.Run; INV-4.
    /// </summary>
    [Test]
    public void CountKmersAsync_InvalidK_ThrowsArgumentOutOfRangeException()
    {
        Assert.Multiple(() =>
        {
            Assert.ThrowsAsync<ArgumentOutOfRangeException>(
                async () => await KmerAnalyzer.CountKmersAsync(Gtagagctgt, 0),
                "k=0 is invalid; the awaited task must throw (INV-4).");
            Assert.ThrowsAsync<ArgumentOutOfRangeException>(
                async () => await KmerAnalyzer.CountKmersAsync(Gtagagctgt, -1),
                "k=-1 is invalid; the awaited task must throw (INV-4).");
        });
    }

    #endregion

    #region Cancellation contract

    /// <summary>
    /// M9: A token already signaled before the call cancels the work; awaiting throws
    /// OperationCanceledException.
    /// Evidence: Task.Run "allows the work to be cancelled if it has not yet started";
    /// Task Cancellation model; INV-3. Deterministic (token pre-canceled, then await).
    /// </summary>
    [Test]
    public void CountKmersAsync_TokenSignaledBeforeCall_ThrowsOperationCanceled()
    {
        using var cts = new CancellationTokenSource();
        cts.Cancel(); // signaled before the operation starts

        // Awaiting a canceled task throws OperationCanceledException; the concrete type is
        // TaskCanceledException (its subclass), per Microsoft Learn — CatchAsync matches both.
        Assert.CatchAsync<OperationCanceledException>(
            async () => await KmerAnalyzer.CountKmersAsync(Gtagagctgt, 3, cts.Token),
            "A pre-signaled token must cancel the task; awaiting throws (INV-3).");
    }

    /// <summary>
    /// M10: Cancellation on a large input (above the internal check interval). The token
    /// is signaled before the call so the operation observes cancellation deterministically
    /// and the awaited task throws OperationCanceledException — no timing dependency.
    /// Evidence: ThrowIfCancellationRequested transitions task to Canceled; INV-3.
    /// </summary>
    [Test]
    public void CountKmersAsync_CancelledLargeInput_ThrowsOperationCanceled()
    {
        // 5000 nt: exceeds the 1000-window internal cancellation check interval.
        var large = string.Concat(Enumerable.Repeat("ACGT", 1250));
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        Assert.CatchAsync<OperationCanceledException>(
            async () => await KmerAnalyzer.CountKmersAsync(large, 4, cts.Token),
            "A signaled token on a large input must cancel; awaiting throws (INV-3).");
    }

    #endregion

    #region SHOULD / COULD — large-input determinism, case, default token, progress, span

    /// <summary>
    /// S1: Large-input determinism — async result equals the synchronous reference on a
    /// fixed 5000 nt sequence (exceeds the check interval).
    /// Evidence: INV-1; the async path targets long inputs.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_LargeInput_EqualsSynchronousReference()
    {
        var large = string.Concat(Enumerable.Repeat("ACGTACGGTA", 500)); // 5000 nt, fixed
        var expected = KmerAnalyzer.CountKmers(large, 6);

        var actual = await KmerAnalyzer.CountKmersAsync(large, 6);

        Assert.Multiple(() =>
        {
            Assert.That(actual, Is.EqualTo(expected),
                "Large-input async result must equal the synchronous reference (INV-1).");
            Assert.That(actual.Values.Sum(), Is.EqualTo(large.Length - 6 + 1),
                "Total k-mers must equal L−k+1 (Wikipedia, INV-2).");
        });
    }

    /// <summary>
    /// S2: Case-insensitive — lowercase input yields the same counts as uppercase.
    /// Evidence: input is uppercased; API consistency with the canonical method.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_LowercaseInput_EqualsUppercaseResult()
    {
        var upper = await KmerAnalyzer.CountKmersAsync(Gtagagctgt, 3);
        var lower = await KmerAnalyzer.CountKmersAsync(Gtagagctgt.ToLowerInvariant(), 3);

        Assert.That(lower, Is.EqualTo(upper),
            "Counting is case-insensitive; lowercase must match uppercase.");
    }

    /// <summary>
    /// S3: Default token path (CancellationToken.None) completes and equals sync.
    /// Evidence: INV-1; default-parameter path.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_DefaultToken_CompletesAndEqualsSync()
    {
        var expected = KmerAnalyzer.CountKmers(Gtagagctgt, 2);

        var actual = await KmerAnalyzer.CountKmersAsync(Gtagagctgt, 2); // no token arg

        Assert.That(actual, Is.EqualTo(expected),
            "Default-token async result must equal the synchronous reference (INV-1).");
    }

    /// <summary>
    /// C1: Progress reporting — the final reported value is 1.0 on completion.
    /// Evidence: synchronous method reports progress 1.0 at completion; preserved by wrapper.
    /// </summary>
    [Test]
    public async Task CountKmersAsync_WithProgress_ReportsCompletion()
    {
        var reports = new List<double>();
        var progress = new SynchronousProgress(reports.Add);

        await KmerAnalyzer.CountKmersAsync(Gtagagctgt, 3, CancellationToken.None, progress);

        Assert.Multiple(() =>
        {
            Assert.That(reports, Is.Not.Empty, "Progress must be reported at least once.");
            Assert.That(reports[^1], Is.EqualTo(1.0).Within(1e-10),
                "Final progress report must be 1.0 at completion.");
        });
    }

    /// <summary>
    /// C2: CountKmersSpan delegate smoke — equals CountKmers (deeply tested in KMER-COUNT-001).
    /// Evidence: API consistency; INV-1 (span variant of the same definition).
    /// </summary>
    [Test]
    public void CountKmersSpan_DelegatesToSpanCount_EqualsCountKmers()
    {
        var expected = KmerAnalyzer.CountKmers(Gtagagctgt, 3);

        var actual = KmerAnalyzer.CountKmersSpan(Gtagagctgt.AsSpan(), 3);

        Assert.That(actual, Is.EqualTo(expected),
            "Span variant must produce the same counts as CountKmers.");
    }

    #endregion

    /// <summary>
    /// Deterministic IProgress that invokes the callback synchronously on Report,
    /// so progress assertions are not subject to SynchronizationContext timing.
    /// </summary>
    private sealed class SynchronousProgress : IProgress<double>
    {
        private readonly Action<double> _handler;
        public SynchronousProgress(Action<double> handler) => _handler = handler;
        public void Report(double value) => _handler(value);
    }
}
