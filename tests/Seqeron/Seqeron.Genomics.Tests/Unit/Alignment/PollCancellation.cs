using System.Diagnostics;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// The deterministic cancellation seam of the B09 tests (audit A5-4): instead of cancelling a running call after a wall-clock delay (which lands in another stage under CPU load), the call runs ONCE with a token that
/// is cancelled AT the n-th poll of <see cref="AlignmentCancellation"/> whose stack contains the named method (the probe of <see cref="AlignmentCancellation.OnPoll"/> runs at the start of every poll of the calling
/// flow, worker threads included). The call must then stop with <see cref="OperationCanceledException"/>: "polled inside X" (the cancelling poll is in the named method), "mid-run" (n &gt;= 2 means rows /
/// iterations were already done) and "no partial result" are asserted without a sleep, a thread or a time bound.
/// </summary>
internal static class PollCancellation
{
    /// <summary>
    /// Runs <paramref name="work"/> with the token cancelled at the <paramref name="nth"/> poll whose stack contains <paramref name="pollInside"/> (any poll when null). Returns the stack of the cancelling poll.
    /// </summary>
    public static string CancelInside(Func<CancellationToken, object?> work, string? pollInside = null, int nth = 3)
    {
        using var cts = new CancellationTokenSource();
        int matching = 0, total = 0;
        string? firing = null;
        OperationCanceledException? caught = null;
        using (AlignmentCancellation.OnPoll(() =>
        {
            Interlocked.Increment(ref total);
            if (cts.IsCancellationRequested)
                return;
            if (pollInside is not null && !new StackTrace(1, false).ToString().Contains(pollInside, StringComparison.Ordinal))
                return;
            if (Interlocked.Increment(ref matching) == nth)
            {
                firing = new StackTrace(1, false).ToString();        // the poll that requests the cancellation: the call must observe it
                cts.Cancel();
            }
        }))
        {
            try
            {
                work(cts.Token);
            }
            catch (OperationCanceledException ex)
            {
                caught = ex;
            }
        }

        Assert.That(caught, Is.Not.Null, pollInside is null ? $"the call finished although the token was cancelled at poll {nth} ({total} polls)" : $"the call finished: only {matching} polls inside {pollInside} ran ({total} polls in all), the token is cancelled at the {nth}th");
        // the named method is the one whose poll requested the cancellation (a Parallel.For worker's exception is re-thrown by the loop, whose stack does not name the polling method)
        return firing ?? string.Empty;
    }

    /// <summary>The call with an already cancelled token: it must throw <see cref="OperationCanceledException"/> before doing the work (no poll of the named method is needed).</summary>
    public static void PreCancelled(Func<CancellationToken, object?> work) =>
        Assert.Catch<OperationCanceledException>(() => work(new CancellationToken(true)));
}
