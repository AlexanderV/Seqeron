using System.Diagnostics;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// The deterministic cancellation seam of the B09 tests (audit A5-4): instead of cancelling a running call after a wall-clock delay (which lands in another stage under CPU load), the call is run once to
/// record every poll of <see cref="AlignmentCancellation"/> (with its stack), and again with a token that is cancelled AT a chosen poll - the one inside the named method, in the middle of its polls. The
/// call must then stop with <see cref="OperationCanceledException"/> raised by that very poll (the stack names the method), having polled exactly that many times: "polled inside X", "mid-run" and
/// "no partial result" are asserted without a sleep, a thread or a time bound.
/// </summary>
internal static class PollCancellation
{
    /// <summary>
    /// Runs <paramref name="work"/> uncancelled (a recording run), then again with the token cancelled at a poll whose stack contains <paramref name="pollInside"/> (the middle one of those; any poll when
    /// <paramref name="pollInside"/> is null). Returns the stack of the exception.
    /// </summary>
    public static string CancelInside(Func<CancellationToken, object?> work, string? pollInside = null, int minimumPolls = 2)
    {
        var stacks = new List<string>();
        using (AlignmentCancellation.OnPoll(() => stacks.Add(new StackTrace(1, false).ToString())))
            work(CancellationToken.None);

        var candidates = Enumerable.Range(0, stacks.Count).Where(i => pollInside is null || stacks[i].Contains(pollInside, StringComparison.Ordinal)).ToList();
        Assert.That(candidates, Has.Count.GreaterThanOrEqualTo(minimumPolls),
            pollInside is null ? $"the call polled {stacks.Count} times (a mid-run cancellation needs at least {minimumPolls})" : $"the poll inside {pollInside} ran {candidates.Count} times in {stacks.Count} polls (need {minimumPolls}); first stack:\n{(stacks.Count > 0 ? stacks[0] : "none")}");
        int target = candidates[candidates.Count / 2] + 1;      // 1-based number of the poll that fires

        using var cts = new CancellationTokenSource();
        int polls = 0;
        OperationCanceledException? caught = null;
        using (AlignmentCancellation.OnPoll(() => { if (Interlocked.Increment(ref polls) == target) cts.Cancel(); }))
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

        Assert.That(caught, Is.Not.Null, $"the call finished although the token was cancelled at poll {target} of {stacks.Count}");
        Assert.That(Volatile.Read(ref polls), Is.EqualTo(target), "the call stopped at the poll that observed the cancellation: no later poll ran");
        if (pollInside is not null)
            Assert.That(caught!.StackTrace, Does.Contain(pollInside), $"the cancellation must be observed by the poll inside {pollInside}");
        return caught!.StackTrace ?? string.Empty;
    }
}
