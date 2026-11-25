using System.Runtime.ExceptionServices;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The ambient cancellation scope of the aligners (audit A2-40). Every long-running dynamic program polls
/// <see cref="AlignmentCancellation.ThrowIfRequested"/> at its row / iteration boundaries; the public entry points that take a
/// <see cref="CancellationToken"/> (the <c>cancellationToken</c> overloads, <c>*Options.CancellationToken</c>,
/// <see cref="SequenceAligner.WithCancellation{T}(CancellationToken, Func{T})"/>) enter a scope for the duration of the call.
/// The scope is an <see cref="AsyncLocal{T}"/>, so it also reaches the worker threads of <c>Parallel.For</c> (the execution context flows),
/// and nested calls inherit the token of the outermost cancellable scope. A cancelled call throws
/// <see cref="OperationCanceledException"/> and returns no partial result (the result objects are built after the last poll).
/// </summary>
public static class AlignmentCancellation
{
    private static readonly AsyncLocal<CancellationToken> Current = new();
    private static readonly AsyncLocal<Action?> PollProbe = new();

    /// <summary>The token of the innermost cancellable scope (<c>default</c> = never cancelled).</summary>
    public static CancellationToken Token => Current.Value;

    /// <summary>Throws <see cref="OperationCanceledException"/> when the scope's token is cancelled (a plain field read when no scope is active).</summary>
    public static void ThrowIfRequested()
    {
        PollProbe.Value?.Invoke();
        Current.Value.ThrowIfCancellationRequested();
    }

    /// <summary>Polls a local token (a <c>cancellationToken</c> parameter of an overload) and the ambient scope token.</summary>
    public static void Poll(CancellationToken local)
    {
        PollProbe.Value?.Invoke();
        local.ThrowIfCancellationRequested();
        Current.Value.ThrowIfCancellationRequested();
    }

    /// <summary>
    /// Test seam (audit A5-4): <paramref name="probe"/> runs at the start of EVERY poll of the calling flow (the execution context flows to the worker threads of <c>Parallel.For</c>, other tests are not affected),
    /// so a test cancels its source at the n-th poll deterministically instead of after a wall-clock delay. Dispose to remove it.
    /// </summary>
    public static IDisposable OnPoll(Action probe)
    {
        var previous = PollProbe.Value;
        PollProbe.Value = probe;
        return new ProbeScope(previous);
    }

    private sealed class ProbeScope : IDisposable
    {
        private readonly Action? _previous;
        public ProbeScope(Action? previous) => _previous = previous;
        public void Dispose() => PollProbe.Value = _previous;
    }

    /// <summary>Enters a scope for <paramref name="token"/> (a token that cannot be cancelled keeps the surrounding scope); dispose to leave it.</summary>
    public static Scope Enter(CancellationToken token)
    {
        token.ThrowIfCancellationRequested();
        var previous = Current.Value;
        if (token.CanBeCanceled)
            Current.Value = token;
        return new Scope(previous);
    }

    public readonly struct Scope : IDisposable
    {
        private readonly CancellationToken _previous;
        internal Scope(CancellationToken previous) => _previous = previous;
        public void Dispose() => Current.Value = _previous;
    }
}

public static partial class SequenceAligner
{
    /// <summary>
    /// Runs <paramref name="aligner"/> (any call of this class, e.g. <c>() =&gt; SequenceAligner.LocalAlignTopN(a, b, matrixScoring, 5, 10)</c>)
    /// so that it observes <paramref name="cancellationToken"/>: the dynamic programs poll it at their row / iteration boundaries and the call
    /// throws <see cref="OperationCanceledException"/> without a partial result. The one mechanism behind all cancellable overloads; use it for the
    /// overloads (matrix scoring, linear space, ...) that have no dedicated <c>cancellationToken</c> parameter.
    /// </summary>
    /// <exception cref="OperationCanceledException">The token was cancelled before or during the call.</exception>
    public static T WithCancellation<T>(CancellationToken cancellationToken, Func<T> aligner)
    {
        ArgumentNullException.ThrowIfNull(aligner);
        using var scope = AlignmentCancellation.Enter(cancellationToken);
        return aligner();
    }

    /// <summary>
    /// Opens a cancellation scope for the calling thread / async flow: every aligner called until the returned object is disposed observes
    /// <paramref name="cancellationToken"/> (see <see cref="WithCancellation{T}(CancellationToken, Func{T})"/>). Use in a <c>using</c> statement; used
    /// by the MCP tools to honour the request's token. A token that cannot be cancelled keeps the surrounding scope.
    /// </summary>
    /// <exception cref="OperationCanceledException">The token is already cancelled.</exception>
    public static IDisposable BeginCancellationScope(CancellationToken cancellationToken) => AlignmentCancellation.Enter(cancellationToken);

    /// <summary>Runs <paramref name="body"/> for every index in [<paramref name="fromInclusive"/>, <paramref name="toExclusive"/>) in parallel and
    /// rethrows the exception of the LOWEST failing index with its original type (never an <see cref="AggregateException"/>; the sequential
    /// path of the same method throws the same exception: <c>ExceptionDispatchInfo</c> keeps the stack). Honours the ambient token and
    /// <paramref name="maxDegreeOfParallelism"/> (-1 = the runtime's choice, 1 = sequential).</summary>
    public static void ParallelFor(int fromInclusive, int toExclusive, Action<int> body, int maxDegreeOfParallelism = -1)
    {
        if (maxDegreeOfParallelism == 0 || maxDegreeOfParallelism < -1)
            throw new ArgumentOutOfRangeException(nameof(maxDegreeOfParallelism), maxDegreeOfParallelism, "MaxDegreeOfParallelism must be -1 (unbounded) or at least 1.");
        var token = AlignmentCancellation.Token;
        token.ThrowIfCancellationRequested();
        Exception? first = null;
        int firstIndex = int.MaxValue;
        var gate = new object();
        Parallel.For(fromInclusive, toExclusive,
            new ParallelOptions { CancellationToken = token, MaxDegreeOfParallelism = maxDegreeOfParallelism },
            (i, state) =>
            {
                try
                {
                    body(i);
                }
                catch (Exception ex)
                {
                    lock (gate)
                    {
                        if (i < firstIndex) { firstIndex = i; first = ex; }
                    }

                    state.Stop();
                }
            });

        if (first is not null)
            ExceptionDispatchInfo.Capture(first).Throw();
    }
}
