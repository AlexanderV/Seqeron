using SuffixTree.Esa;

namespace SuffixTree.EmProto.Harness;

/// <summary>Thrown when the live native blocks of all arenas that share an <see cref="ArenaBudget"/> would exceed its cap (ADR-004 section 5, RAM-cap enforcement from inside).</summary>
internal sealed class ArenaBudgetExceededException(string message) : InvalidOperationException(message);

/// <summary>
/// The shared counter of the in-process RAM cap (ADR-004 section 4.4: "add an ArenaBudget counter that throws when the sum of live blocks exceeds the cap").
/// Thread-safe; every <see cref="BudgetedArena"/> of the builder reserves its blocks here. The peak is kept for the report.
/// </summary>
internal sealed class ArenaBudget(long capBytes)
{
    private long _live;
    private long _peak;

    /// <summary>The cap in bytes (<see cref="long.MaxValue"/> = unlimited).</summary>
    internal long CapBytes { get; } = capBytes;

    internal long LiveBytes => Interlocked.Read(ref _live);

    internal long PeakBytes => Interlocked.Read(ref _peak);

    /// <summary>Reserves <paramref name="bytes"/>; throws (and reserves nothing) when the cap would be exceeded.</summary>
    internal void Reserve(long bytes, string what)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(bytes);
        long now = Interlocked.Add(ref _live, bytes);
        if (now > CapBytes)
        {
            Interlocked.Add(ref _live, -bytes);
            throw new ArenaBudgetExceededException($"arena budget {CapBytes} bytes exceeded by '{what}' ({bytes} bytes; {now - bytes} live)");
        }
        Atomic.RaiseTo(ref _peak, now);
    }

    internal void Release(long bytes) => Interlocked.Add(ref _live, -bytes);
}

/// <summary>
/// A <see cref="NativeArena"/> whose blocks are charged to an <see cref="ArenaBudget"/>. Same contract as the arena (<c>Mark</c> / <c>ReleaseTo</c> release newest first,
/// <c>Dispose</c> frees all and returns the bytes to the budget); not thread-safe, one instance per task like the arena it wraps.
/// </summary>
internal sealed class BudgetedArena : IDisposable
{
    private readonly NativeArena _arena;
    private readonly ArenaBudget _budget;
    private readonly List<long> _sizes = [];
    private bool _disposed;

    internal BudgetedArena(ArenaBudget budget, ArenaOptions? options = null)
    {
        _budget = budget;
        _arena = new NativeArena(options ?? ArenaOptions.Production);
    }

    /// <summary>Payload bytes currently held by this arena.</summary>
    internal long LiveBytes => _sizes.Sum();

    /// <summary>Allocates <paramref name="length"/> elements, charged to the budget first.</summary>
    internal Buf<T> Allocate<T>(int length, bool zero, string name)
        where T : unmanaged
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        long bytes;
        unsafe
        {
            bytes = (long)length * sizeof(T);
        }
        _budget.Reserve(bytes, name);
        try
        {
            Buf<T> block = _arena.Allocate<T>(length, zero, name);
            _sizes.Add(bytes);
            return block;
        }
        catch
        {
            _budget.Release(bytes);
            throw;
        }
    }

    internal int Mark() => _arena.Mark();

    internal void ReleaseTo(int mark)
    {
        long freed = 0;
        for (int i = _sizes.Count - 1; i >= mark; i--)
        {
            freed += _sizes[i];
            _sizes.RemoveAt(i);
        }
        _arena.ReleaseTo(mark);
        _budget.Release(freed);
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        long all = _sizes.Sum();
        _sizes.Clear();
        _arena.Dispose();
        _budget.Release(all);
    }
}
