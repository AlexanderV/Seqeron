namespace SuffixTree.Esa.External;

/// <summary>Thrown when the live native blocks of all arenas that share an <see cref="ArenaBudget"/> would exceed its cap (ADR-004 section 5, RAM-cap enforcement from inside).</summary>
internal sealed class ArenaBudgetExceededException : InvalidOperationException
{
    internal ArenaBudgetExceededException()
    {
    }

    internal ArenaBudgetExceededException(string message)
        : base(message)
    {
    }

    internal ArenaBudgetExceededException(string message, Exception innerException)
        : base(message, innerException)
    {
    }
}

/// <summary>
/// The shared counter of the in-process RAM cap (ADR-004 section 4.4: "an ArenaBudget counter that throws when the sum of live blocks exceeds the cap"). Thread-safe; every <see cref="BudgetedArena"/> of the builder
/// (one per task) reserves its blocks here. The peak is kept for the report. It lives next to, not inside, <see cref="NativeArena"/>: the budget is a policy of the builder that spans many arenas, and the in-memory
/// kernels keep allocating through the arena exactly as they were verified.
/// </summary>
internal sealed class ArenaBudget
{
    private long _live;
    private long _peak;

    /// <param name="capBytes">The cap in bytes (<see cref="long.MaxValue"/> = unlimited).</param>
    internal ArenaBudget(long capBytes)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(capBytes);
        CapBytes = capBytes;
    }

    /// <summary>The cap in bytes.</summary>
    internal long CapBytes { get; }

    /// <summary>Bytes reserved now.</summary>
    internal long LiveBytes => Interlocked.Read(ref _live);

    /// <summary>The largest <see cref="LiveBytes"/> that was ever reached.</summary>
    internal long PeakBytes => Interlocked.Read(ref _peak);

    /// <summary>Reserves <paramref name="bytes"/>; throws <see cref="ArenaBudgetExceededException"/> (and reserves nothing) when the cap would be exceeded. Exactly the cap is allowed.</summary>
    internal void Reserve(long bytes, string what)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(bytes);
        long seen = Interlocked.Read(ref _live);
        while (true)
        {
            // compared as "bytes > cap - live" so that a cap of long.MaxValue cannot wrap around
            if (bytes > CapBytes - seen) throw new ArenaBudgetExceededException($"arena budget {CapBytes} bytes exceeded by '{what}' ({bytes} bytes; {seen} live)");
            long now = seen + bytes;
            long previous = Interlocked.CompareExchange(ref _live, now, seen);
            if (previous == seen)
            {
                Atomic.RaiseTo(ref _peak, now);
                return;
            }
            seen = previous;
        }
    }

    /// <summary>Gives <paramref name="bytes"/> back; more than is live is a bug of the caller and throws.</summary>
    internal void Release(long bytes)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(bytes);
        long seen = Interlocked.Read(ref _live);
        while (true)
        {
            if (bytes > seen) throw new InvalidOperationException($"{bytes} bytes released but only {seen} are live");
            long previous = Interlocked.CompareExchange(ref _live, seen - bytes, seen);
            if (previous == seen) return;
            seen = previous;
        }
    }
}

/// <summary>
/// A <see cref="NativeArena"/> whose blocks are charged to an <see cref="ArenaBudget"/> before they are allocated, so a refused block is never allocated. Same contract as the arena (<c>Mark</c> / <c>ReleaseTo</c>
/// release newest first, <c>Dispose</c> frees all and returns the bytes to the budget); not thread-safe, one instance per task like the arena it wraps. The budget counts payload bytes, not guard bands.
/// </summary>
internal sealed class BudgetedArena : IDisposable
{
    private readonly NativeArena _arena;
    private readonly List<long> _sizes = [];
    private long _live;
    private bool _disposed;

    /// <param name="budget">The shared cap.</param>
    /// <param name="options">Arena options (guard bands in the tests).</param>
    internal BudgetedArena(ArenaBudget budget, ArenaOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(budget);
        Budget = budget;
        _arena = new NativeArena(options ?? ArenaOptions.Production);
    }

    internal ArenaBudget Budget { get; }

    /// <summary>Payload bytes currently held by this arena.</summary>
    internal long LiveBytes => _live;

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
        Budget.Reserve(bytes, name); // a negative length is refused here
        bool allocated = false;
        try
        {
            Buf<T> block = _arena.Allocate<T>(length, zero, name);
            allocated = true;
            _sizes.Add(bytes);
            _live += bytes;
            return block;
        }
        finally
        {
            if (!allocated) Budget.Release(bytes);
        }
    }

    internal int Mark()
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        return _arena.Mark();
    }

    /// <summary>Frees every block allocated since <paramref name="mark"/>, newest first, and gives their bytes back (also when the guard bands of a block were damaged and the arena throws).</summary>
    internal void ReleaseTo(int mark)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        ArgumentOutOfRangeException.ThrowIfNegative(mark);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(mark, _sizes.Count);
        long freed = 0;
        for (int i = _sizes.Count - 1; i >= mark; i--)
        {
            freed += _sizes[i];
            _sizes.RemoveAt(i);
        }
        _live -= freed;
        try
        {
            _arena.ReleaseTo(mark);
        }
        finally
        {
            Budget.Release(freed);
        }
    }

    /// <summary>Throws <see cref="EsaMemoryViolationException"/> when a guard band of a guarded arena was modified; does nothing otherwise.</summary>
    internal void VerifyCanaries() => _arena.VerifyCanaries();

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        long all = _live;
        _live = 0;
        _sizes.Clear();
        try
        {
            _arena.Dispose();
        }
        finally
        {
            Budget.Release(all);
        }
    }
}
