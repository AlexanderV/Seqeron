#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

namespace SuffixTree.Esa.External;

/// <summary>Thrown when the live native blocks of all arenas that share an <see cref="ArenaBudget"/> would exceed its cap.</summary>
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

/// <summary>The shared counter of the in-process RAM cap (ADR-004 section 4.4). Thread-safe.</summary>
internal sealed class ArenaBudget
{
    internal ArenaBudget(long capBytes) => throw new NotImplementedException();

    /// <summary>The cap in bytes (<see cref="long.MaxValue"/> = unlimited).</summary>
    internal long CapBytes => throw new NotImplementedException();

    internal long LiveBytes => throw new NotImplementedException();

    internal long PeakBytes => throw new NotImplementedException();

    /// <summary>Reserves <paramref name="bytes"/>; throws (and reserves nothing) when the cap would be exceeded.</summary>
    internal void Reserve(long bytes, string what) => throw new NotImplementedException();

    /// <summary>Gives <paramref name="bytes"/> back.</summary>
    internal void Release(long bytes) => throw new NotImplementedException();
}

/// <summary>A <see cref="NativeArena"/> whose blocks are charged to an <see cref="ArenaBudget"/>. Not thread-safe, one instance per task.</summary>
internal sealed unsafe class BudgetedArena : IDisposable
{
    internal BudgetedArena(ArenaBudget budget, ArenaOptions? options = null) => throw new NotImplementedException();

    internal ArenaBudget Budget => throw new NotImplementedException();

    /// <summary>Payload bytes currently held by this arena.</summary>
    internal long LiveBytes => throw new NotImplementedException();

    internal Buf<T> Allocate<T>(int length, bool zero, string name)
        where T : unmanaged
        => throw new NotImplementedException();

    internal int Mark() => throw new NotImplementedException();

    internal void ReleaseTo(int mark) => throw new NotImplementedException();

    /// <summary>Throws when a guard band of a guarded arena was modified.</summary>
    internal void VerifyCanaries() => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
