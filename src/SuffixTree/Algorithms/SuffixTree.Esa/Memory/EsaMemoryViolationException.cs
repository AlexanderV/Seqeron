namespace SuffixTree.Esa;

/// <summary>What the safety infrastructure detected (PROCESS.md amendment 2026-11-03, review MAJ-1).</summary>
internal enum MemoryViolationKind
{
    /// <summary>An access through a checked <see cref="Buf{T}"/> outside its logical length (only in a build with <c>ESA_BOUNDS</c>).</summary>
    OutOfBounds,

    /// <summary>A guard band around a <see cref="NativeArena"/> block was modified: some code wrote outside the logical length.</summary>
    CanaryOverwritten,
}

/// <summary>An out-of-bounds access (checked build) or an overwritten canary (<see cref="NativeArena.VerifyCanaries"/>).</summary>
internal sealed class EsaMemoryViolationException : InvalidOperationException
{
    public EsaMemoryViolationException()
    {
    }

    public EsaMemoryViolationException(string message)
        : base(message)
    {
    }

    public EsaMemoryViolationException(string message, Exception innerException)
        : base(message, innerException)
    {
    }

    public EsaMemoryViolationException(MemoryViolationKind kind, string message)
        : base(message)
    {
        Kind = kind;
    }

    /// <summary>The class of the violation.</summary>
    public MemoryViolationKind Kind { get; }
}
