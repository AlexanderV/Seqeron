namespace SuffixTree.Esa;

/// <summary>What the safety infrastructure detected (PROCESS.md amendment 2026-11-03, review MAJ-1).</summary>
public enum MemoryViolationKind
{
    /// <summary>An access through a checked <see cref="Buf{T}"/> outside its logical length (only in a build with <c>ESA_BOUNDS</c>).</summary>
    OutOfBounds,

    /// <summary>A guard band around a <see cref="NativeArena"/> block was modified: some code wrote outside the logical length.</summary>
    CanaryOverwritten,
}

/// <summary>An out-of-bounds access (checked build) or an overwritten canary (<see cref="NativeArena.VerifyCanaries"/>). Deliberately not an <see cref="InvalidOperationException"/>: no handler for ordinary input errors may swallow it.</summary>
public sealed class EsaMemoryViolationException : Exception
{
    /// <summary>Creates the exception without a message.</summary>
    public EsaMemoryViolationException()
    {
    }

    /// <summary>Creates the exception with a message.</summary>
    /// <param name="message">What was detected.</param>
    public EsaMemoryViolationException(string message)
        : base(message)
    {
    }

    /// <summary>Creates the exception with a message and a cause.</summary>
    /// <param name="message">What was detected.</param>
    /// <param name="innerException">The cause.</param>
    public EsaMemoryViolationException(string message, Exception innerException)
        : base(message, innerException)
    {
    }

    /// <summary>Creates the exception for a class of violation.</summary>
    /// <param name="kind">The class of the violation.</param>
    /// <param name="message">What was detected.</param>
    public EsaMemoryViolationException(MemoryViolationKind kind, string message)
        : base(message)
    {
        Kind = kind;
    }

    /// <summary>The class of the violation.</summary>
    public MemoryViolationKind Kind { get; }
}
