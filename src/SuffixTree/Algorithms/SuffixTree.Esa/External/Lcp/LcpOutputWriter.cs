#pragma warning disable S2325 // stub: the members become instance members of the implementation
namespace SuffixTree.Esa.External;

/// <summary>
/// The write side of B3: converts dense LCP values, in strictly ascending rank order, into the base array of the chosen width (little endian, escape value in the narrow kinds) and the overflow list (pairs
/// <c>(rank, value)</c> of two 4-byte fields, ADR-003 section 4 kind 5), both through <see cref="MappedWindowWriter"/> ranges (stream into memory-mapped files). The lengths of both ranges are known before the first
/// entry is written because the width was chosen from the counts of B2. Disposing disposes both writers. Not thread-safe.
/// </summary>
internal sealed class LcpOutputWriter : IDisposable
{
    /// <param name="arena">Owner of the two conversion buffers.</param>
    /// <param name="kind">The width chosen from the counts of B2.</param>
    /// <param name="n">Text length; the base range has <c>n + 1</c> entries of the width of <paramref name="kind"/>.</param>
    /// <param name="overflowCount">Entries of the overflow list (from the counts); the overflow range has <c>8 * overflowCount</c> bytes.</param>
    /// <param name="baseWriter">The range of the base array.</param>
    /// <param name="overflowWriter">The range of the overflow list.</param>
    /// <param name="sliceEntries">Entries converted per step.</param>
    internal LcpOutputWriter(BudgetedArena arena, LcpKind kind, long n, long overflowCount, MappedWindowWriter baseWriter, MappedWindowWriter overflowWriter, int sliceEntries) => throw new NotImplementedException();

    /// <summary>Appends the LCP values of ranks <c>firstRank ..</c>; the ranks must continue where the previous call ended.</summary>
    internal void Write(long firstRank, ReadOnlySpan<uint> values) => throw new NotImplementedException();

    /// <summary>Checks that every entry and exactly the announced number of escapes were written.</summary>
    internal void Complete(long n) => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
