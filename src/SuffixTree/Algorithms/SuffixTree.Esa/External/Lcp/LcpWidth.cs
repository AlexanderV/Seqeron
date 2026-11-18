namespace SuffixTree.Esa.External;

/// <summary>The width rule of ADR-002 section 4 on counts that are known before the first LCP entry is written.</summary>
internal static class LcpWidth
{
    /// <summary>
    /// <c>U8</c> + overflow iff <c>16 * escapes8 &lt;= n</c>, else <c>U16</c> + overflow iff <c>8 * escapes16 &lt;= n</c>, else <c>U32</c>: the rule of <see cref="CompactLcp.SelectKind"/> on 64-bit counts
    /// (that method takes an <c>int</c> n; the tests pin both to each other).
    /// </summary>
    internal static LcpKind Select(long n, long escapes8, long escapes16) => throw new NotImplementedException();

    /// <summary>Bytes per entry of the base array.</summary>
    internal static int Bytes(LcpKind kind) => throw new NotImplementedException();

    /// <summary>The escape value of the base array (a value at or above it is stored in the overflow list); <see cref="uint.MaxValue"/> for U32, which has none.</summary>
    internal static uint Escape(LcpKind kind) => throw new NotImplementedException();

    /// <summary>Entries of the overflow list for the counts of values at or above 255 and 65,535.</summary>
    internal static long OverflowCount(LcpKind kind, long escapes8, long escapes16) => throw new NotImplementedException();
}
