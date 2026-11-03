namespace SuffixTree.Esa;

/// <summary>Size limits of the 32-bit index path (ADR-002 section 2 item 7).</summary>
internal static class EsaLimits
{
    /// <summary>
    /// Guard for a text of <paramref name="textLength"/> symbols: the suffix array has <c>n + 1</c> entries and must fit one
    /// managed array, so <c>n + 1 &lt;= Array.MaxLength</c>. Larger inputs belong to the 64-bit path (stage S2).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">The length is negative.</exception>
    /// <exception cref="NotSupportedException">The suffix array would not fit the 32-bit path.</exception>
    internal static void EnsureTextLength(long textLength)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(textLength);
        // n + 1 > Array.MaxLength, written without the addition so that long.MaxValue cannot overflow
        if (textLength >= Array.MaxLength)
            throw new NotSupportedException(
                $"A text of {textLength} symbols needs {textLength + 1} suffix-array entries, more than the 32-bit path supports ({Array.MaxLength}); use the 64-bit path.");
    }

    /// <summary>The same guard for a symbol reader; only <see cref="ISymbolReader.Length"/> is read, never a symbol.</summary>
    internal static void EnsureReaderFits<TReader>(TReader text)
        where TReader : ISymbolReader, allows ref struct
        => EnsureTextLength(text.Length);
}
