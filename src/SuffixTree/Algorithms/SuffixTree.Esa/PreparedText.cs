namespace SuffixTree.Esa;

/// <summary>
/// A text after the order-preserving dense remap: symbol <c>k</c> is the k-th smallest distinct UTF-16 code unit of the original
/// text, so the suffix array of the remapped sequence equals the suffix array of the original string. The original code units are
/// recoverable through <see cref="Alphabet"/> (the stored text of the container keeps the original units, ADR-003).
/// </summary>
internal sealed class PreparedText
{
    private readonly byte[]? _bytes;
    private readonly ushort[]? _ushorts;

    internal PreparedText(char[] alphabet, byte[]? bytes, ushort[]? ushorts)
    {
        Alphabet = alphabet;
        _bytes = bytes;
        _ushorts = ushorts;
        Width = bytes is not null ? SymbolWidth.Byte : SymbolWidth.UShort;
        Length = (bytes?.Length ?? ushorts?.Length) ?? 0;
    }

    /// <summary>Number of symbols n.</summary>
    internal int Length { get; }

    /// <summary>Number of distinct code units (sigma); symbols are in <c>[0, sigma)</c>.</summary>
    internal int AlphabetSize => Alphabet.Length;

    /// <summary>Byte when sigma &lt;= 256, otherwise UShort.</summary>
    internal SymbolWidth Width { get; }

    /// <summary><c>Alphabet[symbol]</c> is the original code unit; strictly ascending.</summary>
    internal char[] Alphabet { get; }

    /// <summary>The symbols when <see cref="Width"/> is Byte, otherwise empty.</summary>
    internal ReadOnlySpan<byte> Bytes => _bytes;

    /// <summary>The symbols when <see cref="Width"/> is UShort, otherwise empty.</summary>
    internal ReadOnlySpan<ushort> UShorts => _ushorts;

    /// <summary>Reader over <see cref="Bytes"/>.</summary>
    internal ByteSymbols ByteReader() => new(_bytes);

    /// <summary>Reader over <see cref="UShorts"/>.</summary>
    internal UShortSymbols UShortReader() => new(_ushorts);

    /// <summary>The symbol at <paramref name="index"/> whatever the width (convenience for tests and cold paths).</summary>
    internal int SymbolAt(int index) => _bytes is not null ? _bytes[index] : _ushorts![index];
}
