namespace SuffixTree.Esa;

/// <summary>Storage width of a prepared symbol sequence.</summary>
internal enum SymbolWidth : byte
{
    /// <summary>One byte per symbol (alphabet of at most 256 symbols after the dense remap).</summary>
    Byte,

    /// <summary>Two bytes per symbol (alphabet of up to 65,536 UTF-16 code units).</summary>
    UShort,

    /// <summary>Four bytes per symbol (reduced strings of the SA-IS recursion, integer alphabets).</summary>
    Int,
}

/// <summary>
/// Read access to a sequence of non-negative integer symbols. Implementations are (ref) structs so that the generic kernels are
/// specialised by the JIT per width and the call is inlined (ADR-002 section 2 item 2); there is no virtual dispatch.
/// </summary>
internal interface ISymbolReader
{
    /// <summary>Number of symbols.</summary>
    int Length { get; }

    /// <summary>The symbol at <paramref name="index"/> (0-based, bounds-checked by the underlying span).</summary>
    int this[int index] { get; }
}

/// <summary>Symbols stored as bytes.</summary>
internal readonly ref struct ByteSymbols(ReadOnlySpan<byte> data) : ISymbolReader
{
    private readonly ReadOnlySpan<byte> _data = data;

    public int Length => _data.Length;

    public int this[int index] => _data[index];
}

/// <summary>Symbols stored as 16-bit unsigned integers.</summary>
internal readonly ref struct UShortSymbols(ReadOnlySpan<ushort> data) : ISymbolReader
{
    private readonly ReadOnlySpan<ushort> _data = data;

    public int Length => _data.Length;

    public int this[int index] => _data[index];
}

/// <summary>Symbols stored as 32-bit integers.</summary>
internal readonly ref struct IntSymbols(ReadOnlySpan<int> data) : ISymbolReader
{
    private readonly ReadOnlySpan<int> _data = data;

    public int Length => _data.Length;

    public int this[int index] => _data[index];
}
