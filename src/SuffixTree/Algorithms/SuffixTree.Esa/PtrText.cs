namespace SuffixTree.Esa;

/// <summary>Read access to the symbols of a level through <see cref="Buf{T}"/> (a raw load in Release); the pointer passes of the kernels are generic over it, so each width gets its own specialised code.</summary>
internal interface IPtrText
{
    /// <summary>Number of symbols.</summary>
    int Length { get; }

    /// <summary>The symbol at <paramref name="index"/>.</summary>
    int this[int index] { get; }
}

/// <summary>Symbols stored as bytes (alphabet of at most 256 symbols after the dense remap).</summary>
internal readonly struct BytePtrText(Buf<byte> data, int length) : IPtrText
{
    public int Length { get; } = length;

    public int this[int index] => data[index];
}

/// <summary>Symbols stored as 16-bit unsigned integers.</summary>
internal readonly struct UShortPtrText(Buf<ushort> data, int length) : IPtrText
{
    public int Length { get; } = length;

    public int this[int index] => data[index];
}

/// <summary>Symbols stored as 32-bit integers (the reduced strings of the recursion, which live in the suffix-array buffer).</summary>
internal readonly struct IntPtrText(Buf<int> data, int length) : IPtrText
{
    public int Length { get; } = length;

    public int this[int index] => data[index];
}
