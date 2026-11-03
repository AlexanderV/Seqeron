namespace SuffixTree.Esa;

/// <summary>
/// The symbols of one SA-IS level TOGETHER with their S/L types, as the induction, naming and bookkeeping passes read them (iteration 1o.1c). A pass that needs the
/// symbol and the type of the same position (every induction step does) calls <see cref="Get"/>, which is one memory load in the packed representation: the induction
/// passes of a large text are bound by the cache misses of these random reads, and the packed form halves them.
/// </summary>
internal interface ITypedText
{
    /// <summary>The symbol and the type (1 = S-type, 0 = L-type) of position <paramref name="i"/> in <c>[0, n - 1]</c>.</summary>
    void Get(int i, out int symbol, out int isS);

    /// <summary>The symbol at <paramref name="i"/>.</summary>
    int Symbol(int i);

    /// <summary>The type of position <paramref name="i"/> in <c>[0, n - 1]</c>: 1 for S-type, 0 for L-type.</summary>
    int IsS(int i);

    /// <summary>True when positions <paramref name="x"/> and <paramref name="y"/> (both in <c>[0, n - 1]</c>) have the same symbol and the same type.</summary>
    bool Same(int x, int y);

    /// <summary>Hints that position <paramref name="i"/> (any int: the address is not checked) will be read soon.</summary>
    void Prefetch(int i);
}

/// <summary>Where the classification writes: the type array of a split text, or the packed bytes.</summary>
internal interface ITypeSink
{
    /// <summary>The virtual sentinel at position n is S-type (never read by the passes).</summary>
    void SetSentinel(int n);

    /// <summary>Position <paramref name="i"/> has symbol <paramref name="symbol"/> and type <paramref name="isS"/>.</summary>
    void Set(int i, int symbol, int isS);
}

/// <summary>Symbols in their own array (any width) and the types as one byte per position in a second array.</summary>
internal readonly struct SplitTyped<TText>(TText text, Buf<byte> types) : ITypedText
    where TText : struct, IPtrText
{
    public void Get(int i, out int symbol, out int isS)
    {
        symbol = text[i];
        isS = types[i];
    }

    public int Symbol(int i) => text[i];

    public int IsS(int i) => types[i];

    public bool Same(int x, int y) => text[x] == text[y] && types[x] == types[y];

    public void Prefetch(int i)
    {
        text.Prefetch(i);
        types.Prefetch(i);
    }
}

/// <summary>Writes the types of a split text.</summary>
internal readonly struct SplitSink(Buf<byte> types) : ITypeSink
{
    public void SetSentinel(int n) => types[n] = 1;

    public void Set(int i, int symbol, int isS) => types[i] = (byte)isS;
}

/// <summary>
/// Byte symbols below <see cref="Limit"/> with their types in ONE byte per position, <c>(symbol &lt;&lt; 1) | isS</c> (DNA: 3 bits). The original text is not read again after the
/// classification; one load gives symbol and type, and two positions are the same symbol-and-type iff their bytes are equal.
/// </summary>
internal readonly struct PackedBytes(Buf<byte> packed) : ITypedText
{
    /// <summary>Alphabets of fewer symbols than this fit in the packed byte.</summary>
    internal const int Limit = 128;

    public void Get(int i, out int symbol, out int isS)
    {
        int v = packed[i];
        symbol = v >> 1;
        isS = v & 1;
    }

    public int Symbol(int i) => packed[i] >> 1;

    public int IsS(int i) => packed[i] & 1;

    public bool Same(int x, int y) => packed[x] == packed[y];

    public void Prefetch(int i) => packed.Prefetch(i);
}

/// <summary>Writes the packed bytes.</summary>
internal readonly struct PackedSink(Buf<byte> packed) : ITypeSink
{
    public void SetSentinel(int n) => packed[n] = 1;

    public void Set(int i, int symbol, int isS) => packed[i] = (byte)((symbol << 1) | isS);
}
