namespace SuffixTree.Esa.Tests.Large;

/// <summary>The text of a closed form for the long descent (work package S2.4): the dense symbol at a position, computed (no memory at any n).</summary>
internal readonly struct GenSymbols(ClosedFormText form) : IWideText
{
    /// <inheritdoc />
    public int Symbol(long position) => form.CodeAt(position);
}

/// <summary>A factor of a closed-form text as the pattern of the long descent: <c>text[start .. start + length)</c>, any length (a pattern of 2^32 symbols is no span).</summary>
internal readonly struct FactorPattern(ClosedFormText form, long start, long length) : IWidePattern
{
    /// <inheritdoc />
    public long Length => length;

    /// <inheritdoc />
    public int this[long index]
    {
        get
        {
            if ((ulong)index >= (ulong)length) throw new ArgumentOutOfRangeException(nameof(index), index, $"the pattern has {length} symbols");
            return form.CodeAt(start + index);
        }
    }
}

/// <summary>The pattern <c>a^k</c> (every symbol the dense symbol 0), any length.</summary>
internal readonly struct RunPattern(long length) : IWidePattern
{
    /// <inheritdoc />
    public long Length => length;

    /// <inheritdoc />
    public int this[long index]
    {
        get
        {
            if ((ulong)index >= (ulong)length) throw new ArgumentOutOfRangeException(nameof(index), index, $"the pattern has {length} symbols");
            return 0;
        }
    }
}

/// <summary>The long descent over the generated views of a closed-form text, in the two widths of a container (work package S2.4): the readers are <see cref="PosGen{TWidth}"/> and every read is checked.</summary>
internal static class GenDescent
{
    /// <summary>The LCP array of <paramref name="f"/> as U64 (no escape): the kind that holds a depth of 2^32 and more.</summary>
    internal static CompactLcp Lcp64(ClosedFormText f) =>
        CompactLcp.FromReaders(LcpKind.U64, new PosGen<W8>(f, GenSection.LcpBase), new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);

    /// <summary>The interval of a pattern over 5-byte readers and a U64 LCP array (the layout of a container above 2^32 - 4 symbols).</summary>
    internal static bool Width5<TPat>(ClosedFormText f, TPat pat, out long lo, out long hi)
        where TPat : struct, IWidePattern
    {
        var sa = new PosGen<W5>(f, GenSection.Sa);
        var cells = new WideCells<CheckedAccess, PosGen<W5>, PosGen<W8>>(new PosGen<W5>(f, GenSection.ChildPacked), new PosGen<W8>(f, GenSection.LcpBase), Lcp64(f), tagged: false, f.Length);
        return WideSearch.Run<CheckedAccess, PosGen<W5>, WideCells<CheckedAccess, PosGen<W5>, PosGen<W8>>, GenSymbols, TPat>(sa, cells, new GenSymbols(f), f.Length, pat, out lo, out hi);
    }

    /// <summary>The interval of a pattern over 4-byte readers and a U32 LCP array (the layout of a container of at most 2^32 - 4 symbols).</summary>
    internal static bool Width4<TPat>(ClosedFormText f, TPat pat, out long lo, out long hi)
        where TPat : struct, IWidePattern
    {
        var sa = new PosGen<W4>(f, GenSection.Sa);
        CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U32, new PosGen<W4>(f, GenSection.LcpBase), new PosGen<W4>(f, GenSection.LcpPairsOfPowerOfA, long.MaxValue), isChecked: true);
        var cells = new WideCells<CheckedAccess, PosGen<W4>, PosGen<W4>>(new PosGen<W4>(f, GenSection.ChildPacked), new PosGen<W4>(f, GenSection.LcpBase), lcp, tagged: false, f.Length);
        return WideSearch.Run<CheckedAccess, PosGen<W4>, WideCells<CheckedAccess, PosGen<W4>, PosGen<W4>>, GenSymbols, TPat>(sa, cells, new GenSymbols(f), f.Length, pat, out lo, out hi);
    }

    /// <summary>The interval of a pattern of <c>a^n</c> over an 8- or 16-bit LCP array with an escape and the overflow list of pairs (a^n has an escaped value at every rank above the escape).</summary>
    internal static bool EscapedPowerOfA<TPat>(PowerOfA f, int escape, TPat pat, out long lo, out long hi)
        where TPat : struct, IWidePattern
    {
        var sa = new PosGen<W5>(f, GenSection.Sa);
        var childPacked = new PosGen<W5>(f, GenSection.ChildPacked);
        var pairs = new PosGen<W5>(f, GenSection.LcpPairsOfPowerOfA, escape);
        if (escape == 255)
        {
            CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U8, new PosGen<W1>(f, GenSection.LcpBase, escape), pairs, isChecked: true);
            var cells = new WideCells<CheckedAccess, PosGen<W5>, PosGen<W1>>(childPacked, new PosGen<W1>(f, GenSection.LcpBase, escape), lcp, tagged: false, f.Length);
            return WideSearch.Run<CheckedAccess, PosGen<W5>, WideCells<CheckedAccess, PosGen<W5>, PosGen<W1>>, GenSymbols, TPat>(sa, cells, new GenSymbols(f), f.Length, pat, out lo, out hi);
        }
        else
        {
            CompactLcp lcp = CompactLcp.FromReaders(LcpKind.U16, new PosGen<W2>(f, GenSection.LcpBase, escape), pairs, isChecked: true);
            var cells = new WideCells<CheckedAccess, PosGen<W5>, PosGen<W2>>(childPacked, new PosGen<W2>(f, GenSection.LcpBase, escape), lcp, tagged: false, f.Length);
            return WideSearch.Run<CheckedAccess, PosGen<W5>, WideCells<CheckedAccess, PosGen<W5>, PosGen<W2>>, GenSymbols, TPat>(sa, cells, new GenSymbols(f), f.Length, pat, out lo, out hi);
        }
    }
}
