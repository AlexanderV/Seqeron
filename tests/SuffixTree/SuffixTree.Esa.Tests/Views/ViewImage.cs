using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Fuzz;

namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The sections of an index of the S1 in-memory builder re-encoded at a stored width (4 or 5 bytes, ADR-005 section 6.1: five-byte entries followed by three bytes of padding inside the payload; the pairs of the overflow list as two fields
/// of the width) and at an LCP kind (U8, U16, U32 or U64), each section in memory that ends against an inaccessible page (<see cref="GuardedImage"/>): a view that reads one byte past its payload ends the test host. The encoding is written here
/// from the byte layout and shares nothing with the readers; the logical values (<see cref="Sa"/>, <see cref="Lcp"/>, the child cells) are the in-memory index's and are what every view must give back.
/// </summary>
internal sealed unsafe class ViewImage : IDisposable
{
    private readonly GuardedImage _sa;
    private readonly GuardedImage _child;
    private readonly GuardedImage _lcpBase;
    private readonly GuardedImage? _pairs;

    private ViewImage(LcpChildCase c, int posWidth, LcpKind kind, bool tagged)
    {
        Case = c;
        PosWidth = posWidth;
        Kind = kind;
        Tagged = tagged;
        Count = c.N + 1L;
        Sa = [.. c.Sa.Select(v => (long)v)];
        Lcp = [.. c.Lcp.ToArray().Select(v => (long)v)];
        uint[] cells = c.Child(tagged ? ChildCodec.Tagged : ChildCodec.Packed);
        Cells = [.. cells.Select(v => (ulong)v)];
        int padding = posWidth == 5 ? PosReader.Padding40 : 0;
        _sa = new GuardedImage(Encode(Sa.Select(v => (ulong)v), posWidth, padding), flush: true);
        _child = new GuardedImage(Encode(Cells, posWidth, padding), flush: true);
        (byte[] baseBytes, List<(long Rank, long Value)> pairs) = EncodeLcp(Lcp, kind);
        _lcpBase = new GuardedImage(baseBytes, flush: true);
        Pairs = pairs;
        if (pairs.Count > 0)
        {
            var fields = new List<ulong>();
            foreach ((long rank, long value) in pairs)
            {
                fields.Add((ulong)rank);
                fields.Add((ulong)value);
            }
            _pairs = new GuardedImage(Encode(fields, posWidth, padding), flush: true);
        }
    }

    internal LcpChildCase Case { get; }

    internal int PosWidth { get; }

    internal LcpKind Kind { get; }

    internal bool Tagged { get; }

    /// <summary>n + 1.</summary>
    internal long Count { get; }

    internal long[] Sa { get; }

    internal long[] Lcp { get; }

    /// <summary>The raw child cells (the packed or the tagged codec of the S1 index).</summary>
    internal ulong[] Cells { get; }

    /// <summary>The overflow pairs in rank order.</summary>
    internal List<(long Rank, long Value)> Pairs { get; }

    internal byte* SaPointer => _sa.Data;

    internal byte* ChildPointer => _child.Data;

    internal byte* LcpBasePointer => _lcpBase.Data;

    internal byte* PairsPointer => _pairs?.Data;

    internal static ViewImage Of(string text, int posWidth, LcpKind kind, bool tagged = false) => new(LcpChildCase.Of(text), posWidth, kind, tagged);

    internal static ViewImage Of(LcpChildCase c, int posWidth, LcpKind kind, bool tagged = false) => new(c, posWidth, kind, tagged);

    /// <summary>The bytes of <paramref name="values"/> at <paramref name="width"/> bytes per element, little endian, and <paramref name="padding"/> zero bytes at the end.</summary>
    internal static byte[] Encode(IEnumerable<ulong> values, int width, int padding)
    {
        ulong[] v = [.. values];
        var bytes = new byte[(width * v.Length) + padding];
        for (int i = 0; i < v.Length; i++)
        {
            for (int j = 0; j < width; j++) bytes[(width * i) + j] = (byte)(v[i] >> (8 * j));
        }
        return bytes;
    }

    /// <summary>The base array of the LCP at the width of <paramref name="kind"/> and the pairs <c>(rank, value)</c> of its escapes.</summary>
    internal static (byte[] Base, List<(long Rank, long Value)> Pairs) EncodeLcp(long[] values, LcpKind kind)
    {
        int width = kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, LcpKind.U32 => 4, _ => 8 };
        long escape = kind switch { LcpKind.U8 => 255, LcpKind.U16 => 65_535, _ => long.MaxValue };
        var pairs = new List<(long, long)>();
        var stored = new ulong[values.Length];
        for (int i = 0; i < values.Length; i++)
        {
            if (values[i] >= escape)
            {
                stored[i] = (ulong)escape;
                pairs.Add((i, values[i]));
            }
            else
            {
                stored[i] = (ulong)values[i];
            }
        }
        return (Encode(stored, width, 0), pairs);
    }

    public void Dispose()
    {
        _sa.Dispose();
        _child.Dispose();
        _lcpBase.Dispose();
        _pairs?.Dispose();
    }
}
