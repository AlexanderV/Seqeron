using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>
/// The write side of B3: converts dense LCP values, in strictly ascending rank order, into the base array of the chosen width (little endian, escape value in the narrow kinds) and the overflow list (pairs
/// <c>(rank, value)</c> of two 4-byte fields, ADR-003 section 4 kind 5), both through <see cref="MappedWindowWriter"/> ranges (stream into memory-mapped files). The lengths of both ranges are known before the first
/// entry is written because the width was chosen from the counts of B2. Disposing disposes both writers. Not thread-safe.
/// </summary>
internal sealed class LcpOutputWriter : IDisposable
{
    private readonly LcpKind _kind;
    private readonly uint _escape;
    private readonly int _width;
    private readonly long _overflowExpected;
    private readonly int _slice;
    private readonly MappedWindowWriter _base;
    private readonly MappedWindowWriter _overflow;
    private readonly Buf<byte> _bytes;
    private readonly Buf<uint> _pairs;
    private long _nextRank;
    private long _overflowWritten;

    /// <param name="arena">Owner of the two conversion buffers.</param>
    /// <param name="kind">The width chosen from the counts of B2.</param>
    /// <param name="n">Text length; the base range has <c>n + 1</c> entries of the width of <paramref name="kind"/>.</param>
    /// <param name="overflowCount">Entries of the overflow list (from the counts); the overflow range has <c>8 * overflowCount</c> bytes.</param>
    /// <param name="baseWriter">The range of the base array.</param>
    /// <param name="overflowWriter">The range of the overflow list.</param>
    /// <param name="sliceEntries">Entries converted per step.</param>
    internal LcpOutputWriter(BudgetedArena arena, LcpKind kind, long n, long overflowCount, MappedWindowWriter baseWriter, MappedWindowWriter overflowWriter, int sliceEntries)
    {
        ArgumentNullException.ThrowIfNull(arena);
        ArgumentNullException.ThrowIfNull(baseWriter);
        ArgumentNullException.ThrowIfNull(overflowWriter);
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfNegative(overflowCount);
        _width = LcpWidth.Bytes(kind);
        if (baseWriter.Length != (n + 1) * _width) throw new ArgumentException($"the base range has {baseWriter.Length} bytes, {n + 1} entries of {_width} bytes need {(n + 1) * _width}", nameof(baseWriter));
        if (overflowWriter.Length != overflowCount * 8) throw new ArgumentException($"the overflow range has {overflowWriter.Length} bytes, {overflowCount} pairs need {overflowCount * 8}", nameof(overflowWriter));
        _kind = kind;
        _escape = LcpWidth.Escape(kind);
        _slice = Math.Max(1, sliceEntries);
        _overflowExpected = overflowCount;
        _bytes = arena.Allocate<byte>(_slice * _width, zero: false, "lcp.out");
        _pairs = arena.Allocate<uint>(_slice * 2, zero: false, "lcp.pairs");
        _base = baseWriter;
        _overflow = overflowWriter;
    }

    /// <summary>Appends the LCP values of ranks <c>firstRank ..</c>; the ranks must continue where the previous call ended.</summary>
    internal void Write(long firstRank, ReadOnlySpan<uint> values)
    {
        if (firstRank != _nextRank) throw new InvalidOperationException($"LCP values must arrive in rank order: expected rank {_nextRank}, got {firstRank}");
        for (int at = 0; at < values.Length; at += _slice)
        {
            ReadOnlySpan<uint> part = values.Slice(at, Math.Min(_slice, values.Length - at));
            Span<byte> bytes = _bytes.AsSpan(part.Length * _width);
            int pairs = Convert(firstRank + at, part, bytes);
            _base.Append(bytes);
            if (pairs > 0)
            {
                if (_overflowWritten + pairs > _overflowExpected) throw new InvalidOperationException($"more than the {_overflowExpected} announced escapes: the suffix array is not a permutation");
                _overflow.Append(MemoryMarshal.AsBytes(_pairs.AsSpan(pairs * 2)));
                _overflowWritten += pairs;
            }
        }
        _nextRank = firstRank + values.Length;
    }

    /// <summary>Checks that every entry and exactly the announced number of escapes were written.</summary>
    internal void Complete(long n)
    {
        if (_nextRank != n + 1) throw new InvalidOperationException($"{_nextRank} LCP entries written, {n + 1} expected");
        if (_overflowWritten != _overflowExpected) throw new InvalidOperationException($"{_overflowWritten} escapes met, {_overflowExpected} announced from PLCP: the suffix array is not a permutation");
    }

    public void Dispose()
    {
        _base.Dispose();
        _overflow.Dispose();
    }

    /// <summary>Fills the byte buffer with the entries of <paramref name="part"/> and the pair buffer with its escapes; returns the number of escapes.</summary>
    private int Convert(long firstRank, ReadOnlySpan<uint> part, Span<byte> bytes)
    {
        int pairs = 0;
        Span<uint> list = _pairs.AsSpan(part.Length * 2);
        switch (_kind)
        {
            case LcpKind.U8:
                for (int i = 0; i < part.Length; i++)
                {
                    uint v = part[i];
                    if (v >= _escape)
                    {
                        list[2 * pairs] = (uint)(firstRank + i);
                        list[(2 * pairs) + 1] = v;
                        pairs++;
                        bytes[i] = (byte)_escape;
                    }
                    else
                    {
                        bytes[i] = (byte)v;
                    }
                }
                break;
            case LcpKind.U16:
                Span<ushort> u16 = MemoryMarshal.Cast<byte, ushort>(bytes);
                for (int i = 0; i < part.Length; i++)
                {
                    uint v = part[i];
                    if (v >= _escape)
                    {
                        list[2 * pairs] = (uint)(firstRank + i);
                        list[(2 * pairs) + 1] = v;
                        pairs++;
                        u16[i] = (ushort)_escape;
                    }
                    else
                    {
                        u16[i] = (ushort)v;
                    }
                }
                break;
            default:
                part.CopyTo(MemoryMarshal.Cast<byte, uint>(bytes));
                break;
        }
        return pairs;
    }
}
