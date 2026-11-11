using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Ph4;

/// <summary>The width rule of ADR-002 section 4 on counts that are known before the first LCP entry is written.</summary>
internal static class LcpWidth
{
    /// <summary>
    /// <c>U8</c> + overflow iff <c>16 * escapes8 &lt;= n</c>, else <c>U16</c> + overflow iff <c>8 * escapes16 &lt;= n</c>, else <c>U32</c>: the rule of <see cref="CompactLcp.SelectKind"/> on 64-bit counts
    /// (that method takes an <c>int</c> n; the tests pin both to each other).
    /// </summary>
    internal static LcpKind Select(long n, long escapes8, long escapes16)
    {
        if (escapes16 > escapes8) throw new ArgumentException("every value >= 65,535 is also >= 255, so escapes16 cannot exceed escapes8");
        if (16 * escapes8 <= n) return LcpKind.U8;
        if (8 * escapes16 <= n) return LcpKind.U16;
        return LcpKind.U32;
    }

    /// <summary>Bytes per entry of the base array.</summary>
    internal static int Bytes(LcpKind kind) => kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 };

    /// <summary>The escape value of the base array (a value at or above it is stored in the overflow list); <see cref="uint.MaxValue"/> for U32, which has none.</summary>
    internal static uint Escape(LcpKind kind) => kind switch { LcpKind.U8 => byte.MaxValue, LcpKind.U16 => ushort.MaxValue, _ => uint.MaxValue };

    /// <summary>Entries of the overflow list for the counts of values at or above 255 and 65,535.</summary>
    internal static long OverflowCount(LcpKind kind, long escapes8, long escapes16) => kind switch { LcpKind.U8 => escapes8, LcpKind.U16 => escapes16, _ => 0 };
}

/// <summary>
/// The write side of B3: converts dense LCP values, in strictly ascending rank order, into the base array of the chosen width (little endian, escape value in the narrow kinds) and the overflow list
/// (pairs <c>(rank, value)</c> of two 4-byte fields, ADR-003 section 4 kind 5), both through <see cref="MappedWindowWriter"/> (the owner's rule: stream into memory-mapped files). The lengths of both files
/// are known before the first entry is written because the width was chosen from the counts of B2. Not thread-safe.
/// </summary>
internal sealed unsafe class LcpOutputWriter : IDisposable
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
    /// <param name="n">Text length; the base array has <c>n + 1</c> entries.</param>
    /// <param name="overflowCount">Entries of the overflow list (from the counts).</param>
    /// <param name="basePath">The base array file.</param>
    /// <param name="overflowPath">The overflow file (created empty when there is no overflow).</param>
    /// <param name="windowBytes">Size of the sliding windows.</param>
    /// <param name="sliceEntries">Entries converted per step.</param>
    internal LcpOutputWriter(BudgetedArena arena, LcpKind kind, long n, long overflowCount, string basePath, string overflowPath, int windowBytes, int sliceEntries)
    {
        _kind = kind;
        _width = LcpWidth.Bytes(kind);
        _escape = LcpWidth.Escape(kind);
        _slice = Math.Max(1, sliceEntries);
        _overflowExpected = overflowCount;
        _bytes = arena.Allocate<byte>(_slice * _width, zero: false, "lcp.out");
        _pairs = arena.Allocate<uint>(_slice * 2, zero: false, "lcp.pairs");
        _base = new MappedWindowWriter(basePath, (n + 1) * _width, windowBytes);
        _overflow = new MappedWindowWriter(overflowPath, overflowCount * 8, windowBytes);
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
