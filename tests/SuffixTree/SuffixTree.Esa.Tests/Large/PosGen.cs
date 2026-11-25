using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>The stored width that a generated section stands for (the byte count of <see cref="IPosReader.ElementBytes"/>).</summary>
internal interface IGenWidth
{
    static abstract int Bytes { get; }
}

internal readonly struct W1 : IGenWidth
{
    public static int Bytes => 1;
}

internal readonly struct W2 : IGenWidth
{
    public static int Bytes => 2;
}

internal readonly struct W4 : IGenWidth
{
    public static int Bytes => 4;
}

internal readonly struct W5 : IGenWidth
{
    public static int Bytes => 5;
}

internal readonly struct W8 : IGenWidth
{
    public static int Bytes => 8;
}

/// <summary>What a <see cref="PosGen{TWidth}"/> computes.</summary>
internal enum GenSection
{
    /// <summary><c>SA[r]</c> of the closed form (n + 1 elements).</summary>
    Sa,

    /// <summary><c>min(LCP[r], escape)</c>: the base array of the LCP of a kind (the escape is 255 for U8, 65,535 for U16, none above).</summary>
    LcpBase,

    /// <summary>The overflow pairs of the LCP of <c>a^n</c> for an escape value: pair i is <c>(escape + 1 + i, escape + i)</c>, as 2 elements each.</summary>
    LcpPairsOfPowerOfA,

    /// <summary>A hash of the rank below 50: an LCP array that is no text's (the pyramid's structure does not care).</summary>
    Synthetic,

    /// <summary>The cell of CHILD codec 1 (packed).</summary>
    ChildPacked,

    /// <summary>The cell of CHILD codec 2 (tagged: the kind in bits 30-31; texts of at most 2^30 - 1 symbols).</summary>
    ChildTagged,

    /// <summary><c>PSI[r] = ISA[SA[r] + 1]</c> of the closed form, 0 for the empty suffix and for the suffix of the last symbol (work package S2.5).</summary>
    Psi,

    /// <summary>The entries of the RMQ section of the closed form, level 0 first (fanout 32): the closed-form minimum of the LCP values of the block (work package S2.5).</summary>
    Rmq,
}

/// <summary>
/// The test-only reader of work package S2.3 (ADR-005 section 3.7, layer 1): a section of a container that is computed, not stored, from a closed form of S2.1 (<see cref="ClosedFormText"/>), so that the views run at n = 2^31 + 10,
/// 2^32 + 10 and 2^33 and above with no memory. Every read is checked (a test reader has no <c>Trusted</c> mode: a bound that a view forgets is a failure here). The type argument is the stored width that the section stands for.
/// </summary>
internal readonly struct PosGen<TWidth> : IPosReader
    where TWidth : struct, IGenWidth
{
    private readonly ClosedFormText _form;
    private readonly GenSection _section;
    private readonly long _count;
    private readonly long _escape;
    private readonly (long Start, long Count)[]? _levels;

    internal PosGen(ClosedFormText form, GenSection section, long escape = long.MaxValue)
    {
        _form = form;
        _section = section;
        _escape = escape;
        _levels = section == GenSection.Rmq ? [.. SparseContainer.RmqLevels(form.Length + 1)] : null;
        _count = section switch
        {
            GenSection.LcpPairsOfPowerOfA => 2 * Math.Max(0, form.Length - escape),
            GenSection.Rmq => _levels!.Sum(l => l.Count),
            _ => form.Length + 1,
        };
    }

    public static int ElementBytes => TWidth.Bytes;

    public long Count => _count;

    public ulong Get(long index)
    {
        if ((ulong)index >= (ulong)_count) throw new InvalidDataException($"a generated section of {_count} elements is read at element {index}: the file is damaged");
        switch (_section)
        {
            case GenSection.Sa:
                return (ulong)_form.SuffixAt(index);
            case GenSection.LcpBase:
                return (ulong)Math.Min(_form.Lcp(index), _escape);
            case GenSection.LcpPairsOfPowerOfA:
                {
                    long pair = index >> 1;
                    return (index & 1) == 0 ? (ulong)(_escape + 1 + pair) : (ulong)(_escape + pair);
                }
            case GenSection.Synthetic:
                return Hash(index) % 50;
            case GenSection.ChildPacked:
                return (ulong)_form.PackedCell(index);
            case GenSection.Psi:
                return (ulong)SparseContainer.PsiOf(_form, index);
            case GenSection.Rmq:
                return (ulong)RmqEntry(index);
            default:
                return Tagged(_form, index);
        }
    }

    /// <summary>The entry <paramref name="index"/> of the RMQ section: the level that holds it, then the minimum of the LCP values of its block.</summary>
    private long RmqEntry(long index)
    {
        long rows = _form.Length + 1;
        long block = SesaFormat.RmqFanout;
        foreach ((long start, long count) in _levels!)
        {
            if (index < start + count)
            {
                long first = (index - start) * block;
                return _form.LcpMin(first, Math.Min(rows - 1, first + block - 1));
            }
            block *= SesaFormat.RmqFanout;
        }
        throw new InvalidDataException($"the RMQ section of {_count} entries is read at entry {index}");
    }

    /// <summary>The cell of CHILD codec 2 of a closed form: the next l-index (kind 1), else the up field (kind 2), else the down field (kind 3), else 0.</summary>
    internal static ulong Tagged(ClosedFormText f, long rank)
    {
        long next = f.Next(rank);
        if (next != 0) return (1UL << 30) | (ulong)next;
        long up = f.Up(rank);
        if (up != 0) return (2UL << 30) | (ulong)up;
        long down = f.Down(rank);
        return down != 0 ? (3UL << 30) | (ulong)down : 0;
    }

    internal static ulong Hash(long i)
    {
        ulong x = (ulong)i + 0x9E3779B97F4A7C15UL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9UL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBUL;
        return x ^ (x >> 31);
    }
}
