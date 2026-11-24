using System.IO.MemoryMappedFiles;
using System.Runtime.CompilerServices;

// A minimal reader of a SESA container for the descent measurement (ADR-003 section 3: the 128-byte header, the 64-byte section entries at 128 + 64 i).
// Independent of the product reader on purpose (the product types are internal); the product is used only as the oracle (Descent.cs --check).

internal sealed record Section(int Kind, long Offset, long Length, long Count, int ElemWidth, int Codec, long Param0, long Param1);

internal sealed unsafe class Container : IDisposable
{
    private readonly MemoryMappedFile _map;
    private readonly MemoryMappedViewAccessor _view;
    internal byte* Base { get; }
    internal long N { get; }
    internal int PosWidth { get; }
    internal int Sigma { get; }
    internal Dictionary<int, Section> Sections { get; } = new();

    internal Container(string path)
    {
        var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, 1, FileOptions.None);
        _map = MemoryMappedFile.CreateFromFile(fs, null, 0, MemoryMappedFileAccess.Read, HandleInheritability.None, false);
        _view = _map.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read);
        byte* p = null;
        _view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
        Base = p + _view.PointerOffset;
        N = *(long*)(Base + 16);
        int count = *(int*)(Base + 48);
        PosWidth = *(int*)(Base + 72);
        Sigma = *(int*)(Base + 76);
        for (int i = 0; i < count; i++)
        {
            byte* e = Base + 128 + (64 * i);
            var s = new Section(*(int*)e, *(long*)(e + 8), *(long*)(e + 16), *(long*)(e + 24), *(int*)(e + 32), *(int*)(e + 36), *(long*)(e + 40), *(long*)(e + 48));
            Sections[s.Kind] = s;
        }
    }

    internal byte* At(int kind) => Base + Sections[kind].Offset;

    internal bool Has(int kind) => Sections.ContainsKey(kind);

    public void Dispose()
    {
        _view.SafeMemoryMappedViewHandle.ReleasePointer();
        _view.Dispose();
        _map.Dispose();
    }

    /// <summary>The dense letters of the alphabet section (code units, ascending).</summary>
    internal char[] Letters()
    {
        Section a = Sections[2];
        var r = new char[a.Count];
        for (int i = 0; i < r.Length; i++) r[i] = (char)((ushort*)(Base + a.Offset))[i];
        return r;
    }
}

/// <summary>The b-bit packed text of TEXT codec 4 (big-endian bit stream, 64 zero bytes of padding), read as MappedSearch does.</summary>
internal readonly unsafe struct PackedText(byte* text, int bits)
{
    internal int Bits => bits;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal int Symbol(long i)
    {
        long bit = i * bits;
        uint w = System.Buffers.Binary.BinaryPrimitives.ReverseEndianness(Unsafe.ReadUnaligned<ushort>(text + (bit >> 3)));
        return (int)(((w << (int)(bit & 7)) & 0xFFFF) >> (16 - bits));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal ulong Window(long i)
    {
        long bit = i * bits;
        byte* q = text + (bit >> 3);
        int s = (int)(bit & 7);
        ulong w = System.Buffers.Binary.BinaryPrimitives.ReverseEndianness(Unsafe.ReadUnaligned<ulong>(q));
        return (w << s) | ((ulong)q[8] >> (8 - s));
    }

    /// <summary>Whether text symbols [start, start + length) equal the pattern symbols [from, from + length) (the pattern packed in words).</summary>
    internal bool Equal(long start, ulong[] packed, int from, int length)
    {
        int per = 64 / bits;
        for (int done = 0; done < length; done += per)
        {
            ulong x = Window(start + done) ^ PatternWindow(packed, bits, from + done);
            int take = Math.Min(per, length - done);
            if ((x >> (64 - (take * bits))) != 0) return false;
        }
        return true;
    }

    private static ulong PatternWindow(ulong[] packed, int bits, int symbol)
    {
        long bit = (long)symbol * bits;
        int w = (int)(bit >> 6);
        int s = (int)(bit & 63);
        ulong high = packed[w] << s;
        return s == 0 ? high : high | (packed[w + 1] >> (64 - s));
    }

    internal static ulong[] Pack(ReadOnlySpan<byte> pattern, int bits)
    {
        var words = new ulong[(((long)pattern.Length * bits) + 63) / 64 + 2];
        long bit = 0;
        for (int i = 0; i < pattern.Length; i++, bit += bits)
        {
            int w = (int)(bit >> 6);
            int room = 64 - (int)(bit & 63) - bits;
            ulong v = pattern[i];
            if (room >= 0) words[w] |= v << room;
            else
            {
                words[w] |= v >> -room;
                words[w + 1] |= v << (64 + room);
            }
        }
        return words;
    }
}
