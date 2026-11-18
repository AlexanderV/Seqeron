using System.Runtime.InteropServices;
using System.Security.Cryptography;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The oracle of the packed text (work package 1b.7b): everything is computed bit by bit from the definition of the layout (a big-endian bit stream, symbol i in the stream bits [i*b, (i+1)*b), the first bit of a
/// byte is its most significant bit, 64 zero bytes after the data), independently of the production code.
/// </summary>
internal static class TextOracle
{
    /// <summary>The file bytes of the packed text: data, unused bits zero, 64 zero bytes.</summary>
    internal static byte[] Pack(ReadOnlySpan<byte> symbols, int bits)
    {
        var bytes = new byte[(((long)symbols.Length * bits) + 7) / 8 + 64];
        for (int i = 0; i < symbols.Length; i++)
        {
            if (bits < 8 && (symbols[i] >> bits) != 0) throw new ArgumentOutOfRangeException(nameof(symbols), "the symbol does not fit");
            for (int b = 0; b < bits; b++)
            {
                long pos = ((long)i * bits) + b;
                if (((symbols[i] >> (bits - 1 - b)) & 1) != 0) bytes[pos / 8] |= (byte)(0x80 >> (int)(pos % 8));
            }
        }
        return bytes;
    }

    internal static int KeySymbols(int bits) => 64 / bits;

    /// <summary>The first 64 / b symbols at <paramref name="i"/> as one number, the symbols past the end are 0, the unused low bits of the 64 are 0.</summary>
    internal static ulong Key(ReadOnlySpan<byte> symbols, long i, int bits)
    {
        int ks = KeySymbols(bits);
        ulong key = 0;
        for (int k = 0; k < ks; k++)
        {
            long at = i + k;
            ulong s = at < symbols.Length ? symbols[(int)at] : 0UL;
            key = (key << bits) | s;
        }
        return key;
    }

    /// <summary>The 64 stream bits that start at symbol <paramref name="i"/> (the bits past the end are 0).</summary>
    internal static ulong Window(ReadOnlySpan<byte> symbols, long i, int bits)
    {
        ulong w = 0;
        for (int j = 0; j < 64; j++)
        {
            long idx = (i * bits) + j;
            long sym = idx / bits;
            int bit = 0;
            if (sym < symbols.Length) bit = (symbols[(int)sym] >> (bits - 1 - (int)(idx % bits))) & 1;
            w = (w << 1) | (uint)bit;
        }
        return w;
    }

    internal static long Lce(ReadOnlySpan<byte> symbols, long i, long j, long max)
    {
        long cap = Math.Min(max, symbols.Length - Math.Max(i, j));
        long k = 0;
        while (k < cap && symbols[(int)(i + k)] == symbols[(int)(j + k)]) k++;
        return k;
    }

    internal static int CompareSuffixes(ReadOnlySpan<byte> symbols, long i, long j)
    {
        ReadOnlySpan<byte> a = i >= symbols.Length ? [] : symbols[(int)i..];
        ReadOnlySpan<byte> b = j >= symbols.Length ? [] : symbols[(int)j..];
        return Math.Sign(a.SequenceCompareTo(b));
    }

    internal static byte[] Random(int n, int sigma, ulong seed, int repeatEvery = 0)
    {
        var rng = new TestSupport.SplitMix64(seed);
        var s = new byte[n];
        for (int i = 0; i < n; i++) s[i] = (byte)rng.NextInt(sigma);
        if (repeatEvery > 0)
        {
            for (int i = repeatEvery; i < n; i++)
            {
                if ((i / repeatEvery) % 3 != 0) s[i] = s[i - repeatEvery];
            }
        }
        return s;
    }

    /// <summary>Every text of length 0..<paramref name="maxLength"/> over <paramref name="sigma"/> symbols.</summary>
    internal static IEnumerable<byte[]> AllTexts(int sigma, int maxLength)
    {
        for (int length = 0; length <= maxLength; length++)
        {
            var text = new byte[length];
            while (true)
            {
                yield return (byte[])text.Clone();
                int p = length - 1;
                while (p >= 0 && text[p] == sigma - 1) text[p--] = 0;
                if (p < 0) break;
                text[p]++;
            }
        }
    }

    /// <summary>The information a reader checks against, computed from dense symbols (letters = the symbols themselves, counts from the text, the hash of the letter bytes).</summary>
    internal static PackedTextInfo Info(ReadOnlySpan<byte> symbols, int bits, int sigma)
    {
        var counts = new long[sigma];
        foreach (byte s in symbols) counts[s]++;
        var letters = new byte[sigma];
        for (int i = 0; i < sigma; i++) letters[i] = (byte)(i + 'A');
        var input = new byte[symbols.Length];
        for (int i = 0; i < input.Length; i++) input[i] = letters[symbols[i]];
        return new PackedTextInfo(symbols.Length, bits, letters, counts, Convert.ToHexStringLower(SHA256.HashData(input)));
    }

    /// <summary>The dense symbols of an input of bytes (the rank of the byte among the distinct bytes) and the distinct bytes.</summary>
    internal static (byte[] Dense, byte[] Letters) Dense(ReadOnlySpan<byte> input)
    {
        var seen = new bool[256];
        foreach (byte b in input) seen[b] = true;
        var letters = new List<byte>();
        var rank = new byte[256];
        for (int v = 0; v < 256; v++)
        {
            if (!seen[v]) continue;
            rank[v] = (byte)letters.Count;
            letters.Add((byte)v);
        }
        var dense = new byte[input.Length];
        for (int i = 0; i < dense.Length; i++) dense[i] = rank[input[i]];
        return (dense, [.. letters]);
    }

    internal static int BitsOf(int sigma)
    {
        int b = 1;
        while ((1 << b) < sigma) b++;
        return b;
    }

    /// <summary>Writes the packed bytes of <paramref name="symbols"/> to a file and returns its path.</summary>
    internal static string WriteFile(TempDir dir, string name, ReadOnlySpan<byte> symbols, int bits)
    {
        string path = dir.File(name);
        File.WriteAllBytes(path, Pack(symbols, bits));
        return path;
    }

    /// <summary>Sets the symbol at <paramref name="position"/> of a packed file (read-modify-write of the bits; for sparse files of any size).</summary>
    internal static void Plant(string path, long dataOffset, int bits, long position, byte symbol)
    {
        using var h = File.OpenHandle(path, FileMode.Open, FileAccess.ReadWrite);
        var one = new byte[1];
        for (int b = 0; b < bits; b++)
        {
            long bitPos = (position * bits) + b;
            long at = dataOffset + (bitPos / 8);
            RandomAccess.Read(h, one, at);
            byte mask = (byte)(0x80 >> (int)(bitPos % 8));
            if (((symbol >> (bits - 1 - b)) & 1) != 0) one[0] |= mask;
            else one[0] &= (byte)~mask;
            RandomAccess.Write(h, one, at);
        }
    }
}

/// <summary>A byte array pinned for the lifetime of the object, with a reader over it (the in-RAM way to build a reader: no file).</summary>
internal sealed unsafe class PinnedText : IDisposable
{
    private readonly GCHandle _handle;

    internal PinnedText(byte[] fileBytes, long n, int bits)
    {
        FileBytes = fileBytes;
        _handle = GCHandle.Alloc(fileBytes, GCHandleType.Pinned);
        Reader = MappedTextReader.FromMemory((byte*)_handle.AddrOfPinnedObject(), fileBytes.Length, n, bits);
    }

    internal PinnedText(byte[] symbols, int bits)
        : this(TextOracle.Pack(symbols, bits), symbols.Length, bits)
    {
    }

    internal byte[] FileBytes { get; }

    internal MappedTextReader Reader { get; }

    public void Dispose()
    {
        Reader.Dispose();
        _handle.Free();
    }
}
