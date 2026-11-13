using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Pipeline;

/// <summary>What PH0 found: the length, the alphabet (the distinct byte values, increasing; symbol i is <c>Letters[i]</c>), the symbol counts C[], the width b and the SHA-256 of the input bytes.</summary>
internal sealed record PrepareResult(long N, int Bits, byte[] Letters, long[] Counts, string Sha256)
{
    internal int Sigma => Letters.Length;

    /// <summary>The content of the <c>.pk.meta</c> file: the keys that <see cref="Experiments.E3.ReadMeta"/> reads (<c>n</c>, <c>bits</c>) plus the alphabet scan, the counts and the hash.</summary>
    internal string MetaText(string name) => string.Create(CultureInfo.InvariantCulture,
        $"n={N}\nbits={Bits}\nsigma={Sigma}\nletters_hex={Convert.ToHexString(Letters)}\ncounts={string.Join(',', Counts)}\nsha256={Sha256}\nname={name}\n");
}

/// <summary>
/// PH0 of ADR-004 section 4.2: the alphabet scan, the symbol counts, the SHA-256 of the text, and the b-bit packed text streamed into its file. Two sequential passes over the input in blocks of a fixed number of symbols
/// (the width b is known only after the first one): the first pass reads, counts and hashes; the second reads again, maps the bytes to dense symbols (increasing byte order) and packs one block at a time with
/// <see cref="PackedText.PackInto"/> into a <see cref="MappedWindowWriter"/> that creates the file at its final length (data bytes plus the 64 zero bytes of padding). Memory: one block of the input, one of symbols, one of packed bytes.
/// </summary>
internal static class Prepare
{
    /// <summary>Symbols per block that keep the packed blocks byte aligned (8 symbols of b bits are b bytes).</summary>
    internal const int BlockMultiple = 8;

    internal static PrepareResult Run(string inputPath, string packedPath, int blockSymbols, int windowBytes, Action? betweenPasses = null)
    {
        if (blockSymbols < BlockMultiple || blockSymbols % BlockMultiple != 0) throw new ArgumentOutOfRangeException(nameof(blockSymbols), blockSymbols, "a block holds a multiple of 8 symbols, so that every packed block ends on a byte boundary");
        long n = new FileInfo(inputPath).Length;
        var counts = new long[256];
        string sha;
        var block = new byte[blockSymbols];
        using (var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read))
        using (var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256))
        {
            long at = 0;
            while (at < n)
            {
                int take = (int)Math.Min(blockSymbols, n - at);
                input.ReadExactAt(block.AsSpan(0, take), at);
                foreach (byte b in block.AsSpan(0, take)) counts[b]++;
                hash.AppendData(block.AsSpan(0, take));
                at += take;
            }
            sha = Convert.ToHexStringLower(hash.GetHashAndReset());
        }

        betweenPasses?.Invoke();
        var letters = new List<byte>();
        var map = new byte[256];
        for (int v = 0; v < 256; v++)
        {
            if (counts[v] == 0) continue;
            map[v] = (byte)letters.Count;
            letters.Add((byte)v);
        }
        int bits = Math.Max(1, PackedText.BitsFor(letters.Count));
        long dataBytes = PackedText.DataBytes(n, bits);
        var packed = new byte[(blockSymbols * (long)bits / 8) + 1];
        var seen = new long[256];
        unsafe
        {
            using var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read);
            using var writer = new MappedWindowWriter(packedPath, PackedText.TotalBytes(n, bits), windowBytes);
            long at = 0;
            while (at < n)
            {
                int take = (int)Math.Min(blockSymbols, n - at);
                input.ReadExactAt(block.AsSpan(0, take), at);
                for (int i = 0; i < take; i++)
                {
                    byte symbol = map[block[i]];
                    block[i] = symbol;
                    seen[symbol]++;
                }
                int bytes = (int)PackedText.DataBytes(take, bits);
                fixed (byte* dst = packed) PackedText.PackInto(block.AsSpan(0, take), bits, dst);
                writer.Append(packed.AsSpan(0, bytes));
                at += take;
            }
            for (int s = 0; s < letters.Count; s++)
            {
                if (seen[s] != counts[letters[s]]) throw new InvalidOperationException($"symbol {s} occurred {seen[s]} times in the second pass and {counts[letters[s]]} times in the first: the input changed between the two passes");
            }
            if (writer.Position != dataBytes) throw new InvalidOperationException($"the packed text has {writer.Position} data bytes, expected {dataBytes}: the input changed between the two passes");
        }
        return new PrepareResult(n, bits, [.. letters], [.. letters.Select(l => counts[l])], sha);
    }
}
