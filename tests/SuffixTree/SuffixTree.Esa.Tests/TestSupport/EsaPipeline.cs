using System.Buffers.Binary;
using System.Security.Cryptography;

namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>The reference pipeline (text preparation, SA-IS, Phi/PLCP LCP) and the hashing used by the libsais vectors.</summary>
internal static class EsaPipeline
{
    /// <summary>Suffix array (n + 1 layout) and compact LCP of <paramref name="text"/> by <see cref="SaisReference"/> and <see cref="LcpReference"/>.</summary>
    internal static (int[] Sa, CompactLcp Lcp) Build(string text)
    {
        PreparedText prepared = TextPreparation.Prepare(text);
        if (prepared.Width == SymbolWidth.Byte)
        {
            ByteSymbols reader = prepared.ByteReader();
            int[] sa = SaisReference.Build(reader, prepared.AlphabetSize);
            return (sa, LcpReference.Build(reader, sa));
        }
        else
        {
            UShortSymbols reader = prepared.UShortReader();
            int[] sa = SaisReference.Build(reader, prepared.AlphabetSize);
            return (sa, LcpReference.Build(reader, sa));
        }
    }

    /// <summary>
    /// SHA-256 (lower-case hex) of <c>sa[1..n]</c> as n little-endian int32: the libsais SA, which has no entry for the empty suffix.
    /// </summary>
    internal static string HashSa(ReadOnlySpan<int> sa)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        AppendInt32Le(hash, sa.Slice(1));
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }

    /// <summary>
    /// SHA-256 (lower-case hex) of <c>lcp[1..n]</c> as n little-endian int32 (our <c>LCP[1] = 0</c>, <c>LCP[r] = lcp(SA[r-1], SA[r])</c>),
    /// which is the libsais <c>LCP[0..n-1]</c>.
    /// </summary>
    internal static string HashLcp(CompactLcp lcp, int n)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        var block = new int[4096];
        var cursor = lcp.CreateCursor();
        for (int r = 1; r <= n;)
        {
            int k = Math.Min(block.Length, n - r + 1);
            for (int j = 0; j < k; j++) block[j] = cursor.Get(r + j);
            AppendInt32Le(hash, block.AsSpan(0, k));
            r += k;
        }
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }

    private static void AppendInt32Le(IncrementalHash hash, ReadOnlySpan<int> values)
    {
        Span<byte> buffer = stackalloc byte[4 * 1024];
        while (!values.IsEmpty)
        {
            int k = Math.Min(values.Length, 1024);
            for (int j = 0; j < k; j++) BinaryPrimitives.WriteInt32LittleEndian(buffer.Slice(4 * j, 4), values[j]);
            hash.AppendData(buffer.Slice(0, 4 * k));
            values = values.Slice(k);
        }
    }
}
