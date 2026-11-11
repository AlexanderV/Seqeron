using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Tests.Build;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Ph4;

/// <summary>What one PH4 build produced: the stats, the two output files as bytes and the files decoded into the in-repository compact form.</summary>
internal sealed record LcpOutcome(ExtLcpStats Stats, CompactLcp Compact, byte[] BaseBytes, byte[] OverflowBytes);

/// <summary>
/// The oracles of gate G1 of the prototype E4 (ADR-004 section 5): the kernel <see cref="LcpKernel"/> (three passes in RAM, compact form) and the safe <see cref="LcpReference"/>, and the harness that runs
/// the file builder on a suffix array file and decodes its two output files (base array, pairs <c>(rank, value)</c>) into a <see cref="CompactLcp"/> for the comparison.
/// </summary>
internal static class LcpHarness
{
    /// <summary>Options with tiny buffers so that every flush, window and chunk boundary is crossed on texts of a few hundred symbols.</summary>
    internal static ExtLcpOptions Options(string scratch, long r, int threads = 1, LcpMode mode = LcpMode.Distributed, long roundChunks = 0) => new()
    {
        ScratchDir = scratch,
        ChunkPositions = r,
        Threads = threads,
        ReadBufferBytes = 64,
        TupleBufferBytes = 48,
        WindowBytes = 128,
        Mode = mode,
        RoundChunks = roundChunks,
        ArenaGuards = ArenaOptions.Guarded(0xCD), // every buffer that is read before it is written shows as 0xCDCDCDCD
    };

    internal static string WriteSa(uint[] sa, string dir)
    {
        string path = Path.Combine(dir, "sa-" + Guid.NewGuid().ToString("N") + ".bin");
        File.WriteAllBytes(path, MemoryMarshal.AsBytes(sa.AsSpan()).ToArray());
        return path;
    }

    /// <summary>Runs the builder on <paramref name="sa"/> and decodes the result; the scratch directory must be empty afterwards.</summary>
    internal static LcpOutcome Run(PackedText text, uint[] sa, ExtLcpOptions o, string dir, ILcpChunkConsumer? consumer = null)
    {
        string id = Guid.NewGuid().ToString("N");
        string saPath = WriteSa(sa, dir);
        string basePath = Path.Combine(dir, "lcp-" + id + ".base");
        string overflowPath = Path.Combine(dir, "lcp-" + id + ".ovf");
        try
        {
            ExtLcpStats st = ExternalLcp.Build(saPath, text, basePath, overflowPath, o, consumer);
            byte[] b = File.ReadAllBytes(basePath);
            byte[] v = File.ReadAllBytes(overflowPath);
            CompactLcp c = Decode(st.Kind, text.Length, b, v);
            Assert.That(Directory.Exists(o.ScratchDir) ? Directory.GetFileSystemEntries(o.ScratchDir) : [], Is.Empty, "scratch files left behind");
            return new LcpOutcome(st, c, b, v);
        }
        finally
        {
            File.Delete(saPath);
            File.Delete(basePath);
            File.Delete(overflowPath);
        }
    }

    /// <summary>Decodes the base array (little endian, width of the kind) and the overflow pairs into the in-repository representation.</summary>
    internal static CompactLcp Decode(LcpKind kind, long n, byte[] baseBytes, byte[] overflowBytes)
    {
        int width = kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 };
        Assert.That(baseBytes.Length, Is.EqualTo((n + 1) * width), "base file length");
        Assert.That(overflowBytes.Length % 8, Is.Zero, "overflow file holds pairs of 4 bytes");
        ReadOnlySpan<uint> pairs = MemoryMarshal.Cast<byte, uint>(overflowBytes);
        var ranks = new int[pairs.Length / 2];
        var values = new int[pairs.Length / 2];
        for (int i = 0; i < ranks.Length; i++)
        {
            ranks[i] = (int)pairs[2 * i];
            values[i] = (int)pairs[(2 * i) + 1];
        }
        int length = checked((int)n + 1);
        return kind switch
        {
            LcpKind.U8 => CompactLcp.FromParts(kind, length, baseBytes, null, null, ranks, values),
            LcpKind.U16 => CompactLcp.FromParts(kind, length, null, MemoryMarshal.Cast<byte, ushort>(baseBytes).ToArray(), null, ranks, values),
            _ => CompactLcp.FromParts(kind, length, null, null, MemoryMarshal.Cast<byte, int>(baseBytes).ToArray(), ranks, values),
        };
    }

    /// <summary>The kernel of the S1a work (Phi, PLCP, pass 3 in RAM).</summary>
    internal static CompactLcp Kernel(DenseCase c, uint[] sa)
        => LcpKernel.Build(new PreparedText(StringFamilies.Alphabet(Math.Max(1, c.Sigma)), c.Symbols, null), MemoryMarshal.Cast<uint, int>(sa));

    /// <summary>The safe reference (never optimised).</summary>
    internal static CompactLcp Reference(DenseCase c, uint[] sa)
        => LcpReference.Build(new ByteSymbols(c.Symbols), MemoryMarshal.Cast<uint, int>(sa));

    /// <summary>Null when equal in kind, length, base array and overflow lists; otherwise a description of the first difference.</summary>
    internal static string? Diff(CompactLcp expected, CompactLcp actual)
    {
        if (expected.Kind != actual.Kind) return $"kind {actual.Kind} instead of {expected.Kind}";
        if (expected.Length != actual.Length) return $"length {actual.Length} instead of {expected.Length}";
        bool same = expected.Kind switch
        {
            LcpKind.U8 => expected.BaseU8.AsSpan().SequenceEqual(actual.BaseU8),
            LcpKind.U16 => expected.BaseU16.AsSpan().SequenceEqual(actual.BaseU16),
            _ => expected.BaseU32.AsSpan().SequenceEqual(actual.BaseU32),
        };
        if (!same) return "base arrays differ";
        if (!expected.OverflowRanks.AsSpan().SequenceEqual(actual.OverflowRanks)) return "overflow ranks differ";
        if (!expected.OverflowValues.AsSpan().SequenceEqual(actual.OverflowValues)) return "overflow values differ";
        return null;
    }

    /// <summary>Checks one text against both oracles (compact form against the kernel, dense values against the reference) and returns the outcome.</summary>
    internal static LcpOutcome CheckAgainstOracles(DenseCase c, int bits, ExtLcpOptions o, string dir, string label, bool reference = true, bool scalar = false)
    {
        uint[] sa = c.ReferenceSa();
        using PackedText t = PackedText.Pack(c.Symbols, bits);
        t.UseVectors = t.UseVectors && !scalar;
        LcpOutcome got = Run(t, sa, o, dir);
        string? d = Diff(Kernel(c, sa), got.Compact);
        if (d is not null) Assert.Fail($"{label}: compact form differs from LcpKernel: {d}");
        if (reference && !Reference(c, sa).ToArray().AsSpan().SequenceEqual(got.Compact.ToArray()))
            Assert.Fail($"{label}: decoded values differ from LcpReference");
        return got;
    }

    internal static IEnumerable<byte[]> AllStrings(int sigma, int minLength, int maxLength)
    {
        for (int length = minLength; length <= maxLength; length++)
        {
            long total = 1;
            for (int i = 0; i < length; i++) total *= sigma;
            for (long code = 0; code < total; code++)
            {
                var s = new byte[length];
                long c = code;
                for (int i = length - 1; i >= 0; i--)
                {
                    s[i] = (byte)(c % sigma);
                    c /= sigma;
                }
                yield return s;
            }
        }
    }
}
