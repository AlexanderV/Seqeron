using System.Runtime.InteropServices;
using System.Security.Cryptography;
using SuffixTree.EmProto.Pipeline;
using SuffixTree.EmProto.Tests.Ph5;
using SuffixTree.EmProto.Tests.WholePipeline;
using SuffixTree.EmProto.Verify;
using SuffixTree.Esa;

namespace SuffixTree.EmProto.Tests.E8;

/// <summary>Output sets for the verifier tests: written by the pipeline itself (tiny parameters) or from the in-RAM structures of the product (any width of the compact LCP), copies, and byte-level edits.</summary>
internal static class VerifierHarness
{
    internal static readonly string[] Extensions = ["pk", "pk.meta", "sa", "lcp", "ovf", "cld", "stats"];

    /// <summary>The pipeline on <paramref name="text"/> with the tiny parameters of variant <paramref name="variant"/>; returns the output prefix (the input file is <c>in.txt</c> in the same directory).</summary>
    internal static string FromPipeline(string dir, string text, int variant)
    {
        EndToEnd.Run(PipelineChecks.Options(dir, text, variant, 1 + (variant % 3)));
        return Path.Combine(dir, "out");
    }

    /// <summary>All seven output files from the in-RAM structures (SA-IS, LCP kernel, packed child table of the product): the width of the LCP is the kernel's rule.</summary>
    internal static string FromReference(string dir, string text)
    {
        ChildReference r = ChildReference.Of(text);
        (byte[] symbols, char[] letters) = PipelineChecks.Dense(text);
        string prefix = Path.Combine(dir, "out");
        File.WriteAllBytes(Path.Combine(dir, "in.txt"), PipelineChecks.Latin1(text));
        File.WriteAllBytes(prefix + ".pk", PipelineChecks.ReferencePacked(symbols, letters.Length));
        long[] counts = letters.Select(c => (long)text.Count(x => x == c)).ToArray();
        int bits = Math.Max(1, (int)Math.Ceiling(Math.Log2(Math.Max(1, letters.Length))));
        var prep = new PrepareResult(text.Length, bits, [.. letters.Select(c => (byte)c)], counts, Convert.ToHexStringLower(SHA256.HashData(PipelineChecks.Latin1(text))));
        File.WriteAllText(prefix + ".pk.meta", prep.MetaText("in.txt"));
        File.WriteAllBytes(prefix + ".sa", MemoryMarshal.AsBytes(r.Index.SuffixArray.Select(v => (uint)v).ToArray().AsSpan()).ToArray());
        r.WriteLcpFiles(dir, "out");
        File.WriteAllBytes(prefix + ".cld", MemoryMarshal.AsBytes(r.Packed.AsSpan()).ToArray());
        int width = r.Lcp.Kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 };
        long max = 0;
        for (int i = 0; i <= text.Length; i++) max = Math.Max(max, r.Lcp[i]);
        File.WriteAllText(prefix + ".stats", string.Create(System.Globalization.CultureInfo.InvariantCulture,
            $"nodecount={r.Statistics.NodeCount}\ninternal={r.Statistics.InternalCount}\nlrsdepth={r.Statistics.LrsDepth}\nlrsrank={r.Statistics.LrsRank}\nn={text.Length}\nsigma={letters.Length}\nlcpwidth={width}\noverflow={r.Lcp.OverflowCount}\nmaxlcp={max}\n"));
        return prefix;
    }

    /// <summary>Copies the seven outputs (and the input text) of <paramref name="prefix"/> into <paramref name="dir"/>; returns the new prefix.</summary>
    internal static string Clone(string prefix, string dir)
    {
        Directory.CreateDirectory(dir);
        foreach (string e in Extensions) File.Copy(prefix + "." + e, Path.Combine(dir, "out." + e), true);
        string input = Path.Combine(Path.GetDirectoryName(prefix)!, "in.txt");
        if (File.Exists(input)) File.Copy(input, Path.Combine(dir, "in.txt"), true);
        return Path.Combine(dir, "out");
    }

    internal static VerifyOptions Options(string prefix, int threads = 2, int blockRanks = 5, long cap = 65536, bool input = true, bool child = true, bool hash = true, long stackCap = 1L << 26)
        => new()
        {
            Prefix = prefix,
            InputPath = input ? Path.Combine(Path.GetDirectoryName(prefix)!, "in.txt") : null,
            Threads = threads,
            BlockRanks = blockRanks,
            LceCap = cap,
            Child = child,
            Hash = hash,
            StackCapEntries = stackCap,
            ScratchDir = Path.GetDirectoryName(prefix),
        };

    internal static VerifyResult Check(string prefix, int threads = 2, int blockRanks = 5, long cap = 65536, bool input = true) => OutputVerifier.Run(Options(prefix, threads, blockRanks, cap, input));

    internal static void FlipByte(string path, long offset, byte mask)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.ReadWrite);
        fs.Position = offset;
        int b = fs.ReadByte();
        fs.Position = offset;
        fs.WriteByte((byte)(b ^ mask));
    }

    internal static uint ReadUInt(string path, long index)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read);
        fs.Position = index * 4;
        Span<byte> b = stackalloc byte[4];
        fs.ReadExactly(b);
        return BitConverter.ToUInt32(b);
    }

    internal static void WriteUInt(string path, long index, uint value)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.ReadWrite);
        fs.Position = index * 4;
        fs.Write(BitConverter.GetBytes(value));
    }

    internal static void Truncate(string path, long length)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Write);
        fs.SetLength(length);
    }

    /// <summary>Rewrites one <c>key=value</c> line of a statistics or metadata file.</summary>
    internal static void SetKey(string path, string key, string value)
    {
        string[] lines = File.ReadAllLines(path);
        for (int i = 0; i < lines.Length; i++) if (lines[i].StartsWith(key + "=", StringComparison.Ordinal)) lines[i] = key + "=" + value;
        File.WriteAllText(path, string.Join('\n', lines) + "\n");
    }

    internal static string Key(string path, string key) => File.ReadAllLines(path).First(l => l.StartsWith(key + "=", StringComparison.Ordinal))[(key.Length + 1)..];

    internal static string Describe(VerifyResult r) => string.Join("; ", r.Checks.Select(c => $"{c.Name}={c.Status}({c.Detail})"));

    /// <summary>A text whose compact LCP is U8 with an overflow list (a 300-symbol block twice in random DNA).</summary>
    internal static string SmallWithEscapes()
    {
        string block = global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(300, 21);
        return global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(1500, 22) + block + global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(60, 23) + block;
    }
}
