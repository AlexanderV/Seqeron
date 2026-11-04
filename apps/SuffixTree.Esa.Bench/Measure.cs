using System.Buffers.Binary;
using System.Diagnostics;
using System.Globalization;
using System.Security.Cryptography;
using SuffixTree.Esa;

namespace SuffixTree.Esa.Bench;

/// <summary>One build of the reference pipeline with the stopwatch ticks of its three steps.</summary>
internal sealed class BuildResult
{
    internal required int[] Sa { get; init; }

    internal required CompactLcp Lcp { get; init; }

    internal required PreparedText Prepared { get; init; }

    internal long PrepareTicks { get; init; }

    internal long SaisTicks { get; init; }

    internal long LcpTicks { get; init; }

    internal double PrepareSeconds => PrepareTicks / (double)Stopwatch.Frequency;

    internal double SaisSeconds => SaisTicks / (double)Stopwatch.Frequency;

    internal double LcpSeconds => LcpTicks / (double)Stopwatch.Frequency;

    internal double SaSeconds => PrepareSeconds + SaisSeconds;

    internal double TotalSeconds => SaSeconds + LcpSeconds;
}

internal static class Pipeline
{
    /// <summary>Text preparation (dense order-preserving remap), <see cref="SaisReference"/>, <see cref="LcpReference"/>; the optional traces add one timestamp per phase.</summary>
    internal static BuildResult Build(string text, SaisTrace? saisTrace = null, LcpTrace? lcpTrace = null, SaisKernelKind saisKernel = SaisKernelKind.Reference, LcpKernelKind lcpKernel = LcpKernelKind.Reference)
    {
        long t0 = Stopwatch.GetTimestamp();
        PreparedText prepared = TextPreparation.Prepare(text);
        long t1 = Stopwatch.GetTimestamp();
        int[] sa;
        long t2;
        CompactLcp lcp;
        if (prepared.Width == SymbolWidth.Byte)
        {
            ByteSymbols reader = prepared.ByteReader();
            sa = SaisKernel.Build(saisKernel, prepared, saisTrace);
            t2 = Stopwatch.GetTimestamp();
            lcp = lcpKernel == LcpKernelKind.Pointer ? LcpKernel.Build(prepared, sa, lcpTrace) : LcpReference.Build(reader, sa, lcpTrace);
        }
        else
        {
            UShortSymbols reader = prepared.UShortReader();
            sa = SaisKernel.Build(saisKernel, prepared, saisTrace);
            t2 = Stopwatch.GetTimestamp();
            lcp = lcpKernel == LcpKernelKind.Pointer ? LcpKernel.Build(prepared, sa, lcpTrace) : LcpReference.Build(reader, sa, lcpTrace);
        }
        long t3 = Stopwatch.GetTimestamp();
        return new BuildResult { Sa = sa, Lcp = lcp, Prepared = prepared, PrepareTicks = t1 - t0, SaisTicks = t2 - t1, LcpTicks = t3 - t2 };
    }

    /// <summary>Loads a file whose bytes are the text (Latin-1 widening, no intermediate byte[] of the whole file).</summary>
    internal static string LoadText(string path)
    {
        using var stream = File.OpenRead(path);
        return string.Create(checked((int)stream.Length), stream, static (span, s) =>
        {
            var buffer = new byte[1 << 16];
            int pos = 0, read;
            while ((read = s.Read(buffer)) > 0)
                for (int i = 0; i < read; i++) span[pos++] = (char)buffer[i];
        });
    }
}

internal static class Hashing
{
    /// <summary>SHA-256 of <c>sa[1..n]</c> as little-endian int32 (the libsais SA; see analysis/sa_vectors.tsv).</summary>
    internal static string Sa(ReadOnlySpan<int> sa)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        Append(hash, sa.Slice(1));
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }

    /// <summary>SHA-256 of <c>lcp[1..n]</c> as little-endian int32 (the libsais LCP).</summary>
    internal static string Lcp(CompactLcp lcp, int n)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        var block = new int[4096];
        var cursor = lcp.CreateCursor();
        for (int r = 1; r <= n;)
        {
            int k = Math.Min(block.Length, n - r + 1);
            for (int j = 0; j < k; j++) block[j] = cursor.Get(r + j);
            Append(hash, block.AsSpan(0, k));
            r += k;
        }
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }

    private static void Append(IncrementalHash hash, ReadOnlySpan<int> values)
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

    /// <summary>The expected hashes of corpus <paramref name="name"/> in <c>sa_vectors.tsv</c>, or null when it is not listed.</summary>
    internal static (string Sa, string Lcp)? ExpectedVector(string manifestPath, string name)
    {
        foreach (string line in File.ReadLines(manifestPath))
        {
            if (line.StartsWith('#')) continue;
            string[] f = line.Split('\t');
            if (f.Length == 6 && f[0] == name) return (f[4], f[5]);
        }
        return null;
    }
}

internal static class Process
{
    /// <summary>A kB field of <c>/proc/self/status</c> in bytes (VmHWM, VmRSS), 0 when unavailable.</summary>
    internal static long StatusBytes(string key)
    {
        try
        {
            foreach (string line in File.ReadLines("/proc/self/status"))
            {
                if (!line.StartsWith(key + ":", StringComparison.Ordinal)) continue;
                string digits = new(line.Where(char.IsAsciiDigit).ToArray());
                return long.Parse(digits, CultureInfo.InvariantCulture) * 1024;
            }
        }
        catch (IOException)
        {
            // not Linux: the row is then reported as 0
        }
        return 0;
    }
}

internal static class Emit
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    internal static void Row(string op, string metric, double value, string unit, string notes = "") =>
        Console.WriteLine(string.Join('\t', "R", op, metric, value.ToString("R", Inv), unit, notes));

    internal static void Row(string op, string metric, long value, string unit, string notes = "") =>
        Console.WriteLine(string.Join('\t', "R", op, metric, value.ToString(Inv), unit, notes));

    internal static void Row(string op, string metric, string value, string unit, string notes = "") =>
        Console.WriteLine(string.Join('\t', "R", op, metric, value, unit, notes));
}
