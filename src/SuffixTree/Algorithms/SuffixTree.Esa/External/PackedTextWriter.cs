using System.Security.Cryptography;

namespace SuffixTree.Esa.External;

/// <summary>Options of <see cref="PackedTextWriter"/>.</summary>
internal sealed class PackedTextWriterOptions
{
    /// <summary>Symbols per block of the two passes (a positive multiple of <see cref="PackedTextWriter.BlockMultiple"/>, so that every packed block ends on a byte boundary).</summary>
    internal int BlockSymbols { get; init; } = 1 << 20;

    /// <summary>Window of the mapped output (at least 1 MiB by default: the .NET view accessor flushes every window when it is disposed, ADR-004 section 5a, E9).</summary>
    internal int WindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    /// <summary>Smallest width b; the width is <c>max(MinBits, ceil(log2 sigma))</c>. The default is the width of the prototype, 1.</summary>
    internal int MinBits { get; init; } = PackedTextFormat.MinBits;

    /// <summary>Longest input accepted; the default is the design limit of the 32-bit positions, <see cref="PackedTextFormat.MaxSymbols"/>.</summary>
    internal long MaxSymbols { get; init; } = PackedTextFormat.MaxSymbols;

    /// <summary>Counters of the input reads and of the mapped writes (default none).</summary>
    internal IoCounters? Counters { get; init; }

    /// <summary>How the output file is created (reservation of the disk space).</summary>
    internal MappedFileOptions? File { get; init; }

    /// <summary>Whether the finished file is forced to storage before it is renamed over the target (default true: the atomic finish of ADR-003 section 6).</summary>
    internal bool Fsync { get; init; } = true;

    /// <summary>Called between the two passes (tests change the input there).</summary>
    internal Action? BetweenPasses { get; init; }
}

/// <summary>
/// The streaming writer of the packed text, the part of PH0 of ADR-004 section 4.2 that the codec needs: two sequential passes over the input in blocks of a fixed number of symbols (the width b is known only after the first):
/// the first reads, counts and hashes (<see cref="Scan"/>); the second reads again, maps the bytes to dense symbols (increasing byte order), packs one block at a time with <see cref="PackedTextFormat.PackInto"/> and appends
/// it to a <see cref="MappedWindowWriter"/> whose range is the payload (<see cref="Pack"/>, so the container writer can aim it at the TEXT section of its file), and ends with the 64 zero bytes of padding.
/// <see cref="WriteFile"/> does both into a file of its own, <c>&lt;target&gt;.tmp</c>, and renames it over the target when it is complete. Memory: one block of the input, one of symbols, one of packed bytes. The stage
/// journal, the resume and the rest of PH0 belong to the builder (1b.7e).
/// </summary>
internal static class PackedTextWriter
{
    /// <summary>Symbols per block that keep the packed blocks byte aligned (8 symbols of b bits are b bytes).</summary>
    internal const int BlockMultiple = 8;

    /// <summary>The length of the input in bytes, through an open handle (so that a symbolic link to the text counts the text, not the link: <c>FileInfo.Length</c> of a link is the length of its target path).</summary>
    internal static long InputLength(string inputPath)
    {
        using var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read);
        return input.Length;
    }

    /// <summary>Pass 1: the length, the alphabet, the symbol counts, the width and the SHA-256 of the input bytes. Creates nothing.</summary>
    internal static PackedTextInfo Scan(string inputPath, PackedTextWriterOptions? options = null)
    {
        options ??= new PackedTextWriterOptions();
        CheckOptions(options);
        long n = InputLength(inputPath);
        if (n > options.MaxSymbols) throw new NotSupportedException($"the input has {n} symbols, above the limit of {options.MaxSymbols} (the positions of the builder are 32-bit: ADR-004 section 4.3)");
        var counts = new long[256];
        var block = new byte[options.BlockSymbols];
        string sha;
        using (var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read, options.Counters))
        using (var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256))
        {
            for (long at = 0; at < n; at += block.Length)
            {
                int take = (int)Math.Min(block.Length, n - at);
                input.ReadExactAt(block.AsSpan(0, take), at);
                foreach (byte b in block.AsSpan(0, take)) counts[b]++;
                hash.AppendData(block.AsSpan(0, take));
            }
            sha = Convert.ToHexStringLower(hash.GetHashAndReset());
        }
        var letters = new List<byte>();
        for (int v = 0; v < 256; v++)
        {
            if (counts[v] != 0) letters.Add((byte)v);
        }
        int bits = Math.Max(options.MinBits, PackedTextFormat.BitsFor(letters.Count));
        return new PackedTextInfo(n, bits, [.. letters], [.. letters.Select(l => counts[l])], sha);
    }

    /// <summary>
    /// Pass 2: packs the input into <paramref name="target"/>, which must be at its start and have exactly <see cref="PackedTextInfo.TotalBytes"/> bytes to go. Every byte of the range is written, the padding included, so a range
    /// that held something else (a section of an older file) is overwritten. An input that changed since <see cref="Scan"/> (length, alphabet or counts) is an <see cref="InvalidOperationException"/>.
    /// </summary>
    internal static void Pack(string inputPath, PackedTextInfo info, MappedWindowWriter target, PackedTextWriterOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(info);
        ArgumentNullException.ThrowIfNull(target);
        options ??= new PackedTextWriterOptions();
        CheckOptions(options);
        long room = target.Length - target.Position;
        if (room != info.TotalBytes) throw new ArgumentException($"the target has {room} bytes to go but {info.N} symbols of {info.Bits} bits need {info.TotalBytes}", nameof(target));
        long n = InputLength(inputPath);
        if (n != info.N) throw new InvalidOperationException($"the input has {n} bytes, {info.N} were scanned: the input changed between the two passes");
        var map = new int[256];
        Array.Fill(map, -1);
        for (int s = 0; s < info.Letters.Length; s++) map[info.Letters[s]] = s;
        int bits = info.Bits;
        var block = new byte[options.BlockSymbols];
        var packed = new byte[(options.BlockSymbols * (long)bits / 8) + 1];
        var seen = new long[info.Letters.Length];
        long start = target.Position;
        using (var input = new CountedFile(inputPath, FileMode.Open, FileAccess.Read, options.Counters))
        {
            for (long at = 0; at < n; at += block.Length)
            {
                int take = (int)Math.Min(block.Length, n - at);
                input.ReadExactAt(block.AsSpan(0, take), at);
                for (int i = 0; i < take; i++)
                {
                    int symbol = map[block[i]];
                    if (symbol < 0) throw new InvalidOperationException($"the byte {block[i]} at position {at + i} was not in the alphabet of the scan: the input changed between the two passes");
                    block[i] = (byte)symbol;
                    seen[symbol]++;
                }
                int bytes = (int)PackedTextFormat.DataBytes(take, bits);
                PackedTextFormat.PackInto(block.AsSpan(0, take), bits, packed);
                target.Append(packed.AsSpan(0, bytes));
            }
        }
        for (int s = 0; s < seen.Length; s++)
        {
            if (seen[s] != info.Counts[s]) throw new InvalidOperationException($"symbol {s} occurred {seen[s]} times in the second pass and {info.Counts[s]} times in the first: the input changed between the two passes");
        }
        if (target.Position - start != info.DataBytes) throw new InvalidOperationException($"the packed text has {target.Position - start} data bytes, expected {info.DataBytes}: the input changed between the two passes");
        target.Append(new byte[PackedTextFormat.PadBytes]);
    }

    /// <summary>Writes the packed text of <paramref name="inputPath"/> to <paramref name="packedPath"/>: the temporary file <c>&lt;packedPath&gt;.tmp</c> is made, filled, forced to storage and renamed over the target; a failure removes it and leaves the target as it was.</summary>
    internal static PackedTextInfo WriteFile(string inputPath, string packedPath, PackedTextWriterOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(inputPath);
        ArgumentNullException.ThrowIfNull(packedPath);
        options ??= new PackedTextWriterOptions();
        string? directory = Path.GetDirectoryName(Path.GetFullPath(packedPath));
        if (directory is not null && !Directory.Exists(directory)) throw new DirectoryNotFoundException($"the directory of {packedPath} does not exist");
        PackedTextInfo info = Scan(inputPath, options);
        options.BetweenPasses?.Invoke();
        string temporary = packedPath + ".tmp";
        MappedOutputFile? file = null;
        try
        {
            file = MappedOutputFile.Create(temporary, info.TotalBytes, options.File);
            using (MappedWindowWriter writer = file.CreateWriter(options.WindowBytes, options.Counters)) Pack(inputPath, info, writer, options);
            file.Complete(options.Fsync);
            file.Dispose();
            file = null;
            File.Move(temporary, packedPath, overwrite: true);
        }
        catch
        {
            file?.Dispose();
            DeleteQuietly(temporary);
            throw;
        }
        return info;
    }

    private static void CheckOptions(PackedTextWriterOptions o)
    {
        if (o.BlockSymbols < BlockMultiple || o.BlockSymbols % BlockMultiple != 0) throw new ArgumentOutOfRangeException(nameof(o), o.BlockSymbols, "a block holds a positive multiple of 8 symbols, so that every packed block ends on a byte boundary");
        ArgumentOutOfRangeException.ThrowIfLessThan(o.MinBits, PackedTextFormat.MinBits);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(o.MinBits, PackedTextFormat.MaxBits);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.WindowBytes, 1);
    }

    private static void DeleteQuietly(string path)
    {
        try
        {
            File.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // the original failure is what the caller needs
        }
    }
}
