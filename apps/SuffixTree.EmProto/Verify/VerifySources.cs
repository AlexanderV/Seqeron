using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using SuffixTree.Esa;
using SuffixTree.Esa.External;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Verify;

/// <summary>
/// Where <see cref="OutputVerifier"/> gets the seven outputs of an end-to-end run, by extension (<c>pk</c>, <c>pk.meta</c>, <c>sa</c>, <c>lcp</c>, <c>ovf</c>, <c>cld</c>, <c>stats</c>): separate files that carry a prefix
/// (the outputs of <c>e6</c>, <see cref="FileVerifySource"/>) or the sections of a SESA container, read in place from its mapping (<see cref="ContainerVerifySource"/>, work package 1b.7f): the checks are the same eleven.
/// </summary>
internal interface IVerifySource
{
    /// <summary>A short name for messages and for the scratch file of the permutation bitmap (the prefix, or the container's file name).</summary>
    string Label { get; }

    /// <summary>The directory that holds the scratch file by default.</summary>
    string DefaultScratchDirectory { get; }

    /// <summary>Whether the output exists.</summary>
    bool Exists(string extension);

    /// <summary>Bytes of the output.</summary>
    long Length(string extension);

    /// <summary>The output as a block of memory (a mapped file, or a section of the container's mapping: nothing is copied).</summary>
    IMemoryBlock Open(string extension);

    /// <summary>The text of <c>pk.meta</c> or <c>stats</c>.</summary>
    string ReadText(string extension);

    /// <summary>The SHA-256 of the output (lower-case hexadecimal).</summary>
    string Hash(string extension, int bufferBytes);

    /// <summary>What the SHA-256 of the extension is the hash of, for the message of the row.</summary>
    string Describe(string extension);

    /// <summary>Whether the child table holds tagged cells (CHILD codec 2) when the source knows it (a container does); null = the option of the run decides. STUB of the red commit of work package 1b.9a.</summary>
    bool? ChildTagged => null;
}

/// <summary>The files <c>&lt;prefix&gt;.pk</c>, <c>.pk.meta</c>, <c>.sa</c>, <c>.lcp</c>, <c>.ovf</c>, <c>.cld</c>, <c>.stats</c> of an <c>e6</c> run.</summary>
internal sealed class FileVerifySource(string prefix) : IVerifySource
{
    private string P(string extension) => prefix + "." + extension;

    public string Label => prefix;

    public string DefaultScratchDirectory => Path.GetDirectoryName(Path.GetFullPath(prefix))!;

    public bool Exists(string extension) => File.Exists(P(extension));

    public long Length(string extension) => new FileInfo(P(extension)).Length;

    public IMemoryBlock Open(string extension) => new MappedFile(P(extension), writable: false);

    public string ReadText(string extension) => MetaFiles.ReadText(P(extension));

    public string Hash(string extension, int bufferBytes) => OutputVerifier.HashFile(P(extension), bufferBytes);

    public string Describe(string extension) => $"{Path.GetFileName(prefix)}.{extension}";
}

/// <summary>
/// A SESA container as the seven outputs: <c>pk</c> is the TEXT section, <c>sa</c> the SA, <c>lcp</c> the LCP, <c>ovf</c> the LCP_OVERFLOW list (empty when there is none), <c>cld</c> the CHILD section; <c>pk.meta</c> and
/// <c>stats</c> are written out as text from the header, the ALPHABET and the COUNTS (the format of the prototype's files), so that the same checks run on them. The container is mapped once, as one view, by
/// <see cref="SesaFile"/> at its header level; every section is a pointer into that mapping and nothing is extracted. The header level has checked the table, so the sections lie inside the file (the CRC-32C of the payloads
/// is <c>SuffixIndexFactory.Open</c> at the <c>Sections</c> level's business; this verifier compares the content itself).
/// </summary>
internal sealed unsafe class ContainerVerifySource : IVerifySource, IDisposable
{
    private readonly SesaFile _file;
    private readonly string _path;

    internal ContainerVerifySource(string path)
    {
        _path = path;
        _file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Header });
    }

    /// <summary>The mapped container.</summary>
    internal SesaFile File => _file;

    public string Label => Path.GetFileName(_path);

    public string DefaultScratchDirectory => Path.GetDirectoryName(Path.GetFullPath(_path))!;

    public bool Exists(string extension) => true;

    private SesaSection? Section(string extension) => extension switch
    {
        "pk" => _file.Parsed.Text,
        "sa" => _file.Parsed.Sa,
        "lcp" => _file.Parsed.Lcp,
        "ovf" => _file.Parsed.Overflow,
        "cld" => _file.Parsed.Child,
        _ => null,
    };

    public long Length(string extension) => extension switch
    {
        "pk.meta" => Encoding.Latin1.GetByteCount(ReadText("pk.meta")),
        "stats" => Encoding.Latin1.GetByteCount(ReadText("stats")),
        _ => Section(extension)?.Length ?? 0,
    };

    public IMemoryBlock Open(string extension)
    {
        SesaSection section = Section(extension) ?? throw new InvalidDataException($"the container has no section for '.{extension}'");
        return new SectionBlock(_file.Data + section.Offset, section.Length);
    }

    public string ReadText(string extension)
    {
        SesaParsed p = _file.Parsed;
        SesaHeader h = p.Header;
        switch (extension)
        {
            case "pk.meta":
                {
                    string letters = Convert.ToHexString(p.Letters.Select(c => (byte)c).ToArray());
                    long[] counts = p.Counts64?.Select(c => (long)c).ToArray() ?? RecountSymbols(p);
                    string sha = Convert.ToHexStringLower(h.TextSha256);
                    return string.Create(CultureInfo.InvariantCulture, $"n={h.TextLength}\nbits={p.Bits}\nsigma={h.Sigma}\nletters_hex={letters}\ncounts={string.Join(',', counts)}\nsha256={sha}\nname={Label}\n");
                }
            case "stats":
                {
                    long n = h.TextLength;
                    long lrsRank = h.LrsDepth == 0 ? 0 : h.LrsRank;
                    // the container has no MAX_LCP field: the longest repeat is the deepest internal node, so the two are the same number (LRS_DEPTH)
                    var stats = new StatsInfo(h.NodeCount, h.NodeCount - 1 - (n + 1), h.LrsDepth, lrsRank, n, h.Sigma, p.Lcp.ElemWidth, p.Lcp.Param0, h.LrsDepth);
                    return MetaFiles.StatsText(stats);
                }
            default:
                throw new ArgumentOutOfRangeException(nameof(extension), extension, "not a text output");
        }
    }

    public string Hash(string extension, int bufferBytes)
    {
        using var hash = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        if (extension is "pk.meta" or "stats")
        {
            hash.AppendData(Encoding.Latin1.GetBytes(ReadText(extension)));
        }
        else if (Section(extension) is { } section)
        {
            const long piece = 1L << 30;
            for (long done = 0; done < section.Length; done += piece)
                hash.AppendData(new ReadOnlySpan<byte>(_file.Data + section.Offset + done, (int)Math.Min(piece, section.Length - done)));
        }
        return Convert.ToHexStringLower(hash.GetHashAndReset());
    }

    public string Describe(string extension) => extension is "pk.meta" or "stats" ? $"the {extension} text derived from the header of {Label}" : $"the payload of the section of {Label} that plays '.{extension}'";

    public void Dispose() => _file.Dispose();

    /// <summary>A container without a COUNTS section: the counts are a recount of the text (so the check of the counts against the metadata can only find damage to the text itself).</summary>
    private long[] RecountSymbols(SesaParsed p)
    {
        var counts = new long[p.Header.Sigma];
        MappedTextReader reader = _file.TextReader;
        var block = new byte[1 << 16];
        for (long at = 0; at < p.Header.TextLength; at += block.Length)
        {
            int take = (int)Math.Min(block.Length, p.Header.TextLength - at);
            reader.ReadSymbols(at, block.AsSpan(0, take));
            for (int k = 0; k < take; k++) if (block[k] < counts.Length) counts[block[k]]++;
        }
        return counts;
    }

    /// <summary>A section of the container's mapping as a block; it owns nothing.</summary>
    private sealed class SectionBlock(byte* pointer, long bytes) : IMemoryBlock
    {
        public byte* Pointer => pointer;

        public long Bytes => bytes;

        public void Dispose()
        {
        }
    }
}
