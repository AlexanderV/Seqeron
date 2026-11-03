using System.Security.Cryptography;

namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// Reader for <c>sha256sum</c>-style manifests (<c>MANIFEST.sha256</c> of the benchmark inputs; the same format is used for the
/// libsais vector manifest of 1a.6): each line is <c>&lt;64 hex digits&gt;  &lt;relative path&gt;</c>, a <c>*</c> before the path
/// (binary mode) is accepted, blank lines and <c>#</c> comment lines are ignored.
/// </summary>
public sealed class CorpusManifest
{
    /// <summary>One manifest line.</summary>
    public readonly record struct Entry(string Sha256, string Path);

    private readonly Dictionary<string, Entry> _byPath;

    private CorpusManifest(List<Entry> entries)
    {
        Entries = entries;
        _byPath = entries.ToDictionary(e => e.Path, StringComparer.Ordinal);
    }

    /// <summary>The entries in file order.</summary>
    public IReadOnlyList<Entry> Entries { get; }

    /// <summary>Parses manifest text; throws <see cref="FormatException"/> on a malformed line or a duplicate path.</summary>
    public static CorpusManifest Parse(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        var entries = new List<Entry>();
        var seen = new HashSet<string>(StringComparer.Ordinal);
        int lineNo = 0;
        foreach (string raw in text.Split('\n'))
        {
            lineNo++;
            string line = raw.TrimEnd('\r');
            if (line.Trim().Length == 0 || line.StartsWith('#')) continue;
            if (line.Length < 66 || !IsHex(line.AsSpan(0, 64)) || line[64] != ' ')
                throw new FormatException($"manifest line {lineNo}: expected '<64 hex digits>  <path>'");
            string path = line.Substring(64).TrimStart(' ');
            if (path.StartsWith('*')) path = path.Substring(1);
            if (path.Length == 0)
                throw new FormatException($"manifest line {lineNo}: missing path");
            if (!seen.Add(path))
                throw new FormatException($"manifest line {lineNo}: duplicate path '{path}'");
            entries.Add(new Entry(line.Substring(0, 64).ToLowerInvariant(), path));
        }
        return new CorpusManifest(entries);
    }

    /// <summary>Reads and parses a manifest file.</summary>
    public static CorpusManifest Load(string filePath) => Parse(File.ReadAllText(filePath));

    /// <summary>The repository root: the nearest ancestor of the test binaries that contains <c>Seqeron.sln</c>.</summary>
    public static string RepositoryRoot()
    {
        for (var dir = new DirectoryInfo(AppContext.BaseDirectory); dir is not null; dir = dir.Parent)
            if (File.Exists(System.IO.Path.Combine(dir.FullName, "Seqeron.sln"))) return dir.FullName;
        throw new DirectoryNotFoundException("Seqeron.sln not found above " + AppContext.BaseDirectory);
    }

    /// <summary>Directory of the committed benchmark manifest and (after regeneration) its <c>data/</c> and <c>queries/</c> folders.</summary>
    public static string BenchDirectory() =>
        System.IO.Path.Combine(RepositoryRoot(), "docs", "refactoring", "suffix-tree-upgrade", "analysis", "bench");

    /// <summary>Loads <c>analysis/bench/MANIFEST.sha256</c> from the repository.</summary>
    public static CorpusManifest LoadBenchManifest() => Load(System.IO.Path.Combine(BenchDirectory(), "MANIFEST.sha256"));

    /// <summary>Looks an entry up by its relative path.</summary>
    public bool TryGet(string path, out Entry entry) => _byPath.TryGetValue(path, out entry);

    /// <summary>
    /// True when the file <c>baseDirectory/path</c> has the listed SHA-256. Throws <see cref="KeyNotFoundException"/> for an
    /// unlisted path and <see cref="FileNotFoundException"/> for a missing file (a missing input is a failure, never a skip).
    /// </summary>
    public bool Verify(string path, string baseDirectory)
    {
        if (!_byPath.TryGetValue(path, out var entry)) throw new KeyNotFoundException($"'{path}' is not listed in the manifest");
        string full = System.IO.Path.Combine(baseDirectory, path);
        if (!File.Exists(full)) throw new FileNotFoundException("corpus file missing", full);
        return string.Equals(ComputeSha256Hex(full), entry.Sha256, StringComparison.Ordinal);
    }

    /// <summary>Lower-case hex SHA-256 of a file.</summary>
    public static string ComputeSha256Hex(string filePath)
    {
        using var stream = File.OpenRead(filePath);
        return Convert.ToHexStringLower(SHA256.HashData(stream));
    }

    private static bool IsHex(ReadOnlySpan<char> s)
    {
        foreach (char c in s)
            if (!char.IsAsciiHexDigit(c)) return false;
        return true;
    }
}
