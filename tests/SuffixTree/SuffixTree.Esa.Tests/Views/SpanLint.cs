using System.Text.RegularExpressions;

namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The source scan of risk R8 (ADR-005 section 10, work package S2.3): <c>Span&lt;T&gt;</c> and <c>ReadOnlySpan&lt;T&gt;</c> take an <see cref="int"/> length, so a view of a section of more than <c>int.MaxValue</c> bytes that is made from a pointer with
/// a narrowed length is the classic failure of a 64-bit port (a wrapped length reads the wrong bytes, or none). The scan finds every <c>new Span&lt;</c>, <c>new ReadOnlySpan&lt;</c>, <c>AsSpan(start, length)</c> and
/// <c>MemoryMarshal.CreateSpan</c> / <c>CreateReadOnlySpan</c> whose LENGTH argument is a cast to <c>int</c> or <c>uint</c>, a <c>checked(...)</c> expression or an expression that names a section; each such site must be on the allow-list with the reason why
/// its length is bounded (a buffer of a few KB, a chunk of a pass) or why it refuses a length that does not fit; a new site fails, and so does an allow-list entry that no longer matches (the list cannot rot).
/// </summary>
internal static partial class SpanLint
{
    /// <summary>One call that creates a span: where, what and its length argument (whitespace normalised).</summary>
    internal sealed record Site(string File, int Line, string Call, string Length)
    {
        /// <summary>The key of the allow-list: the file, and the call with its length expression.</summary>
        internal string Key => $"{File}|{Call}|{Length}";
    }

    [GeneratedRegex(@"new\s+(?:ReadOnlySpan|Span)\s*<[^>(]*>\s*\(|\.AsSpan\s*\(|CreateReadOnlySpan\s*\(|CreateSpan\s*\(")]
    private static partial Regex CreatesASpan();

    [GeneratedRegex(@"\(\s*u?int\s*\)|\bchecked\s*\(|\b[Ss]ection\b|\bSection\.")]
    private static partial Regex Narrows();

    /// <summary>Every span creation of a source text with its length argument, in file order.</summary>
    internal static IEnumerable<Site> Find(string file, string source)
    {
        foreach (Match m in CreatesASpan().Matches(source))
        {
            int open = m.Index + m.Length - 1;
            List<string>? args = Arguments(source, open);
            if (args is null) continue;
            string call = m.Value.TrimEnd('(').Trim();
            call = Regex.Replace(call, @"\s+", " ");
            string? length = LengthOf(call, args);
            if (length is null) continue;
            int line = 1 + source.AsSpan(0, m.Index).Count('\n');
            yield return new Site(file, line, call, Regex.Replace(length, @"\s+", " ").Trim());
        }
    }

    /// <summary>The creations whose length argument narrows to an integer or comes from a section: the ones that must be accounted for.</summary>
    internal static IEnumerable<Site> Hazards(string file, string source) => Find(file, source).Where(s => Narrows().IsMatch(s.Length));

    /// <summary>The length argument of a call, or null when the call has none (a span of a whole array, <c>AsSpan(start)</c>).</summary>
    private static string? LengthOf(string call, List<string> args)
    {
        if (call.StartsWith("new ", StringComparison.Ordinal)) return args.Count switch { 2 => args[1], 3 => args[2], _ => null };
        return args.Count >= 2 ? args[1] : null; // AsSpan(start, length), CreateSpan(ref, length)
    }

    /// <summary>The top-level arguments of the call whose opening parenthesis is at <paramref name="open"/>, or null when it is not closed.</summary>
    private static List<string>? Arguments(string s, int open)
    {
        var args = new List<string>();
        int depth = 0, start = open + 1, i = open;
        while (i < s.Length)
        {
            char c = s[i];
            if (c == '"')
            {
                i = SkipString(s, i) + 1;
                continue;
            }
            if (c is '(' or '[' or '{') depth++;
            else if (c is ')' or ']' or '}')
            {
                depth--;
                if (depth == 0)
                {
                    string last = s[start..i];
                    if (last.Trim().Length > 0 || args.Count > 0) args.Add(last);
                    return args;
                }
            }
            else if (c == ',' && depth == 1)
            {
                args.Add(s[start..i]);
                start = i + 1;
            }
            i++;
        }
        return null;
    }

    private static int SkipString(string s, int at)
    {
        int i = at + 1;
        while (i < s.Length)
        {
            if (s[i] == '\\') i += 2;
            else if (s[i] == '"') return i;
            else i++;
        }
        return s.Length - 1;
    }

    /// <summary>The root of the repository (the directory of <c>Seqeron.sln</c>), found from the test assembly.</summary>
    internal static string RepositoryRoot()
    {
        string? dir = AppContext.BaseDirectory;
        while (dir is not null && !File.Exists(Path.Combine(dir, "Seqeron.sln"))) dir = Path.GetDirectoryName(dir);
        return dir ?? throw new InvalidOperationException("Seqeron.sln was not found above " + AppContext.BaseDirectory);
    }

    /// <summary>The hazards of every source file below <paramref name="directory"/> (bin and obj excluded), the file names relative to it.</summary>
    internal static List<Site> ScanDirectory(string directory)
    {
        var sites = new List<Site>();
        foreach (string path in Directory.EnumerateFiles(directory, "*.cs", SearchOption.AllDirectories).Order(StringComparer.Ordinal))
        {
            string relative = Path.GetRelativePath(directory, path).Replace('\\', '/');
            if (relative.StartsWith("bin/", StringComparison.Ordinal) || relative.StartsWith("obj/", StringComparison.Ordinal)) continue;
            sites.AddRange(Hazards(relative, File.ReadAllText(path)));
        }
        return sites;
    }
}
