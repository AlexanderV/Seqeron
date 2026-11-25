using System.Security.Cryptography;
using System.Text;
using System.Text.RegularExpressions;

namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The source scan of work package S2.4 (ADR-005 section 3.3: "long in registers"; the rule of the package: no new unchecked narrowing anywhere): in the files of the search side, a cast to <c>int</c> or <c>uint</c> that is not inside a
/// <c>checked(...)</c> expression is a place where a 64-bit rank, count, position or length of a text above 2^31 symbols can wrap silently. Every such cast must be on the allow-list (<see cref="NarrowingLintTests"/>) with the reason why
/// it cannot wrap (an index into a span or an array, a value below an <see cref="int"/> that the line above proves, an <see cref="int"/> that the 32-bit loop of S1 reads); a new cast fails the test, and so does an allow-list entry that
/// no longer matches (the list cannot rot). The 32-bit loop of <c>MappedSearch</c> (the descent of S1, whose casts are those of the code that was gated) is not scanned: it is pinned by a hash of its text instead
/// (<see cref="PinnedLoop"/>), so that "the S1 hot path is untouched" is a test and not a sentence.
/// </summary>
internal static partial class NarrowingLint
{
    /// <summary>One cast: where, to what type, of what.</summary>
    internal sealed record Cast(string File, int Line, string Type, string Operand)
    {
        /// <summary>The key of the allow-list: the file and the cast with its operand (whitespace normalised).</summary>
        internal string Key => $"{File}|({Type}){Operand}";
    }

    [GeneratedRegex(@"\(\s*(int|uint)\s*\)")]
    private static partial Regex CastRegex();

    [GeneratedRegex(@"\bchecked\s*\(")]
    private static partial Regex CheckedRegex();

    /// <summary>The source with comments, string and character literals and <c>checked(...)</c> groups blanked (newlines kept, so that line numbers stay).</summary>
    internal static string Strip(string source)
    {
        var sb = new StringBuilder(source.Length);
        int i = 0;
        while (i < source.Length)
        {
            char c = source[i];
            if (c == '/' && i + 1 < source.Length && source[i + 1] == '/')
            {
                while (i < source.Length && source[i] != '\n') i++;
            }
            else if (c == '/' && i + 1 < source.Length && source[i + 1] == '*')
            {
                int end = source.IndexOf("*/", i + 2, StringComparison.Ordinal);
                int stop = end < 0 ? source.Length : end + 2;
                for (; i < stop; i++) sb.Append(source[i] == '\n' ? '\n' : ' ');
            }
            else if (c is '"' or '\'')
            {
                bool verbatim = i > 0 && source[i - 1] == '@';
                int j = i + 1;
                while (j < source.Length)
                {
                    if (!verbatim && source[j] == '\\') j += 2;
                    else if (source[j] == c)
                    {
                        if (verbatim && j + 1 < source.Length && source[j + 1] == c) j += 2;
                        else break;
                    }
                    else j++;
                }
                int stop = Math.Min(j + 1, source.Length);
                for (; i < stop; i++) sb.Append(source[i] == '\n' ? '\n' : ' ');
            }
            else
            {
                sb.Append(c);
                i++;
            }
        }
        string text = sb.ToString();
        // the checked(...) groups: a narrowing inside one of them throws on overflow, it never wraps
        var result = new StringBuilder(text);
        foreach (Match m in CheckedRegex().Matches(text))
        {
            int open = m.Index + m.Length - 1;
            int depth = 0;
            for (int k = open; k < text.Length; k++)
            {
                if (text[k] == '(') depth++;
                else if (text[k] == ')' && --depth == 0)
                {
                    for (int z = m.Index; z <= k; z++) if (result[z] != '\n') result[z] = ' ';
                    break;
                }
            }
        }
        return result.ToString();
    }

    /// <summary>Every cast to <c>int</c> or <c>uint</c> of a source text (not a <c>typeof(int)</c>), in file order.</summary>
    internal static IEnumerable<Cast> Find(string file, string source)
    {
        string text = Strip(source);
        foreach (Match m in CastRegex().Matches(text))
        {
            if (m.Index > 0 && (char.IsLetterOrDigit(text[m.Index - 1]) || text[m.Index - 1] is '_' or '>')) continue; // typeof(int), sizeof(int), default(int), a generic call
            string operand = Operand(text, m.Index + m.Length);
            if (operand.Length == 0) continue;
            int line = 1 + text.AsSpan(0, m.Index).Count('\n');
            yield return new Cast(file, line, m.Groups[1].Value, operand);
        }
    }

    /// <summary>The operand of the cast that ends at <paramref name="at"/>: a unary operator, then a parenthesised group or a name with its member accesses, calls and indices.</summary>
    private static string Operand(string s, int at)
    {
        int i = at;
        while (i < s.Length && char.IsWhiteSpace(s[i])) i++;
        int start = i;
        while (i < s.Length && s[i] is '-' or '!' or '~') i++;
        if (i < s.Length && s[i] == '(')
        {
            i = Group(s, i);
        }
        else
        {
            while (i < s.Length && (char.IsLetterOrDigit(s[i]) || s[i] is '_' or '.')) i++;
        }
        while (i < s.Length)
        {
            if (s[i] is '(' or '[') i = Group(s, i);
            else if (s[i] == '.' || (s[i] == '?' && i + 1 < s.Length && s[i + 1] == '.'))
            {
                i += s[i] == '?' ? 2 : 1;
                while (i < s.Length && (char.IsLetterOrDigit(s[i]) || s[i] is '_')) i++;
            }
            else break;
        }
        return Regex.Replace(s[start..i], @"\s+", " ").Trim();
    }

    private static int Group(string s, int open)
    {
        int depth = 0;
        for (int k = open; k < s.Length; k++)
        {
            if (s[k] is '(' or '[') depth++;
            else if (s[k] is ')' or ']' && --depth == 0) return k + 1;
        }
        return s.Length;
    }

    /// <summary>The root of the repository (the directory of <c>Seqeron.sln</c>).</summary>
    internal static string RepositoryRoot() => SpanLint.RepositoryRoot();

    /// <summary>The production directory that holds the files of the search side.</summary>
    internal static string EsaDirectory() => Path.Combine(RepositoryRoot(), "src", "SuffixTree", "Algorithms", "SuffixTree.Esa");

    /// <summary>The files of the search side (S2.4) and of the navigation side (S2.5) that the lint scans, relative to <see cref="EsaDirectory"/>.</summary>
    internal static readonly string[] Files =
    [
        "WideSearch.cs", "MappedSearch.cs", "EsaIndex.cs", "EsaIndex.Search.cs", "EsaIndex.Analysis.cs", "EsaIndex.Traverse.cs", "EsaIndex.Matching.cs", "EsaIndex.Navigation.cs",
        "LeftmostOccurrence.cs", "LcpIntervalScan.cs", "ChildTable.cs", "EsaSuffixIndex.cs", "EsaSuffixIndex.Arrays.cs", "EsaSuffixIndex.Tree.cs", "External/Sesa/SesaFile.cs", "IntBridge.cs",
        // the navigation side (S2.5): the navigator and its handles, the pack, the pyramid, the views that they read
        "EsaNavigator.cs", "EsaNavigationPack.cs", "EsaNode.cs", "LcpPyramid.cs", "SaView.cs", "CompactLcp.cs", "ChildLayouts.cs", "EsaText.cs", "EsaSuffixTree.cs", "PosReaders.cs",
    ];

    // ------------------------------------------------------------------------------------ the 32-bit loop of MappedSearch

    private const string ReadersStart = "/// <summary>One LCP entry of a container by rank";
    private const string ReadersEnd = "/// <summary>\n/// The pattern descent over the sections of a container (work package 1b.9a)";
    private const string LoopStart = "    // No AggressiveOptimization (work package 1b.9g)";
    private const string SymbolDoc = "    /// <summary>\n    /// The dense symbol at <paramref name=\"i\"/> (inside the text).";
    private const string EqualDoc = "    /// <summary>\n    /// Whether the text symbols <c>[start, start + length)</c> (inside the text) equal the pattern symbols";
    private const string ThrowDoc = "    [MethodImpl(MethodImplOptions.NoInlining)]\n    internal static void ThrowDamaged()";

    /// <summary>
    /// The text of <c>MappedSearch.cs</c> split into the 32-bit loop of S1 (the readers and cells of the codecs, <c>Run</c> with <c>Window</c>, and <c>Equal</c> with the packed-pattern helpers: the code that the A/B gates measured) and the rest
    /// (the constructor, the overloads of <c>TryMatch</c>, <c>Symbol</c>, <c>ThrowDamaged</c>).
    /// </summary>
    internal static (string Pinned, string Other) SplitMappedSearch(string source)
    {
        string text = source.Replace("\r\n", "\n", StringComparison.Ordinal);
        int a = Find(text, ReadersStart), b = Find(text, ReadersEnd), c = Find(text, LoopStart), d = Find(text, SymbolDoc), e = Find(text, EqualDoc), f = Find(text, ThrowDoc);
        if (!(a < b && b < c && c < d && d < e && e < f)) throw new InvalidOperationException("the anchors of MappedSearch.cs are not in the expected order");
        string pinned = text[a..b] + "\n" + text[c..d] + "\n" + text[e..f];
        string rest = text[..a] + text[b..c] + text[d..e] + text[f..];
        return (pinned, rest);

        static int Find(string haystack, string needle)
        {
            int at = haystack.IndexOf(needle, StringComparison.Ordinal);
            return at >= 0 ? at : throw new InvalidOperationException($"the anchor '{needle[..Math.Min(40, needle.Length)]}' is not in MappedSearch.cs");
        }
    }

    /// <summary>The SHA-256 (lower-case hex) of a pinned text.</summary>
    internal static string Hash(string pinned) => Convert.ToHexStringLower(SHA256.HashData(Encoding.UTF8.GetBytes(pinned)));

    /// <summary>The hash of the 32-bit loop of S1 as it is at the head of S2.3 (<c>0bfb2114</c>, unchanged since S1): <see cref="SplitMappedSearch"/> of <c>MappedSearch.cs</c> there.</summary>
    internal const string PinnedLoop = "18daf01fb88dffd9a0b460e316e93c75dd9b51d02c0e3aae287ccb1a4a2b8e73";
}
