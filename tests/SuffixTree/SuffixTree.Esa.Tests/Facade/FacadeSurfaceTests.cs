using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// Every member of <see cref="ISuffixTreeSearch"/>, <see cref="ISuffixTreeAnalysis"/>, <see cref="ISuffixTreeDiagnostics"/> and <see cref="ISuffixTree"/> through the backend factory
/// (work package 1b.6), with positive, absent, empty, single-character, longer-than-text and boundary arguments and the argument exceptions. The fixture runs for the Ukkonen tree and the
/// v6 persistent tree as well: an assertion that one of them fails is not backend-agnostic and does not belong here. Expected values are hand-derived or come from the definitions.
/// </summary>
[TestFixtureSource(typeof(SuffixTreeBackends), nameof(SuffixTreeBackends.AllFixtures))]
public class FacadeSurfaceTests(SuffixTreeBackend backend) : PortedBackendTestBase(backend)
{
    private static readonly string[] Corpus =
    [
        "", "a", "aa", "ab", "ba", "aba", "abab", "banana", "mississippi", "abracadabra", "aaaaaaaaaa", "abababababab", "abcabxabcd", "xabxac", "GATTACAGATTACA",
        "ACGTACGTTACGT", "\uD800a\uD800", "\0￿\0", "a$b$", "🧬αβγ🧪$", "the quick brown fox jumps over the lazy dog",
    ];

    // ------------------------------------------------------------------ header values

    [Test]
    public void Banana_HandDerivedHeader()
    {
        ISuffixTree tree = Tree("banana");
        Assert.Multiple(() =>
        {
            Assert.That(tree.Text.ToString(), Is.EqualTo("banana"));
            Assert.That(tree.Text.Length, Is.EqualTo(6));
            Assert.That(tree.LeafCount, Is.EqualTo(6));
            Assert.That(tree.MaxDepth, Is.EqualTo(6));
            Assert.That(tree.IsEmpty, Is.False);
            // root, 7 leaves (six suffixes and the terminator-only one), internal nodes a, ana, na: 1 + 7 + 3
            Assert.That(tree.NodeCount, Is.EqualTo(11));
        });
    }

    [Test]
    public void EmptyText_HeaderAndEverySearchMember()
    {
        ISuffixTree tree = Tree("");
        Assert.Multiple(() =>
        {
            Assert.That(tree.Text.Length, Is.EqualTo(0));
            Assert.That(tree.IsEmpty, Is.True);
            Assert.That(tree.NodeCount, Is.EqualTo(1));
            Assert.That(tree.LeafCount, Is.EqualTo(0));
            Assert.That(tree.MaxDepth, Is.EqualTo(0));
            Assert.That(tree.Contains(""), Is.True);
            Assert.That(tree.Contains("a"), Is.False);
            Assert.That(tree.CountOccurrences(""), Is.EqualTo(0));
            Assert.That(tree.CountOccurrences("a"), Is.EqualTo(0));
            Assert.That(tree.FindAllOccurrences(""), Is.Empty);
            Assert.That(tree.FindAllOccurrences("a"), Is.Empty);
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(""));
            Assert.That(tree.LongestRepeatedSubstringMemory().Length, Is.EqualTo(0));
            Assert.That(tree.GetAllSuffixes(), Is.Empty);
            Assert.That(tree.EnumerateSuffixes(), Is.Empty);
            Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(0));
            Assert.That(tree.CountDistinctSubstringsByLength(3), Is.EqualTo(new long[] { 1 }));
            Assert.That(tree.FindAllLongestRepeatedSubstrings(), Is.Empty);
            Assert.That(tree.FindMaximalRepeatedPairs(1), Is.Empty);
            Assert.That(tree.LongestCommonSubstring("abc"), Is.EqualTo(""));
            Assert.That(tree.LongestCommonSubstringInfo("abc"), Is.EqualTo(("", -1, -1)));
            Assert.That(tree.FindAllLongestCommonSubstrings("abc").PositionsInText, Is.Empty);
            Assert.That(tree.FindAllDistinctLongestCommonSubstrings("abc"), Is.Empty);
            Assert.That(tree.FindExactMatchAnchors("abc", 1), Is.Empty);
            Assert.That(tree.FindMaximalExactMatches("abc", 1), Is.Empty);
            Assert.That(tree.FindMaximalUniqueMatches("abc", 1), Is.Empty);
        });
    }

    [Test]
    public void SingleCharacter_HeaderAndEverySearchMember()
    {
        ISuffixTree tree = Tree("x");
        Assert.Multiple(() =>
        {
            Assert.That(tree.NodeCount, Is.EqualTo(3), "root, the leaf of x, the terminator-only leaf");
            Assert.That(tree.LeafCount, Is.EqualTo(1));
            Assert.That(tree.MaxDepth, Is.EqualTo(1));
            Assert.That(tree.Contains("x"), Is.True);
            Assert.That(tree.Contains("xx"), Is.False);
            Assert.That(tree.CountOccurrences("x"), Is.EqualTo(1));
            Assert.That(tree.CountOccurrences(""), Is.EqualTo(1));
            Assert.That(tree.FindAllOccurrences("x"), Is.EqualTo(new[] { 0 }));
            Assert.That(tree.FindAllOccurrences(""), Is.EqualTo(new[] { 0 }));
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(""));
            Assert.That(tree.GetAllSuffixes(), Is.EqualTo(new[] { "x" }));
            Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(1));
            Assert.That(tree.CountDistinctSubstringsByLength(5), Is.EqualTo(new long[] { 1, 1 }));
            Assert.That(tree.LongestCommonSubstring("axb"), Is.EqualTo("x"));
            Assert.That(tree.LongestCommonSubstringInfo("axb"), Is.EqualTo(("x", 0, 1)));
            Assert.That(tree.FindMaximalExactMatches("axb", 1), Is.EqualTo(new[] { (0, 1, 1) }));
        });
    }

    // ------------------------------------------------------------------ search members, positive / absent / empty / boundary

    [Test]
    public void SearchMembers_AgainstTheDefinition_OnTheCorpus()
    {
        foreach (string text in Corpus)
        {
            ISuffixTree tree = Tree(text);
            FacadeChecks.CheckSearch(tree, text, FacadeChecks.Patterns(text, 8), Backend.Name);
        }
    }

    [Test]
    public void Mississippi_HandDerivedCounts()
    {
        ISuffixTree tree = Tree("mississippi");
        Assert.Multiple(() =>
        {
            Assert.That(tree.CountOccurrences("i"), Is.EqualTo(4));
            Assert.That(tree.CountOccurrences("s"), Is.EqualTo(4));
            Assert.That(tree.CountOccurrences("ss"), Is.EqualTo(2));
            Assert.That(tree.CountOccurrences("issi"), Is.EqualTo(2));
            Assert.That(tree.CountOccurrences("ississippi"), Is.EqualTo(1));
            Assert.That(tree.CountOccurrences("mississippi"), Is.EqualTo(1));
            Assert.That(tree.CountOccurrences("mississippix"), Is.EqualTo(0));
            Assert.That(tree.FindAllOccurrences("issi").Order(), Is.EqualTo(new[] { 1, 4 }));
            Assert.That(tree.FindAllOccurrences("ppi"), Is.EqualTo(new[] { 8 }));
            Assert.That(tree.FindAllOccurrences("i").Order(), Is.EqualTo(new[] { 1, 4, 7, 10 }));
            Assert.That(tree.FindAllOccurrences("").Order(), Is.EqualTo(Enumerable.Range(0, 11)));
        });
    }

    [Test]
    public void PatternLongerThanText_IsAbsent_IncludingEveryPrefixPlusOneSymbol()
    {
        foreach (string text in new[] { "a", "ab", "banana", "aaaa" })
        {
            ISuffixTree tree = Tree(text);
            Assert.Multiple(() =>
            {
                Assert.That(tree.Contains(text + "a"), Is.False, text);
                Assert.That(tree.Contains(text + text), Is.False, text);
                Assert.That(tree.CountOccurrences(text + "a"), Is.EqualTo(0), text);
                Assert.That(tree.FindAllOccurrences(text + text), Is.Empty, text);
                Assert.That(tree.Contains(new string('z', text.Length + 1)), Is.False, text);
            });
        }
    }

    [Test]
    public void SpanOverloads_EqualTheStringOverloads_AndAcceptSlices()
    {
        ISuffixTree tree = Tree("abracadabra");
        string padded = "xxabraxx";
        Assert.Multiple(() =>
        {
            Assert.That(tree.Contains(padded.AsSpan(2, 4)), Is.True);
            Assert.That(tree.CountOccurrences(padded.AsSpan(2, 4)), Is.EqualTo(2));
            Assert.That(tree.FindAllOccurrences(padded.AsSpan(2, 4)).Order(), Is.EqualTo(new[] { 0, 7 }));
            Assert.That(tree.Contains(padded.AsSpan(0, 4)), Is.False);
            Assert.That(tree.Contains(ReadOnlySpan<char>.Empty), Is.True);
            Assert.That(tree.CountOccurrences(ReadOnlySpan<char>.Empty), Is.EqualTo(11));
            Assert.That(tree.FindAllOccurrences(ReadOnlySpan<char>.Empty), Has.Count.EqualTo(11));
            Assert.That(tree.LongestCommonSubstring("cadabx".AsSpan(0, 5)), Is.EqualTo("cadab"));
        });
    }

    // ------------------------------------------------------------------ argument exceptions, equal for every backend

    [Test]
    public void NullArguments_ThrowArgumentNullException_WithTheInterfaceParameterNames()
    {
        ISuffixTree tree = Tree("banana");
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.Contains((string)null!))!.ParamName, Is.EqualTo("value"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.FindAllOccurrences((string)null!))!.ParamName, Is.EqualTo("pattern"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.CountOccurrences((string)null!))!.ParamName, Is.EqualTo("pattern"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.LongestCommonSubstring((string)null!))!.ParamName, Is.EqualTo("other"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.LongestCommonSubstringInfo(null!))!.ParamName, Is.EqualTo("other"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.FindAllLongestCommonSubstrings(null!))!.ParamName, Is.EqualTo("other"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.FindAllDistinctLongestCommonSubstrings(null!))!.ParamName, Is.EqualTo("other"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.FindExactMatchAnchors(null!, 1))!.ParamName, Is.EqualTo("query"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.FindMaximalExactMatches(null!, 1))!.ParamName, Is.EqualTo("query"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => tree.FindMaximalUniqueMatches(null!, 1))!.ParamName, Is.EqualTo("query"));
            Assert.Throws<ArgumentNullException>(() => tree.Traverse(null!));
        });
    }

    [Test]
    public void OutOfRangeArguments_ThrowArgumentOutOfRangeException()
    {
        ISuffixTree tree = Tree("banana");
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => tree.CountDistinctSubstringsByLength(0))!.ParamName, Is.EqualTo("maxLength"));
            Assert.Throws<ArgumentOutOfRangeException>(() => tree.CountDistinctSubstringsByLength(-1));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => tree.FindMaximalExactMatches("an", 0))!.ParamName, Is.EqualTo("minLength"));
            Assert.Throws<ArgumentOutOfRangeException>(() => tree.FindMaximalExactMatches("an", -3));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => tree.FindMaximalUniqueMatches("an", 0))!.ParamName, Is.EqualTo("minLength"));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => tree.FindMaximalUniqueMatches("an", 1, (MumUniqueness)99))!.ParamName, Is.EqualTo("uniqueness"));
            Assert.Throws<ArgumentOutOfRangeException>(() => tree.FindMaximalRepeatedPairs(0));
        });
    }

    [Test]
    public void MinLengthOfAnchors_ZeroOrNegative_IsAnEmptyList_NotAnException()
    {
        ISuffixTree tree = Tree("banana");
        Assert.Multiple(() =>
        {
            Assert.That(tree.FindExactMatchAnchors("banana", 0), Is.Empty);
            Assert.That(tree.FindExactMatchAnchors("banana", -5), Is.Empty);
            Assert.That(tree.FindExactMatchAnchors("", 1), Is.Empty);
        });
    }

    // ------------------------------------------------------------------ analysis against the definitions

    [Test]
    public void AnalysisMembers_AgainstTheDefinitions_OnTheCorpus()
    {
        foreach (string text in Corpus)
            FacadeChecks.CheckAnalysis(Tree(text), text, Backend.Name);
    }

    [Test]
    public void LongestRepeatedSubstring_HandDerived()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Tree("banana").LongestRepeatedSubstring(), Is.EqualTo("ana"));
            Assert.That(Tree("mississippi").LongestRepeatedSubstring(), Is.EqualTo("issi"));
            Assert.That(Tree("aaaaa").LongestRepeatedSubstring(), Is.EqualTo("aaaa"));
            Assert.That(Tree("abcdef").LongestRepeatedSubstring(), Is.EqualTo(""));
            Assert.That(Tree("abcab").LongestRepeatedSubstring(), Is.EqualTo("ab"));
            Assert.That(Tree("ababab").LongestRepeatedSubstring(), Is.EqualTo("abab"));
        });
    }

    [Test]
    public void LongestRepeatedSubstring_TiesAreComparedByProperty()
    {
        const string text = "xabyzcdwcd";   // "cd" is the only repeat of length 2
        const string tie = "abxcdyabzcd";    // "ab" and "cd" both occur twice, nothing longer repeats
        var ties = AnalysisReference.AllLongestRepeated(tie);
        Assert.That(ties.Select(t => t.Substring), Is.EqualTo(new[] { "ab", "cd" }));
        string lrs = Tree(tie).LongestRepeatedSubstring();
        Assert.That(lrs, Is.AnyOf("ab", "cd"), "any of the tied substrings is a valid representative (ISuffixTreeAnalysis remarks)");
        Assert.That(Tree(text).LongestRepeatedSubstring(), Is.EqualTo("cd"));
        var all = Tree(tie).FindAllLongestRepeatedSubstrings();
        Assert.That(all.Select(a => a.Substring), Is.EqualTo(new[] { "ab", "cd" }));
        Assert.That(all[0].Positions, Is.EqualTo(new[] { 0, 6 }));
        Assert.That(all[1].Positions, Is.EqualTo(new[] { 3, 9 }));
    }

    [Test]
    public void SuffixEnumeration_IsInOrdinalOrder_AndLazy()
    {
        ISuffixTree tree = Tree("banana");
        string[] expected = ["a", "ana", "anana", "banana", "na", "nana"];
        Assert.That(tree.GetAllSuffixes(), Is.EqualTo(expected));
        Assert.That(tree.EnumerateSuffixes(), Is.EqualTo(expected));
        // the enumerable can be enumerated twice
        IEnumerable<string> lazy = tree.EnumerateSuffixes();
        Assert.That(lazy.Count(), Is.EqualTo(6));
        Assert.That(lazy.First(), Is.EqualTo("a"));
    }

    [Test]
    public void DistinctSubstringCounts_HandDerived()
    {
        ISuffixTree tree = Tree("banana");
        Assert.Multiple(() =>
        {
            Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(15));
            Assert.That(tree.CountDistinctSubstringsByLength(6), Is.EqualTo(new long[] { 1, 3, 3, 3, 3, 2, 1 }));
            Assert.That(tree.CountDistinctSubstringsByLength(2), Is.EqualTo(new long[] { 1, 3, 3 }));
            Assert.That(tree.CountDistinctSubstringsByLength(100), Is.EqualTo(new long[] { 1, 3, 3, 3, 3, 2, 1 }), "clamped to the text length");
        });
    }

    // ------------------------------------------------------------------ matching members

    [Test]
    public void MatchingMembers_AgainstTheDefinitions_OnTheCorpus()
    {
        var rng = new SplitMix64(11);
        foreach (string text in Corpus)
        {
            ISuffixTree tree = Tree(text);
            foreach (string query in FacadeChecks.Queries(text, rng))
                FacadeChecks.CheckMatching(tree, null, text, query, [0, 1, 2, 4], Backend.Name);
        }
    }

    [Test]
    public void MumUmmerCases_HandDerived()
    {
        ISuffixTree tree = Tree("GATTACAGATTACA");
        Assert.Multiple(() =>
        {
            Assert.That(tree.FindMaximalExactMatches("ATTACAGGATTACAT", 3), Is.EqualTo(new[] { (1, 0, 7), (8, 0, 6), (0, 7, 7), (7, 7, 7) }));
            Assert.That(tree.FindMaximalUniqueMatches("ATTACAGGATTACAT", 3), Is.EqualTo(new[] { (1, 0, 7) }));
            Assert.That(tree.FindMaximalUniqueMatches("ATTACAGGATTACAT", 3, MumUniqueness.Reference), Is.EqualTo(new[] { (1, 0, 7) }));
            Assert.That(tree.FindExactMatchAnchors("ATTACAGGATTACAT", 3).Select(a => (a.PositionInQuery, a.Length)), Is.EqualTo(new[] { (0, 7), (7, 7) }));
        });
    }

    [Test]
    public void LongestCommonSubstring_HandDerived()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Tree("xabxac").LongestCommonSubstring("abcabxabcd"), Is.EqualTo("abxa"));
            Assert.That(Tree("banana").LongestCommonSubstring("bandana"), Is.EqualTo("ban"));
            Assert.That(Tree("abc").LongestCommonSubstring("xyz"), Is.EqualTo(""));
            Assert.That(Tree("abc").LongestCommonSubstring(""), Is.EqualTo(""));
            Assert.That(Tree("abc").LongestCommonSubstringInfo("zzbczz"), Is.EqualTo(("bc", 1, 2)));
        });
    }

    [Test]
    public void FindAllLongestCommonSubstrings_ReturnsAllPositions_AndNoMatchIsEmpty()
    {
        var (substring, positionsInText, positionsInOther) = Tree("abcabc").FindAllLongestCommonSubstrings("xabcy");
        Assert.Multiple(() =>
        {
            Assert.That(substring, Is.EqualTo("abc"));
            Assert.That(positionsInText, Is.EqualTo(new[] { 0, 3 }));
            Assert.That(positionsInOther, Is.EqualTo(new[] { 1 }));
        });
        var none = Tree("abc").FindAllLongestCommonSubstrings("xyz");
        Assert.Multiple(() =>
        {
            Assert.That(none.Substring, Is.EqualTo(string.Empty));
            Assert.That(none.PositionsInText, Is.Empty);
            Assert.That(none.PositionsInOther, Is.Empty);
        });
    }

    [Test]
    public void LongestCommonSubstringInfo_SpanOverload_EqualsTheStringOverload()
    {
        // the span overload is a member of the concrete classes (PersistentSuffixTree, EsaSuffixTree), not of ISuffixTree
        foreach (var (text, other) in new[] { ("banana", "bandana"), ("mississippi", "ssi"), ("abc", "def"), ("🧬αβγ🧪", "xxαβγyy"), ("A\u0000B\uFFFFC", "\u0000B\uFFFF"), ("", "abc"), ("abc", "") })
        {
            ISuffixTree tree = Tree(text);
            (string Substring, int PositionInText, int PositionInOther) viaString = tree.LongestCommonSubstringInfo(other);
            (string Substring, int PositionInText, int PositionInOther)? viaSpan = tree switch
            {
                EsaSuffixTree esa => esa.LongestCommonSubstringInfo(other.AsSpan()),
                PersistentSuffixTree persistent => persistent.LongestCommonSubstringInfo(other.AsSpan()),
                _ => null,
            };
            if (viaSpan is null) continue;
            Assert.That(viaSpan.Value, Is.EqualTo(viaString), $"[{Backend.Name}] '{text}' / '{other}'");
            if (viaString.Substring.Length == 0) Assert.That(viaString, Is.EqualTo(("", -1, -1)));
        }
    }

    // ------------------------------------------------------------------ diagnostics

    [Test]
    public void Traverse_BananaEventStream_HandDerived()
    {
        // root: (start, end, leafCount, childCount, depth); children in key order, -1 (the terminator) first
        ISuffixTree tree = Tree("ab");
        List<string> events = FacadeChecks.Events(tree);
        Assert.That(events[0], Is.EqualTo("N(0,0,3,3,0)"), "root: n + 1 leaves, three children (terminator, a, b)");
        Assert.That(events.Count(e => e.StartsWith("E(", StringComparison.Ordinal)), Is.EqualTo(3));
        Assert.That(events.Count(e => e == "X"), Is.EqualTo(3));
        Assert.That(events.Where(e => e.StartsWith("E(", StringComparison.Ordinal)), Is.EqualTo(new[] { "E(-1)", "E(97)", "E(98)" }));
    }

    [Test]
    public void Traverse_VisitsEachNodeOnce_AndBalancesBranches()
    {
        foreach (string text in Corpus)
        {
            ISuffixTree tree = Tree(text);
            List<string> events = FacadeChecks.Events(tree);
            int nodes = events.Count(e => e.StartsWith("N(", StringComparison.Ordinal));
            int enter = events.Count(e => e.StartsWith("E(", StringComparison.Ordinal));
            int exit = events.Count(e => e == "X");
            Assert.That(nodes, Is.EqualTo(tree.NodeCount), text);
            Assert.That(enter, Is.EqualTo(exit), text);
            Assert.That(enter, Is.EqualTo(tree.NodeCount - 1), text);
        }
    }

    [Test]
    public void PrintTree_HeaderRootAndShape()
    {
        ISuffixTree tree = Tree("banana");
        string[] lines = tree.PrintTree().Split('\n', StringSplitOptions.RemoveEmptyEntries).Select(l => l.TrimEnd('\r')).ToArray();
        Assert.Multiple(() =>
        {
            Assert.That(lines[0], Is.EqualTo("Content length: 6"));
            Assert.That(lines[1], Is.EqualTo("0:ROOT"));
            Assert.That(lines.Length, Is.EqualTo(2 + tree.NodeCount - 1), "header, root, one line per other node");
            Assert.That(lines.Count(l => l.EndsWith(" (Leaf)", StringComparison.Ordinal)), Is.EqualTo(7), "six suffixes and the terminator-only leaf");
            Assert.That(lines.Skip(2).Where(l => l.StartsWith("  1:", StringComparison.Ordinal)).Select(l => l.Split(':')[1].TrimStart()[0]), Is.EqualTo(new[] { '#', 'a', 'b', 'n' }));
        });
    }

    [Test]
    public void PrintTree_EmptyTree_IsHeaderAndRoot()
    {
        string[] lines = Tree("").PrintTree().Split('\n', StringSplitOptions.RemoveEmptyEntries).Select(l => l.TrimEnd('\r')).ToArray();
        Assert.That(lines, Is.EqualTo(new[] { "Content length: 0", "0:ROOT" }));
    }
}
