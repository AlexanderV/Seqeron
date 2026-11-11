using System.Text;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Facade;

using Triple = (int PositionInText, int PositionInQuery, int Length);

/// <summary>
/// The checks of every member of <see cref="ISuffixTree"/> (search, analysis, diagnostics) on one candidate tree, written once and used by the default tests, the explicit Exhaustive and
/// Fuzz classes. The oracles are the definitions (<see cref="BruteForce"/>, <see cref="AnalysisReference"/>, the <c>...ByDefinition</c> implementations) and the trees of the other backends
/// that are passed in, never the code under test. What the contract leaves open is checked by property (ADR-001 section 7): the representative of the longest repeated substring on ties,
/// the text position of <c>LongestCommonSubstringInfo</c> and of the anchors, the order of <c>FindAllOccurrences</c>.
/// </summary>
internal static class FacadeChecks
{
    internal static string Short(string s) => s.Length > 40 ? s[..40] + $"...({s.Length})" : s;

    internal static string Render(IEnumerable<int> values) => string.Join(",", values);

    private static void Fail(string context, string text, string what) => Assert.Fail($"{context} text '{Short(text)}': {what}");

    // ---------------------------------------------------------------- search

    /// <summary>Header values, then the search members on <paramref name="patterns"/> against the definition.</summary>
    internal static void CheckSearch(ISuffixTree tree, string text, IEnumerable<string> patterns, string context)
    {
        if (tree.Text.ToString() != text || tree.Text.Length != text.Length) Fail(context, text, "Text is not the indexed text");
        if (tree.IsEmpty != (text.Length == 0)) Fail(context, text, "IsEmpty");
        if (tree.LeafCount != text.Length) Fail(context, text, $"LeafCount {tree.LeafCount}");
        int n = text.Length;
        if (tree.MaxDepth != n) Fail(context, text, $"MaxDepth {tree.MaxDepth}");
        if (n == 0)
        {
            if (tree.NodeCount != 1) Fail(context, text, $"NodeCount of the empty text is {tree.NodeCount}");
        }
        else if (tree.NodeCount < n + 1 || tree.NodeCount > 2 * n + 1)
        {
            Fail(context, text, $"NodeCount {tree.NodeCount} outside [n + 1, 2n + 1]");
        }

        foreach (string pattern in patterns)
        {
            List<int> expected = BruteForce.Occurrences(text, pattern);
            string at = $"pattern '{Short(pattern)}'";
            bool contains = pattern.Length == 0 || expected.Count > 0;
            if (tree.Contains(pattern) != contains || tree.Contains(pattern.AsSpan()) != contains) Fail(context, text, $"{at}: Contains");
            if (tree.CountOccurrences(pattern) != expected.Count || tree.CountOccurrences(pattern.AsSpan()) != expected.Count) Fail(context, text, $"{at}: CountOccurrences {tree.CountOccurrences(pattern)} against {expected.Count}");
            IReadOnlyList<int> found = tree.FindAllOccurrences(pattern);
            IReadOnlyList<int> foundSpan = tree.FindAllOccurrences(pattern.AsSpan());
            if (found.Count != expected.Count || foundSpan.Count != expected.Count) Fail(context, text, $"{at}: FindAllOccurrences count {found.Count} against {expected.Count}");
            if (!found.Order().SequenceEqual(expected) || !foundSpan.Order().SequenceEqual(expected)) Fail(context, text, $"{at}: positions [{Render(found.Order())}] against [{Render(expected)}]");
        }
    }

    // ---------------------------------------------------------------- analysis (no query)

    /// <summary>The query-free analysis members against the definitions.</summary>
    internal static void CheckAnalysis(ISuffixTree tree, string text, string context, bool definitions = true)
    {
        int n = text.Length;
        var ties = AnalysisReference.AllLongestRepeated(text);
        string lrs = tree.LongestRepeatedSubstring();
        int lrsLength = ties.Count == 0 ? 0 : ties[0].Substring.Length;
        if (lrs.Length != lrsLength) Fail(context, text, $"LongestRepeatedSubstring '{Short(lrs)}' has length {lrs.Length}, the definition says {lrsLength}");
        if (ties.Count > 0 && !ties.Any(t => t.Substring == lrs)) Fail(context, text, $"LongestRepeatedSubstring '{Short(lrs)}' is not one of the {ties.Count} longest repeats");
        if (ties.Count == 1 && lrs != ties[0].Substring) Fail(context, text, "LongestRepeatedSubstring on a unique maximum");
        if (tree.LongestRepeatedSubstringMemory().ToString() != lrs) Fail(context, text, "LongestRepeatedSubstringMemory differs from LongestRepeatedSubstring");

        var all = tree.FindAllLongestRepeatedSubstrings();
        if (all.Count != ties.Count) Fail(context, text, $"FindAllLongestRepeatedSubstrings has {all.Count} entries, the definition {ties.Count}");
        for (int i = 0; i < ties.Count; i++)
            if (all[i].Substring != ties[i].Substring || !all[i].Positions.SequenceEqual(ties[i].Positions))
                Fail(context, text, $"FindAllLongestRepeatedSubstrings entry {i}: '{all[i].Substring}' [{Render(all[i].Positions)}] against '{ties[i].Substring}' [{Render(ties[i].Positions)}]");

        if (tree.CountDistinctSubstrings() != AnalysisReference.DistinctTotal(text)) Fail(context, text, $"CountDistinctSubstrings {tree.CountDistinctSubstrings()}");
        foreach (int m in new[] { 1, 2, 3, Math.Max(1, n - 1), Math.Max(1, n), n + 1, n + 7 }.Distinct())
        {
            long[] counts = tree.CountDistinctSubstringsByLength(m);
            long[] expected = AnalysisReference.DistinctByLength(text, m);
            if (!counts.SequenceEqual(expected)) Fail(context, text, $"CountDistinctSubstringsByLength({m}) [{string.Join(",", counts)}] against [{string.Join(",", expected)}]");
        }

        string[] suffixes = AnalysisReference.SortedSuffixes(text);
        var list = tree.GetAllSuffixes();
        if (!list.SequenceEqual(suffixes)) Fail(context, text, "GetAllSuffixes");
        if (!tree.EnumerateSuffixes().SequenceEqual(suffixes)) Fail(context, text, "EnumerateSuffixes");

        foreach (int minLength in new[] { 1, 2, 3, 5 })
        {
            foreach (Func<char, bool>? unique in new Func<char, bool>?[] { null, c => c == 'b' || c == '\0' })
            {
                var pairs = tree.FindMaximalRepeatedPairs(minLength, unique);
                var expected = AnalysisReference.MaximalRepeatedPairs(text, minLength, unique);
                if (!pairs.SequenceEqual(expected)) Fail(context, text, $"FindMaximalRepeatedPairs({minLength}, {(unique is null ? "none" : "unique b/NUL")}) {pairs.Count} pairs against {expected.Count}");
            }
        }
    }

    // ---------------------------------------------------------------- matching (query)

    /// <summary>The query members: LCS family, anchors, MEM, MUM against the definitions and, when given, the Ukkonen tree.</summary>
    internal static void CheckMatching(ISuffixTree tree, ISuffixTree? ukkonen, string text, string query, int[] minLengths, string context)
    {
        string where = $"{context} text '{Short(text)}' query '{Short(query)}'";
        var source = new StringTextSource(text);

        // LCS family
        var expectedDistinct = SuffixTreeAlgorithms.FindAllDistinctLongestCommonSubstringsByDefinition(source, query);
        var distinct = tree.FindAllDistinctLongestCommonSubstrings(query);
        if (distinct.Count != expectedDistinct.Count) Assert.Fail($"{where}: {distinct.Count} distinct longest common substrings, expected {expectedDistinct.Count}");
        for (int i = 0; i < expectedDistinct.Count; i++)
        {
            if (distinct[i].Substring != expectedDistinct[i].Substring
                || !distinct[i].PositionsInText.SequenceEqual(expectedDistinct[i].PositionsInText)
                || !distinct[i].PositionsInOther.SequenceEqual(expectedDistinct[i].PositionsInOther))
                Assert.Fail($"{where}: distinct entry {i} '{distinct[i].Substring}' differs from the definition '{expectedDistinct[i].Substring}'");
        }

        var allLcs = tree.FindAllLongestCommonSubstrings(query);
        var info = tree.LongestCommonSubstringInfo(query);
        string lcs = tree.LongestCommonSubstring(query);
        string lcsSpan = tree.LongestCommonSubstring(query.AsSpan());
        if (expectedDistinct.Count == 0)
        {
            if (allLcs.Substring.Length != 0 || allLcs.PositionsInText.Count != 0 || allLcs.PositionsInOther.Count != 0) Assert.Fail($"{where}: FindAllLongestCommonSubstrings of no match is not empty");
            if (info != ("", -1, -1)) Assert.Fail($"{where}: info of no match is {info}");
            if (lcs.Length != 0 || lcsSpan.Length != 0) Assert.Fail($"{where}: LongestCommonSubstring of no match");
        }
        else
        {
            var first = expectedDistinct[0];
            if (allLcs.Substring != first.Substring || !allLcs.PositionsInText.SequenceEqual(first.PositionsInText) || !allLcs.PositionsInOther.SequenceEqual(first.PositionsInOther))
                Assert.Fail($"{where}: FindAllLongestCommonSubstrings differs from the first distinct entry");
            if (info.Substring != first.Substring || info.PositionInOther != first.PositionsInOther[0] || !first.PositionsInText.Contains(info.PositionInText))
                Assert.Fail($"{where}: info {info}: substring and position in the other string are fixed, the text position must be one of [{Render(first.PositionsInText)}]");
            if (lcs != first.Substring || lcsSpan != first.Substring) Assert.Fail($"{where}: LongestCommonSubstring '{lcs}' / '{lcsSpan}' against '{first.Substring}'");
        }

        if (ukkonen is not null)
        {
            var uDistinct = ukkonen.FindAllDistinctLongestCommonSubstrings(query);
            if (uDistinct.Count != distinct.Count || uDistinct.Where((u, i) => u.Substring != distinct[i].Substring || !u.PositionsInText.SequenceEqual(distinct[i].PositionsInText) || !u.PositionsInOther.SequenceEqual(distinct[i].PositionsInOther)).Any())
                Assert.Fail($"{where}: distinct longest common substrings differ from the Ukkonen tree's");
        }

        foreach (int minLength in minLengths)
        {
            string at = $"{where} L {minLength}";

            // anchors: (query position, length) fixed, the text position is a valid occurrence
            var anchors = tree.FindExactMatchAnchors(query, minLength);
            var shapes = AnchorShapes(text, query, minLength);
            if (!anchors.Select(a => (a.PositionInQuery, a.Length)).SequenceEqual(shapes)) Assert.Fail($"{at}: anchors {string.Join(" ", anchors)} against the definition {string.Join(" ", shapes)}");
            foreach (Triple a in anchors)
                if (string.CompareOrdinal(text, a.PositionInText, query, a.PositionInQuery, a.Length) != 0) Assert.Fail($"{at}: anchor {a} is not an occurrence");
            if (ukkonen is not null && !anchors.Select(a => (a.PositionInQuery, a.Length)).SequenceEqual(ukkonen.FindExactMatchAnchors(query, minLength).Select(a => (a.PositionInQuery, a.Length))))
                Assert.Fail($"{at}: anchors differ from the Ukkonen tree's");

            if (minLength < 1) continue;
            var mems = tree.FindMaximalExactMatches(query, minLength);
            var memsExpected = SuffixTreeAlgorithms.FindMaximalExactMatchesByDefinition(source, query, minLength);
            if (!mems.SequenceEqual(memsExpected)) Assert.Fail($"{at}: MEM {string.Join(" ", mems)} against the definition {string.Join(" ", memsExpected)}");
            if (ukkonen is not null && !mems.SequenceEqual(ukkonen.FindMaximalExactMatches(query, minLength))) Assert.Fail($"{at}: MEM differs from the Ukkonen tree's");
            foreach (MumUniqueness mode in new[] { MumUniqueness.Both, MumUniqueness.Reference })
            {
                var mums = tree.FindMaximalUniqueMatches(query, minLength, mode);
                var mumsExpected = SuffixTreeAlgorithms.FindMaximalUniqueMatchesByDefinition(source, query, minLength, mode);
                if (!mums.SequenceEqual(mumsExpected)) Assert.Fail($"{at} {mode}: MUM {string.Join(" ", mums)} against the definition {string.Join(" ", mumsExpected)}");
                if (ukkonen is not null && !mums.SequenceEqual(ukkonen.FindMaximalUniqueMatches(query, minLength, mode))) Assert.Fail($"{at} {mode}: MUM differs from the Ukkonen tree's");
            }
            if (!tree.FindMaximalUniqueMatches(query, minLength).SequenceEqual(tree.FindMaximalUniqueMatches(query, minLength, MumUniqueness.Both))) Assert.Fail($"{at}: the default uniqueness is not Both");
        }
    }

    /// <summary>The (query position, length) pairs of the anchors by their definition: one per maximal run of matching statistics at least minLength, its first peak.</summary>
    internal static List<(int PositionInQuery, int Length)> AnchorShapes(string text, string query, int minLength)
    {
        var result = new List<(int, int)>();
        if (query.Length == 0 || text.Length == 0 || minLength <= 0) return result;
        var ms = new int[query.Length];
        for (int i = 0; i < query.Length; i++)
        {
            int best = 0;
            for (int length = 1; length <= i + 1; length++)
            {
                if (text.AsSpan().IndexOf(query.AsSpan(i - length + 1, length), StringComparison.Ordinal) < 0) break;
                best = length;
            }
            ms[i] = best;
        }
        int peakLength = 0, peakEnd = -1;
        for (int i = 0; i < ms.Length; i++)
        {
            if (ms[i] >= minLength)
            {
                if (ms[i] > peakLength)
                {
                    peakLength = ms[i];
                    peakEnd = i;
                }
            }
            else if (peakLength >= minLength)
            {
                result.Add((peakEnd - peakLength + 1, peakLength));
                peakLength = 0;
                peakEnd = -1;
            }
        }
        if (peakLength >= minLength) result.Add((peakEnd - peakLength + 1, peakLength));
        return result;
    }

    // ---------------------------------------------------------------- diagnostics

    internal sealed class EventCollector : ISuffixTreeVisitor
    {
        internal readonly List<string> Events = [];

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Events.Add($"N({startIndex},{endIndex},{leafCount},{childCount},{depth})");

        public void EnterBranch(int key) => Events.Add($"E({key})");

        public void ExitBranch() => Events.Add("X");
    }

    /// <summary>Shape of the event stream with the occurrence-dependent start/end indices of internal edges removed: what must be identical between backends.</summary>
    internal static List<string> Events(ISuffixTree tree)
    {
        var c = new EventCollector();
        tree.Traverse(c);
        return c.Events;
    }

    /// <summary>The Traverse event stream, the node counts and (when given) the logical hash equal those of the reference trees.</summary>
    internal static void CheckAgainstTree(ISuffixTree tree, ISuffixTree reference, string referenceName, string text, string context)
    {
        if (tree.NodeCount != reference.NodeCount) Fail(context, text, $"NodeCount {tree.NodeCount} against {referenceName} {reference.NodeCount}");
        if (tree.LeafCount != reference.LeafCount) Fail(context, text, $"LeafCount against {referenceName}");
        if (tree.MaxDepth != reference.MaxDepth) Fail(context, text, $"MaxDepth against {referenceName}");
        List<string> mine = Events(tree);
        List<string> theirs = Events(reference);
        if (!mine.SequenceEqual(theirs))
        {
            int k = 0;
            while (k < mine.Count && k < theirs.Count && mine[k] == theirs[k]) k++;
            Fail(context, text, $"Traverse event stream differs from {referenceName} at event {k}: '{(k < mine.Count ? mine[k] : "<end>")}' against '{(k < theirs.Count ? theirs[k] : "<end>")}'");
        }
    }

    /// <summary>
    /// PrintTree: the header, the root line and the node lines (depth, indentation, label, leaf mark) equal the reference's once the link marks are removed (the marks are
    /// backend specific: the Ukkonen tree prints the code of the first symbol, the v6 tree and the ESA print the character).
    /// </summary>
    internal static void CheckPrintTree(ISuffixTree tree, ISuffixTree reference, string text, string context)
    {
        string mine = StripLinkMarks(tree.PrintTree());
        string theirs = StripLinkMarks(reference.PrintTree());
        if (mine != theirs) Fail(context, text, "PrintTree differs from the reference's once the link marks are removed:\n" + mine + "\n---\n" + theirs);
    }

    internal static string StripLinkMarks(string printed)
    {
        var sb = new StringBuilder();
        foreach (string line in printed.Split('\n'))
        {
            int k = line.IndexOf(" -> [", StringComparison.Ordinal);
            sb.Append(k >= 0 ? line[..k] : line).Append('\n');
        }
        return sb.ToString();
    }

    // ---------------------------------------------------------------- query recipes

    /// <summary>Probe patterns: the empty one, every distinct substring up to <paramref name="maxLength"/>, the whole text, absent patterns of every kind.</summary>
    internal static IEnumerable<string> Patterns(string text, int maxLength = 6)
    {
        yield return string.Empty;
        var seen = new HashSet<string>(StringComparer.Ordinal) { string.Empty };
        for (int i = 0; i < text.Length; i++)
            for (int l = 1; l <= Math.Min(maxLength, text.Length - i); l++)
                if (seen.Add(text.Substring(i, l))) yield return text.Substring(i, l);
        foreach (string p in new[] { text, text + "q", "q" + text, "q", text + text, "\u0001", "ab\u0001" })
            if (seen.Add(p)) yield return p;
        if (text.Length > 0)
        {
            foreach (string p in new[] { text + text[0], new string(text[^1], text.Length + 1), string.Concat(text.AsSpan(0, text.Length - 1), "\u0001"), text[1..] + text[0], text[..^1] })
                if (seen.Add(p)) yield return p;
            string suffix = text[(text.Length / 2)..];
            if (seen.Add(suffix)) yield return suffix;
        }
    }

    /// <summary>Probe queries for the matching members.</summary>
    internal static IEnumerable<string> Queries(string text, SplitMix64 rng)
    {
        yield return string.Empty;
        yield return text;
        yield return "q";
        yield return text + "q" + text;
        if (text.Length > 0)
        {
            int start = rng.NextInt(text.Length);
            string window = text.Substring(start, 1 + rng.NextInt(text.Length - start));
            yield return window;
            yield return "x" + window + "y";
            yield return new string(window.Reverse().ToArray());
            var mutated = text.ToCharArray();
            for (int k = 0; k < Math.Max(1, mutated.Length / 10); k++) mutated[rng.NextInt(mutated.Length)] = (char)('a' + rng.NextInt(3));
            yield return new string(mutated);
            yield return text[(text.Length / 2)..] + text[..(text.Length / 2)];
        }
        yield return new string(Enumerable.Range(0, 1 + rng.NextInt(12)).Select(_ => (char)('a' + rng.NextInt(4))).ToArray());
    }

    /// <summary>Everything: search, analysis, matching, diagnostics, against the definitions and the given reference trees.</summary>
    internal static void CheckEverything(ISuffixTree tree, ISuffixTree? ukkonen, ISuffixTree? v6, string text, ulong seed, string context, int maxPatternLength = 6, bool printTree = true)
    {
        CheckSearch(tree, text, Patterns(text, maxPatternLength), context);
        CheckAnalysis(tree, text, context);
        var rng = new SplitMix64(seed ^ (ulong)text.Length);
        foreach (string query in Queries(text, rng))
            CheckMatching(tree, ukkonen, text, query, [0, 1, 2, 3, 5], context);
        if (ukkonen is not null)
        {
            CheckAgainstTree(tree, ukkonen, "ukkonen", text, context);
            if (printTree) CheckPrintTree(tree, ukkonen, text, context);
        }
        if (v6 is not null)
        {
            CheckAgainstTree(tree, v6, "v6", text, context);
            if (printTree) CheckPrintTree(tree, v6, text, context);
        }
    }
}
