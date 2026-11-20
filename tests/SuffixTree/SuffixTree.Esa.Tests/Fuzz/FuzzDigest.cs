using System.Text;
using SuffixTree.Esa.Tests.Facade;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// The probes of one text: patterns for the search members, queries for the matching members, the minimum lengths. A pure function of the text and a seed, so the same script runs on every engine of a case, on every
/// thread of the concurrency test and again on the shrunk texts of the minimizer.
/// </summary>
internal sealed class FuzzScript
{
    /// <summary>Texts of at most this length get the definition-level patterns of <see cref="FacadeChecks.Patterns"/> (every distinct substring up to length 5).</summary>
    internal const int ExhaustivePatternsUpTo = 120;

    /// <summary>The queries of a text above this length are at most 150 symbols (up to it at most 400, up to <see cref="ExhaustivePatternsUpTo"/> the queries of <see cref="FacadeChecks.Queries"/>).</summary>
    internal const int BoundedUpTo = 2000;

    private FuzzScript(string text, List<string> patterns, List<string> queries)
    {
        Text = text;
        Patterns = patterns;
        Queries = queries;
    }

    /// <summary>The text.</summary>
    internal string Text { get; }

    /// <summary>The patterns of the search members (the empty pattern, present, absent, longer than the text).</summary>
    internal IReadOnlyList<string> Patterns { get; }

    /// <summary>The queries of the matching members.</summary>
    internal IReadOnlyList<string> Queries { get; }

    /// <summary>The minimum lengths of anchors, MEM and MUM (0 only has a meaning for the anchors).</summary>
    internal int[] MinLengths { get; } = [0, 1, 2, 3, 5, 12];

    /// <summary>Builds the script of <paramref name="text"/>.</summary>
    internal static FuzzScript Create(string text, ulong seed)
    {
        int n = text.Length;
        var rng = new SplitMix64(seed ^ ((ulong)n * 0x9E3779B97F4A7C15UL));
        var patterns = new List<string>();
        if (n <= ExhaustivePatternsUpTo)
        {
            patterns.AddRange(FacadeChecks.Patterns(text, 5));
        }
        else
        {
            patterns.Add(string.Empty);
            patterns.Add(text[..1]);
            patterns.Add(text[^1..]);
            patterns.Add(text[..Math.Min(n, 40)]);
            patterns.Add(text[Math.Max(0, n - 40)..]);
            if (n <= 5000) patterns.Add(text);
            patterns.Add(text[..Math.Min(n, 30)] + "q");
            patterns.Add("q" + text[..Math.Min(n, 30)]);
            patterns.Add(new string(text[^1], Math.Min(n + 1, 60)));
            for (int k = 0; k < 36; k++)
            {
                int length = 1 + (k < 12 ? rng.NextInt(6) : rng.NextInt(30));
                int start = rng.NextInt(n - Math.Min(length, n) + 1);
                string p = text.Substring(start, Math.Min(length, n - start));
                patterns.Add(p);
                if (k % 6 == 0) patterns.Add(p.Length > 1 ? string.Concat(p.AsSpan(0, p.Length - 1), "\u0001") : p + "\u0001");
            }
        }

        // patterns that deviate from a present pattern at one position (the first mismatch inside an edge, not at a key): the absent patterns that a descent must refuse after having matched a prefix
        char[] alphabet = [.. text.Take(2000).Distinct().Append('q').Append('\u0001')];
        int present = patterns.Count;
        for (int k = 0; k < 14 && present > 0; k++)
        {
            string basis = patterns[rng.NextInt(present)];
            if (basis.Length < 2) continue;
            int at = rng.NextInt(basis.Length);
            char symbol = alphabet[rng.NextInt(alphabet.Length)];
            if (symbol != basis[at]) patterns.Add(string.Concat(basis.AsSpan(0, at), symbol.ToString(), basis.AsSpan(at + 1)));
        }

        var queries = new List<string>();
        if (n <= ExhaustivePatternsUpTo)
        {
            queries.AddRange(FacadeChecks.Queries(text, rng));
        }
        else
        {
            // the output of the matching members is as long as the number of (query position, text position) pairs of the matches, which on a repetitive text is quadratic: the queries are short
            int cap = n <= BoundedUpTo ? 400 : 150;
            queries.Add(string.Empty);
            queries.Add("q");
            int start = rng.NextInt(n);
            string window = text.Substring(start, Math.Min(1 + rng.NextInt(cap), n - start));
            queries.Add(window);
            queries.Add("x" + window + "y");
            queries.Add(new string(window.Reverse().ToArray()));
            var mutated = window.ToCharArray();
            for (int k = 0; k < Math.Max(1, mutated.Length / 10); k++) mutated[rng.NextInt(mutated.Length)] = (char)('a' + rng.NextInt(3));
            queries.Add(new string(mutated));
            int half = Math.Min(n, cap / 2);
            queries.Add(string.Concat(text.AsSpan(n - half), "q", text.AsSpan(0, half)));
            queries.Add(new string(Enumerable.Range(0, 1 + rng.NextInt(40)).Select(_ => (char)('a' + rng.NextInt(4))).ToArray()));
        }
        return new FuzzScript(text, patterns, queries);
    }
}

/// <summary>Limits that keep the cost of a digest bounded: the members whose cost is quadratic in the text are run on short texts only.</summary>
/// <param name="SuffixesUpTo">Longest text whose <c>GetAllSuffixes</c> is compared in full (above it, the first 50 of <c>EnumerateSuffixes</c>).</param>
/// <param name="PairsUpTo">Longest text for <c>FindMaximalRepeatedPairs</c> (a^n has n^2 / 2 pairs).</param>
/// <param name="PrintUpTo">Longest text whose <c>PrintTree</c> is compared.</param>
internal sealed record DigestLimits(int SuffixesUpTo = 3000, int PairsUpTo = 400, int PrintUpTo = 5000)
{
    internal static DigestLimits Default { get; } = new();
}

/// <summary>
/// The answers of a tree to a <see cref="FuzzScript"/>, member by member, as a list of <c>label=value</c> lines (work package 1b.8): what the contract fixes is compared as it is, what the contract leaves open is
/// reduced to its property (ADR-001 section 7): the order of <c>FindAllOccurrences</c> (sorted), the representative of the longest repeat on a tie (its length, and that it occurs twice), the text position of
/// <c>LongestCommonSubstringInfo</c> and of the anchors (that it is an occurrence of the match). The digest of the in-memory tree, of the file-backed tree, of the Ukkonen tree and of the v6 tree of one text are
/// the same lines, so the engines are compared against each other, and the checksum of the lines is the single-thread answer of the concurrency test.
/// </summary>
internal static class FuzzDigest
{
    internal static ulong Fnv(ulong h, long v)
    {
        unchecked
        {
            for (int i = 0; i < 8; i++)
            {
                h ^= (byte)(v >> (8 * i));
                h *= 0x100000001B3UL;
            }
            return h;
        }
    }

    internal static ulong Fnv(ulong h, string s)
    {
        unchecked
        {
            h = Fnv(h, s.Length);
            foreach (char c in s)
            {
                h ^= c;
                h *= 0x100000001B3UL;
            }
            return h;
        }
    }

    internal const ulong Basis = 0xCBF29CE484222325UL;

    private static string Hex(ulong h) => h.ToString("x16", System.Globalization.CultureInfo.InvariantCulture);

    private static ulong HashInts(IEnumerable<int> values)
    {
        ulong h = Basis;
        foreach (int v in values) h = Fnv(h, v);
        return h;
    }

    private static ulong HashSorted(IReadOnlyList<int> values)
    {
        int[] copy = [.. values];
        Array.Sort(copy);
        return HashInts(copy);
    }

    private static bool Occurs(string text, string match, int position, int at) =>
        position >= 0 && at >= 0 && position <= text.Length - match.Length && string.CompareOrdinal(text, position, match, at, match.Length) == 0;

    private static bool OccursTwice(string text, string s)
    {
        if (s.Length == 0) return true;
        int first = text.IndexOf(s, StringComparison.Ordinal);
        return first >= 0 && text.IndexOf(s, first + 1, StringComparison.Ordinal) >= 0;
    }

    private sealed class HashingVisitor : ISuffixTreeVisitor
    {
        internal ulong Hash = Basis;
        internal long Nodes;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            Nodes++;
            Hash = Fnv(Fnv(Fnv(Fnv(Fnv(Hash, 1), startIndex), endIndex), leafCount), (long)childCount * 1_000_003 + depth);
        }

        public void EnterBranch(int key) => Hash = Fnv(Fnv(Hash, 2), key);

        public void ExitBranch() => Hash = Fnv(Hash, 3);
    }

    /// <summary>The lines of <paramref name="tree"/> for <paramref name="script"/>.</summary>
    internal static List<string> Compute(ISuffixTree tree, string text, FuzzScript script, DigestLimits limits)
    {
        var lines = new List<string>();
        int n = text.Length;
        void Add(string label, string value) => lines.Add(label + "=" + value);

        Add("header", $"len={tree.Text.Length} empty={tree.IsEmpty} leaves={tree.LeafCount} depth={tree.MaxDepth} nodes={tree.NodeCount}");
        AddTextAccess(tree, text, Add);

        // search
        for (int i = 0; i < script.Patterns.Count; i++)
        {
            string p = script.Patterns[i];
            bool contains = tree.Contains(p);
            int count = tree.CountOccurrences(p);
            IReadOnlyList<int> occurrences = tree.FindAllOccurrences(p);
            ulong positions = HashSorted(occurrences);
            bool span = tree.Contains(p.AsSpan()) == contains && tree.CountOccurrences(p.AsSpan()) == count && HashSorted(tree.FindAllOccurrences(p.AsSpan())) == positions;
            Add($"search[{i}]", $"{contains}/{count}/{occurrences.Count}/{Hex(positions)}{(span ? string.Empty : " SPAN-OVERLOAD-DIFFERS")}");
        }

        // analysis without a query
        string lrs = tree.LongestRepeatedSubstring();
        Add("lrs", $"len={lrs.Length} occurs-twice={OccursTwice(text, lrs)} memory-equal={tree.LongestRepeatedSubstringMemory().ToString() == lrs} cached={ReferenceEquals(lrs, tree.LongestRepeatedSubstring()) || lrs.Length == 0}");
        var all = tree.FindAllLongestRepeatedSubstrings();
        ulong h = Basis;
        foreach (var (substring, positions) in all)
        {
            h = Fnv(h, substring);
            foreach (int q in positions) h = Fnv(h, q);
        }
        Add("lrs-all", $"{all.Count}/{Hex(h)}");
        Add("distinct", tree.CountDistinctSubstrings().ToString(System.Globalization.CultureInfo.InvariantCulture));
        foreach (int m in new[] { 1, 2, 3, 20, Math.Max(1, n - 1), Math.Max(1, n), n + 1, n + 7 }.Distinct())
            Add($"distinct-by-length[{m}]", Hex(HashLongs(tree.CountDistinctSubstringsByLength(m))));
        if (n <= limits.SuffixesUpTo)
        {
            Add("suffixes", Hex(HashStrings(tree.GetAllSuffixes())));
            Add("suffixes-enumerated", Hex(HashStrings(tree.EnumerateSuffixes())));
        }
        else
        {
            Add("suffixes-first-50", Hex(HashStrings(tree.EnumerateSuffixes().Take(50))));
        }
        if (n <= limits.PairsUpTo)
        {
            foreach (int minLength in new[] { 1, 2, 3, 5 })
                foreach (Func<char, bool>? unique in new Func<char, bool>?[] { null, c => c == 'b' || c == '\0' })
                {
                    var pairs = tree.FindMaximalRepeatedPairs(minLength, unique);
                    ulong ph = Basis;
                    foreach (var (a, b, l) in pairs) ph = Fnv(Fnv(Fnv(ph, a), b), l);
                    Add($"pairs[{minLength},{(unique is null ? "none" : "b/NUL")}]", $"{pairs.Count}/{Hex(ph)}");
                }
        }

        // matching
        for (int j = 0; j < script.Queries.Count; j++)
            AddMatching(tree, text, script, script.Queries[j], j, Add);

        AddErrors(tree, script, Add);

        // diagnostics
        var visitor = new HashingVisitor();
        tree.Traverse(visitor);
        Add("traverse", $"{visitor.Nodes}/{Hex(visitor.Hash)}");
        if (n <= limits.PrintUpTo) Add("print-tree", Hex(Fnv(Basis, FacadeChecks.StripLinkMarks(tree.PrintTree()))));
        return lines;
    }

    /// <summary>The exceptions of the argument checks: a null argument, a minimum length below 1, an undefined uniqueness, a length of zero; the type and the parameter name must be those of the other backends.</summary>
    private static void AddErrors(ISuffixTree tree, FuzzScript script, Action<string, string> add)
    {
        string q = script.Queries.Count > 2 ? script.Queries[2] : "ab";
        void Probe(string label, Action call)
        {
            string outcome;
            try
            {
                call();
                outcome = "no exception";
            }
            catch (ArgumentException e)
            {
                outcome = $"{e.GetType().Name}({e.ParamName})";
            }
            catch (Exception e)
            {
                outcome = e.GetType().Name;
            }
            add($"error[{label}]", outcome);
        }

        Probe("Contains(null)", () => tree.Contains((string)null!));
        Probe("FindAllOccurrences(null)", () => tree.FindAllOccurrences((string)null!));
        Probe("CountOccurrences(null)", () => tree.CountOccurrences((string)null!));
        Probe("LongestCommonSubstring(null)", () => tree.LongestCommonSubstring((string)null!));
        Probe("LongestCommonSubstringInfo(null)", () => tree.LongestCommonSubstringInfo(null!));
        Probe("FindAllLongestCommonSubstrings(null)", () => tree.FindAllLongestCommonSubstrings(null!));
        Probe("FindAllDistinctLongestCommonSubstrings(null)", () => tree.FindAllDistinctLongestCommonSubstrings(null!));
        Probe("FindExactMatchAnchors(null)", () => tree.FindExactMatchAnchors(null!, 3));
        Probe("FindMaximalExactMatches(null)", () => tree.FindMaximalExactMatches(null!, 3));
        Probe("FindMaximalUniqueMatches(null)", () => tree.FindMaximalUniqueMatches(null!, 3));
        Probe("FindMaximalExactMatches(minLength 0)", () => tree.FindMaximalExactMatches(q, 0));
        Probe("FindMaximalExactMatches(minLength -1)", () => tree.FindMaximalExactMatches(q, -1));
        Probe("FindMaximalUniqueMatches(minLength 0)", () => tree.FindMaximalUniqueMatches(q, 0));
        Probe("FindMaximalUniqueMatches(undefined uniqueness)", () => tree.FindMaximalUniqueMatches(q, 3, (MumUniqueness)99));
        Probe("CountDistinctSubstringsByLength(0)", () => tree.CountDistinctSubstringsByLength(0));
        Probe("CountDistinctSubstringsByLength(-1)", () => tree.CountDistinctSubstringsByLength(-1));
        Probe("Traverse(null)", () => tree.Traverse(null!));
    }

    private static void AddTextAccess(ISuffixTree tree, string text, Action<string, string> add)
    {
        ITextSource source = tree.Text;
        int n = text.Length;
        var sb = new StringBuilder();
        if (n > 0)
        {
            foreach (int i in new[] { 0, n / 2, n - 1 }.Distinct())
                if (source[i] != text[i]) sb.Append($" item[{i}]");
            int length = Math.Min(n, 17);
            int start = n / 3;
            length = Math.Min(length, n - start);
            if (source.Substring(start, length) != text.Substring(start, length)) sb.Append(" substring");
            if (source.Slice(start, length).ToString() != text.Substring(start, length)) sb.Append(" slice");
        }
        if (source.ToString() != text) sb.Append(" to-string");
        add("text", sb.Length == 0 ? "ok" : "BAD" + sb);
    }

    private static ulong HashLongs(long[] values)
    {
        ulong h = Basis;
        foreach (long v in values) h = Fnv(h, v);
        return h;
    }

    private static ulong HashStrings(IEnumerable<string> values)
    {
        ulong h = Basis;
        foreach (string s in values) h = Fnv(h, s);
        return h;
    }

    private static void AddMatching(ISuffixTree tree, string text, FuzzScript script, string query, int j, Action<string, string> add)
    {
        string label = $"query[{j}]";
        string lcs = tree.LongestCommonSubstring(query);
        string lcsSpan = tree.LongestCommonSubstring(query.AsSpan());
        var info = tree.LongestCommonSubstringInfo(query);
        bool infoValid = info.Substring.Length == 0 ? info.PositionInText == -1 && info.PositionInOther == -1 : Occurs(text, info.Substring, info.PositionInText, 0) && Occurs(query, info.Substring, info.PositionInOther, 0);
        add($"{label}.lcs", $"{Hex(Fnv(Basis, lcs))}/{lcs == lcsSpan}/{Hex(Fnv(Basis, info.Substring))}/{info.PositionInOther}/{infoValid}");
        var allLcs = tree.FindAllLongestCommonSubstrings(query);
        add($"{label}.lcs-all", Hex(HashLcs(allLcs.Substring, allLcs.PositionsInText, allLcs.PositionsInOther)));
        var distinct = tree.FindAllDistinctLongestCommonSubstrings(query);
        ulong d = Basis;
        foreach (var (substring, inText, inOther) in distinct) d = Fnv(d, unchecked((long)HashLcs(substring, inText, inOther)));
        add($"{label}.lcs-distinct", $"{distinct.Count}/{Hex(d)}");

        foreach (int minLength in script.MinLengths)
        {
            string at = $"{label}.L{minLength}";
            var anchors = tree.FindExactMatchAnchors(query, minLength);
            ulong ah = Basis;
            bool valid = true;
            foreach (var (positionInText, positionInQuery, length) in anchors)
            {
                ah = Fnv(Fnv(ah, positionInQuery), length);
                valid &= length >= minLength && Occurs(text, query.Substring(positionInQuery, length), positionInText, 0);
            }
            add($"{at}.anchors", $"{anchors.Count}/{Hex(ah)}/{valid}");
            if (minLength < 1) continue;
            add($"{at}.mem", HashTriples(tree.FindMaximalExactMatches(query, minLength)));
            var both = tree.FindMaximalUniqueMatches(query, minLength, MumUniqueness.Both);
            add($"{at}.mum-both", HashTriples(both));
            add($"{at}.mum-reference", HashTriples(tree.FindMaximalUniqueMatches(query, minLength, MumUniqueness.Reference)));
            add($"{at}.mum-default", HashTriples(tree.FindMaximalUniqueMatches(query, minLength)) == HashTriples(both) ? "= both" : "DIFFERS FROM Both");
        }
    }

    private static ulong HashLcs(string substring, IReadOnlyList<int> inText, IReadOnlyList<int> inOther)
    {
        ulong h = Fnv(Basis, substring);
        foreach (int p in inText) h = Fnv(h, p);
        h = Fnv(h, -1);
        foreach (int p in inOther) h = Fnv(h, p);
        return h;
    }

    private static string HashTriples(IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> triples)
    {
        ulong h = Basis;
        foreach (var (a, b, c) in triples) h = Fnv(Fnv(Fnv(h, a), b), c);
        return $"{triples.Count}/{Hex(h)}";
    }

    /// <summary>The checksum of the lines (what the concurrency test compares between threads).</summary>
    internal static ulong Checksum(IEnumerable<string> lines)
    {
        ulong h = Basis;
        foreach (string line in lines) h = Fnv(h, line);
        return h;
    }

    /// <summary>The first line in which two digests differ, as <c>label: a against b</c>; null when they are equal.</summary>
    internal static string? Diff(IReadOnlyList<string> a, IReadOnlyList<string> b, string aName, string bName)
    {
        int k = Math.Min(a.Count, b.Count);
        for (int i = 0; i < k; i++)
            if (a[i] != b[i]) return $"{LabelOf(a[i])}: {aName} {ValueOf(a[i])} against {bName} {ValueOf(b[i])}";
        return a.Count == b.Count ? null : $"digest length: {aName} {a.Count} lines against {bName} {b.Count}";
    }

    internal static string LabelOf(string line) => line[..line.IndexOf('=', StringComparison.Ordinal)];

    /// <summary>
    /// The rule of the ESA engines for the representative of the longest repeated substring on a tie (ADR-001 section 7 leaves it implementation-specific, so Ukkonen and v6 are not held to it): the lexicographically smallest of the
    /// tied substrings that <see cref="ISuffixTree.FindAllLongestRepeatedSubstrings"/> lists. Null when the tree keeps it, else the evidence.
    /// </summary>
    internal static string? LrsRuleViolation(ISuffixTree tree)
    {
        string lrs = tree.LongestRepeatedSubstring();
        var all = tree.FindAllLongestRepeatedSubstrings();
        if (all.Count == 0) return lrs.Length == 0 ? null : $"LongestRepeatedSubstring is '{Shorten(lrs)}' but no substring repeats";
        string smallest = all[0].Substring;
        foreach (var (substring, _) in all)
            if (string.CompareOrdinal(substring, smallest) < 0) smallest = substring;
        return lrs == smallest ? null : $"LongestRepeatedSubstring is '{Shorten(lrs)}', the smallest of the {all.Count} tied substrings is '{Shorten(smallest)}'";
    }

    private static string Shorten(string s) => s.Length <= 60 ? s : string.Concat(s.AsSpan(0, 57), "...");

    private static string ValueOf(string line) => line[(line.IndexOf('=', StringComparison.Ordinal) + 1)..];
}

/// <summary>
/// The answers of the contract surface of a tree (<see cref="ISuffixIndex"/>, <see cref="ISuffixTreeNavigation"/>, <see cref="ISuffixArrayAccess"/>, <see cref="IIndexText"/>) to a fixed walk: the root and the first
/// levels of the tree, each node with its count, parent, suffix link, children, edge symbols, extensions and located positions, and a sample of ranks and rank pairs. The same lines come from the in-memory and the
/// file-backed tree, and from every thread of the concurrency test (the support level of the operations that depend on the navigation pack is not part of the lines: it is the one answer that changes with the history of the tree).
/// </summary>
internal static class FuzzContractDigest
{
    private static readonly IndexOperation[] StableSupport =
    [
        IndexOperation.Count, IndexOperation.ExtendRight, IndexOperation.ExtendLeft, IndexOperation.LazyLocate, IndexOperation.TextAccess,
        IndexOperation.Child, IndexOperation.StringDepth, IndexOperation.SuffixArrayAccess, IndexOperation.LcpAccess,
    ];

    private static string Render(IndexLocus l) => $"{l.Lo},{l.Hi},{l.Length},{l.Aux}";

    /// <summary>The lines of <paramref name="tree"/>.</summary>
    internal static List<string> Compute<T>(T tree, string text, ulong seed, int nodeBudget = 120)
        where T : ISuffixTreeNavigation, ISuffixArrayAccess
    {
        var lines = new List<string>();
        void Add(string label, string value) => lines.Add(label + "=" + value);
        int n = text.Length;
        var rng = new SplitMix64(seed ^ 0xC0DEC0DEUL ^ (ulong)n);

        Add("meta", $"{tree.Engine}/{tree.TextLength}/{Render(tree.Root)}/{tree.Count(tree.Root)}");
        foreach (IndexOperation op in Enum.GetValues<IndexOperation>())
        {
            IndexSupport support = tree.Support(op);
            string? note = tree.SupportNote(op);
            if (Array.IndexOf(StableSupport, op) >= 0) Add($"support[{op}]", $"{support}/{note}");
        }

        IIndexText indexText = ((ISuffixIndex)tree).Text;
        var sb = new StringBuilder();
        foreach (long p in new long[] { 0, n / 2, Math.Max(0, n - 1), n, n + 3 }.Distinct())
            sb.Append(indexText.SymbolAt(p)).Append(',');
        var buffer = new char[19];
        int copied = indexText.CopyTo(n / 3, buffer);
        Add("index-text", $"{indexText.Length}/{sb}/{copied}/{new string(buffer, 0, copied)}");

        // ranks
        var ranks = new List<long> { 0, 1, n / 2, Math.Max(0, n - 1), n };
        for (int i = 0; i < 8; i++) ranks.Add(rng.NextInt(n + 1));
        foreach (long r in ranks.Distinct())
        {
            if (r < 0 || r > n) continue;
            long suffix = tree.SuffixAt(r);
            Add($"rank[{r}]", $"{suffix}/{tree.RankOf(suffix)}/{tree.LcpAt(r)}");
        }
        for (int i = 0; i < 10; i++)
        {
            long a = rng.NextInt(n + 1), b = rng.NextInt(n + 1);
            Add($"pair[{a},{b}]", $"{tree.Lcp(a, b)}/{Render(tree.Lca(a, b))}");
        }

        // the first levels of the tree, breadth first
        var queue = new Queue<IndexLocus>();
        queue.Enqueue(tree.Root);
        var positions = new long[16];
        int visited = 0;
        while (queue.Count > 0 && visited < nodeBudget)
        {
            IndexLocus node = queue.Dequeue();
            visited++;
            string at = $"node[{Render(node)}]";
            bool explicitNode = tree.IsExplicitNode(node);
            bool hasParent = tree.TryParent(node, out IndexLocus parent);
            bool hasLink = tree.TrySuffixLink(node, out IndexLocus link);
            int count = tree.GetChildren(node, []);
            var children = new IndexLocus[Math.Min(count, 1_000_000)];
            if (count > 0 && tree.GetChildren(node, children) != count) throw new InvalidOperationException($"GetChildren returned another count at the second call for the node {Render(node)}");
            int located = tree.Locate(node, 0, positions);
            ulong ph = FuzzDigest.Basis;
            for (int i = 0; i < located; i++) ph = FuzzDigest.Fnv(ph, positions[i]);
            int skipped = tree.Locate(node, 3, positions);
            Add(at, $"{tree.Count(node)}/{explicitNode}/{(hasParent ? Render(parent) : "-")}/{(hasLink ? Render(link) : "-")}/{count}/{located}:{ph:x16}/{skipped}");
            int shown = Math.Min(count, 6);
            for (int i = 0; i < shown; i++)
            {
                IndexLocus child = children[i];
                int first = tree.EdgeSymbol(child, 0);
                bool found = tree.TryChild(node, first, out IndexLocus viaChild);
                bool extended = tree.TryExtendRight(node, first, out IndexLocus viaExtend);
                Add($"{at}.child[{i}]", $"{Render(child)}/first={first}/{found && viaChild == child}/{extended}:{(extended ? Render(viaExtend) : "-")}/left={LeftOf(tree, node, first)}");
                if (queue.Count < nodeBudget) queue.Enqueue(child);
            }
            Add($"{at}.absent", $"{tree.TryExtendRight(node, 0xFFFF, out _)}/{tree.TryChild(node, 0x0100, out _)}/{LeftOf(tree, node, 'z')}");
        }
        return lines;
    }

    /// <summary>
    /// Loci that are not loci, with the answer fixed by the contract (not compared with another engine, which shares the code): a point at the depth of the parent of a node but with the interval of the node lies on the parent, not
    /// inside the edge of the node, and is an <see cref="ArgumentException"/>; so is a locus with an interval that is not a node of its depth. Null when every probe was refused.
    /// </summary>
    internal static string? NegativeLoci<T>(T tree, string text, int nodeBudget = 60)
        where T : ISuffixTreeNavigation
    {
        var queue = new Queue<IndexLocus>();
        queue.Enqueue(tree.Root);
        int visited = 0;
        while (queue.Count > 0 && visited < nodeBudget)
        {
            IndexLocus node = queue.Dequeue();
            visited++;
            if (node.Length > 0 && tree.TryParent(node, out IndexLocus parent) && parent.Length > 0 && parent.Length < node.Length)
            {
                var onTheParent = new IndexLocus(node.Lo, node.Hi, parent.Length, node.Aux);
                foreach ((string what, Action call) in new (string, Action)[] { ("IsExplicitNode", () => tree.IsExplicitNode(onTheParent)), ("TryParent", () => tree.TryParent(onTheParent, out _)) })
                {
                    try
                    {
                        call();
                        return $"{what} of {Render(onTheParent)} (a point at the depth of the parent {Render(parent)} with the interval of the node {Render(node)}) did not throw";
                    }
                    catch (ArgumentException)
                    {
                        // refused, as the contract says
                    }
                }
            }
            int count = tree.GetChildren(node, []);
            var children = new IndexLocus[Math.Min(count, 1_000_000)];
            if (count > 0 && tree.GetChildren(node, children) != count) return $"GetChildren of {Render(node)} returned another count at the second call";
            for (int i = 0; i < Math.Min(count, 4); i++) queue.Enqueue(children[i]);
        }
        return null;
    }

    private static string LeftOf<T>(T tree, IndexLocus node, int symbol)
        where T : ISuffixTreeNavigation
    {
        bool ok = tree.TryExtendLeft(node, symbol, out IndexLocus left);
        return ok ? Render(left) : "-";
    }
}
