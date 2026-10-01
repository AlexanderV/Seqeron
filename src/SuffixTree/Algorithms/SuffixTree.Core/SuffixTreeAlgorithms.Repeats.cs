using System;
using System.Collections.Generic;

namespace SuffixTree;

public static partial class SuffixTreeAlgorithms
{
    /// <summary>Left-character class of a position with no matchable left neighbour (start of text or a unique symbol).</summary>
    private const int UniqueLeftClass = -1;

    /// <summary>
    /// Finds every maximal repeated pair of the tree's text with length ≥ <paramref name="minLength"/>
    /// (Gusfield 1997 §7.12): every triple (FirstPosition i, SecondPosition j, Length L), 0-based,
    /// i &lt; j, with text[i..i+L) = text[j..j+L) that is right-maximal (j + L = n or
    /// text[i+L] ≠ text[j+L]) and left-maximal (i = 0 or text[i−1] ≠ text[j−1]). Copies may overlap
    /// (j &lt; i + L, e.g. tandem repeats). Forward strand only — MUMmer 3 <c>repeat-match -f -n minLength</c>
    /// (1-based there).
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Unique symbols.</b> When <paramref name="isUniqueSymbol"/> is given, every character for which
    /// it returns true is treated as matching nothing, not even itself (separators, IUPAC <c>N</c>, …):
    /// a match never extends through such a character (it counts as a mismatch on the right, and a
    /// preceding unique character makes the pair left-maximal). This equals running the definition on
    /// the text with each unique-symbol occurrence replaced by a distinct fresh symbol, as in
    /// <c>RepeatFinder.FindDirectRepeats</c> (non-ACGT unique) and Vmatch/REPuter separator handling.
    /// Without it (default) every character matches itself — <c>repeat-match</c> behaviour, which also
    /// matches <c>N</c> with <c>N</c>.
    /// </para>
    /// <para>
    /// <b>Algorithm</b> (Gusfield 1997 §7.12.3): a bottom-up traversal keeps, for every node, its leaves
    /// in one linked list per left character (one class for "no left character"). When a child is merged
    /// into its parent v of string depth d ≥ minLength, every leaf p of the child is paired with every
    /// leaf q already merged into v whose left character differs (or is undefined): LCP(p, q) = d exactly
    /// (right-maximal) and the pair is left-maximal. Lists are then concatenated in O(1) per class
    /// (smaller class map into the larger). With unique symbols a subtree whose edge holds the first
    /// unique character at string depth u is "exploded" at depth u: all its leaves are paired with each
    /// other at length u and nothing deeper is emitted. Enumeration O(n log σ' + z) for z pairs
    /// (σ' = distinct left characters), plus O(z log z) for the final sort. Every maximal pair is
    /// emitted exactly once (at the lowest common ancestor, or at the explosion point).
    /// </para>
    /// <para>
    /// <b>Output order</b> (deterministic, identical for every tree implementation): ascending by
    /// FirstPosition, then SecondPosition (a pair (i, j) has exactly one maximal length).
    /// </para>
    /// </remarks>
    /// <param name="tree">Any suffix tree (traversed via <see cref="ISuffixTreeDiagnostics.Traverse"/>).</param>
    /// <param name="minLength">Minimum pair length (≥ 1).</param>
    /// <param name="isUniqueSymbol">Optional predicate: characters that never match (null = none).</param>
    /// <exception cref="ArgumentNullException"><paramref name="tree"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLength"/> &lt; 1.</exception>
    public static IReadOnlyList<(int FirstPosition, int SecondPosition, int Length)> FindMaximalRepeatedPairs(
        ISuffixTree tree, int minLength, Func<char, bool>? isUniqueSymbol = null)
    {
        ArgumentNullException.ThrowIfNull(tree);
        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1);
        ITextSource text = tree.Text;
        if (text.Length <= minLength)
            return Array.Empty<(int, int, int)>();

        var collector = new MaximalPairCollector(text, minLength, isUniqueSymbol);
        tree.Traverse(collector);
        var pairs = collector.Pairs;
        pairs.Sort(static (a, b) => a.FirstPosition != b.FirstPosition
            ? a.FirstPosition.CompareTo(b.FirstPosition)
            : a.SecondPosition.CompareTo(b.SecondPosition));
        return pairs;
    }

    /// <summary>
    /// Bottom-up maximal-pair enumeration over <see cref="ISuffixTreeVisitor"/> events
    /// (<c>VisitNode</c> on entry, <c>ExitBranch</c> after the subtree; leaves are exited at once).
    /// </summary>
    private sealed class MaximalPairCollector : ISuffixTreeVisitor
    {
        private sealed class Frame
        {
            public int StartDepth;
            public int EndDepth;
            public int Representative = -1;
            public Dictionary<int, (int Head, int Tail)> Classes = new();
        }

        private readonly int _n;
        private readonly int _minLength;
        private readonly int[] _cut;       // u(p): length of the longest prefix of suffix p free of unique symbols
        private readonly int[] _leftClass; // text[p-1], or UniqueLeftClass
        private readonly int[] _next;      // linked lists of leaf positions
        private readonly List<Frame> _stack = new();
        private bool _rootSeen;
        private bool _leafPending;

        public MaximalPairCollector(ITextSource text, int minLength, Func<char, bool>? isUniqueSymbol)
        {
            _n = text.Length;
            _minLength = minLength;
            _cut = new int[_n];
            _leftClass = new int[_n];
            _next = new int[_n];

            int nextUnique = _n;
            for (int p = _n - 1; p >= 0; p--)
            {
                if (isUniqueSymbol != null && isUniqueSymbol(text[p]))
                    nextUnique = p;
                _cut[p] = nextUnique - p;
            }
            for (int p = 0; p < _n; p++)
            {
                _leftClass[p] = p == 0 || _cut[p - 1] == 0 ? UniqueLeftClass : text[p - 1];
            }
        }

        public List<(int FirstPosition, int SecondPosition, int Length)> Pairs { get; } = new();

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (!_rootSeen)
            {
                _rootSeen = true;
                _stack.Add(new Frame());
                return;
            }

            if (childCount > 0)
            {
                _stack.Add(new Frame { StartDepth = depth, EndDepth = depth + (endIndex - startIndex) });
                return;
            }

            _leafPending = true;
            int p = startIndex - depth;
            if (p < 0 || p >= _n)
                return; // empty-suffix (terminator-only) leaf

            Frame parent = _stack[^1];
            int u = _cut[p];
            // A single exploded leaf forms no pair; merge it into the parent unless the parent's label is cut.
            bool emit = u >= parent.EndDepth && parent.EndDepth >= _minLength;
            int a = _leftClass[p];
            if (emit)
            {
                foreach (var (b, list) in parent.Classes)
                {
                    if (a != b || a == UniqueLeftClass)
                        EmitLeafWithList(p, list.Head, parent.EndDepth);
                }
            }

            _next[p] = -1;
            if (parent.Classes.TryGetValue(a, out var existing))
            {
                _next[existing.Tail] = p;
                parent.Classes[a] = (existing.Head, p);
            }
            else
            {
                parent.Classes[a] = (p, p);
            }
            if (parent.Representative < 0)
                parent.Representative = p;
        }

        public void EnterBranch(int key) { }

        public void ExitBranch()
        {
            if (_leafPending)
            {
                _leafPending = false;
                return;
            }

            Frame child = _stack[^1];
            _stack.RemoveAt(_stack.Count - 1);
            if (child.Representative < 0)
                return;

            Frame parent = _stack[^1];
            int u = _cut[child.Representative];
            if (u >= child.StartDepth && u < child.EndDepth && u >= _minLength)
                EmitAllPairs(child, u);

            bool emit = u >= parent.EndDepth && parent.EndDepth >= _minLength;
            MergeInto(parent, child, emit);
        }

        private void MergeInto(Frame parent, Frame child, bool emit)
        {
            var large = parent.Classes;
            var small = child.Classes;
            if (small.Count > large.Count)
                (large, small) = (small, large);

            if (emit)
            {
                int length = parent.EndDepth;
                foreach (var (a, la) in small)
                {
                    foreach (var (b, lb) in large)
                    {
                        if (a != b || a == UniqueLeftClass)
                            EmitCross(la.Head, lb.Head, length);
                    }
                }
            }

            foreach (var (a, la) in small)
            {
                if (large.TryGetValue(a, out var lb))
                {
                    _next[lb.Tail] = la.Head;
                    large[a] = (lb.Head, la.Tail);
                }
                else
                {
                    large[a] = la;
                }
            }

            parent.Classes = large;
            if (parent.Representative < 0)
                parent.Representative = child.Representative;
        }

        /// <summary>All left-maximal pairs inside one subtree whose leaves all mismatch at depth <paramref name="length"/>.</summary>
        private void EmitAllPairs(Frame frame, int length)
        {
            var lists = new List<(int Class, int Head)>(frame.Classes.Count);
            foreach (var (c, l) in frame.Classes)
                lists.Add((c, l.Head));

            for (int x = 0; x < lists.Count; x++)
            {
                for (int y = x + 1; y < lists.Count; y++)
                    EmitCross(lists[x].Head, lists[y].Head, length);

                if (lists[x].Class == UniqueLeftClass)
                {
                    for (int p = lists[x].Head; p >= 0; p = _next[p])
                    {
                        for (int q = _next[p]; q >= 0; q = _next[q])
                            Add(p, q, length);
                    }
                }
            }
        }

        private void EmitCross(int headA, int headB, int length)
        {
            for (int p = headA; p >= 0; p = _next[p])
                EmitLeafWithList(p, headB, length);
        }

        private void EmitLeafWithList(int p, int head, int length)
        {
            for (int q = head; q >= 0; q = _next[q])
                Add(p, q, length);
        }

        private void Add(int p, int q, int length)
            => Pairs.Add(p < q ? (p, q, length) : (q, p, length));
    }

    /// <summary>
    /// Every distinct longest common substring of the tree's text and <paramref name="other"/> (all
    /// length ties), each with all 0-based start positions in the text and in <paramref name="other"/>,
    /// both ascending and duplicate-free; ordered by first occurrence in <paramref name="other"/> (the
    /// first entry is the canonical substring of
    /// <see cref="ISuffixTreeAnalysis.FindAllLongestCommonSubstrings"/>). Empty when no character is shared.
    /// </summary>
    /// <remarks>
    /// Matching statistics (Chang &amp; Lawler 1994; Gusfield 1997 §7.8): ms(i) is the length of the
    /// longest suffix of other[0..i] occurring in the text. With L = max ms(i), the occurrences in
    /// <paramref name="other"/> of every length-L common substring end exactly at the i with ms(i) = L
    /// (ms(i) ≥ L there and ms ≤ L everywhere), so grouping those ends by their substring yields all
    /// distinct LCS strings with all their query positions; text positions are the leaves below each
    /// substring's locus. O(|other| + Σ occurrences · log).
    /// </remarks>
    public static IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)>
        FindAllDistinctLongestCommonSubstrings<TNode, TNav>(ref TNav nav, string other)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        ArgumentNullException.ThrowIfNull(other);
        var (maxLen, matches) = CollectMaximalMatchingStatistics<TNode, TNav>(ref nav, other, firstOnly: false);
        var results = new List<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)>();
        if (maxLen == 0)
            return results;

        var index = new Dictionary<string, int>(StringComparer.Ordinal);
        foreach (var match in matches)
        {
            int start = match.MatchEndInOther - maxLen + 1;
            string substring = other.Substring(start, maxLen);
            if (index.TryGetValue(substring, out int slot))
            {
                ((List<int>)results[slot].PositionsInOther).Add(start);
                continue;
            }

            var inText = new List<int>();
            nav.CollectLeaves(match.Node, match.DepthFromRoot, inText);
            SortAndDeduplicateInPlace(inText);
            index[substring] = results.Count;
            results.Add((substring, inText, new List<int> { start }));
        }

        return results;
    }

    /// <summary>
    /// Definition-level counterpart of <see cref="FindAllDistinctLongestCommonSubstrings{TNode, TNav}"/>
    /// (same output and order), used as the default interface implementation: O(|text|·|other|) dynamic
    /// programming for the length, then a scan of the text and <paramref name="other"/> for each distinct
    /// length-L substring.
    /// </summary>
    public static IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)>
        FindAllDistinctLongestCommonSubstringsByDefinition(ITextSource text, string other)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(other);
        var results = new List<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)>();
        string reference = text.Substring(0, text.Length);
        int n = reference.Length, m = other.Length;
        if (n == 0 || m == 0)
            return results;

        // prev[j] / cur[j]: length of the common suffix of reference[..i] and other[..j].
        var prev = new int[m + 1];
        var cur = new int[m + 1];
        int maxLen = 0;
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= m; j++)
            {
                cur[j] = reference[i - 1] == other[j - 1] ? prev[j - 1] + 1 : 0;
                if (cur[j] > maxLen)
                    maxLen = cur[j];
            }
            (prev, cur) = (cur, prev);
        }
        if (maxLen == 0)
            return results;

        var inText = new HashSet<string>(StringComparer.Ordinal);
        for (int r = 0; r + maxLen <= n; r++)
            inText.Add(reference.Substring(r, maxLen));

        var index = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int q = 0; q + maxLen <= m; q++)
        {
            string s = other.Substring(q, maxLen);
            if (!inText.Contains(s))
                continue;
            if (index.TryGetValue(s, out int slot))
            {
                ((List<int>)results[slot].PositionsInOther).Add(q);
                continue;
            }

            var textPositions = new List<int>();
            for (int r = 0; r + maxLen <= n; r++)
            {
                if (string.CompareOrdinal(reference, r, s, 0, maxLen) == 0)
                    textPositions.Add(r);
            }
            index[s] = results.Count;
            results.Add((s, textPositions, new List<int> { q }));
        }
        return results;
    }

    /// <summary>
    /// Streams <paramref name="other"/> against the tree (matching statistics) and returns the maximal
    /// ms value with the loci of every end position attaining it (only the first when
    /// <paramref name="firstOnly"/>), in increasing end position. O(|other|) amortised.
    /// </summary>
    private static (int MaxLength, List<(TNode Node, int MatchEndInOther, int DepthFromRoot)> Matches)
        CollectMaximalMatchingStatistics<TNode, TNav>(ref TNav nav, string other, bool firstOnly)
        where TNav : struct, ISuffixTreeNavigator<TNode>
    {
        ReadOnlySpan<char> otherSpan = other.AsSpan();
        var bestMatches = new List<(TNode Node, int MatchEndInOther, int DepthFromRoot)>();
        if (otherSpan.Length == 0 || nav.Text.Length == 0)
            return (0, bestMatches);

        int maxLen = 0;
        TNode currentNode = nav.Root;
        TNode currentEdge = nav.NullNode;
        int edgeOffset = 0;
        int currentMatchLen = 0;

        // Depth to END of currentNode's edge (= GetNodeDepth(currentNode)).
        // For root this is 0. Updated on suffix link follow (-1), edge
        // completion (+edgeLen), and rescan (+edgeLen per full edge).
        int currentNodeDepth = 0;

        for (int i = 0; i < otherSpan.Length; i++)
        {
            int c = otherSpan[i];

            while (true)
            {
                if (TryConsumeSymbol(ref nav, c, ref currentNode, ref currentEdge, ref edgeOffset, ref currentMatchLen, ref currentNodeDepth))
                    break;

                // Cannot extend — follow suffix link
                if (currentMatchLen == 0) break;
                FollowSuffixLinkAndRescan(ref nav, otherSpan, i, ref currentNode, ref currentEdge, ref edgeOffset, ref currentMatchLen, ref currentNodeDepth);
            }

            bool better = currentMatchLen > maxLen;
            if (better || (currentMatchLen == maxLen && maxLen > 0 && !firstOnly))
            {
                if (better)
                {
                    maxLen = currentMatchLen;
                    bestMatches.Clear();
                }
                TNode matchNode = nav.IsNull(currentEdge) ? currentNode : currentEdge;
                // DepthFromRoot of matchNode:
                // - currentNode: currentNodeDepth - LengthOf(currentNode)
                // - currentEdge: currentNodeDepth (= depth to END of parent = depth to START of child)
                int matchDFR = nav.IsNull(currentEdge)
                    ? currentNodeDepth - nav.LengthOf(currentNode)
                    : currentNodeDepth;
                bestMatches.Add((matchNode, i, matchDFR));
            }
        }

        return (maxLen, bestMatches);
    }
}
