using System;
using System.Collections.Generic;

namespace SuffixTree;

public static partial class SuffixTreeAlgorithms
{
    /// <summary>
    /// Every distinct longest substring that occurs in at least <paramref name="minSupport"/> of the
    /// <paramref name="texts"/> (Gusfield 1997 §7.6, the k-common substring problem; with
    /// minSupport = texts.Count the longest common substring of all strings — Rosalind LCSM), sorted
    /// ordinally. Empty when no non-empty substring reaches the support (e.g. an empty input text with
    /// minSupport = texts.Count).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The texts are concatenated as t₁ $₁ t₂ $₂ … t_k $_k with k distinct separator characters that
    /// occur in no text (<see cref="BuildGeneralizedText"/>), and a suffix tree of the concatenation is
    /// built with <paramref name="buildTree"/> (generalized suffix tree). Each separator occurs once, so
    /// no internal node's path label contains one: the path label of an internal node v is a substring
    /// of the inputs, and the number C(v) of distinct texts owning a leaf below v (Hui 1992 colour-set
    /// size: leaves counted in depth-first order minus one per pair of consecutive same-text leaves at
    /// their lowest common ancestor) is the number of texts containing it. l(q) = the largest string
    /// depth of an internal node with C(v) ≥ q; a length-l(q) substring with support ≥ q always ends at
    /// such a node (otherwise the node below it would be deeper with the same support), so the answer is
    /// the set of path labels of internal nodes of depth l(q) with C(v) ≥ q. For minSupport = 1 the
    /// answer is the longest input text(s).
    /// O(N log N) for N = total length (the tree traversal; ancestor lookup by binary search).
    /// </para>
    /// </remarks>
    /// <param name="texts">The input strings (k ≥ 1; elements non-null).</param>
    /// <param name="minSupport">Minimum number of texts that must contain the substring (1 ≤ q ≤ k).</param>
    /// <param name="buildTree">Builds a suffix tree over a string (e.g. <c>SuffixTree.Build</c>); the tree is
    /// disposed after use when it implements <see cref="IDisposable"/>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="texts"/> or <paramref name="buildTree"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="texts"/> is empty or contains null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minSupport"/> outside [1, k].</exception>
    public static IReadOnlyList<string> FindLongestCommonSubstrings(
        IReadOnlyList<string> texts, int minSupport, Func<string, ISuffixTree> buildTree)
    {
        ValidateTexts(texts);
        ArgumentNullException.ThrowIfNull(buildTree);
        ArgumentOutOfRangeException.ThrowIfLessThan(minSupport, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(minSupport, texts.Count);

        if (minSupport == 1)
        {
            int longest = 0;
            foreach (string t in texts)
                longest = Math.Max(longest, t.Length);
            var whole = new SortedSet<string>(StringComparer.Ordinal);
            if (longest > 0)
            {
                foreach (string t in texts)
                {
                    if (t.Length == longest)
                        whole.Add(t);
                }
            }
            return new List<string>(whole);
        }

        var collector = RunColorCounting(texts, buildTree, out string generalized);
        int best = 0;
        foreach (var node in collector.Nodes)
        {
            if (node.Support >= minSupport && node.Depth > best)
                best = node.Depth;
        }

        var result = new SortedSet<string>(StringComparer.Ordinal);
        if (best > 0)
        {
            foreach (var node in collector.Nodes)
            {
                if (node.Support >= minSupport && node.Depth == best)
                    result.Add(generalized.Substring(node.EdgeEnd - best, best));
            }
        }
        return new List<string>(result);
    }

    /// <summary>
    /// The k-common substring lengths of Gusfield 1997 §7.6: element q (1 ≤ q ≤ k) is l(q), the length
    /// of the longest substring occurring in at least q of the <paramref name="texts"/>; element 0 is
    /// unused (equals l(1)). Non-increasing in q. l(1) is the longest text length.
    /// </summary>
    /// <param name="texts">The input strings (k ≥ 1; elements non-null).</param>
    /// <param name="buildTree">Builds a suffix tree over a string (e.g. <c>SuffixTree.Build</c>).</param>
    public static int[] LongestCommonSubstringLengthsBySupport(
        IReadOnlyList<string> texts, Func<string, ISuffixTree> buildTree)
    {
        ValidateTexts(texts);
        ArgumentNullException.ThrowIfNull(buildTree);
        int k = texts.Count;
        var lengths = new int[k + 1];
        foreach (string t in texts)
            lengths[1] = Math.Max(lengths[1], t.Length);

        if (k >= 2)
        {
            var collector = RunColorCounting(texts, buildTree, out _);
            foreach (var node in collector.Nodes)
            {
                int c = Math.Min(node.Support, k);
                if (c >= 2 && node.Depth > lengths[c])
                    lengths[c] = node.Depth;
            }
            for (int q = k - 1; q >= 2; q--)
                lengths[q] = Math.Max(lengths[q], lengths[q + 1]);
        }
        lengths[0] = lengths[1];
        return lengths;
    }

    /// <summary>
    /// Concatenates <paramref name="texts"/> as t₁ $₁ t₂ $₂ … t_k $_k with k distinct separator characters
    /// that occur in none of the texts (taken from the Private Use Area U+E000 upward, then any other
    /// unused UTF-16 code unit). <c>Owner[p]</c> is the index of the text containing position p, −1 for a
    /// separator; <c>Starts[i]</c> is the offset of text i.
    /// </summary>
    /// <exception cref="ArgumentException">More texts than unused code units.</exception>
    public static (string Text, int[] Owner, int[] Starts) BuildGeneralizedText(IReadOnlyList<string> texts)
    {
        ValidateTexts(texts);
        var used = new bool[char.MaxValue + 1];
        long total = 0;
        foreach (string t in texts)
        {
            foreach (char c in t)
                used[c] = true;
            total += t.Length + 1L;
        }
        if (total > int.MaxValue)
            throw new ArgumentException("The concatenated texts exceed the maximum string length.", nameof(texts));

        var separators = new char[texts.Count];
        int found = 0;
        for (int c = 0xE000; c <= char.MaxValue && found < separators.Length; c++)
        {
            if (!used[c]) separators[found++] = (char)c;
        }
        for (int c = 1; c < 0xE000 && found < separators.Length; c++)
        {
            if (!used[c]) separators[found++] = (char)c;
        }
        if (found < separators.Length)
            throw new ArgumentException("Not enough unused characters to separate the texts.", nameof(texts));

        var chars = new char[total];
        var owner = new int[total];
        var starts = new int[texts.Count];
        int pos = 0;
        for (int i = 0; i < texts.Count; i++)
        {
            starts[i] = pos;
            string t = texts[i];
            t.CopyTo(0, chars, pos, t.Length);
            Array.Fill(owner, i, pos, t.Length);
            pos += t.Length;
            chars[pos] = separators[i];
            owner[pos] = -1;
            pos++;
        }
        return (new string(chars), owner, starts);
    }

    private static void ValidateTexts(IReadOnlyList<string> texts)
    {
        ArgumentNullException.ThrowIfNull(texts);
        if (texts.Count == 0)
            throw new ArgumentException("At least one text is required.", nameof(texts));
        for (int i = 0; i < texts.Count; i++)
        {
            if (texts[i] is null)
                throw new ArgumentException($"Text {i} is null.", nameof(texts));
        }
    }

    /// <summary>Builds the generalized suffix tree, runs the colour counting and disposes the tree if disposable.</summary>
    private static ColorCountingCollector RunColorCounting(
        IReadOnlyList<string> texts, Func<string, ISuffixTree> buildTree, out string generalized)
    {
        var (text, owner, _) = BuildGeneralizedText(texts);
        generalized = text;
        ISuffixTree tree = buildTree(text) ?? throw new InvalidOperationException("buildTree returned null.");
        try
        {
            if (tree.Text.Length != text.Length)
                throw new InvalidOperationException("buildTree returned a tree over a different text.");
            var collector = new ColorCountingCollector(owner, texts.Count);
            tree.Traverse(collector);
            return collector;
        }
        finally
        {
            (tree as IDisposable)?.Dispose();
        }
    }

    /// <summary>
    /// Hui (1992) colour-set sizes for every internal node over <see cref="ISuffixTreeVisitor"/> events:
    /// C(v) = coloured leaves below v − duplicates, where a duplicate is charged to the lowest common
    /// ancestor of two consecutive (in depth-first order) leaves of the same text.
    /// </summary>
    private sealed class ColorCountingCollector : ISuffixTreeVisitor
    {
        private sealed class Frame
        {
            public int EntryTime;
            public int EndDepth;
            public int EdgeEnd;
            public int Leaves;
            public int Duplicates;
        }

        private readonly int[] _owner;
        private readonly int[] _lastLeafTime;
        private readonly List<Frame> _stack = new();
        private int _clock;
        private bool _rootSeen;
        private bool _leafPending;

        public ColorCountingCollector(int[] owner, int textCount)
        {
            _owner = owner;
            _lastLeafTime = new int[textCount];
            Array.Fill(_lastLeafTime, -1);
        }

        /// <summary>Every internal non-root node: string depth, colour-set size, end offset of its edge.</summary>
        public List<(int Depth, int Support, int EdgeEnd)> Nodes { get; } = new();

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (!_rootSeen)
            {
                _rootSeen = true;
                _stack.Add(new Frame { EntryTime = _clock });
                return;
            }

            if (childCount > 0)
            {
                _stack.Add(new Frame
                {
                    EntryTime = _clock,
                    EndDepth = depth + (endIndex - startIndex),
                    EdgeEnd = endIndex,
                });
                return;
            }

            _leafPending = true;
            int p = startIndex - depth;
            if (p < 0 || p >= _owner.Length)
                return;
            int color = _owner[p];
            if (color < 0)
                return;

            int time = ++_clock;
            _stack[^1].Leaves++;
            int previous = _lastLeafTime[color];
            if (previous >= 0)
            {
                // Deepest open ancestor entered before the previous same-colour leaf = their LCA.
                int lo = 0, hi = _stack.Count - 1;
                while (lo < hi)
                {
                    int mid = (lo + hi + 1) >> 1;
                    if (_stack[mid].EntryTime < previous) lo = mid;
                    else hi = mid - 1;
                }
                _stack[lo].Duplicates++;
            }
            _lastLeafTime[color] = time;
        }

        public void EnterBranch(int key) { }

        public void ExitBranch()
        {
            if (_leafPending)
            {
                _leafPending = false;
                return;
            }

            Frame frame = _stack[^1];
            _stack.RemoveAt(_stack.Count - 1);
            Nodes.Add((frame.EndDepth, frame.Leaves - frame.Duplicates, frame.EdgeEnd));
            Frame parent = _stack[^1];
            parent.Leaves += frame.Leaves;
            parent.Duplicates += frame.Duplicates;
        }
    }
}
