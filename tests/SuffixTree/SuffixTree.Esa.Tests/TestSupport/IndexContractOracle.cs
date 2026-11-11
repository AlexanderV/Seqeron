using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.TestSupport;

/// <summary>
/// Definitional oracle of the index contract for the enhanced suffix array (test oracle L0, work package 1b.5). It is written from the
/// definitions and shares no code with the engine: the suffix array and the LCP array come from <see cref="NaiveSuffixArray"/> and
/// <see cref="NaiveLcp"/>, everything about the tree is derived from the occurrence sets of strings in <c>text$</c>.
/// <list type="bullet">
/// <item>A locus string is an <c>int[]</c> over the code units of the text; the terminator <c>-1</c> may only be its last symbol.</item>
/// <item><c>occ(w)</c> = the positions <c>i in 0..n</c> where <c>text$[i..i+|w|)</c> equals <c>w</c> (position <c>n</c> is the empty suffix).</item>
/// <item>The node below <c>w</c> is <c>w</c> extended while all occurrences continue with the same symbol, and stopped after a terminator;
/// <c>w</c> is an explicit node when nothing was added (the root <c>""</c> always is, a string ending in <c>$</c> always is). An absent string has no node.</item>
/// <item>The locus of <c>w</c> is <c>(min rank, max rank, |w|, depth of the node below)</c> over the ranks of the suffixes in <c>occ(w)</c>.</item>
/// <item>Count = the number of positions of <c>occ(w)</c> below <c>n</c> (the empty suffix is no text position); the parent of <c>w</c> is its longest
/// proper prefix that is an explicit node; the children of a node are the nodes below <c>w + s</c> for every symbol <c>s</c> that continues some occurrence
/// (terminator first); the suffix link of an internal node <c>aY</c> is <c>Y</c>.</item>
/// </list>
/// The empty text has only the root (the Ukkonen tree reports one node for it).
/// </summary>
public sealed class IndexContractOracle
{
    private readonly int[] _t;

    /// <summary>Builds the oracle of <paramref name="text"/>.</summary>
    public IndexContractOracle(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        Text = text;
        N = text.Length;
        Sa = NaiveSuffixArray.Build(text);
        Lcp = NaiveLcp.Build(text, Sa);
        Isa = new int[N + 1];
        for (int r = 0; r <= N; r++) Isa[Sa[r]] = r;
        _t = new int[N + 1];
        for (int i = 0; i < N; i++) _t[i] = text[i];
        _t[N] = -1;
    }

    /// <summary>The text.</summary>
    public string Text { get; }

    /// <summary>Its length.</summary>
    public int N { get; }

    /// <summary>The naive suffix array, <c>n + 1</c> entries, <c>SA[0] = n</c>.</summary>
    public int[] Sa { get; }

    /// <summary>The inverse of <see cref="Sa"/>.</summary>
    public int[] Isa { get; }

    /// <summary>The naive LCP array, <c>n + 1</c> entries.</summary>
    public int[] Lcp { get; }

    /// <summary>The symbol of <c>text$</c> at <paramref name="position"/> (<c>-1</c> at <c>n</c>).</summary>
    public int SymbolAt(int position) => _t[position];

    /// <summary>The root string.</summary>
    public static int[] RootString => [];

    /// <summary>The string of <paramref name="s"/> as symbols.</summary>
    public static int[] Symbols(string s) => s.Select(c => (int)c).ToArray();

    /// <summary>The string with a terminator appended.</summary>
    public static int[] WithTerminator(int[] w) => [.. w, -1];

    /// <summary>Printable form (terminator as <c>$</c>).</summary>
    public static string Render(int[] w) => w.Length == 0 ? "\"\"" : "\"" + string.Concat(w.Select(s => s < 0 ? "$" : s is >= 32 and < 127 ? ((char)s).ToString() : $"\\u{s:X4}")) + "\"";

    /// <summary>The key of a string (for dictionaries and sets).</summary>
    public static string Key(int[] w) => string.Join(",", w);

    /// <summary>The positions of <c>text$</c> where <paramref name="w"/> occurs.</summary>
    public List<int> Occ(ReadOnlySpan<int> w)
    {
        var result = new List<int>();
        for (int i = 0; i <= N && i + w.Length <= N + 1; i++)
        {
            bool ok = true;
            for (int k = 0; k < w.Length && ok; k++) ok = _t[i + k] == w[k];
            if (ok) result.Add(i);
        }
        return result;
    }

    /// <summary>True when <paramref name="w"/> occurs in <c>text$</c>.</summary>
    public bool Exists(int[] w) => Occ(w).Count > 0;

    /// <summary>The number of text positions at which <paramref name="w"/> starts (the empty suffix is not one).</summary>
    public long Count(int[] w) => Occ(w).Count(i => i < N);

    /// <summary>The text positions of <paramref name="w"/> in ascending order.</summary>
    public List<int> Positions(int[] w) => Occ(w).Where(i => i < N).ToList();

    /// <summary>The depth of the node below <paramref name="w"/> (see the class remarks); <paramref name="w"/> must occur.</summary>
    public int NodeBelowDepth(int[] w)
    {
        List<int> occ = Occ(w);
        if (occ.Count == 0) throw new InvalidOperationException($"{Render(w)} does not occur in '{Text}'");
        if (N == 0) return 0;
        int d = w.Length;
        while (true)
        {
            if (d > 0 && _t[occ[0] + d - 1] == -1) return d;
            int s = _t[occ[0] + d];
            foreach (int i in occ)
                if (_t[i + d] != s) return d;
            d++;
        }
    }

    /// <summary>The string of the node below <paramref name="w"/>.</summary>
    public int[] NodeBelow(int[] w)
    {
        int d = NodeBelowDepth(w);
        int i0 = Occ(w)[0];
        return _t.AsSpan(i0, d).ToArray();
    }

    /// <summary>True when <paramref name="w"/> is an explicit node (root, internal node or leaf).</summary>
    public bool IsExplicit(int[] w) => w.Length == 0 || NodeBelowDepth(w) == w.Length;

    /// <summary>The rank interval of the suffixes that start with <paramref name="w"/>; the ranks must be contiguous (self-check).</summary>
    public (int Lo, int Hi) RankInterval(int[] w)
    {
        List<int> occ = Occ(w);
        if (occ.Count == 0) throw new InvalidOperationException($"{Render(w)} does not occur in '{Text}'");
        int lo = occ.Min(i => Isa[i]);
        int hi = occ.Max(i => Isa[i]);
        if (hi - lo + 1 != occ.Count) throw new InvalidOperationException($"the ranks of {Render(w)} are not contiguous");
        return (lo, hi);
    }

    /// <summary>The locus of <paramref name="w"/> in the representation of the enhanced suffix array (documented on <see cref="IndexLocus"/>).</summary>
    public IndexLocus LocusOf(int[] w)
    {
        if (N == 0)
        {
            if (w.Length != 0) throw new InvalidOperationException("the empty text has only the root");
            return new IndexLocus(0, 0, 0, 0);
        }
        var (lo, hi) = RankInterval(w);
        return new IndexLocus(lo, hi, w.Length, NodeBelowDepth(w));
    }

    /// <summary>The leaf of the suffix of rank <paramref name="rank"/>: the suffix with the terminator (the root for the empty text).</summary>
    public int[] LeafString(int rank) => N == 0 ? [] : _t.AsSpan(Sa[rank], N - Sa[rank] + 1).ToArray();

    /// <summary>The nearest explicit ancestor above <paramref name="w"/>, null for the root.</summary>
    public int[]? Parent(int[] w)
    {
        if (w.Length == 0) return null;
        for (int len = w.Length - 1; len >= 1; len--)
        {
            int[] prefix = w[..len];
            if (IsExplicit(prefix)) return prefix;
        }
        return [];
    }

    /// <summary>The children of the explicit node <paramref name="w"/>: first symbol (<c>-1</c> first) and node string, ascending by symbol.</summary>
    public List<(int Key, int[] Child)> Children(int[] w)
    {
        var result = new List<(int, int[])>();
        if (N == 0 || (w.Length > 0 && w[^1] == -1)) return result;
        List<int> occ = Occ(w);
        foreach (int s in occ.Select(i => _t[i + w.Length]).Distinct().OrderBy(x => x))
            result.Add((s, NodeBelow([.. w, s])));
        return result;
    }

    /// <summary>The suffix link of the internal node <paramref name="w"/> (the string without its first symbol); null for the root and for a leaf.</summary>
    public static int[]? SuffixLink(int[] w)
    {
        if (w.Length == 0 || w[^1] == -1) return null;
        return w[1..];
    }

    /// <summary>The string <c>c + w</c> when it occurs, else null.</summary>
    public int[]? ExtendLeft(int c, int[] w)
    {
        int[] extended = [c, .. w];
        return Exists(extended) ? extended : null;
    }

    /// <summary>The string <c>w + c</c> when it occurs, else null.</summary>
    public int[]? ExtendRight(int[] w, int c)
    {
        int[] extended = [.. w, c];
        return Exists(extended) ? extended : null;
    }

    /// <summary>The longest common prefix of the suffixes of the two ranks (text only, no terminator); equal ranks give the length of the suffix.</summary>
    public int LcpOfRanks(int a, int b)
    {
        int p = Sa[a], q = Sa[b];
        if (a == b) return N - p;
        int k = 0;
        while (p + k < N && q + k < N && Text[p + k] == Text[q + k]) k++;
        return k;
    }

    /// <summary>The lowest common ancestor of the leaves of two ranks: the longest common prefix of the two terminated suffixes (the leaf itself for equal ranks).</summary>
    public int[] LcaString(int a, int b)
    {
        if (N == 0) return [];
        if (a == b) return LeafString(a);
        int p = Sa[a], q = Sa[b];
        int k = 0;
        while (p + k <= N && q + k <= N && _t[p + k] == _t[q + k]) k++;
        return _t.AsSpan(p, k).ToArray();
    }

    /// <summary>Every explicit node of the tree in preorder (the root first, children by ascending key, terminator first).</summary>
    public List<int[]> Nodes()
    {
        var result = new List<int[]>();
        var stack = new Stack<int[]>();
        stack.Push([]);
        while (stack.Count > 0)
        {
            int[] w = stack.Pop();
            result.Add(w);
            var kids = Children(w);
            for (int i = kids.Count - 1; i >= 0; i--) stack.Push(kids[i].Child);
        }
        return result;
    }

    /// <summary>The distinct alphabet of the text, ascending.</summary>
    public int[] Alphabet() => Text.Select(c => (int)c).Distinct().OrderBy(x => x).ToArray();
}
