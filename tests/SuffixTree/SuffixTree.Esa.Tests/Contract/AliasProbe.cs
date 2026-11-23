using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>What one call of a locus-returning member produced: its answer and the locus it left in the result variable.</summary>
internal readonly record struct Outcome(bool Ok, IndexLocus Locus)
{
    /// <inheritdoc />
    public override string ToString() => $"{Ok} {Locus}";
}

/// <summary>What one call of <c>GetChildren</c> produced: the returned number and the whole destination.</summary>
internal readonly record struct ChildrenOutcome(int Count, IndexLocus[] Buffer)
{
    /// <inheritdoc />
    public override string ToString() => $"{Count} in [{string.Join("; ", Buffer)}]";
}

/// <summary>
/// The two ways to call every member of the index contract that takes a locus by <c>in</c> and hands a locus back (work package 1b.12, finding R01 of the review of 1b.10): <c>...Apart</c> passes a variable for
/// the locus and another one for the result, <c>...Alias</c> passes ONE variable as both (<c>tree.TryExtendRight(l, c, out l)</c>, the natural loop that advances a locus; C# allows it and the callee then sees one
/// storage location under two names), and <see cref="ChildrenAlias"/> puts the node into the very array that receives the children. The contract says the two forms give the same answer.
/// </summary>
internal static class AliasProbe
{
    /// <summary>The content of a destination slot that a call must leave alone.</summary>
    internal static readonly IndexLocus Filler = new(7, 7, 7, 7);

    internal static Outcome ExtendRightApart(ISuffixIndex index, IndexLocus locus, int symbol)
    {
        bool ok = index.TryExtendRight(in locus, symbol, out IndexLocus result);
        return new Outcome(ok, result);
    }

    internal static Outcome ExtendRightAlias(ISuffixIndex index, IndexLocus locus, int symbol)
    {
        bool ok = index.TryExtendRight(in locus, symbol, out locus);
        return new Outcome(ok, locus);
    }

    internal static Outcome ExtendRightApart(ISuffixIndex index, IndexLocus locus, string pattern)
    {
        bool ok = index.TryExtendRight(in locus, pattern.AsSpan(), out IndexLocus result);
        return new Outcome(ok, result);
    }

    internal static Outcome ExtendRightAlias(ISuffixIndex index, IndexLocus locus, string pattern)
    {
        bool ok = index.TryExtendRight(in locus, pattern.AsSpan(), out locus);
        return new Outcome(ok, locus);
    }

    internal static Outcome ExtendLeftApart(ISuffixIndex index, IndexLocus locus, int symbol)
    {
        bool ok = index.TryExtendLeft(in locus, symbol, out IndexLocus result);
        return new Outcome(ok, result);
    }

    internal static Outcome ExtendLeftAlias(ISuffixIndex index, IndexLocus locus, int symbol)
    {
        bool ok = index.TryExtendLeft(in locus, symbol, out locus);
        return new Outcome(ok, locus);
    }

    internal static Outcome ParentApart(ISuffixTreeNavigation nav, IndexLocus locus)
    {
        bool ok = nav.TryParent(in locus, out IndexLocus result);
        return new Outcome(ok, result);
    }

    internal static Outcome ParentAlias(ISuffixTreeNavigation nav, IndexLocus locus)
    {
        bool ok = nav.TryParent(in locus, out locus);
        return new Outcome(ok, locus);
    }

    internal static Outcome ChildApart(ISuffixTreeNavigation nav, IndexLocus node, int key)
    {
        bool ok = nav.TryChild(in node, key, out IndexLocus result);
        return new Outcome(ok, result);
    }

    internal static Outcome ChildAlias(ISuffixTreeNavigation nav, IndexLocus node, int key)
    {
        bool ok = nav.TryChild(in node, key, out node);
        return new Outcome(ok, node);
    }

    internal static Outcome SuffixLinkApart(ISuffixTreeNavigation nav, IndexLocus node)
    {
        bool ok = nav.TrySuffixLink(in node, out IndexLocus result);
        return new Outcome(ok, result);
    }

    internal static Outcome SuffixLinkAlias(ISuffixTreeNavigation nav, IndexLocus node)
    {
        bool ok = nav.TrySuffixLink(in node, out node);
        return new Outcome(ok, node);
    }

    /// <summary><c>GetChildren</c> into a fresh destination of <paramref name="size"/> slots filled with <see cref="Filler"/>.</summary>
    internal static ChildrenOutcome ChildrenApart(ISuffixTreeNavigation nav, IndexLocus node, int size)
    {
        IndexLocus[] buffer = Fresh(size);
        int count = nav.GetChildren(in node, buffer);
        return new ChildrenOutcome(count, buffer);
    }

    /// <summary><c>GetChildren</c> with the node stored in slot <paramref name="slot"/> of the destination and passed from there: the children overwrite the very variable that holds the node.</summary>
    internal static ChildrenOutcome ChildrenAlias(ISuffixTreeNavigation nav, IndexLocus node, int size, int slot)
    {
        IndexLocus[] buffer = Fresh(size);
        buffer[slot] = node;
        int count = nav.GetChildren(in buffer[slot], buffer);
        return new ChildrenOutcome(count, buffer);
    }

    /// <summary>The destination that a correct <c>GetChildren</c> leaves: the children in front when they fit, nothing written otherwise, every other slot as it was (the node in <paramref name="slot"/>, <see cref="Filler"/> elsewhere).</summary>
    internal static IndexLocus[] ExpectedBuffer(IReadOnlyList<IndexLocus> children, int size, IndexLocus node, int slot)
    {
        IndexLocus[] buffer = Fresh(size);
        if (slot >= 0) buffer[slot] = node;
        if (children.Count <= size)
            for (int i = 0; i < children.Count; i++) buffer[i] = children[i];
        return buffer;
    }

    private static IndexLocus[] Fresh(int size)
    {
        var buffer = new IndexLocus[size];
        Array.Fill(buffer, Filler);
        return buffer;
    }
}

/// <summary>
/// The loci of a small text for the member-by-member tests: the root, every distinct non-empty prefix of every terminated suffix of the text (the points inside edges, the explicit nodes, the leaves with and without their terminator
/// and the leaf of the empty suffix), as strings of the definitional oracle.
/// </summary>
internal static class AliasCatalog
{
    /// <summary>Every locus string of the text, the root first.</summary>
    internal static List<int[]> Loci(IndexContractOracle oracle)
    {
        var seen = new HashSet<string>(StringComparer.Ordinal) { IndexContractOracle.Key([]) };
        List<int[]> all = [[]];
        if (oracle.N == 0) return all; // the empty text has only the root
        for (int i = 0; i <= oracle.N; i++)
        {
            for (int length = 1; i + length <= oracle.N + 1; length++)
            {
                int[] w = new int[length];
                for (int k = 0; k < length; k++) w[k] = oracle.SymbolAt(i + k);
                if (seen.Add(IndexContractOracle.Key(w))) all.Add(w);
            }
        }
        return all;
    }

    /// <summary>The symbols to try next to a locus: the alphabet of the text, a symbol that is not in it, the code units at both ends of the range, and the values that are no code unit.</summary>
    internal static int[] Symbols(IndexContractOracle oracle)
        => [.. oracle.Alphabet(), 'z', 0, 0xFFFF, -1, -2, 65536, int.MaxValue, int.MinValue];

    /// <summary>The patterns to append to a locus: every substring of the text up to <paramref name="maxLength"/>, a pattern that is absent, one that occurs but not behind most loci, and the empty one.</summary>
    internal static List<string> Patterns(string text, int maxLength)
    {
        var seen = new HashSet<string>(StringComparer.Ordinal) { string.Empty };
        var all = new List<string> { string.Empty };
        for (int i = 0; i < text.Length; i++)
            for (int length = 1; length <= Math.Min(maxLength, text.Length - i); length++)
                if (seen.Add(text.Substring(i, length))) all.Add(text.Substring(i, length));
        foreach (string absent in new[] { "q", "zz", text + "q", "\u0001" })
            if (seen.Add(absent)) all.Add(absent);
        return all;
    }

    /// <summary>The expected answer of a member for the string <paramref name="expected"/> (null: it must answer false and leave the default locus).</summary>
    internal static Outcome Expect(IndexContractOracle oracle, int[]? expected)
        => expected is null ? new Outcome(false, default) : new Outcome(true, oracle.LocusOf(expected));
}
