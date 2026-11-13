using System.Reflection;
using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Ph5;

namespace SuffixTree.EmProto.Tests.Ph5;

/// <summary>
/// The in-RAM oracle of PH5 for one text (gate G1 of the prototype E5, ADR-004 section 5): the <see cref="EsaIndex"/> (SA-IS, LCP kernel, child-table scan of the product), the raw arrays of its single-array child
/// layouts (read through reflection, because the fields are private), the three plain arrays <c>next / up / down</c> of the reference layout, the previous / next smaller values of the <see cref="LcpPyramid"/>
/// and the header statistics. Everything PH5 writes is compared with these.
/// </summary>
internal sealed class ChildReference
{
    private ChildReference(string text, EsaIndex index)
    {
        Text = text;
        Index = index;
        N = text.Length;
        Tagged = RawArray(index.WithLayout(ChildLayoutKind.Tagged).Child.Tagged, "_cld");
        Packed = RawArray(index.WithLayout(ChildLayoutKind.Packed).Child.Packed, "_cld");
        Reference3Child reference3 = index.WithLayout(ChildLayoutKind.Reference3).Child.Reference3;
        Next = Ints(reference3, "_next");
        Up = Ints(reference3, "_up");
        Down = Ints(reference3, "_down");
        using LcpPyramid pyramid = LcpPyramid.Build(index.Lcp, 4);
        Psv = new int[N + 1];
        Nsv = new int[N + 1];
        for (int p = 0; p <= N; p++)
        {
            Psv[p] = pyramid.PreviousSmaller(p);
            Nsv[p] = pyramid.NextSmaller(p);
        }
    }

    internal string Text { get; }

    internal EsaIndex Index { get; }

    internal int N { get; }

    internal CompactLcp Lcp => Index.Lcp;

    internal EsaStatistics Statistics => Index.Statistics;

    /// <summary>The raw field array of the in-RAM tagged layout.</summary>
    internal uint[] Tagged { get; }

    /// <summary>The raw field array of the in-RAM packed layout.</summary>
    internal uint[] Packed { get; }

    internal int[] Next { get; }

    internal int[] Up { get; }

    internal int[] Down { get; }

    /// <summary>Previous smaller value (strict), -1 for none.</summary>
    internal int[] Psv { get; }

    /// <summary>Next smaller value (strict), n + 1 for none.</summary>
    internal int[] Nsv { get; }

    internal static ChildReference Of(string text) => new(text, EsaIndex.Build(text, ChildLayoutKind.Tagged));

    private static uint[] RawArray<TLayout>(TLayout layout, string field)
        where TLayout : struct
        => (uint[])typeof(TLayout).GetField(field, BindingFlags.NonPublic | BindingFlags.Instance)!.GetValue(layout)!;

    private static int[] Ints<TLayout>(TLayout layout, string field)
        where TLayout : struct
        => (int[])typeof(TLayout).GetField(field, BindingFlags.NonPublic | BindingFlags.Instance)!.GetValue(layout)!;

    /// <summary>Writes the compact LCP in the two-file form of PH4 (base array, pairs <c>(rank, value)</c>) and returns the two paths.</summary>
    internal (string Base, string Overflow) WriteLcpFiles(string dir, string name)
    {
        string basePath = Path.Combine(dir, name + ".lcp");
        string overflowPath = Path.Combine(dir, name + ".ovf");
        CompactLcp c = Lcp;
        byte[] bytes = c.Kind switch
        {
            LcpKind.U8 => c.BaseU8!,
            LcpKind.U16 => MemoryMarshal.AsBytes(c.BaseU16.AsSpan()).ToArray(),
            _ => MemoryMarshal.AsBytes(c.BaseU32.AsSpan()).ToArray(),
        };
        File.WriteAllBytes(basePath, bytes);
        var pairs = new uint[c.OverflowCount * 2];
        for (int i = 0; i < c.OverflowCount; i++)
        {
            pairs[2 * i] = (uint)c.OverflowRanks[i];
            pairs[(2 * i) + 1] = (uint)c.OverflowValues[i];
        }
        File.WriteAllBytes(overflowPath, MemoryMarshal.AsBytes(pairs.AsSpan()).ToArray());
        return (basePath, overflowPath);
    }
}

/// <summary>Texts with a known shape of the compact LCP.</summary>
internal static class ChildSamples
{
    /// <summary>Random DNA with a block of 600 symbols twice: about 700 values at or above 255, so the base array stays U8 (16 x escapes &lt;= n) and the overflow list is not empty.</summary>
    internal static string U8WithEscapes()
    {
        string block = global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(600, 11);
        return global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(20_000, 5) + block + global::SuffixTree.Esa.Tests.TestSupport.StringFamilies.RandomDna(100, 6) + block;
    }
}

/// <summary>Decoders of the files PH5 writes.</summary>
internal static class ChildFiles
{
    internal static uint[] ReadUInts(string path) => MemoryMarshal.Cast<byte, uint>(File.ReadAllBytes(path)).ToArray();

    /// <summary>A PSV or NSV file as ints: the stored 0xFFFFFFFF (no smaller value to the left) is -1.</summary>
    internal static int[] ReadSmaller(string path) => ReadUInts(path).Select(v => v == uint.MaxValue ? -1 : (int)v).ToArray();

    internal static int ValueOf(uint tagged) => (int)(tagged & ((1u << 30) - 1));

    /// <summary>The tagged single array as the three plain arrays of the reference layout (0 = undefined).</summary>
    internal static (int[] Next, int[] Up, int[] Down) TaggedToThree(uint[] cld)
    {
        var next = new int[cld.Length];
        var up = new int[cld.Length];
        var down = new int[cld.Length];
        for (int k = 0; k < cld.Length; k++)
        {
            switch (cld[k] >> 30)
            {
                case 1: next[k] = ValueOf(cld[k]); break;
                case 2: up[k] = ValueOf(cld[k]); break;
                case 3: down[k] = ValueOf(cld[k]); break;
                default: break;
            }
        }
        return (next, up, down);
    }

    /// <summary>The dense LCP values of a compact LCP as uints (rank 0 included).</summary>
    internal static uint[] Dense(CompactLcp lcp) => lcp.ToArray().Select(v => (uint)v).ToArray();
}

/// <summary>Read side of a child structure, as <see cref="IChildTable"/> but over decoded arrays (the test decoders of layout A and layout F).</summary>
internal interface IChildView
{
    int FirstLIndex(int lb, int rb, bool lastChild);

    int NextLIndex(int k, int rb);
}

/// <summary>Layout A, tagged codec: the field kind is in bits 30-31.</summary>
internal sealed class TaggedView(uint[] cld) : IChildView
{
    public int FirstLIndex(int lb, int rb, bool lastChild) => ChildFiles.ValueOf(cld[lastChild ? lb : rb]);

    public int NextLIndex(int k, int rb) => (cld[k] >> 30) == 1 ? ChildFiles.ValueOf(cld[k]) : -1;
}

/// <summary>Layout A, packed codec: next and down are told apart by comparing LCP values (the read rule of the in-RAM packed layout).</summary>
internal sealed class PackedView(uint[] cld, CompactLcp lcp) : IChildView
{
    public int FirstLIndex(int lb, int rb, bool lastChild) => (int)cld[lastChild ? lb : rb];

    public int NextLIndex(int k, int rb)
    {
        if (k >= rb) return -1;
        int lk = lcp[k];
        if (lcp[k + 1] == lk) return k + 1;
        int v = (int)cld[k];
        return lcp[v] == lk ? v : -1;
    }
}

/// <summary>
/// Layout F: children from the next-smaller-value array and the LCP values. The first l-index of <c>[lb..rb]</c> is the first position of the minimum of <c>LCP(lb..rb]</c>, which the chain of
/// next smaller values from <c>lb + 1</c> reaches (it visits the strict prefix minima); the l-index after <c>k</c> is the first position of the minimum of <c>LCP(k..rb]</c> when that minimum equals <c>LCP[k]</c>.
/// </summary>
internal sealed class NsvView(int[] nsv, CompactLcp lcp) : IChildView
{
    public int FirstLIndex(int lb, int rb, bool lastChild) => FirstMinimum(lb + 1, rb);

    public int NextLIndex(int k, int rb)
    {
        if (k >= rb) return -1;
        int x = FirstMinimum(k + 1, rb);
        return lcp[x] == lcp[k] ? x : -1;
    }

    private int FirstMinimum(int from, int rb)
    {
        int x = from;
        while (nsv[x] <= rb) x = nsv[x];
        return x;
    }
}

/// <summary>The walk of <c>EsaIndex.Traverse</c> over any <see cref="IChildView"/> (same events, same order), so that the structures PH5 wrote can be compared with the in-RAM tree and, through the logical hash, with the Ukkonen tree.</summary>
internal static class ChildTraversal
{
    internal static void Traverse(string text, int[] sa, CompactLcp lcp, IChildView child, ISuffixTreeVisitor visitor)
    {
        int n = text.Length;
        if (n == 0)
        {
            visitor.VisitNode(0, 0, 0, 0, 0);
            return;
        }
        int[] occurrence = LeftmostOccurrence.Compute(sa, lcp);
        int firstOfRoot = child.FirstLIndex(0, n, true);
        visitor.VisitNode(0, 0, n + 1, CountChildren(child, firstOfRoot, n), 0);
        var frames = new List<(int A, int NextBoundary, int Rb, int Ell)> { (0, firstOfRoot, n, 0) };
        while (frames.Count > 0)
        {
            int top = frames.Count - 1;
            (int fa, int fnext, int frb, int fell) = frames[top];
            if (fa > frb)
            {
                frames.RemoveAt(top);
                if (top > 0) visitor.ExitBranch();
                continue;
            }
            bool lastChild = fnext < 0;
            int a = fa;
            int b = lastChild ? frb : fnext - 1;
            int next = lastChild ? -1 : child.NextLIndex(fnext, frb);
            frames[top] = (b + 1, next, frb, fell);
            int p = sa[a] + fell;
            visitor.EnterBranch(p < n ? text[p] : -1);
            if (a == b)
            {
                visitor.VisitNode(p, -1, 1, 0, fell);
                visitor.ExitBranch();
                continue;
            }
            int first = child.FirstLIndex(a, b, lastChild);
            int ell = lcp[first];
            int occ = occurrence[first];
            visitor.VisitNode(occ + fell, occ + ell, b - a + 1, CountChildren(child, first, b), fell);
            frames.Add((a, first, b, ell));
        }
    }

    private static int CountChildren(IChildView child, int first, int rb)
    {
        int count = 1;
        for (int k = first; k >= 0; k = child.NextLIndex(k, rb)) count++;
        return count;
    }
}

/// <summary>Records the events of a traversal as text lines.</summary>
internal sealed class EventRecorder : ISuffixTreeVisitor
{
    internal List<string> Events { get; } = [];

    public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Events.Add($"N {startIndex} {endIndex} {leafCount} {childCount} {depth}");

    public void EnterBranch(int key) => Events.Add($"E {key}");

    public void ExitBranch() => Events.Add("X");
}

/// <summary>An <see cref="ISuffixTree"/> that has only a text and a traversal (what <c>CalculateLogicalHash</c> reads), the traversal being any delegate.</summary>
internal sealed class WalkTree(string text, Action<ISuffixTreeVisitor> traverse) : ISuffixTree
{
    public ITextSource Text { get; } = new StringTextSource(text);

    public int NodeCount => throw Unsupported();

    public int LeafCount => throw Unsupported();

    public int MaxDepth => throw Unsupported();

    public bool IsEmpty => text.Length == 0;

    public void Traverse(ISuffixTreeVisitor visitor) => traverse(visitor);

    public string PrintTree() => throw Unsupported();

    public bool Contains(string value) => throw Unsupported();

    public bool Contains(ReadOnlySpan<char> value) => throw Unsupported();

    public IReadOnlyList<int> FindAllOccurrences(string pattern) => throw Unsupported();

    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern) => throw Unsupported();

    public int CountOccurrences(string pattern) => throw Unsupported();

    public int CountOccurrences(ReadOnlySpan<char> pattern) => throw Unsupported();

    public string LongestRepeatedSubstring() => throw Unsupported();

    public IReadOnlyList<string> GetAllSuffixes() => throw Unsupported();

    public IEnumerable<string> EnumerateSuffixes() => throw Unsupported();

    public string LongestCommonSubstring(string other) => throw Unsupported();

    public string LongestCommonSubstring(ReadOnlySpan<char> other) => throw Unsupported();

    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other) => throw Unsupported();

    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other) => throw Unsupported();

    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(string query, int minLength) => throw Unsupported();

    private static NotSupportedException Unsupported() => new("the walk tree has only a text and a traversal");
}

/// <summary>Runs PH5 on a reference text and compares everything it writes (G1 of the prototype E5).</summary>
internal static class ChildChecks
{
    /// <summary>Options with tiny buffers so that every spill, refill and window path is crossed on texts of a few hundred symbols.</summary>
    internal static ChildBuildOptions Tiny(string scratch, int blockEntries = 3, int windowBytes = 32, int slots = 2, ChildCodec codec = ChildCodec.Tagged, bool flush = true) => new()
    {
        ScratchDir = scratch,
        BlockEntries = blockEntries,
        WindowBytes = windowBytes,
        WindowSlots = slots,
        FlushWindows = flush,
        SpillFlushBytes = 64,
        ReadBufferBytes = 64,
        Codec = codec,
        ArenaGuards = ArenaOptions.Guarded(0xCD),
    };

    /// <summary>
    /// Builds layout A (with the codec of <paramref name="o"/>), the PSV and, when <paramref name="both"/>, the NSV of <paramref name="r"/> from its LCP files and fails with a description of the first difference
    /// from the in-RAM structures. Returns the run's statistics.
    /// </summary>
    internal static ChildBuildStats RunAndCompare(ChildReference r, ChildBuildOptions o, string dir, string label, bool both = true)
    {
        (string lcpBase, string lcpOverflow) = r.WriteLcpFiles(dir, "lcp");
        var outputs = new ChildBuildOutputs
        {
            Child = Path.Combine(dir, "child.bin"),
            Psv = both ? Path.Combine(dir, "psv.bin") : null,
            Nsv = both ? Path.Combine(dir, "nsv.bin") : null,
        };
        ChildBuildStats st = ExternalChild.BuildFromFiles(r.N, lcpBase, lcpOverflow, outputs, o);
        CompareChild(r, outputs.Child, o.Codec, label);
        if (both)
        {
            CompareArrays(r.Psv, ChildFiles.ReadSmaller(outputs.Psv!), "PSV", label);
            CompareArrays(r.Nsv, ChildFiles.ReadSmaller(outputs.Nsv!), "NSV", label);
        }
        CompareStats(r, st.Forward, label);
        Assert.That(Directory.Exists(o.ScratchDir) ? Directory.GetFileSystemEntries(o.ScratchDir) : [], Is.Empty, $"{label}: scratch files left behind");
        return st;
    }

    internal static void CompareChild(ChildReference r, string childPath, ChildCodec codec, string label)
    {
        uint[] got = ChildFiles.ReadUInts(childPath);
        uint[] expected = codec == ChildCodec.Tagged ? r.Tagged : r.Packed;
        if (got.Length != expected.Length) Assert.Fail($"{label}: CHILD has {got.Length} entries, expected {expected.Length}");
        for (int k = 0; k < got.Length; k++)
        {
            if (got[k] != expected[k]) Assert.Fail($"{label}: CHILD[{k}] = 0x{got[k]:X8}, expected 0x{expected[k]:X8} ({codec})");
        }
        if (codec == ChildCodec.Tagged)
        {
            (int[] next, int[] up, int[] down) = ChildFiles.TaggedToThree(got);
            CompareArrays(r.Next, next, "next", label);
            CompareArrays(r.Up, up, "up", label);
            CompareArrays(r.Down, down, "down", label);
        }
    }

    internal static void CompareArrays(int[] expected, int[] got, string what, string label)
    {
        if (expected.Length != got.Length) Assert.Fail($"{label}: {what} has {got.Length} entries, expected {expected.Length}");
        for (int k = 0; k < expected.Length; k++)
        {
            if (expected[k] != got[k]) Assert.Fail($"{label}: {what}[{k}] = {got[k]}, expected {expected[k]}");
        }
    }

    internal static void CompareStats(ChildReference r, ChildScanResult s, string label)
    {
        Assert.That(s.NodeCount, Is.EqualTo(r.Statistics.NodeCount), $"{label}: NODE_COUNT");
        Assert.That(s.InternalCount, Is.EqualTo(r.Statistics.InternalCount), $"{label}: internal nodes");
        Assert.That((int)s.LrsDepth, Is.EqualTo(r.Statistics.LrsDepth), $"{label}: LRS_DEPTH");
        Assert.That((int)s.LrsRank, Is.EqualTo(r.Statistics.LrsRank), $"{label}: LRS_RANK");
    }

    /// <summary>The event stream of the in-RAM index, of layout A (both codecs) and of layout F, and the logical hashes of all of them and of the Ukkonen tree.</summary>
    internal static void CompareTraversals(ChildReference r, string childTagged, string childPacked, string nsvPath, string label, bool ukkonen)
    {
        var expected = new EventRecorder();
        r.Index.Traverse(expected);
        CompactLcp lcp = r.Lcp;
        int[] sa = r.Index.SuffixArray;
        var views = new (string Name, IChildView View)[]
        {
            ("A tagged", new TaggedView(ChildFiles.ReadUInts(childTagged))),
            ("A packed", new PackedView(ChildFiles.ReadUInts(childPacked), lcp)),
            ("F", new NsvView(ChildFiles.ReadSmaller(nsvPath), lcp)),
        };
        foreach ((string name, IChildView view) in views)
        {
            var got = new EventRecorder();
            ChildTraversal.Traverse(r.Text, sa, lcp, view, got);
            if (!expected.Events.SequenceEqual(got.Events))
            {
                int i = 0;
                while (i < expected.Events.Count && i < got.Events.Count && expected.Events[i] == got.Events[i]) i++;
                Assert.Fail($"{label}: traversal through {name} differs from the in-RAM index at event {i}: '{(i < got.Events.Count ? got.Events[i] : "<end>")}' instead of '{(i < expected.Events.Count ? expected.Events[i] : "<end>")}'");
            }
        }
        if (!ukkonen) return;
        byte[] hashUkkonen = Persistent.SuffixTreeSerializer.CalculateLogicalHash(global::SuffixTree.SuffixTree.Build(r.Text));
        foreach ((string name, IChildView view) in views)
        {
            byte[] hash = Persistent.SuffixTreeSerializer.CalculateLogicalHash(new WalkTree(r.Text, v => ChildTraversal.Traverse(r.Text, sa, lcp, view, v)));
            Assert.That(Convert.ToHexString(hash), Is.EqualTo(Convert.ToHexString(hashUkkonen)), $"{label}: logical hash through {name} vs Ukkonen");
        }
    }
}
