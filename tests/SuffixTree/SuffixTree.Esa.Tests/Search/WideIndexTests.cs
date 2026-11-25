using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Search;

/// <summary>
/// Work package S2.4: every operation of the index at 64 bits over sections that the 32-bit engine of S1 cannot read (5-byte entries, a 64-bit LCP array, the wide forms of the 4-byte layouts; <see cref="WideImage"/>) against the same
/// operation of the index that was built in memory: search (<c>Contains</c>, the counts of both widths, the lists, the ranks), the pages of <c>Locate</c>, the analysis (the longest repeats, the distinct substrings, the suffixes), the
/// events of <c>Traverse</c> (the 64-bit visitor and the 32-bit one) and the contract (<see cref="EsaSuffixIndex"/>: loci, counts, pages, children, edge symbols, the arrays). Texts of at most a few thousand symbols, every layout, both modes.
/// </summary>
[TestFixture]
public class WideIndexTests
{
    private static readonly EsaAccessMode[] Modes = [EsaAccessMode.Checked, EsaAccessMode.Trusted];

    private static string RandomText(SplitMix64 rng, int n, int sigma)
    {
        var chars = new char[n];
        for (int i = 0; i < n; i++) chars[i] = (char)('A' + rng.NextInt(sigma));
        return new string(chars);
    }

    public static IEnumerable<string> Texts()
    {
        var rng = new SplitMix64(0x5EED2400);
        yield return "mississippi";
        yield return "abracadabra";
        yield return "aaaaaaaaaaaaaaaaaaaa";
        yield return "abababababababababab";
        yield return "banana";
        yield return "a";
        yield return RandomText(rng, 400, 4);
        yield return RandomText(rng, 300, 2);
        yield return RandomText(rng, 250, 20);
        string block = RandomText(rng, 60, 4);
        yield return RandomText(rng, 80, 4) + block + RandomText(rng, 80, 4) + block + new string('A', 300);
    }

    /// <summary>Every layout and mode: the memory index, the index over the sections, the contract views of both.</summary>
    private static void ForEachLayout(string text, Action<string, EsaIndex, EsaIndex> check)
    {
        EsaIndex memory = EsaIndex.Build(text);
        LcpChildCase c = LcpChildCase.Of(text);
        char[] alphabet = [.. text.Distinct().Order()];
        foreach ((int width, LcpKind kind, bool tagged) in WideImage.Layouts)
        {
            using var image = new WideImage(c, width, kind, tagged);
            foreach (EsaAccessMode mode in Modes)
            {
                EsaIndex wide = image.Index(mode == EsaAccessMode.Checked, alphabet);
                check($"width {width} {kind}{(tagged ? " tagged" : string.Empty)} {mode} '{Show(text)}'", memory, wide);
            }
        }
    }

    private static string Show(string s) => s.Length <= 24 ? s : s[..24] + "...";

    private static IEnumerable<string> Patterns(string text, SplitMix64 rng)
    {
        for (int i = 0; i < Math.Min(text.Length, 40); i++)
        {
            for (int len = 1; len <= Math.Min(8, text.Length - i); len++) yield return text.Substring(i, len);
        }
        for (int i = 0; i < 60; i++)
        {
            int len = 1 + rng.NextInt(Math.Min(60, text.Length));
            yield return text.Substring(rng.NextInt(text.Length - len + 1), len);
        }
        yield return text;
        yield return text + text[0];
        yield return "ZZ";
        yield return new string(text[0], 12);
        yield return ((char)(text.Max() + 1)).ToString();   // the letter just above the alphabet of the text: no letter, no match
        yield return text[0] + ((char)(text.Max() + 1)).ToString();
    }

    // ------------------------------------------------------------------------------------------ search

    [Test]
    public void TheSearch_ContainsCountsAndRanks_EqualTheMemoryIndex([ValueSource(nameof(Texts))] string text)
    {
        var rng = new SplitMix64(0xFACE5);
        string[] patterns = [.. Patterns(text, rng)];
        ForEachLayout(text, (label, memory, wide) =>
        {
            Assert.That(wide.TextLength, Is.EqualTo(text.Length), label);
            Assert.That(wide.Length, Is.EqualTo(text.Length), label);
            foreach (string p in patterns)
            {
                Assert.That(wide.Contains(p), Is.EqualTo(memory.Contains(p)), $"{label} Contains '{Show(p)}'");
                Assert.That(wide.LongCount(p), Is.EqualTo((long)memory.CountOccurrences(p)), $"{label} LongCount '{Show(p)}'");
                Assert.That(wide.CountOccurrences(p), Is.EqualTo(memory.CountOccurrences(p)), $"{label} CountOccurrences '{Show(p)}'");
                bool expected = memory.TryMatch(p.AsSpan(), out int elo, out int ehi);
                bool found = wide.TryMatch(p.AsSpan(), out long lo, out long hi);
                bool found32 = wide.TryMatch(p.AsSpan(), out int lo32, out int hi32);
                Assert.That((found, lo, hi), Is.EqualTo((expected, (long)elo, (long)ehi)), $"{label} TryMatch '{Show(p)}'");
                Assert.That((found32, lo32, hi32), Is.EqualTo((expected, elo, ehi)), $"{label} TryMatch (32-bit) '{Show(p)}'");
                Assert.That(wide.FindAllOccurrences(p), Is.EqualTo(memory.FindAllOccurrences(p)), $"{label} FindAllOccurrences '{Show(p)}'");
            }
            Assert.That(wide.LongCount(string.Empty), Is.EqualTo((long)text.Length), label);
            Assert.That(wide.CountOccurrences(string.Empty), Is.EqualTo(text.Length), label);
            Assert.That(wide.FindAllOccurrences(string.Empty), Is.EqualTo(memory.FindAllOccurrences(string.Empty)), label);
            Assert.That(wide.Contains(string.Empty), Is.True, label);
            Assert.That(() => wide.TryMatch(string.Empty.AsSpan(), out long _, out long _), Throws.TypeOf<ArgumentException>(), label);
        });
    }

    // ------------------------------------------------------------------------------------------ the pages of Locate

    [Test]
    public void LocatePages_OfEveryMatch_EqualTheMemoryIndex_AtBothPositionWidths([ValueSource(nameof(Texts))] string text)
    {
        var rng = new SplitMix64(0xBA6E5);
        string[] patterns = [.. Patterns(text, rng).Take(60)];
        ForEachLayout(text, (label, memory, wide) =>
        {
            foreach (string p in patterns)
            {
                if (!memory.TryMatch(p.AsSpan(), out int lo, out int hi)) continue;
                long count = (long)hi - lo + 1;
                foreach (int capacity in new[] { 1, 3, 7, 64 })
                {
                    foreach (long skip in new long[] { 0, 1, count / 2, count - 1, count, count + 5 })
                    {
                        var expectedLong = new long[capacity];
                        var actualLong = new long[capacity];
                        int e = memory.Locate(lo, hi, skip, expectedLong);
                        int a = wide.Locate(lo, hi, skip, actualLong);
                        Assert.That((a, actualLong), Is.EqualTo((e, expectedLong)), $"{label} pattern '{Show(p)}' skip {skip} capacity {capacity} (long)");
                        var expectedInt = new int[capacity];
                        var actualInt = new int[capacity];
                        Assert.That(wide.Locate(lo, hi, skip, actualInt), Is.EqualTo(memory.Locate(lo, hi, skip, expectedInt)), label);
                        Assert.That(actualInt, Is.EqualTo(expectedInt), $"{label} pattern '{Show(p)}' skip {skip} capacity {capacity} (int)");
                    }
                }
            }
            Assert.That(() => wide.Locate(1, 1, -1, new long[2]), Throws.TypeOf<ArgumentOutOfRangeException>(), label);
            Assert.That(() => wide.Locate(3, 2, 0, new long[2]), Throws.TypeOf<ArgumentOutOfRangeException>(), label);
            Assert.That(() => wide.Locate(1, text.Length + 1, 0, new long[2]), Throws.TypeOf<ArgumentOutOfRangeException>(), label);
            Assert.That(() => wide.Locate(-1, 1, 0, new long[2]), Throws.TypeOf<ArgumentOutOfRangeException>(), label);
        });
    }

    [TestCase(1L << 33, 5L, 1L << 32, 0L, 100, 5L, 100)]
    [TestCase(1L << 33, 5L, 1L << 32, 3L, 100, 8L, 100)]
    [TestCase(1L << 33, (1L << 31) - 5, (1L << 32) + 7, 0L, 16, (1L << 31) - 5, 16)]
    [TestCase(1L << 33, (1L << 31) - 5, (1L << 32) + 7, 4L, 16, (1L << 31) - 1, 16)]
    [TestCase(1L << 33, (1L << 31) - 5, (1L << 31) + 4, 0L, 100, (1L << 31) - 5, 10)]
    [TestCase(1L << 33, (1L << 31) - 5, (1L << 31) + 4, 9L, 100, (1L << 31) + 4, 1)]
    [TestCase(1L << 33, (1L << 31) - 5, (1L << 31) + 4, 10L, 100, (1L << 31) - 5, 0)]
    [TestCase(1L << 33, 1L << 32, 1L << 33, (1L << 32) - 1, 3, (1L << 33) - 1, 2)]
    [TestCase(1L << 33, 1L << 32, 1L << 33, long.MaxValue, 3, 1L << 32, 0)]
    [TestCase(1L << 33, 7L, 7L, 0L, 0, 7L, 0)]
    [TestCase(1L << 33, 7L, 7L, 0L, 5, 7L, 1)]
    [TestCase(10L, 1L, 10L, 0L, int.MaxValue, 1L, 10)]
    public void PageOf_IsTheArithmeticOfLocate_AtRanksAbove2To31(long length, long lo, long hi, long skip, int capacity, long expectedStart, int expectedCount)
    {
        int count = EsaIndex.PageOf(length, lo, hi, skip, capacity, out long start);
        Assert.That(count, Is.EqualTo(expectedCount));
        if (expectedCount > 0) Assert.That(start, Is.EqualTo(expectedStart));
    }

    [Test]
    public void PageOf_RefusesANegativeSkipAndAnIntervalThatIsNoRankInterval()
    {
        Assert.That(() => EsaIndex.PageOf(100, 1, 5, -1, 2, out _), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => EsaIndex.PageOf(100, 6, 5, 0, 2, out _), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => EsaIndex.PageOf(100, 1, 101, 0, 2, out _), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => EsaIndex.PageOf(100, -1, 5, 0, 2, out _), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(EsaIndex.PageOf(100, 1, 100, 0, 2, out long start), Is.EqualTo(2));
        Assert.That(start, Is.EqualTo(1));
    }

    // ------------------------------------------------------------------------------------------ analysis

    [Test]
    public void TheAnalysis_Statistics_Repeats_DistinctSubstringsAndSuffixes_EqualTheMemoryIndex([ValueSource(nameof(Texts))] string text)
    {
        ForEachLayout(text, (label, memory, wide) =>
        {
            Assert.That(wide.LongestRepeatedSubstring(), Is.EqualTo(memory.LongestRepeatedSubstring()), label);
            Assert.That(wide.LongestRepeatedSubstringMemory().ToString(), Is.EqualTo(memory.LongestRepeatedSubstringMemory().ToString()), label);
            Assert.That(wide.LrsLength, Is.EqualTo((long)memory.LongestRepeatedSubstring().Length), label);
            Assert.That(wide.LrsPosition, Is.EqualTo(memory.Statistics.LrsDepth == 0 ? -1 : (long)memory.SuffixArray[(int)memory.Statistics.LrsRank]), label);
            Assert.That(wide.CountDistinctSubstrings(), Is.EqualTo(memory.CountDistinctSubstrings()), label);
            foreach (int max in new[] { 1, 2, 5, text.Length, text.Length + 7 })
                Assert.That(wide.CountDistinctSubstringsByLength(max), Is.EqualTo(memory.CountDistinctSubstringsByLength(max)), $"{label} by length {max}");
            Assert.That(wide.GetAllSuffixes(), Is.EqualTo(memory.GetAllSuffixes()), label);
            Assert.That(wide.EnumerateSuffixes().ToArray(), Is.EqualTo(memory.EnumerateSuffixes().ToArray()), label);
            var expectedTies = memory.FindAllLongestRepeatedSubstrings();
            var actualTies = wide.FindAllLongestRepeatedSubstrings();
            Assert.That(actualTies.Select(t => t.Substring), Is.EqualTo(expectedTies.Select(t => t.Substring)), label);
            Assert.That(actualTies.Select(t => t.Positions), Is.EqualTo(expectedTies.Select(t => t.Positions)), label);
            Assert.That(wide.Statistics, Is.EqualTo(memory.Statistics), label);
        });
    }

    // ------------------------------------------------------------------------------------------ Traverse

    private sealed class Recorder : ISuffixTreeVisitor, IWideTreeVisitor
    {
        internal List<string> Events { get; } = [];

        void ISuffixTreeVisitor.VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Events.Add($"V {startIndex} {endIndex} {leafCount} {childCount} {depth}");

        void IWideTreeVisitor.VisitNode(long startIndex, long endIndex, long leafCount, long childCount, long depth) => Events.Add($"V {startIndex} {endIndex} {leafCount} {childCount} {depth}");

        public void EnterBranch(int key) => Events.Add($"E {key}");

        public void ExitBranch() => Events.Add("X");
    }

    [Test]
    public void TheEventsOfTraverse_AtBothWidthsOfTheVisitor_EqualTheMemoryIndex_AndCountTheNodesOfTheHeader([ValueSource(nameof(Texts))] string text)
    {
        ForEachLayout(text, (label, memory, wide) =>
        {
            var expected = new Recorder();
            memory.Traverse((ISuffixTreeVisitor)expected);
            var narrow = new Recorder();
            wide.Traverse((ISuffixTreeVisitor)narrow);
            var wideEvents = new Recorder();
            wide.Traverse<IWideTreeVisitor>(wideEvents);
            Assert.That(narrow.Events, Is.EqualTo(expected.Events), $"{label}: the 32-bit visitor");
            Assert.That(wideEvents.Events, Is.EqualTo(expected.Events), $"{label}: the 64-bit visitor");
            Assert.That(wideEvents.Events.Count(e => e[0] == 'V'), Is.EqualTo(wide.Statistics.NodeCount), $"{label}: one event per node");
            Assert.That(wideEvents.Events[0], Is.EqualTo($"V 0 0 {text.Length + 1} {expected.Events[0].Split(' ')[4]} 0"), $"{label}: the root carries n + 1 leaves");
        });
    }

    [Test]
    public void TheTraverseOfAnEmptyText_IsOneEvent_AtBothWidths()
    {
        EsaIndex memory = EsaIndex.Build(string.Empty);
        var events = new Recorder();
        memory.Traverse<IWideTreeVisitor>(events);
        Assert.That(events.Events, Is.EqualTo(new[] { "V 0 0 0 0 0" }));
    }

    [Test]
    public void ANodeThatTheThirtyTwoBitVisitorCannotBeTold_IsNotSupported_NamingTheMember()
    {
        var narrow = new NarrowingVisitor(new Recorder());
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => narrow.VisitNode(0, 0, int.MaxValue + 1L, 2, 0))!;
        Assert.That(e.Message, Does.Contain("EsaIndex.Traverse").And.Contain("IWideTreeVisitor"));
        Assert.That(() => narrow.VisitNode(0, -1, 1, 0, 3), Throws.Nothing, "the end of a leaf is -1: boundless");
        Assert.That(() => narrow.VisitNode(0, int.MaxValue + 1L, 1, 0, 3), Throws.TypeOf<NotSupportedException>());
        Assert.That(() => narrow.VisitNode(0, 0, 1, 0, int.MaxValue + 2L), Throws.TypeOf<NotSupportedException>());
    }

    // ------------------------------------------------------------------------------------------ the contract

    [Test]
    public void TheContractIndex_Loci_Counts_Pages_Children_AndTheArrays_EqualTheMemoryIndex([ValueSource(nameof(Texts))] string text)
    {
        var rng = new SplitMix64(0xC0FFEE);
        string[] patterns = [.. Patterns(text, rng).Take(80)];
        ForEachLayout(text, (label, memory, wide) =>
        {
            var m = new EsaSuffixIndex(memory);
            var w = new EsaSuffixIndex(wide);
            Assert.That(w.TextLength, Is.EqualTo(m.TextLength), label);
            Assert.That(w.Root, Is.EqualTo(m.Root), label);
            foreach (string p in patterns)
            {
                bool expected = m.TryExtendRight(m.Root, p.AsSpan(), out IndexLocus el);
                bool actual = w.TryExtendRight(w.Root, p.AsSpan(), out IndexLocus al);
                Assert.That((actual, al), Is.EqualTo((expected, el)), $"{label} ExtendRight '{Show(p)}'");
                if (!expected) continue;
                Assert.That(w.Count(al), Is.EqualTo(m.Count(el)), $"{label} Count '{Show(p)}'");
                Assert.That(w.IsExplicitNode(al), Is.EqualTo(m.IsExplicitNode(el)), $"{label} IsExplicitNode '{Show(p)}'");
                var ep = new long[5];
                var ap = new long[5];
                foreach (long skip in new long[] { 0, 1, 3 })
                {
                    Assert.That(w.Locate(al, skip, ap), Is.EqualTo(m.Locate(el, skip, ep)), $"{label} Locate '{Show(p)}'");
                    Assert.That(ap, Is.EqualTo(ep), $"{label} Locate '{Show(p)}' skip {skip}");
                }
                if (m.IsExplicitNode(el))
                {
                    var ec = new IndexLocus[text.Length + 2];
                    var ac = new IndexLocus[text.Length + 2];
                    int n = m.GetChildren(el, ec);
                    Assert.That(w.GetChildren(al, ac), Is.EqualTo(n), $"{label} GetChildren '{Show(p)}'");
                    Assert.That(ac.Take(n), Is.EqualTo(ec.Take(n)), $"{label} children of '{Show(p)}'");
                    foreach (int symbol in new[] { 'A', 'B', 'a', 'b', 'i', 's' })
                    {
                        bool e2 = m.TryChild(el, symbol, out IndexLocus ech);
                        bool a2 = w.TryChild(al, symbol, out IndexLocus ach);
                        Assert.That((a2, ach), Is.EqualTo((e2, ech)), $"{label} TryChild '{Show(p)}' {symbol}");
                    }
                    foreach (long offset in new long[] { 0, 1, 5, text.Length - 1, text.Length })
                        Assert.That(w.EdgeSymbol(al, offset), Is.EqualTo(m.EdgeSymbol(el, offset)), $"{label} EdgeSymbol '{Show(p)}' {offset}");
                }
            }
            foreach (long r in new long[] { 0, 1, 2, text.Length / 2, text.Length - 1, text.Length }.Where(x => x <= text.Length).Distinct())
            {
                Assert.That(w.SuffixAt(r), Is.EqualTo(m.SuffixAt(r)), $"{label} SA[{r}]");
                Assert.That(w.LcpAt(r), Is.EqualTo(m.LcpAt(r)), $"{label} LCP[{r}]");
            }
            Assert.That(() => w.SuffixAt(text.Length + 1), Throws.TypeOf<ArgumentOutOfRangeException>(), label);
            Assert.That(() => w.LcpAt(-1), Throws.TypeOf<ArgumentOutOfRangeException>(), label);
        });
    }

    [Test]
    public void TheContractIndex_NavigationOperations_OverEveryLayout_EqualTheMemoryIndex()
    {
        // work package S2.5: the navigation pack reads the views at any width (until S2.4 the 5-byte layouts were refused by name here: the 32-bit pack could not read them); the full comparison is Navigation.WideNavigationTests
        string text = "mississippimississippi";
        ForEachLayout(text, (label, memory, wide) =>
        {
            var m = new EsaSuffixIndex(memory);
            var w = new EsaSuffixIndex(wide);
            Navigation.WideNavigationTests.CompareAll(label, m, w, text);
            foreach (string p in new[] { "i", "ss", "issi", "ppi", "mississippi", "s", "p" })
            {
                m.TryExtendRight(m.Root, p.AsSpan(), out IndexLocus el);
                w.TryExtendRight(w.Root, p.AsSpan(), out IndexLocus al);
                Assert.That(al, Is.EqualTo(el), $"{label} locus of '{p}'");
            }
        });
    }

    // ------------------------------------------------------------------------------------------ the 32-bit surface refuses by name

    [Test]
    public void TheCountOfAnIntervalAbove2To31_OnTheThirtyTwoBitSurface_IsNotSupported_NamingTheLongMember()
    {
        // no text of 2^31 symbols is built here: the bridge is tested on the numbers
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => IntBridge.Length((1L << 31) + 3, "EsaIndex.CountOccurrences (the 64-bit member is LongCount)"))!;
        Assert.That(e.Message, Does.Contain("LongCount"));
    }
}
