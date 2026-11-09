using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Analysis;

/// <summary>
/// <c>GetAllSuffixes</c> and <c>EnumerateSuffixes</c> of the in-memory ESA (work package 1b.3): the suffixes in lexicographic (ordinal) order, a plain walk of the suffix array.
/// Oracles: the definition (sort of all suffixes, L0), the Ukkonen tree (its own native order, which must be the same list), and the leaf order of the <c>Traverse</c> stream of the ESA
/// (the leaves of the event stream are in suffix order; <see cref="TraverseSuffixes"/>).
/// </summary>
[TestFixture]
public class EsaSuffixEnumerationTests
{
    /// <summary>The suffixes in the order of the leaves of the <c>Traverse</c> stream: a leaf with start <c>s</c> and parent depth <c>d</c> is the suffix at <c>s - d</c>; the empty-suffix leaf is skipped.</summary>
    internal static List<string> TraverseSuffixes(EsaIndex esa)
    {
        var collector = new LeafCollector(esa.Text);
        esa.Traverse(collector);
        return collector.Suffixes;
    }

    private sealed class LeafCollector(string text) : ISuffixTreeVisitor
    {
        private bool _rootSeen;

        internal List<string> Suffixes { get; } = [];

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (!_rootSeen)
            {
                _rootSeen = true;
                return;
            }
            if (childCount != 0) return;
            int p = startIndex - depth;
            if (p < text.Length) Suffixes.Add(text.Substring(p));
        }

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }

    private static void Check(string text, string context)
    {
        EsaIndex esa = EsaIndex.Build(text);
        string[] expected = AnalysisReference.SortedSuffixes(text);
        Assert.That(esa.GetAllSuffixes(), Is.EqualTo(expected), $"{context}: GetAllSuffixes against the definition");
        Assert.That(esa.EnumerateSuffixes(), Is.EqualTo(expected), $"{context}: EnumerateSuffixes against the definition");
        UkkonenTree ukkonen = UkkonenTree.Build(text);
        Assert.That(esa.GetAllSuffixes(), Is.EqualTo(ukkonen.GetAllSuffixes()), $"{context}: GetAllSuffixes against Ukkonen");
        Assert.That(esa.EnumerateSuffixes(), Is.EqualTo(ukkonen.EnumerateSuffixes()), $"{context}: EnumerateSuffixes against Ukkonen");
        Assert.That(esa.GetAllSuffixes(), Is.EqualTo(TraverseSuffixes(esa)), $"{context}: against the leaf order of Traverse");
    }

    [Test]
    public void Banana_AllSuffixesInOrder()
    {
        Assert.That(EsaIndex.Build("banana").GetAllSuffixes(), Is.EqualTo(new[] { "a", "ana", "anana", "banana", "na", "nana" }));
        Assert.That(EsaIndex.Build("banana").EnumerateSuffixes(), Is.EqualTo(new[] { "a", "ana", "anana", "banana", "na", "nana" }));
    }

    [Test]
    public void EmptyAndSingle()
    {
        EsaIndex empty = EsaIndex.Build("");
        Assert.That(empty.GetAllSuffixes(), Is.Empty);
        Assert.That(empty.EnumerateSuffixes(), Is.Empty);
        EsaIndex one = EsaIndex.Build("a");
        Assert.That(one.GetAllSuffixes(), Is.EqualTo(new[] { "a" }));
        Assert.That(one.EnumerateSuffixes(), Is.EqualTo(new[] { "a" }));
    }

    [Test]
    public void AllEqualText_ShorterPrefixesFirst()
    {
        EsaIndex esa = EsaIndex.Build("aaaa");
        Assert.That(esa.GetAllSuffixes(), Is.EqualTo(new[] { "a", "aa", "aaa", "aaaa" }));
    }

    [Test]
    public void EverySuffixOnce_CountIsTheTextLength_AndUnicodeOrderIsOrdinal()
    {
        foreach (string text in new[] { "abc", "banana", "mississippi", "￿\0a￿a", "é中😀é中", "𐀀\uD800" })
        {
            IReadOnlyList<string> all = EsaIndex.Build(text).GetAllSuffixes();
            Assert.That(all, Has.Count.EqualTo(text.Length), $"'{text}'");
            Assert.That(all.Distinct().Count(), Is.EqualTo(text.Length), $"'{text}'");
            for (int i = 1; i < all.Count; i++)
                Assert.That(string.CompareOrdinal(all[i - 1], all[i]), Is.LessThan(0), $"'{text}' strictly ascending at {i}");
        }
    }

    [Test]
    public void AllSmallStrings_AgreeWithTheDefinitionUkkonenAndTheTraverseLeafOrder()
    {
        foreach (var (name, text) in AnalysisCorpus.Definitional()) Check(text, name);
    }

    [Test]
    public void LargerTexts_AreStrictlyAscending_EverySuffixOnce()
    {
        var rng = new SplitMix64(0x1B3C);
        for (int i = 0; i < 6; i++)
        {
            int n = 1000 + rng.NextInt(3001);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 5 == 0))
            {
                EsaIndex esa = EsaIndex.Build(text);
                string? previous = null;
                int count = 0;
                long lengthSum = 0;
                foreach (string suffix in esa.EnumerateSuffixes())
                {
                    if (previous is not null && string.CompareOrdinal(previous, suffix) >= 0) Assert.Fail($"{name} n {n}: suffix {count} is not above its predecessor");
                    if (!text.EndsWith(suffix, StringComparison.Ordinal)) Assert.Fail($"{name} n {n}: suffix {count} is not a suffix of the text");
                    previous = suffix;
                    count++;
                    lengthSum += suffix.Length;
                }
                Assert.That(count, Is.EqualTo(n), $"{name} n {n}");
                Assert.That(lengthSum, Is.EqualTo((long)n * (n + 1) / 2), $"{name} n {n}: the lengths are 1..n, each once (strict order and suffix-ness make the suffixes distinct)");
            }
        }
    }

    [Test]
    public void LargeText_FirstAndLastSuffixesEqualUkkonensLazyEnumeration()
    {
        string text = StringFamilies.ThueMorse(8000);
        EsaIndex esa = EsaIndex.Build(text);
        UkkonenTree ukkonen = UkkonenTree.Build(text);
        Assert.That(esa.EnumerateSuffixes().Take(40), Is.EqualTo(ukkonen.EnumerateSuffixes().Take(40)));
        Assert.That(esa.EnumerateSuffixes().Skip(7960), Is.EqualTo(ukkonen.EnumerateSuffixes().Skip(7960)));
    }

    [Test]
    public void EnumerateSuffixes_IsLazy_NothingIsMaterialisedBeyondWhatIsTaken()
    {
        string text = new('a', 20000);        // the full list would be 2e8 characters = 400 MB
        EsaIndex esa = EsaIndex.Build(text);
        _ = esa.EnumerateSuffixes().Take(2).ToList();   // JIT
        long before = GC.GetAllocatedBytesForCurrentThread();
        IEnumerable<string> lazy = esa.EnumerateSuffixes();
        List<string> first = lazy.Take(3).ToList();
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(first, Is.EqualTo(new[] { "a", "aa", "aaa" }));
        Assert.That(allocated, Is.LessThan(64 * 1024), "the call itself allocates the enumerator and three short strings, not the list");
    }

    [Test]
    public void EnumerateSuffixes_CanBeEnumeratedTwice_AndGetAllSuffixesReturnsAFreshList()
    {
        EsaIndex esa = EsaIndex.Build("mississippi");
        IEnumerable<string> e = esa.EnumerateSuffixes();
        List<string> firstPass = e.ToList();
        List<string> secondPass = e.ToList();
        Assert.That(secondPass, Is.EqualTo(firstPass));
        IReadOnlyList<string> a = esa.GetAllSuffixes();
        IReadOnlyList<string> b = esa.GetAllSuffixes();
        Assert.That(a, Is.Not.SameAs(b));
        Assert.That(a, Is.EqualTo(b));
    }

    [Test]
    public void ResultIsIndependentOfTheChildTableLayout()
    {
        foreach (string text in new[] { "banana", "mississippi", "aabbaabb", new string('a', 30), StringFamilies.Fibonacci(60) })
            foreach (ChildLayoutKind layout in Enum.GetValues<ChildLayoutKind>())
                Assert.That(EsaIndex.Build(text, layout).GetAllSuffixes(), Is.EqualTo(AnalysisReference.SortedSuffixes(text)), $"layout {layout} '{text}'");
    }

    [Test]
    public void ConcurrentReaders_GetTheSameAnswers()
    {
        string text = StringFamilies.RandomDna(1500, 9);
        EsaIndex esa = EsaIndex.Build(text);
        string[] expected = AnalysisReference.SortedSuffixes(text);
        Parallel.For(0, 4, new ParallelOptions { MaxDegreeOfParallelism = 4 }, _ =>
        {
            for (int k = 0; k < 5; k++)
            {
                Assert.That(esa.GetAllSuffixes(), Is.EqualTo(expected));
                Assert.That(esa.EnumerateSuffixes(), Is.EqualTo(expected));
            }
        });
    }
}
