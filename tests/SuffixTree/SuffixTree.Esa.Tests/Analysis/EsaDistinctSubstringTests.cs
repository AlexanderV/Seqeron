using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Analysis;

/// <summary>
/// <c>CountDistinctSubstrings</c> and <c>CountDistinctSubstringsByLength</c> of the in-memory ESA (work package 1b.3), computed from the LCP array alone. Oracles: the definition
/// (sets of substrings, L0), the Ukkonen tree (L2: its own native walk), and the <c>Traverse</c>-default path of the generic algorithm over the ESA (<see cref="EsaTreeShim"/>).
/// </summary>
[TestFixture]
public class EsaDistinctSubstringTests
{
    private static int[] MaxLengths(int n) => new[] { 1, 2, 3, 5, 8, n / 2, n - 1, n, n + 1, n + 7, int.MaxValue }.Where(m => m >= 1).Distinct().ToArray();

    private static void CheckTotals(EsaIndex esa, UkkonenTree ukkonen, string text, string context)
    {
        long total = esa.CountDistinctSubstrings();
        Assert.That(total, Is.EqualTo(AnalysisReference.DistinctTotal(text)), $"{context}: total against the definition");
        Assert.That(total, Is.EqualTo(ukkonen.CountDistinctSubstrings()), $"{context}: total against Ukkonen");
        Assert.That(total, Is.EqualTo(SuffixTreeAlgorithms.CountDistinctSubstrings(new EsaTreeShim(esa))), $"{context}: total against the Traverse-default path");
    }

    private static void CheckByLength(EsaIndex esa, UkkonenTree ukkonen, string text, int maxLength, string context)
    {
        long[] actual = esa.CountDistinctSubstringsByLength(maxLength);
        Assert.That(actual, Is.EqualTo(AnalysisReference.DistinctByLength(text, maxLength)), $"{context} m {maxLength}: against the definition");
        Assert.That(actual, Is.EqualTo(ukkonen.CountDistinctSubstringsByLength(maxLength)), $"{context} m {maxLength}: against Ukkonen");
        Assert.That(actual, Is.EqualTo(SuffixTreeAlgorithms.CountDistinctSubstringsByLength(new EsaTreeShim(esa), maxLength)), $"{context} m {maxLength}: against the Traverse-default path");
    }

    [Test]
    public void Banana_LiteratureValues()
    {
        // b,a,n | ba,an,na | ban,ana,nan | bana,anan,nana | banan,anana | banana
        EsaIndex esa = EsaIndex.Build("banana");
        Assert.That(esa.CountDistinctSubstringsByLength(6), Is.EqualTo(new long[] { 1, 3, 3, 3, 3, 2, 1 }));
        Assert.That(esa.CountDistinctSubstrings(), Is.EqualTo(15));
    }

    [TestCase("abcd", 10)]      // all n(n+1)/2 substrings are distinct
    [TestCase("aaaa", 4)]       // one substring per length
    [TestCase("a", 1)]
    [TestCase("ATTTGGATT", 35)]
    [TestCase("", 0)]
    public void KnownTotals(string text, long expected) => Assert.That(EsaIndex.Build(text).CountDistinctSubstrings(), Is.EqualTo(expected));

    [Test]
    public void EmptyText_HasNoSubstrings_AndTheEmptyOne()
    {
        EsaIndex esa = EsaIndex.Build("");
        Assert.That(esa.CountDistinctSubstrings(), Is.EqualTo(0));
        foreach (int m in new[] { 1, 5, int.MaxValue }) Assert.That(esa.CountDistinctSubstringsByLength(m), Is.EqualTo(new long[] { 1 }), $"m {m}");
    }

    [Test]
    public void SingleSymbol_AndAllEqualText()
    {
        Assert.That(EsaIndex.Build("a").CountDistinctSubstringsByLength(1), Is.EqualTo(new long[] { 1, 1 }));
        Assert.That(EsaIndex.Build("a").CountDistinctSubstringsByLength(9), Is.EqualTo(new long[] { 1, 1 }));
        foreach (int n in new[] { 2, 3, 10, 257, 5000 })
        {
            EsaIndex esa = EsaIndex.Build(new string('b', n));
            Assert.That(esa.CountDistinctSubstrings(), Is.EqualTo(n), $"b^{n}");
            Assert.That(esa.CountDistinctSubstringsByLength(int.MaxValue), Is.All.EqualTo(1L), $"b^{n}");
            Assert.That(esa.CountDistinctSubstringsByLength(int.MaxValue), Has.Length.EqualTo(n + 1), $"b^{n}");
        }
    }

    [Test]
    public void AllDistinctSymbols_EverySubstringIsDistinct()
    {
        string text = StringFamilies.AllDistinct(300, 5);
        EsaIndex esa = EsaIndex.Build(text);
        Assert.That(esa.CountDistinctSubstrings(), Is.EqualTo(300L * 301 / 2));
        long[] counts = esa.CountDistinctSubstringsByLength(300);
        for (int i = 1; i <= 300; i++) Assert.That(counts[i], Is.EqualTo(300 - i + 1), $"length {i}");
    }

    [Test]
    public void MaxLengthBelowOne_Throws_AsTheOtherTrees()
    {
        EsaIndex esa = EsaIndex.Build("acgt");
        UkkonenTree ukkonen = UkkonenTree.Build("acgt");
        foreach (int m in new[] { 0, -1, int.MinValue })
        {
            ArgumentOutOfRangeException? expected = Assert.Throws<ArgumentOutOfRangeException>(() => ukkonen.CountDistinctSubstringsByLength(m));
            ArgumentOutOfRangeException? actual = Assert.Throws<ArgumentOutOfRangeException>(() => esa.CountDistinctSubstringsByLength(m), $"m {m}");
            Assert.That(actual!.ParamName, Is.EqualTo(expected!.ParamName), $"m {m}");
        }
        // checked before the text is looked at: the empty text throws too
        Assert.Throws<ArgumentOutOfRangeException>(() => EsaIndex.Build("").CountDistinctSubstringsByLength(0));
    }

    [Test]
    public void MaxLengthBeyondTheText_IsClampedToTheTextLength()
    {
        EsaIndex esa = EsaIndex.Build("acgtac");
        Assert.That(esa.CountDistinctSubstringsByLength(int.MaxValue), Is.EqualTo(AnalysisReference.DistinctByLength("acgtac", 6)));
        Assert.That(esa.CountDistinctSubstringsByLength(7), Has.Length.EqualTo(7));
        Assert.That(esa.CountDistinctSubstringsByLength(6), Has.Length.EqualTo(7));
        Assert.That(esa.CountDistinctSubstringsByLength(5), Has.Length.EqualTo(6));
        Assert.That(esa.CountDistinctSubstringsByLength(1), Has.Length.EqualTo(2));
    }

    [Test]
    public void LengthsAddUpToTheTotal_AndEachCallReturnsAFreshArray()
    {
        EsaIndex esa = EsaIndex.Build(StringFamilies.Fibonacci(500));
        long[] counts = esa.CountDistinctSubstringsByLength(500);
        Assert.That(counts.Skip(1).Sum(), Is.EqualTo(esa.CountDistinctSubstrings()));
        counts[1] = -1;
        Assert.That(esa.CountDistinctSubstringsByLength(500)[1], Is.EqualTo(2), "the Fibonacci word has two distinct symbols");
    }

    [Test]
    public void AllSmallStrings_AgreeWithTheDefinitionUkkonenAndTheTraverseDefaultPath()
    {
        foreach (var (name, text) in AnalysisCorpus.Definitional())
        {
            EsaIndex esa = EsaIndex.Build(text);
            UkkonenTree ukkonen = UkkonenTree.Build(text);
            CheckTotals(esa, ukkonen, text, name);
            foreach (int m in MaxLengths(text.Length)) CheckByLength(esa, ukkonen, text, m, name);
        }
    }

    [Test]
    public void LargerTexts_AgreeWithUkkonenAndTheTraverseDefaultPath()
    {
        var rng = new SplitMix64(0x1B3B);
        for (int i = 0; i < 10; i++)
        {
            int n = 2000 + rng.NextInt(18001);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 6 == 0))
            {
                EsaIndex esa = EsaIndex.Build(text);
                UkkonenTree ukkonen = UkkonenTree.Build(text);
                var shim = new EsaTreeShim(esa);
                Assert.That(esa.CountDistinctSubstrings(), Is.EqualTo(ukkonen.CountDistinctSubstrings()), $"{name} n {n}: total");
                Assert.That(esa.CountDistinctSubstrings(), Is.EqualTo(SuffixTreeAlgorithms.CountDistinctSubstrings(shim)), $"{name} n {n}: total (Traverse)");
                foreach (int m in new[] { 1, 7, 100, n / 2, n, n + 1 })
                {
                    long[] actual = esa.CountDistinctSubstringsByLength(m);
                    Assert.That(actual, Is.EqualTo(ukkonen.CountDistinctSubstringsByLength(m)), $"{name} n {n} m {m}: Ukkonen");
                    Assert.That(actual, Is.EqualTo(SuffixTreeAlgorithms.CountDistinctSubstringsByLength(shim, m)), $"{name} n {n} m {m}: Traverse");
                }
            }
        }
    }

    [Test]
    public void TotalIsComputedInLong_WhenItExceedsInt()
    {
        // n (n + 1) / 2 must be evaluated in long: a text of 120,001 symbols (almost all substrings distinct) already has about 7.2e9 of them, more than int.MaxValue
        string text = string.Concat(Enumerable.Range(0, 3).Select(k => StringFamilies.AllDistinct(40000, (ulong)k + 1))) + "ሴ";
        EsaIndex esa = EsaIndex.Build(text);
        long total = esa.CountDistinctSubstrings();
        Assert.That(total, Is.GreaterThan(int.MaxValue));
        Assert.That(total, Is.EqualTo(SuffixTreeAlgorithms.CountDistinctSubstrings(new EsaTreeShim(esa))));
    }

    [Test]
    public void NoAllocationProportionalToTheTextForTheTotal_AndOnlyTheResultForTheLengths()
    {
        string text = StringFamilies.RandomDna(200_000, 11);
        EsaIndex esa = EsaIndex.Build(text);
        esa.CountDistinctSubstrings();                           // JIT
        esa.CountDistinctSubstringsByLength(10);
        long before = GC.GetAllocatedBytesForCurrentThread();
        esa.CountDistinctSubstrings();
        long total = GC.GetAllocatedBytesForCurrentThread() - before;
        before = GC.GetAllocatedBytesForCurrentThread();
        esa.CountDistinctSubstringsByLength(10);
        long small = GC.GetAllocatedBytesForCurrentThread() - before;
        before = GC.GetAllocatedBytesForCurrentThread();
        long[] full = esa.CountDistinctSubstringsByLength(text.Length);
        long large = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(total, Is.LessThan(1024), "CountDistinctSubstrings allocates nothing proportional to n");
        Assert.That(small, Is.LessThan(1024), "ByLength(10) allocates the 11-element result only");
        Assert.That(large, Is.LessThan(8L * (text.Length + 1) + 1024), "ByLength(n) allocates the result array only: 8 bytes per length, no second array");
        Assert.That(full, Has.Length.EqualTo(text.Length + 1));
    }

    [Test]
    public void ResultIsIndependentOfTheChildTableLayout()
    {
        foreach (string text in new[] { "banana", "mississippi", "aabbaabb", "abcabcxyzxyz", new string('a', 50), StringFamilies.Fibonacci(100) })
            foreach (ChildLayoutKind layout in Enum.GetValues<ChildLayoutKind>())
            {
                EsaIndex esa = EsaIndex.Build(text, layout);
                Assert.That(esa.CountDistinctSubstrings(), Is.EqualTo(AnalysisReference.DistinctTotal(text)), $"layout {layout} '{text}'");
                Assert.That(esa.CountDistinctSubstringsByLength(text.Length), Is.EqualTo(AnalysisReference.DistinctByLength(text, text.Length)), $"layout {layout} '{text}'");
            }
    }

    [Test]
    public void ConcurrentReaders_GetTheSameAnswers()
    {
        EsaIndex esa = EsaIndex.Build(StringFamilies.RandomProtein(20000, 3));
        long total = esa.CountDistinctSubstrings();
        long[] byLength = esa.CountDistinctSubstringsByLength(50);
        Parallel.For(0, 4, new ParallelOptions { MaxDegreeOfParallelism = 4 }, _ =>
        {
            for (int k = 0; k < 20; k++)
            {
                Assert.That(esa.CountDistinctSubstrings(), Is.EqualTo(total));
                Assert.That(esa.CountDistinctSubstringsByLength(50), Is.EqualTo(byLength));
            }
        });
    }
}
