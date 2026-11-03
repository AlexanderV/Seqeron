using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Sais;

/// <summary>1a.3: the safe SA-IS reference (int32) against the definitional oracle and the linear verifier.</summary>
[TestFixture]
public class SaisReferenceTests
{
    private static void AssertEqualsNaive(string text, string context)
    {
        int[] expected = NaiveSuffixArray.Build(text);
        int[] actual = SaisReference.BuildSuffixArray(text);
        if (!expected.AsSpan().SequenceEqual(actual))
            Assert.Fail($"{context}: SA differs from the naive oracle (n = {text.Length}): text '{Preview(text)}'");
    }

    private static string Preview(string s) => s.Length <= 60 ? s : string.Concat(s.AsSpan(0, 60), "...");

    [Test]
    public void KnownExamples()
    {
        Assert.That(SaisReference.BuildSuffixArray(""), Is.EqualTo(new[] { 0 }));
        Assert.That(SaisReference.BuildSuffixArray("x"), Is.EqualTo(new[] { 1, 0 }));
        Assert.That(SaisReference.BuildSuffixArray("ab"), Is.EqualTo(new[] { 2, 0, 1 }));
        Assert.That(SaisReference.BuildSuffixArray("ba"), Is.EqualTo(new[] { 2, 1, 0 }));
        Assert.That(SaisReference.BuildSuffixArray("banana"), Is.EqualTo(new[] { 6, 5, 3, 1, 0, 4, 2 }));
        Assert.That(SaisReference.BuildSuffixArray("mississippi"), Is.EqualTo(NaiveSuffixArray.Build("mississippi")));
        Assert.That(SaisReference.BuildSuffixArray("aaaa"), Is.EqualTo(new[] { 4, 3, 2, 1, 0 }));
    }

    [Test]
    public void Verifier_AcceptsTheNaiveSa_AndRejectsPerturbedOnes()
    {
        foreach (string s in ExhaustiveStrings.All(2, 8))
        {
            int[] sa = NaiveSuffixArray.Build(s);
            Assert.That(SuffixArrayVerifier.Check(s, sa), Is.Null, s);
            for (int a = 1; a < sa.Length; a++)
                for (int b = a + 1; b < sa.Length; b++)
                {
                    int[] bad = (int[])sa.Clone();
                    (bad[a], bad[b]) = (bad[b], bad[a]);
                    Assert.That(SuffixArrayVerifier.Check(s, bad), Is.Not.Null, $"{s} swap {a},{b}");
                }
        }
        Assert.That(SuffixArrayVerifier.Check("ab", new[] { 2, 0 }), Is.Not.Null);
        Assert.That(SuffixArrayVerifier.Check("ab", new[] { 2, 0, 0 }), Is.Not.Null);
        Assert.That(SuffixArrayVerifier.Check("ab", new[] { 0, 0, 1 }), Is.Not.Null);
    }

    [Test]
    public void Exhaustive_AllStrings_Sigma2UpTo14_Sigma3UpTo9_Sigma4UpTo7()
    {
        var trace = new SaisTrace();
        long count = 0;
        int maxDepthSeen = 0;
        foreach (var (sigma, maxLen) in new[] { (2, 14), (3, 9), (4, 7) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
            {
                int[] expected = NaiveSuffixArray.Build(s);
                trace.Levels.Clear();
                int[] actual = SaisReference.BuildSuffixArray(s, trace);
                if (!expected.AsSpan().SequenceEqual(actual))
                    Assert.Fail($"sigma {sigma}: SA differs from the naive oracle on '{s}'");
                count++;
                maxDepthSeen = Math.Max(maxDepthSeen, trace.Depth);
            }
        Assert.That(count, Is.EqualTo(32767 + 29524 + 21845));
        TestContext.Out.WriteLine($"exhaustive default run: {count} strings, maximum SA-IS level count {maxDepthSeen}");
        Assert.That(maxDepthSeen, Is.GreaterThanOrEqualTo(2));
    }

    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    [TestCase(4)]
    [TestCase(20)]
    [TestCase(256)]
    [TestCase(1000)]
    public void RandomTexts_UpTo2000_MatchNaive(int sigma)
    {
        var rng = new SplitMix64(1000UL + (ulong)sigma);
        for (int i = 0; i < 300; i++)
        {
            int n = rng.NextInt(2001);
            ulong seed = rng.Next();
            AssertEqualsNaive(StringFamilies.Random(n, sigma, seed), $"random sigma {sigma} seed {seed} n {n}");
        }
    }

    [Test]
    public void AdversarialFamilies_UpTo2000_MatchNaive_AndReachDeepRecursion()
    {
        var trace = new SaisTrace();
        int maxDepth = 0;
        foreach (int n in new[] { 0, 1, 2, 3, 4, 5, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 127, 128, 129, 255, 256, 257, 500, 1000, 1023, 1024, 1025, 2000 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 3))
            {
                trace.Levels.Clear();
                int[] actual = SaisReference.BuildSuffixArray(text, trace);
                if (!NaiveSuffixArray.Build(text).AsSpan().SequenceEqual(actual))
                    Assert.Fail($"{name} n {n}: SA differs from the naive oracle");
                maxDepth = Math.Max(maxDepth, trace.Depth);
            }
        // Fibonacci words shrink by about 2.6 per level: n = 2000 needs 8 levels, so the third and deeper recursion levels are checked
        // against the naive oracle here (binary strings up to length 14 reach only 2 levels, up to 22 only 3: review MIN-10 a).
        Assert.That(maxDepth, Is.GreaterThanOrEqualTo(6));
    }

    [Test]
    public void LmsMaximisingTexts_MatchNaive()
    {
        foreach (int n in new[] { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 31, 32, 33, 1999, 2000 })
            AssertEqualsNaive(StringFamilies.AbRepeat(n), $"abab n {n}");
    }

    [Test]
    public void LmsSubstringEquality_MatchesTheDefinition_OnAllSmallStrings()
    {
        long pairs = 0, equalPairs = 0;
        foreach (var (sigma, maxLen) in new[] { (2, 14), (3, 9), (4, 7) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
            {
                var p = TextPreparation.Prepare(s);
                var reader = p.ByteReader();
                bool[] isS = SaisReference.ClassifyTypes(reader);
                int n = s.Length;
                var lms = new List<int>();
                for (int i = 1; i < n; i++)
                    if (SaisReference.IsLms(isS, i)) lms.Add(i);
                lms.Add(n); // the sentinel is the last LMS position
                for (int x = 0; x + 1 < lms.Count; x++)
                    for (int y = 0; y + 1 < lms.Count; y++)
                    {
                        if (x == y) continue;
                        bool expected = DefinitionalLmsEqual(s, isS, lms, x, y);
                        pairs++;
                        if (expected) equalPairs++;
                        Assert.That(SaisReference.LmsSubstringsEqual(reader, isS, lms[x], lms[y]), Is.EqualTo(expected),
                            $"'{s}': LMS substrings at {lms[x]} and {lms[y]}");
                    }
            }
        Assert.That(equalPairs, Is.GreaterThan(10_000), "the enumeration must contain many equal pairs");
        TestContext.Out.WriteLine($"{pairs} ordered LMS pairs checked, {equalPairs} equal");
    }

    /// <summary>Same length, same symbols and same types, where the substring ending at the sentinel equals nothing.</summary>
    private static bool DefinitionalLmsEqual(string s, bool[] isS, List<int> lms, int x, int y)
    {
        int ax = lms[x], bx = lms[x + 1], ay = lms[y], by = lms[y + 1]; // inclusive end positions
        if (bx == s.Length || by == s.Length) return false;
        if (bx - ax != by - ay) return false;
        for (int d = 0; d <= bx - ax; d++)
            if (s[ax + d] != s[ay + d] || isS[ax + d] != isS[ay + d]) return false;
        return true;
    }

    [Test]
    public void Utf16Texts_AreSortedByCodeUnit()
    {
        AssertEqualsNaive("￿\uDC00\uD800\0a😀􏿿", "mixed units");
        foreach (ulong seed in new ulong[] { 1, 2, 3 })
        {
            AssertEqualsNaive(StringFamilies.Bmp(1500, seed), "bmp");
            AssertEqualsNaive(StringFamilies.SurrogatePairs(1500, seed), "pairs");
            AssertEqualsNaive(StringFamilies.LoneSurrogates(1500, seed), "lone");
            AssertEqualsNaive(StringFamilies.NulAndFfff(1500, seed), "nul/ffff");
        }
    }

    [Test]
    public void AllDistinctAlphabet_UpTo65536_IsValid()
    {
        foreach (int n in new[] { 255, 256, 257, 1000, 65535, 65536 })
        {
            string s = StringFamilies.AllDistinct(n, 6);
            Assert.That(SuffixArrayVerifier.Check(s, SaisReference.BuildSuffixArray(s)), Is.Null, $"n {n}");
        }
    }

    [Test]
    public void BoundaryLengths_Around2To16_PassTheLinearVerifier()
    {
        foreach (int n in new[] { 65535, 65536, 65537 })
            foreach (var (name, text) in new[]
            {
                ("a^n", StringFamilies.AllSame(n)),
                ("abab", StringFamilies.AbRepeat(n)),
                ("fibonacci", StringFamilies.Fibonacci(n)),
                ("random-dna", StringFamilies.RandomDna(n, 4)),
                ("random-binary", StringFamilies.Random(n, 2, 5)),
                ("thue-morse", StringFamilies.ThueMorse(n)),
            })
            {
                string? problem = SuffixArrayVerifier.Check(text, SaisReference.BuildSuffixArray(text));
                Assert.That(problem, Is.Null, $"{name} n {n}: {problem}");
            }
    }

    [Test]
    public void RecursionDepth_AtMostLogN_AndEveryLevelAtLeastHalves()
    {
        foreach (int n in new[] { 1000, 100_000 })
            foreach (var (name, text) in new[]
            {
                ("a^n", StringFamilies.AllSame(n)),
                ("fibonacci", StringFamilies.Fibonacci(n)),
                ("thue-morse", StringFamilies.ThueMorse(n)),
                ("abab", StringFamilies.AbRepeat(n)),
                ("random-dna", StringFamilies.RandomDna(n, 9)),
                ("periodic-aab-p5", StringFamilies.PeriodicAaB(n, 5)),
            })
            {
                var trace = new SaisTrace();
                int[] sa = SaisReference.BuildSuffixArray(text, trace);
                Assert.That(SuffixArrayVerifier.Check(text, sa), Is.Null, name);
                AssertDepthBounds(trace, n, name);
            }
    }

    internal static void AssertDepthBounds(SaisTrace trace, int n, string name)
    {
        int bound = (int)Math.Ceiling(Math.Log2(Math.Max(n, 1))) + 1;
        Assert.That(trace.Depth, Is.InRange(1, bound), $"{name}: levels {trace.Depth} > ceil(log2 n) + 1 = {bound}");
        Assert.That(trace.Levels[0].Length, Is.EqualTo(n), name);
        for (int k = 0; k + 1 < trace.Depth; k++)
        {
            var cur = trace.Levels[k];
            var next = trace.Levels[k + 1];
            Assert.That(next.Length, Is.EqualTo(cur.LmsCount), $"{name} level {k}: the reduced string has one symbol per LMS position");
            Assert.That(cur.LmsCount, Is.LessThanOrEqualTo(cur.Length / 2), $"{name} level {k}: m <= n/2");
            Assert.That(next.Length, Is.LessThan(cur.Length), $"{name} level {k}");
        }
    }

    [Test]
    public void Trace_RecordsLevelsOfTheKnownExample()
    {
        var trace = new SaisTrace();
        SaisReference.BuildSuffixArray("mmiissiissiippii", trace);
        Assert.That(trace.Depth, Is.GreaterThanOrEqualTo(1));
        Assert.That(trace.Levels[0].Length, Is.EqualTo(16));
        Assert.That(trace.Levels[0].AlphabetSize, Is.EqualTo(4)); // i, m, p, s
        Assert.That(trace.Levels[0].LmsCount, Is.GreaterThan(0));
    }

    [Test]
    public void SizeGuard_FakeReadersAboveTheLimit_ThrowBeforeAnyAllocationOrRead()
    {
        Assert.Throws<NotSupportedException>(() => SaisReference.Build(new FakeReader(int.MaxValue), 4));
        Assert.Throws<NotSupportedException>(() => SaisReference.Build(new FakeReader(Array.MaxLength), 4));
        var span = new int[4];
        Assert.Throws<NotSupportedException>(() => SaisReference.Build(new FakeReader(Array.MaxLength), 4, span));
    }

    [Test]
    public void InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentException>(() => SaisReference.Build(new IntSymbols(new[] { 0, 1, 2 }), 2));           // symbol 2 >= sigma
        Assert.Throws<ArgumentException>(() => SaisReference.Build(new IntSymbols(new[] { 0, -1 }), 2));             // negative symbol
        Assert.Throws<ArgumentException>(() => SaisReference.Build(new IntSymbols(new[] { 0, 1 }), 2, new int[2]));   // sa must have n + 1 entries
        Assert.Throws<ArgumentOutOfRangeException>(() => SaisReference.Build(new IntSymbols(new[] { 0, 1 }), -1));
    }

    [Test]
    public void IntegerReader_GivesTheSameSaAsTheStringPath()
    {
        var rng = new SplitMix64(5);
        for (int i = 0; i < 200; i++)
        {
            int n = rng.NextInt(300);
            string s = StringFamilies.Random(n, 1 + rng.NextInt(30), rng.Next());
            var p = TextPreparation.Prepare(s);
            int[] ints = Enumerable.Range(0, n).Select(p.SymbolAt).ToArray();
            Assert.That(SaisReference.Build(new IntSymbols(ints), p.AlphabetSize), Is.EqualTo(SaisReference.BuildSuffixArray(s)));
        }
    }

    private readonly struct FakeReader(int length) : ISymbolReader
    {
        public int Length { get; } = length;

        public int this[int index] => throw new InvalidOperationException("the size guard must fire before any symbol is read");
    }
}
