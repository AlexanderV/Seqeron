using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Sais;

/// <summary>
/// Heavy checks of the SA-IS reference. Category <c>Exhaustive</c>, <c>[Explicit]</c>: not part of the default run, run with
/// <c>dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter Category=Exhaustive</c>.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log)")]
public class SaisExhaustiveTests
{
    [Test]
    public void AllBinaryStrings_UpTo22_MatchNaive_AndReachThreeLevels()
    {
        var (count, maxDepth, example) = RunExhaustive(2, 22);
        TestContext.Out.WriteLine($"binary <= 22: {count} strings, max level count {maxDepth}, first string at that depth: '{example}'");
        Assert.That(count, Is.EqualTo((1L << 23) - 1));
        // Review MIN-10 (a) asked for depth >= 4 somewhere in the exhaustive run: binary strings up to length 22 reach 3 levels only
        // (first at length 18); the deeper levels are exercised against the naive oracle by the Fibonacci / Thue-Morse families
        // (SaisReferenceTests.AdversarialFamilies_UpTo2000_MatchNaive requires >= 6 levels) and at n = 10^6 below.
        Assert.That(maxDepth, Is.GreaterThanOrEqualTo(3));
    }

    [Test]
    public void AllTernaryStrings_UpTo13_MatchNaive()
    {
        var (count, maxDepth, example) = RunExhaustive(3, 13);
        TestContext.Out.WriteLine($"ternary <= 13: {count} strings, max level count {maxDepth}, first string at that depth: '{example}'");
        Assert.That(count, Is.EqualTo((long)((Math.Pow(3, 14) - 1) / 2)));
    }

    [Test]
    public void AllQuaternaryStrings_UpTo11_MatchNaive()
    {
        var (count, maxDepth, example) = RunExhaustive(4, 11);
        TestContext.Out.WriteLine($"4-ary <= 11: {count} strings, max level count {maxDepth}, first string at that depth: '{example}'");
        Assert.That(count, Is.EqualTo((long)((Math.Pow(4, 12) - 1) / 3)));
    }

    [Test]
    public void RecursionDepth_At1e6_IsBelowTheLogBound_AndHalvesOnEveryLevel()
    {
        const int n = 1_000_000;
        foreach (var (name, text) in new[]
        {
            ("a^n", StringFamilies.AllSame(n)),
            ("fibonacci", StringFamilies.Fibonacci(n)),
            ("thue-morse", StringFamilies.ThueMorse(n)),
            ("abab", StringFamilies.AbRepeat(n)),
            ("random-dna", StringFamilies.RandomDna(n, 21)),
            ("random-binary", StringFamilies.Random(n, 2, 22)),
            ("periodic-aab-p7", StringFamilies.PeriodicAaB(n, 7)),
            ("debruijn-2", StringFamilies.DeBruijn(2, 20, n)),
            ("debruijn-4", StringFamilies.DeBruijn(4, 10, n)),
        })
        {
            var trace = new SaisTrace();
            int[] sa = SaisReference.BuildSuffixArray(text, trace);
            string? problem = SuffixArrayVerifier.Check(text, sa);
            Assert.That(problem, Is.Null, $"{name}: {problem}");
            SaisReferenceTests.AssertDepthBounds(trace, n, name);
            TestContext.Out.WriteLine($"{name}: n = {n}, levels = {trace.Depth} (bound {(int)Math.Ceiling(Math.Log2(n)) + 1}), level sizes {string.Join(",", trace.Levels.Select(l => l.Length))}");
        }
    }

    private static (long Count, int MaxDepth, string Example) RunExhaustive(int sigma, int maxLen)
    {
        var trace = new SaisTrace();
        long count = 0;
        int maxDepth = 0;
        string example = "";
        foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
        {
            trace.Levels.Clear();
            int[] actual = SaisReference.BuildSuffixArray(s, trace);
            if (!NaiveSuffixArray.Build(s).AsSpan().SequenceEqual(actual))
                Assert.Fail($"sigma {sigma}: SA differs from the naive oracle on '{s}'");
            count++;
            if (trace.Depth > maxDepth)
            {
                maxDepth = trace.Depth;
                example = s;
            }
        }
        return (count, maxDepth, example);
    }
}
