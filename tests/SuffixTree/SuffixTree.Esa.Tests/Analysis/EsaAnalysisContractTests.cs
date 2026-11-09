using System.Reflection;
using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Analysis;

/// <summary>
/// Cross-cutting checks of the native analysis operations (work package 1b.3): the signatures are exactly those of the interfaces that the facade of work package 1b.6 forwards to,
/// and the answers do not depend on how the LCP array is stored (the three representations of <see cref="CompactLcp"/>, escapes included).
/// </summary>
[TestFixture]
public class EsaAnalysisContractTests
{
    [TestCase(typeof(ISuffixTreeAnalysis), "LongestRepeatedSubstring")]
    [TestCase(typeof(ISuffixTreeAnalysis), "LongestRepeatedSubstringMemory")]
    [TestCase(typeof(ISuffixTreeAnalysis), "GetAllSuffixes")]
    [TestCase(typeof(ISuffixTreeAnalysis), "EnumerateSuffixes")]
    [TestCase(typeof(ISuffixTree), "FindAllLongestRepeatedSubstrings")]
    [TestCase(typeof(ISuffixTree), "CountDistinctSubstrings")]
    [TestCase(typeof(ISuffixTree), "CountDistinctSubstringsByLength")]
    public void Signature_IsThatOfTheInterfaceMember(Type contract, string name)
    {
        MethodInfo expected = contract.GetMethod(name, BindingFlags.Public | BindingFlags.Instance)!;
        MethodInfo? actual = typeof(EsaIndex).GetMethod(name, BindingFlags.NonPublic | BindingFlags.Public | BindingFlags.Instance);
        Assert.That(actual, Is.Not.Null, name);
        Assert.That(actual!.ReturnType, Is.EqualTo(expected.ReturnType), $"{name}: return type");
        Assert.That(actual.GetParameters().Select(p => p.ParameterType), Is.EqualTo(expected.GetParameters().Select(p => p.ParameterType)), $"{name}: parameters");
    }

    private static EsaIndex WithLcpKind(EsaIndex source, LcpKind kind) =>
        EsaIndex.FromArrays(source.Text, source.Prepared, source.SuffixArray, CompactLcp.FromDense(source.Lcp.ToArray(), kind), EsaIndex.DefaultLayoutFor(source.Length));

    private static void AssertSameAnswers(EsaIndex expected, EsaIndex actual, string context)
    {
        Assert.That(actual.CountDistinctSubstrings(), Is.EqualTo(expected.CountDistinctSubstrings()), $"{context}: total");
        foreach (int m in new[] { 1, 3, expected.Length, expected.Length + 1 }.Where(m => m >= 1))
            Assert.That(actual.CountDistinctSubstringsByLength(m), Is.EqualTo(expected.CountDistinctSubstringsByLength(m)), $"{context}: ByLength({m})");
        Assert.That(actual.LongestRepeatedSubstring(), Is.EqualTo(expected.LongestRepeatedSubstring()), $"{context}: LRS");
        var a = actual.FindAllLongestRepeatedSubstrings();
        var e = expected.FindAllLongestRepeatedSubstrings();
        Assert.That(a.Select(x => x.Substring), Is.EqualTo(e.Select(x => x.Substring)), $"{context}: ties");
        Assert.That(a.Select(x => x.Positions), Is.EqualTo(e.Select(x => x.Positions)), $"{context}: positions");
    }

    [Test]
    public void AnswersDoNotDependOnTheLcpRepresentation()
    {
        foreach (string text in new[] { "banana", "mississippi", "abcabcxyzxyz", new string('a', 300), StringFamilies.Fibonacci(700), StringFamilies.RandomDna(400, 3) + StringFamilies.RandomDna(400, 3) })
        {
            EsaIndex reference = EsaIndex.Build(text);
            foreach (LcpKind kind in Enum.GetValues<LcpKind>())
                AssertSameAnswers(reference, WithLcpKind(reference, kind), $"{kind} n {text.Length}");
        }
    }

    [Test]
    public void EscapedLcpValues_AboveSixtyFiveThousand_AreResolved()
    {
        // a^70000: LCP values up to 69,999, so even the U16 representation has escapes; every answer has a closed form
        const int n = 70_000;
        EsaIndex reference = EsaIndex.Build(new string('a', n));
        Assert.That(reference.Lcp.Kind, Is.Not.EqualTo(LcpKind.U8).Or.Property("OverflowCount").GreaterThan(0));
        foreach (LcpKind kind in Enum.GetValues<LcpKind>())
        {
            EsaIndex index = WithLcpKind(reference, kind);
            Assert.That(index.CountDistinctSubstrings(), Is.EqualTo(n), kind.ToString());
            long[] byLength = index.CountDistinctSubstringsByLength(n);
            Assert.That(byLength, Has.Length.EqualTo(n + 1), kind.ToString());
            Assert.That(byLength, Is.All.EqualTo(1L), kind.ToString());
            Assert.That(index.LongestRepeatedSubstring(), Has.Length.EqualTo(n - 1), kind.ToString());
            var all = index.FindAllLongestRepeatedSubstrings();
            Assert.That(all, Has.Count.EqualTo(1), kind.ToString());
            Assert.That(all[0].Positions, Is.EqualTo(new[] { 0, 1 }), kind.ToString());
        }
    }

    [Test]
    public void LongRepeatInsideRandomText_IsFoundAtTheEnds()
    {
        // the repeat starts at rank-order extremes: all A's block sorts first, a long copy at the end
        string block = StringFamilies.RandomDna(3000, 21);
        string text = block + "T" + block;
        var all = EsaIndex.Build(text).FindAllLongestRepeatedSubstrings();
        UkkonenTree ukkonen = UkkonenTree.Build(text);
        Assert.That(all.Select(e => (e.Substring, e.Positions)), Is.EqualTo(ukkonen.FindAllLongestRepeatedSubstrings().Select(e => (e.Substring, e.Positions))));
        Assert.That(all[0].Substring, Is.EqualTo(block));
        Assert.That(all[0].Positions, Is.EqualTo(new[] { 0, 3001 }));
    }

    [Test]
    public void ManyTies_AreAllListed_AndTheirPositionsAreBoundedByTheTextLength()
    {
        // 600 distinct 8-mers, each twice, every copy followed by its own separator code unit: no repeat contains a separator, so the ties are exactly the 600 8-mers
        var rng = new SplitMix64(0x1B3D);
        var seen = new HashSet<string>();
        var parts = new List<string>();
        while (parts.Count < 600)
        {
            string part = string.Concat(Enumerable.Range(0, 8).Select(_ => "ACGT"[rng.NextInt(4)]));
            if (seen.Add(part)) parts.Add(part);
        }
        string text = string.Concat(parts.Select((part, i) => part + (char)(0x100 + i))) + string.Concat(parts.Select((part, i) => part + (char)(0x1000 + i)));
        var all = EsaIndex.Build(text).FindAllLongestRepeatedSubstrings();
        var expected = UkkonenTree.Build(text).FindAllLongestRepeatedSubstrings();
        Assert.That(all.Select(e => e.Substring), Is.EqualTo(expected.Select(e => e.Substring)));
        Assert.That(all.Select(e => e.Positions), Is.EqualTo(expected.Select(e => e.Positions)));
        Assert.That(all.Sum(e => e.Positions.Count), Is.LessThanOrEqualTo(text.Length));
        Assert.That(all, Has.Count.EqualTo(600));
    }
}
