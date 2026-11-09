using System.Runtime.InteropServices;
using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Analysis;

/// <summary>
/// <c>LongestRepeatedSubstring</c>, <c>LongestRepeatedSubstringMemory</c> and <c>FindAllLongestRepeatedSubstrings</c> of the in-memory ESA (work package 1b.3). Oracles: the
/// definition (<see cref="AnalysisReference"/>, L0), the Ukkonen tree (L2) and the <c>Traverse</c>-default path of the generic algorithm over the ESA (<see cref="EsaTreeShim"/>).
/// The list operation must be identical (same entries, same order); the single representative is implementation-specific on a length tie (<c>ISuffixTreeAnalysis.cs:13</c>),
/// so it is compared by property: maximal length, one of the tied substrings, and the ESA rule (the lexicographically smallest tied substring).
/// </summary>
[TestFixture]
public class EsaLongestRepeatTests
{
    private static void AssertSameList(IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> expected, IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> actual, string what, string context)
    {
        string Show(IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> l) =>
            string.Join(" | ", l.Take(6).Select(e => $"'{(e.Substring.Length > 40 ? e.Substring[..40] + "..." : e.Substring)}'{e.Substring.Length}@[{string.Join(",", e.Positions.Take(8))}]"));
        if (expected.Count != actual.Count)
            Assert.Fail($"{context}: {what}: {actual.Count} entries, expected {expected.Count}: actual {Show(actual)} vs expected {Show(expected)}");
        for (int i = 0; i < expected.Count; i++)
        {
            if (!string.Equals(expected[i].Substring, actual[i].Substring, StringComparison.Ordinal) || !expected[i].Positions.SequenceEqual(actual[i].Positions))
                Assert.Fail($"{context}: {what}: entry {i} differs: actual {Show(actual.Skip(i).ToList())} vs expected {Show(expected.Skip(i).ToList())}");
        }
    }

    private static IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> AsList(List<(string Substring, int[] Positions)> reference) =>
        reference.Select(e => (e.Substring, (IReadOnlyList<int>)e.Positions)).ToList();

    private static void AssertRepresentative(string text, string lrs, IReadOnlyList<(string Substring, IReadOnlyList<int> Positions)> all, UkkonenTree ukkonen, string context)
    {
        string fromUkkonen = ukkonen.LongestRepeatedSubstring();
        Assert.That(lrs.Length, Is.EqualTo(fromUkkonen.Length), $"{context}: length of LongestRepeatedSubstring");
        if (all.Count == 0)
        {
            Assert.That(lrs, Is.Empty, $"{context}: no repeat");
            return;
        }
        Assert.That(all.Select(e => e.Substring), Does.Contain(lrs), $"{context}: one of the tied substrings");
        Assert.That(lrs, Is.EqualTo(all.Select(e => e.Substring).Min(StringComparer.Ordinal)), $"{context}: the ESA tie rule (lexicographically smallest)");
        Assert.That(all.Select(e => e.Substring), Does.Contain(fromUkkonen), $"{context}: the Ukkonen representative is a tied substring too");
        if (all.Count == 1) Assert.That(lrs, Is.EqualTo(fromUkkonen), $"{context}: a unique longest repeat is the same string in every implementation");
    }

    private static void Check(string text, string context)
    {
        EsaIndex esa = EsaIndex.Build(text);
        UkkonenTree ukkonen = UkkonenTree.Build(text);
        var all = esa.FindAllLongestRepeatedSubstrings();
        AssertSameList(AsList(AnalysisReference.AllLongestRepeated(text)), all, "against the definition", context);
        AssertSameList(ukkonen.FindAllLongestRepeatedSubstrings(), all, "against Ukkonen", context);
        AssertSameList(SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings(new EsaTreeShim(esa)), all, "against the Traverse-default path", context);

        string lrs = esa.LongestRepeatedSubstring();
        AssertRepresentative(text, lrs, all, ukkonen, context);
        ReadOnlyMemory<char> memory = esa.LongestRepeatedSubstringMemory();
        Assert.That(memory.Length, Is.EqualTo(lrs.Length), $"{context}: Memory length");
        Assert.That(memory.Span.SequenceEqual(lrs), Is.True, $"{context}: Memory content");
    }

    [Test]
    public void KnownTexts_LiteratureValues()
    {
        // the expected values of AllLongestRepeatedSubstringsTests (Python brute force, Gusfield 1997 section 7.1)
        var literature = new (string Text, (string, int[])[] Expected)[]
        {
            ("banana", [("ana", [1, 3])]),
            ("mississippi", [("issi", [1, 4])]),
            ("abracadabra", [("abra", [0, 7])]),
            ("ATATAT", [("ATAT", [0, 2])]),
            ("aaaa", [("aaa", [0, 1])]),
            ("aabb", [("a", [0, 1]), ("b", [2, 3])]),
            ("abcabcxyzxyz", [("abc", [0, 3]), ("xyz", [6, 9])]),
            ("GATTACAGATTTACA", [("TTACA", [2, 10])]),
            ("abcxbcaxcab", [("ab", [0, 9]), ("bc", [1, 4]), ("ca", [5, 8])]),
            ("xyzzyxxz", [("x", [0, 5, 6]), ("y", [1, 4]), ("z", [2, 3, 7])]),
        };
        foreach (var (text, expected) in literature)
            AssertSameList(expected.Select(e => (e.Item1, (IReadOnlyList<int>)e.Item2)).ToList(), EsaIndex.Build(text).FindAllLongestRepeatedSubstrings(), "literature", $"'{text}'");
    }

    [Test]
    public void EmptyOneAndNoRepeat_HaveNoRepeatedSubstring()
    {
        foreach (string text in new[] { "", "a", "ab", "abc", "abcd", "\0￿", "zyx" })
        {
            EsaIndex esa = EsaIndex.Build(text);
            Assert.That(esa.LongestRepeatedSubstring(), Is.Empty, $"'{text}'");
            Assert.That(esa.LongestRepeatedSubstringMemory().Length, Is.EqualTo(0), $"'{text}'");
            Assert.That(esa.FindAllLongestRepeatedSubstrings(), Is.Empty, $"'{text}'");
            Check(text, $"'{text}'");
        }
    }

    [Test]
    public void AllEqualText_RepeatIsTheTextWithoutItsLastSymbol()
    {
        foreach (int n in new[] { 2, 3, 4, 10, 100, 1000 })
        {
            string text = new('a', n);
            EsaIndex esa = EsaIndex.Build(text);
            Assert.That(esa.LongestRepeatedSubstring(), Is.EqualTo(new string('a', n - 1)), $"a^{n}");
            var all = esa.FindAllLongestRepeatedSubstrings();
            Assert.That(all, Has.Count.EqualTo(1), $"a^{n}");
            Assert.That(all[0].Substring, Is.EqualTo(new string('a', n - 1)), $"a^{n}");
            Assert.That(all[0].Positions, Is.EqualTo(new[] { 0, 1 }), $"a^{n}");
        }
    }

    [Test]
    public void Ties_AreListedByFirstOccurrence_AndTheRepresentativeIsTheSmallest()
    {
        // dc at 0 and 2, ba at 4 and 6, xy at 8 and 10: no repeat of length 3; the list is ordered by first occurrence, the representative is the lexicographically smallest
        EsaIndex ties = EsaIndex.Build("dcdcbabaxyxy");
        AssertSameList(
            [("dc", new[] { 0, 2 }), ("ba", new[] { 4, 6 }), ("xy", new[] { 8, 10 })],
            ties.FindAllLongestRepeatedSubstrings(), "ties", "'dcdcbabaxyxy'");
        Assert.That(ties.LongestRepeatedSubstring(), Is.EqualTo("ba"));
    }

    [Test]
    public void Memory_IsASliceOfTheOriginalText()
    {
        EsaIndex esa = EsaIndex.Build("xabcyabcz");
        ReadOnlyMemory<char> memory = esa.LongestRepeatedSubstringMemory();
        Assert.That(MemoryMarshal.TryGetString(memory, out string? backing, out int start, out int length), Is.True, "backed by a string, not a copy");
        Assert.That(backing, Is.SameAs(esa.Text));
        Assert.That(length, Is.EqualTo(3));
        Assert.That(esa.Text.Substring(start, length), Is.EqualTo("abc"));
    }

    [Test]
    public void AllSmallStrings_AgreeWithTheDefinitionUkkonenAndTheTraverseDefaultPath()
    {
        foreach (var (name, text) in AnalysisCorpus.Definitional()) Check(text, name);
    }

    [Test]
    public void PositionsAreAscendingAndEveryOccurrenceIsListed_OnLargerTexts()
    {
        // above the reach of the quadratic definition: Ukkonen and the Traverse-default path as oracles, plus the properties of the list
        var rng = new SplitMix64(0x1B3A);
        for (int i = 0; i < 10; i++)
        {
            int n = 2000 + rng.NextInt(18001);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 6 == 0))
            {
                EsaIndex esa = EsaIndex.Build(text);
                UkkonenTree ukkonen = UkkonenTree.Build(text);
                var all = esa.FindAllLongestRepeatedSubstrings();
                AssertSameList(ukkonen.FindAllLongestRepeatedSubstrings(), all, "against Ukkonen", $"{name} n {n}");
                AssertSameList(SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings(new EsaTreeShim(esa)), all, "against the Traverse-default path", $"{name} n {n}");
                AssertRepresentative(text, esa.LongestRepeatedSubstring(), all, ukkonen, $"{name} n {n}");
                foreach (var (substring, positions) in all)
                {
                    Assert.That(positions, Is.Ordered.Ascending, $"{name} n {n}");
                    Assert.That(positions, Is.EqualTo(BruteForce.Occurrences(text, substring)), $"{name} n {n}: every occurrence is listed");
                }
            }
        }
    }

    [Test]
    public void ResultIsIndependentOfTheChildTableLayout()
    {
        foreach (string text in new[] { "banana", "mississippi", "aabbaabb", "abcabcxyzxyz", new string('a', 50), StringFamilies.Fibonacci(100) })
        {
            var expected = EsaIndex.Build(text, ChildLayoutKind.Tagged);
            foreach (ChildLayoutKind layout in Enum.GetValues<ChildLayoutKind>())
            {
                EsaIndex esa = EsaIndex.Build(text, layout);
                AssertSameList(expected.FindAllLongestRepeatedSubstrings(), esa.FindAllLongestRepeatedSubstrings(), $"layout {layout}", $"'{text}'");
                Assert.That(esa.LongestRepeatedSubstring(), Is.EqualTo(expected.LongestRepeatedSubstring()), $"layout {layout} '{text}'");
            }
        }
    }

    [Test]
    public void ConcurrentReaders_GetTheSameAnswers()
    {
        string text = StringFamilies.RandomDna(5000, 7) + StringFamilies.RandomDna(5000, 7);   // one long repeat
        EsaIndex esa = EsaIndex.Build(text);
        string expectedLrs = esa.LongestRepeatedSubstring();
        var expectedAll = esa.FindAllLongestRepeatedSubstrings();
        Parallel.For(0, 4, new ParallelOptions { MaxDegreeOfParallelism = 4 }, _ =>
        {
            for (int k = 0; k < 20; k++)
            {
                Assert.That(esa.LongestRepeatedSubstring(), Is.EqualTo(expectedLrs));
                AssertSameList(expectedAll, esa.FindAllLongestRepeatedSubstrings(), "concurrent", "dna");
            }
        });
        Assert.That(expectedLrs.Length, Is.EqualTo(5000));
    }
}
