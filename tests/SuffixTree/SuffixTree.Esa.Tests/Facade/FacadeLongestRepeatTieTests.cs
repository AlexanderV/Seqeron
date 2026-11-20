using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The tie-break of two longest repeated substrings of equal length (the gap of the audit of 1b.8 that the harness did not see, because it compared only the length, the double occurrence and the identity of the cached string).
/// ADR-001 section 7 leaves the representative implementation-specific <i>by contract</i> (<c>ISuffixTreeAnalysis.LongestRepeatedSubstring</c>: "on a length tie the returned representative is implementation-specific; its
/// length is always maximal"), and the other engines may differ. Both ESA engines, the in-memory tree and the file-backed tree (whose answer comes from the header of the container), keep one rule and document it: <b>the
/// lexicographically smallest of the tied substrings</b> (the deepest interval that closes first in the scan, i.e. the lowest rank, since two tied intervals are disjoint). The harness compares the two engines with each other and
/// with the rule (<c>FullFuzzCheck</c>); these tests pin the rule on known answers, on every combination of codec, access mode and validation level, and on random texts.
/// </summary>
[TestFixture]
public class FacadeLongestRepeatTieTests
{
    private static readonly string[] Ties =
    [
        "efghXefghYabcdZabcd",           // the smallest tied substring ("abcd") is the later one in the text
        "abcdXabcdYefghZefgh",           // ... and the earlier one
        "bbbXbbbYaaaZaaa",
        "cherry pie and apple pie and cherry tart and apple tart",
        "mississippi",                   // "issi" occurs at 1 and 4 (overlapping) and is the only longest repeat
        "abcabc defdef",                 // "abc", "def" (and "bc"...): ties of 3
        "xyzxyz abcabc mnomno",          // three tied substrings of length 3, the smallest in the middle of the text
        "aaaa",
        "ab",
    ];

    /// <summary>The longest substring that occurs at least twice (overlaps allowed), the smallest by ordinal order among those of that length; empty when no symbol repeats.</summary>
    private static string ByDefinition(string text)
    {
        for (int length = text.Length - 1; length >= 1; length--)
        {
            string? best = null;
            var seen = new HashSet<string>(StringComparer.Ordinal);
            for (int i = 0; i + length <= text.Length; i++)
            {
                string s = text.Substring(i, length);
                if (!seen.Add(s) && (best is null || string.CompareOrdinal(s, best) < 0)) best = s;
            }
            if (best is not null) return best;
        }
        return string.Empty;
    }

    [TestCaseSource(nameof(Ties))]
    public void TheRepresentative_OfEqualLongestRepeats_IsTheLexicographicallySmallest_InTheInMemoryTree(string text)
    {
        string expected = ByDefinition(text);
        using EsaSuffixTree tree = EsaSuffixTree.Build(text);
        Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(expected));
        Assert.That(tree.LongestRepeatedSubstringMemory().ToString(), Is.EqualTo(expected));
        string[] tied = [.. tree.FindAllLongestRepeatedSubstrings().Select(t => t.Substring)];
        if (expected.Length == 0) Assert.That(tied, Is.Empty);
        else Assert.That(tied.Min(StringComparer.Ordinal), Is.EqualTo(expected));
    }

    [TestCaseSource(nameof(Ties))]
    public void TheRepresentative_OfEqualLongestRepeats_IsTheLexicographicallySmallest_InTheFileBackedTree_EveryCodecAccessAndLevel(string text)
    {
        string expected = ByDefinition(text);
        foreach (ChildCodec codec in new[] { ChildCodec.Packed, ChildCodec.Tagged })
        {
            foreach (EsaAccessMode access in new[] { EsaAccessMode.Checked, EsaAccessMode.Trusted })
            {
                foreach (EsaValidationLevel level in new[] { EsaValidationLevel.Header, EsaValidationLevel.Structure })
                {
                    using EsaSuffixTree tree = SesaTestFile.FileTree(text, access: access, level: level, codec: codec);
                    Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(expected), $"{codec} {access} {level}");
                    Assert.That(tree.LongestRepeatedSubstringMemory().ToString(), Is.EqualTo(expected), $"{codec} {access} {level}");
                }
            }
        }
    }

    [TestCaseSource(nameof(Ties))]
    public void TheRepresentative_OfEqualLongestRepeats_IsTheLexicographicallySmallest_InAContainerOfTheEndToEndBuilder(string text)
    {
        // the header of a container that the production builder wrote (its own scan keeps the first deepest interval), at the Structure level, which recomputes it from the lcp-interval tree
        string expected = ByDefinition(text);
        using var dir = new TempDir();
        SesaRun.Outcome built = SesaRun.Build(dir, text, SesaRun.Tiny(text.Length, 2, psi: true), "-t");
        foreach (EsaValidationLevel level in new[] { EsaValidationLevel.Header, EsaValidationLevel.Structure })
        {
            using EsaSuffixTree tree = SesaTestFile.Open(built.Target, level);
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(expected), $"{level}");
        }
    }

    [Test]
    public void TheRepresentative_OnRandomTexts_IsTheLexicographicallySmallest_InBothEngines_AndTheEnginesAgree()
    {
        var rng = new SplitMix64(20261120);
        int ties = 0;
        for (int i = 0; i < 400; i++)
        {
            int sigma = 2 + rng.NextInt(4);
            string text = StringFamilies.Random(8 + rng.NextInt(40), sigma, (ulong)i + 1);
            string expected = ByDefinition(text);
            using EsaSuffixTree memory = EsaSuffixTree.Build(text);
            using EsaSuffixTree file = SesaTestFile.FileTree(text, codec: i % 2 == 0 ? ChildCodec.Packed : ChildCodec.Tagged);
            Assert.That(memory.LongestRepeatedSubstring(), Is.EqualTo(expected), $"memory, '{text}'");
            Assert.That(file.LongestRepeatedSubstring(), Is.EqualTo(expected), $"file, '{text}'");
            if (memory.FindAllLongestRepeatedSubstrings().Count > 1) ties++;
        }
        Assert.That(ties, Is.GreaterThan(20), "the random texts include ties");
    }
}
