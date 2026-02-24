using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests.Validation;

[TestFixture]
[Category("Validation")]
public class LcsAndEnumerationContractTests
{
    private string? _tempDir;

    [SetUp]
    public void SetUp()
    {
        _tempDir = Path.Combine(Path.GetTempPath(), "SuffixTree_LcsContracts_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_tempDir);
    }

    [TearDown]
    public void TearDown()
    {
        if (_tempDir != null && Directory.Exists(_tempDir))
        {
            try { Directory.Delete(_tempDir, true); } catch { }
        }
    }

    [TestCase(false, TestName = "LcsInfoAndFindAll_MatchReference_Heap")]
    [TestCase(true, TestName = "LcsInfoAndFindAll_MatchReference_Mmf")]
    public void LcsInfoAndFindAll_MatchReference(bool useMappedStorage)
    {
        var cases = new (string Text, string Other)[]
        {
            ("banana", "bandana"),
            ("mississippi", "ssi"),
            ("abracadabra", "cadabrax"),
            ("abcxyzabc", "xyzabc"),
            ("aaaaab", "baaaa"),
            ("abc", "def"),
            ("🧬αβγ🧪", "xxαβγyy"),
            ("A\u0000B\uFFFFC", "\u0000B\uFFFF"),
        };

        foreach (var (text, other) in cases)
        {
            using var persistent = CreateTree(text, useMappedStorage);
            var reference = global::SuffixTree.SuffixTree.Build(text);

            string expectedLcs = reference.LongestCommonSubstring(other);
            var expectedInfo = reference.LongestCommonSubstringInfo(other);
            var expectedAll = reference.FindAllLongestCommonSubstrings(other);

            string actualLcs = persistent.LongestCommonSubstring(other);
            string actualLcsSpan = persistent.LongestCommonSubstring(other.AsSpan());
            var actualInfo = persistent.LongestCommonSubstringInfo(other);
            var actualInfoSpan = persistent.LongestCommonSubstringInfo(other.AsSpan());
            var actualAll = persistent.FindAllLongestCommonSubstrings(other);

            Assert.Multiple(() =>
            {
                Assert.That(actualLcs, Is.EqualTo(expectedLcs), $"text=\"{text}\", other=\"{other}\": LCS mismatch");
                Assert.That(actualLcsSpan, Is.EqualTo(actualLcs), $"text=\"{text}\", other=\"{other}\": span LCS mismatch");

                Assert.That(actualInfoSpan, Is.EqualTo(actualInfo), $"text=\"{text}\", other=\"{other}\": span LCSInfo mismatch");
                Assert.That(actualInfo.Substring, Is.EqualTo(actualLcs), $"text=\"{text}\", other=\"{other}\": info substring mismatch");
                Assert.That(actualInfo.Substring, Is.EqualTo(expectedInfo.Substring),
                    $"text=\"{text}\", other=\"{other}\": LCSInfo substring mismatch");

                Assert.That(actualAll.Substring, Is.EqualTo(expectedAll.Substring), $"text=\"{text}\", other=\"{other}\": FindAll substring mismatch");
                Assert.That(actualAll.Substring, Is.EqualTo(actualLcs), $"text=\"{text}\", other=\"{other}\": FindAll substring != LCS");

                Assert.That(actualAll.PositionsInText.OrderBy(x => x).ToList(),
                    Is.EqualTo(expectedAll.PositionsInText.OrderBy(x => x).ToList()),
                    $"text=\"{text}\", other=\"{other}\": PositionsInText mismatch");
                Assert.That(actualAll.PositionsInOther.OrderBy(x => x).ToList(),
                    Is.EqualTo(expectedAll.PositionsInOther.OrderBy(x => x).ToList()),
                    $"text=\"{text}\", other=\"{other}\": PositionsInOther mismatch");

                if (actualInfo.Substring.Length == 0)
                {
                    Assert.That(actualInfo.PositionInText, Is.EqualTo(-1), $"text=\"{text}\", other=\"{other}\": info.PositionInText");
                    Assert.That(actualInfo.PositionInOther, Is.EqualTo(-1), $"text=\"{text}\", other=\"{other}\": info.PositionInOther");
                }
                else
                {
                    Assert.That(expectedAll.PositionsInText.Contains(actualInfo.PositionInText), Is.True,
                        $"text=\"{text}\", other=\"{other}\": info.PositionInText should be one of FindAll positions");
                    Assert.That(expectedAll.PositionsInOther.Contains(actualInfo.PositionInOther), Is.True,
                        $"text=\"{text}\", other=\"{other}\": info.PositionInOther should be one of FindAll positions");
                }
            });

            AssertAllPositionsAreValid(text, other, actualAll);
            AssertInfoPositionBelongsToAll(actualInfo, actualAll);
        }
    }

    [TestCase(false, TestName = "EnumerateSuffixes_MatchesGetAllAndReference_Heap")]
    [TestCase(true, TestName = "EnumerateSuffixes_MatchesGetAllAndReference_Mmf")]
    public void EnumerateSuffixes_MatchesGetAllAndReference(bool useMappedStorage)
    {
        var texts = new List<string>
        {
            string.Empty,
            "a",
            "banana",
            "mississippi",
            "abcabcabc",
            "🧬αβγ🧪$",
        };

        for (int seed = 0; seed < 10; seed++)
        {
            var rng = new Random(seed + 5000);
            texts.Add(GenerateRandom(rng, rng.Next(1, 80), "abc"));
        }

        foreach (string text in texts)
        {
            using var persistent = CreateTree(text, useMappedStorage);
            var reference = global::SuffixTree.SuffixTree.Build(text);

            var enumerated = persistent.EnumerateSuffixes().ToList();
            var fromGetAll = persistent.GetAllSuffixes().ToList();
            var expected = reference.GetAllSuffixes().ToList();

            Assert.Multiple(() =>
            {
                Assert.That(enumerated, Is.EqualTo(fromGetAll), $"text=\"{text}\": Enumerate vs GetAll mismatch");
                Assert.That(enumerated, Is.EqualTo(expected), $"text=\"{text}\": persistent suffix enumeration mismatch");
            });
        }
    }

    [TestCase(false, TestName = "LongestRepeatedSubstringMemory_MatchesStringResult_Heap")]
    [TestCase(true, TestName = "LongestRepeatedSubstringMemory_MatchesStringResult_Mmf")]
    public void LongestRepeatedSubstringMemory_MatchesStringResult(bool useMappedStorage)
    {
        string[] texts =
        {
            "",
            "banana",
            "mississippi",
            "abcdefg",
            "aaaaaa",
            "abcabcxyzabcabc",
            "🧬αβγ🧪🧬αβγ🧪",
        };

        foreach (string text in texts)
        {
            using var tree = CreateTree(text, useMappedStorage);
            ISuffixTree contract = tree;

            string lrs = tree.LongestRepeatedSubstring();
            string memoryLrs = contract.LongestRepeatedSubstringMemory().ToString();

            Assert.That(memoryLrs, Is.EqualTo(lrs), $"text=\"{text}\": memory/string LRS mismatch");
            if (lrs.Length > 0)
            {
                Assert.That(tree.CountOccurrences(lrs.AsSpan()), Is.GreaterThanOrEqualTo(2),
                    $"text=\"{text}\": LRS from memory must repeat");
            }
        }
    }

    private PersistentSuffixTree CreateTree(string text, bool useMappedStorage)
    {
        if (!useMappedStorage)
            return (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(new StringTextSource(text));

        string filePath = Path.Combine(_tempDir!, $"tree_{Guid.NewGuid():N}.st");
        return (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(new StringTextSource(text), filePath);
    }

    private static void AssertAllPositionsAreValid(
        string text,
        string other,
        (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) all)
    {
        if (all.Substring.Length == 0)
        {
            Assert.That(all.PositionsInText, Is.Empty);
            Assert.That(all.PositionsInOther, Is.Empty);
            return;
        }

        foreach (int pos in all.PositionsInText)
        {
            Assert.That(pos, Is.GreaterThanOrEqualTo(0), "PositionsInText contains negative index");
            Assert.That(pos + all.Substring.Length, Is.LessThanOrEqualTo(text.Length), "PositionsInText exceeds text length");
            Assert.That(text.Substring(pos, all.Substring.Length), Is.EqualTo(all.Substring),
                "PositionsInText does not point to the returned substring");
        }

        foreach (int pos in all.PositionsInOther)
        {
            Assert.That(pos, Is.GreaterThanOrEqualTo(0), "PositionsInOther contains negative index");
            Assert.That(pos + all.Substring.Length, Is.LessThanOrEqualTo(other.Length), "PositionsInOther exceeds other length");
            Assert.That(other.Substring(pos, all.Substring.Length), Is.EqualTo(all.Substring),
                "PositionsInOther does not point to the returned substring");
        }
    }

    private static void AssertInfoPositionBelongsToAll(
        (string Substring, int PositionInText, int PositionInOther) info,
        (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) all)
    {
        if (info.Substring.Length == 0)
        {
            Assert.That(info.PositionInText, Is.EqualTo(-1));
            Assert.That(info.PositionInOther, Is.EqualTo(-1));
            return;
        }

        Assert.That(all.PositionsInText.Contains(info.PositionInText), Is.True,
            "LCSInfo.PositionInText should be one of FindAll positions");
        Assert.That(all.PositionsInOther.Contains(info.PositionInOther), Is.True,
            "LCSInfo.PositionInOther should be one of FindAll positions");
    }

    private static string GenerateRandom(Random rng, int length, string alphabet)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }
}
