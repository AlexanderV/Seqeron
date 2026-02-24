using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests.Validation
{
    [TestFixture]
    [Category("Validation")]
    public class SpanOverloadContractTests
    {
        private string? _tempDir;

        [SetUp]
        public void SetUp()
        {
            _tempDir = Path.Combine(Path.GetTempPath(), "SuffixTree_Span_" + Guid.NewGuid().ToString("N"));
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

        [TestCase(false, TestName = "SpanOverloads_ContainsCountFindAll_MatchString_Heap")]
        [TestCase(true, TestName = "SpanOverloads_ContainsCountFindAll_MatchString_Mmf")]
        public void SpanOverloads_ContainsCountFindAll_MatchString(bool useMappedStorage)
        {
            string[] texts =
            {
                "",
                "banana",
                "mississippi",
                "aaaaa",
                "abcabcabc",
                "🧬αβγ🧪$",
                "A\u0000B\uFFFFC"
            };

            foreach (string text in texts)
            {
                using var treeDisposable = (IDisposable)CreateTree(text, useMappedStorage);
                var tree = (ISuffixTree)treeDisposable;

                var patterns = BuildPatternSet(text, seed: text.Length + 17);
                foreach (string pattern in patterns)
                {
                    bool containsString = tree.Contains(pattern);
                    bool containsSpan = tree.Contains(pattern.AsSpan());
                    Assert.That(containsSpan, Is.EqualTo(containsString),
                        $"Contains mismatch: text=\"{text}\", pattern=\"{pattern}\"");

                    int countString = tree.CountOccurrences(pattern);
                    int countSpan = tree.CountOccurrences(pattern.AsSpan());
                    Assert.That(countSpan, Is.EqualTo(countString),
                        $"Count mismatch: text=\"{text}\", pattern=\"{pattern}\"");

                    var posString = tree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                    var posSpan = tree.FindAllOccurrences(pattern.AsSpan()).OrderBy(x => x).ToList();
                    Assert.That(posSpan, Is.EqualTo(posString),
                        $"FindAll mismatch: text=\"{text}\", pattern=\"{pattern}\"");
                }
            }
        }

        [TestCase(false, TestName = "SpanOverloads_FindAll_MatchesBruteForce_Heap")]
        [TestCase(true, TestName = "SpanOverloads_FindAll_MatchesBruteForce_Mmf")]
        public void SpanOverloads_FindAll_MatchesBruteForce(bool useMappedStorage)
        {
            for (int seed = 0; seed < 20; seed++)
            {
                var rng = new Random(seed);
                string text = GenerateRandom(rng, rng.Next(1, 120), "abcd");

                using var treeDisposable = (IDisposable)CreateTree(text, useMappedStorage);
                var tree = (ISuffixTree)treeDisposable;

                for (int q = 0; q < 12; q++)
                {
                    string pattern;
                    if (q == 0)
                    {
                        pattern = string.Empty;
                    }
                    else if (q % 3 == 0)
                    {
                        pattern = BuildMissingPattern(text);
                    }
                    else
                    {
                        int start = rng.Next(text.Length);
                        int len = rng.Next(1, Math.Min(16, text.Length - start + 1));
                        pattern = text.Substring(start, len);
                    }

                    var expected = BruteForcePositions(text, pattern);
                    var actual = tree.FindAllOccurrences(pattern.AsSpan()).OrderBy(x => x).ToList();

                    Assert.That(actual, Is.EqualTo(expected),
                        $"Seed={seed}, pattern=\"{pattern}\"");
                    Assert.That(tree.CountOccurrences(pattern.AsSpan()), Is.EqualTo(expected.Count),
                        $"Seed={seed}, pattern=\"{pattern}\": Count mismatch");
                    Assert.That(tree.Contains(pattern.AsSpan()), Is.EqualTo(expected.Count > 0 || pattern.Length == 0),
                        $"Seed={seed}, pattern=\"{pattern}\": Contains mismatch");
                }
            }
        }

        [TestCase(false, TestName = "SpanOverload_Lcs_MatchesStringAndOracle_Heap")]
        [TestCase(true, TestName = "SpanOverload_Lcs_MatchesStringAndOracle_Mmf")]
        public void SpanOverload_Lcs_MatchesStringAndOracle(bool useMappedStorage)
        {
            var fixedPairs = new (string Text, string Other)[]
            {
                ("banana", "bandana"),
                ("mississippi", "ssi"),
                ("abc", "def"),
                ("abcdxyz", "xyzabcd"),
                ("🧬αβγ🧪", "xxαβγyy"),
                ("aaaaab", "baaaa"),
            };

            foreach (var (text, other) in fixedPairs)
            {
                using var treeDisposable = (IDisposable)CreateTree(text, useMappedStorage);
                var tree = (ISuffixTree)treeDisposable;

                AssertLcsContract(tree, text, other, context: "fixed");
            }

            for (int seed = 0; seed < 20; seed++)
            {
                var rng = new Random(seed + 1000);
                string text = GenerateRandom(rng, rng.Next(5, 70), "abc");
                string other = GenerateRandom(rng, rng.Next(5, 70), "abc");

                using var treeDisposable = (IDisposable)CreateTree(text, useMappedStorage);
                var tree = (ISuffixTree)treeDisposable;

                AssertLcsContract(tree, text, other, context: $"seed={seed}");
            }
        }

        private ISuffixTree CreateTree(string text, bool useMappedStorage)
        {
            if (!useMappedStorage)
                return PersistentSuffixTreeFactory.Create(new StringTextSource(text));

            string filePath = Path.Combine(_tempDir!, $"span_{Guid.NewGuid():N}.st");
            return PersistentSuffixTreeFactory.Create(new StringTextSource(text), filePath);
        }

        private static List<string> BuildPatternSet(string text, int seed)
        {
            var patterns = new List<string> { string.Empty, BuildMissingPattern(text) };

            if (text.Length == 0)
            {
                patterns.Add("a");
                return patterns;
            }

            patterns.Add(text);
            patterns.Add(text.Substring(0, 1));
            patterns.Add(text.Substring(text.Length - 1, 1));
            patterns.Add(text.Substring(0, Math.Min(8, text.Length)));
            patterns.Add(text.Substring(Math.Max(0, text.Length - Math.Min(8, text.Length))));

            var rng = new Random(seed);
            for (int i = 0; i < 20; i++)
            {
                int start = rng.Next(text.Length);
                int len = rng.Next(1, text.Length - start + 1);
                patterns.Add(text.Substring(start, len));
            }

            return patterns.Distinct(StringComparer.Ordinal).ToList();
        }

        private static List<int> BruteForcePositions(string text, string pattern)
        {
            var result = new List<int>();
            if (pattern.Length == 0)
            {
                result.AddRange(Enumerable.Range(0, text.Length));
                return result;
            }

            for (int i = 0; i <= text.Length - pattern.Length; i++)
            {
                if (text.AsSpan(i, pattern.Length).SequenceEqual(pattern.AsSpan()))
                    result.Add(i);
            }

            return result;
        }

        private static void AssertLcsContract(ISuffixTree tree, string text, string other, string context)
        {
            string lcsString = tree.LongestCommonSubstring(other);
            string lcsSpan = tree.LongestCommonSubstring(other.AsSpan());

            Assert.That(lcsSpan, Is.EqualTo(lcsString), $"{context}: span/string LCS mismatch");

            if (lcsString.Length > 0)
            {
                Assert.That(text.Contains(lcsString, StringComparison.Ordinal), Is.True,
                    $"{context}: LCS not in text");
                Assert.That(other.Contains(lcsString, StringComparison.Ordinal), Is.True,
                    $"{context}: LCS not in other");
            }

            int expectedLength = LongestCommonSubstringLength(text, other);
            Assert.That(lcsString.Length, Is.EqualTo(expectedLength),
                $"{context}: LCS length mismatch");
        }

        private static int LongestCommonSubstringLength(string a, string b)
        {
            if (a.Length == 0 || b.Length == 0)
                return 0;

            var dp = new int[b.Length + 1];
            int best = 0;
            for (int i = 1; i <= a.Length; i++)
            {
                for (int j = b.Length; j >= 1; j--)
                {
                    if (a[i - 1] == b[j - 1])
                    {
                        dp[j] = dp[j - 1] + 1;
                        if (dp[j] > best)
                            best = dp[j];
                    }
                    else
                    {
                        dp[j] = 0;
                    }
                }
            }

            return best;
        }

        private static string BuildMissingPattern(string text)
        {
            string[] candidates = { "\u0000", "\uFFFF", "§§§", text + "\u0000", text + "\uFFFF", "not-present" };
            foreach (var candidate in candidates)
            {
                if (!text.Contains(candidate, StringComparison.Ordinal))
                    return candidate;
            }

            return text + "|missing|";
        }

        private static string GenerateRandom(Random rng, int length, string alphabet)
        {
            var chars = new char[length];
            for (int i = 0; i < length; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            return new string(chars);
        }
    }
}
