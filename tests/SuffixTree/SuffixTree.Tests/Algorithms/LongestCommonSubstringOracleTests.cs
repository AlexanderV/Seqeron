using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    [TestFixture]
    [Category("Algorithms")]
    public class LongestCommonSubstringOracleTests
    {
        [Test]
        public void LcsInfoAndFindAll_AgreeWithBruteForce_OnRandomInputs()
        {
            for (int seed = 0; seed < 30; seed++)
            {
                var rng = new Random(seed);
                string text = GenerateRandom(rng, rng.Next(1, 90), "abcd");
                string other = GenerateRandom(rng, rng.Next(1, 90), "abcd");

                var tree = SuffixTree.Build(text);

                string lcs = tree.LongestCommonSubstring(other);
                string lcsSpan = tree.LongestCommonSubstring(other.AsSpan());
                var info = tree.LongestCommonSubstringInfo(other);
                var all = tree.FindAllLongestCommonSubstrings(other);

                int expectedLength = LongestCommonSubstringLength(text, other);

                Assert.Multiple(() =>
                {
                    Assert.That(lcsSpan, Is.EqualTo(lcs), $"Seed={seed}: span/string LCS mismatch");
                    Assert.That(lcs.Length, Is.EqualTo(expectedLength), $"Seed={seed}: LCS length mismatch");

                    Assert.That(info.Substring, Is.EqualTo(lcs), $"Seed={seed}: LCSInfo substring mismatch");
                    Assert.That(all.Substring, Is.EqualTo(lcs), $"Seed={seed}: FindAllLCS substring mismatch");
                });

                if (lcs.Length == 0)
                {
                    Assert.Multiple(() =>
                    {
                        Assert.That(info.PositionInText, Is.EqualTo(-1), $"Seed={seed}: info.PositionInText");
                        Assert.That(info.PositionInOther, Is.EqualTo(-1), $"Seed={seed}: info.PositionInOther");
                        Assert.That(all.PositionsInText, Is.Empty, $"Seed={seed}: all.PositionsInText");
                        Assert.That(all.PositionsInOther, Is.Empty, $"Seed={seed}: all.PositionsInOther");
                    });
                    continue;
                }

                Assert.That(info.PositionInText, Is.GreaterThanOrEqualTo(0), $"Seed={seed}: info.PositionInText");
                Assert.That(info.PositionInOther, Is.GreaterThanOrEqualTo(0), $"Seed={seed}: info.PositionInOther");
                Assert.That(info.PositionInText + lcs.Length, Is.LessThanOrEqualTo(text.Length), $"Seed={seed}: info text bounds");
                Assert.That(info.PositionInOther + lcs.Length, Is.LessThanOrEqualTo(other.Length), $"Seed={seed}: info other bounds");
                Assert.That(text.Substring(info.PositionInText, lcs.Length), Is.EqualTo(lcs), $"Seed={seed}: info text match");
                Assert.That(other.Substring(info.PositionInOther, lcs.Length), Is.EqualTo(lcs), $"Seed={seed}: info other match");

                var expectedTextPositions = BruteForcePositions(text, lcs).ToHashSet();
                var expectedOtherPositions = BruteForcePositions(other, lcs).ToHashSet();

                var actualTextPositions = all.PositionsInText.ToHashSet();
                var actualOtherPositions = all.PositionsInOther.ToHashSet();

                Assert.Multiple(() =>
                {
                    Assert.That(actualTextPositions.SetEquals(expectedTextPositions), Is.True,
                        $"Seed={seed}: PositionsInText mismatch for \"{lcs}\"");
                    Assert.That(actualOtherPositions.SetEquals(expectedOtherPositions), Is.True,
                        $"Seed={seed}: PositionsInOther mismatch for \"{lcs}\"");
                    Assert.That(actualTextPositions.Contains(info.PositionInText), Is.True,
                        $"Seed={seed}: info.PositionInText is not among FindAll positions");
                    Assert.That(actualOtherPositions.Contains(info.PositionInOther), Is.True,
                        $"Seed={seed}: info.PositionInOther is not among FindAll positions");
                });
            }
        }

        private static string GenerateRandom(Random rng, int length, string alphabet)
        {
            var chars = new char[length];
            for (int i = 0; i < length; i++)
            {
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            }

            return new string(chars);
        }

        private static List<int> BruteForcePositions(string text, string pattern)
        {
            var result = new List<int>();
            if (pattern.Length == 0)
            {
                return result;
            }

            for (int i = 0; i <= text.Length - pattern.Length; i++)
            {
                if (text.AsSpan(i, pattern.Length).SequenceEqual(pattern.AsSpan()))
                {
                    result.Add(i);
                }
            }

            return result;
        }

        private static int LongestCommonSubstringLength(string a, string b)
        {
            if (a.Length == 0 || b.Length == 0)
            {
                return 0;
            }

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
                        {
                            best = dp[j];
                        }
                    }
                    else
                    {
                        dp[j] = 0;
                    }
                }
            }

            return best;
        }
    }
}
