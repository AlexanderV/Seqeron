using System;
using System.Linq;
using System.Text;
using NUnit.Framework;

namespace SuffixTree.Tests.Robustness
{
    /// <summary>
    /// Stress tests for performance and memory under load.
    /// </summary>
    [TestFixture]
    public class StressTests
    {
        #region Large Text Construction

        [Test]
        [Category("Stress")]
        public void Build_100KCharacters_Succeeds()
        {
            var text = GenerateRandomString(100_000);

            var st = SuffixTree.Build(text);

            Assert.That(st.Text.Length, Is.EqualTo(100_000));
        }

        [Test]
        [Category("Stress")]
        public void Build_HighlyRepetitive_Succeeds()
        {
            var text = new string('a', 100_000);

            var st = SuffixTree.Build(text);

            Assert.Multiple(() =>
            {
                Assert.That(st.Text.Length, Is.EqualTo(100_000));
                Assert.That(st.Contains(new string('a', 50_000)), Is.True);
            });
        }

        #endregion

        #region Many Queries

        [Test]
        [Category("Stress")]
        public void ManyContainsQueries_Succeeds()
        {
            var text = GenerateRandomString(10_000);
            var st = SuffixTree.Build(text);
            var random = new Random(12345);

            for (int i = 0; i < 10_000; i++)
            {
                var start = random.Next(text.Length);
                var len = random.Next(1, Math.Min(100, text.Length - start + 1));
                var pattern = text.Substring(start, len);

                Assert.That(st.Contains(pattern), Is.True);
            }
        }

        #endregion

        #region LCS Stress

        [Test]
        [Category("Stress")]
        public void LCS_LargeStrings_Succeeds()
        {
            var common = GenerateRandomString(1000);
            var text1 = GenerateRandomString(5000) + common + GenerateRandomString(5000);
            var text2 = GenerateRandomString(3000) + common + GenerateRandomString(3000);

            var st = SuffixTree.Build(text1);
            var lcs = st.LongestCommonSubstring(text2);

            Assert.That(lcs.Length, Is.GreaterThanOrEqualTo(1000));
        }

        #endregion

        #region LRS Stress

        [Test]
        [Category("Stress")]
        public void LRS_LargeRepetitiveText_Succeeds()
        {
            var pattern = "abcdefghij";
            var text = string.Concat(Enumerable.Repeat(pattern, 1000));

            var st = SuffixTree.Build(text);
            var lrs = st.LongestRepeatedSubstring();

            // The repeated pattern should be found
            Assert.That(lrs.Length, Is.GreaterThanOrEqualTo(pattern.Length));
        }

        #endregion

        #region GetAllSuffixes Stress

        [Test]
        [Category("Stress")]
        public void GetAllSuffixes_LargeText_Succeeds()
        {
            var text = GenerateRandomString(10_000);
            var st = SuffixTree.Build(text);

            var count = st.GetAllSuffixes().Count();

            Assert.That(count, Is.EqualTo(10_000));
        }

        #endregion

        #region Helper Methods

        private static string GenerateRandomString(int length, int seed = 42)
        {
            var random = new Random(seed);
            var sb = new StringBuilder(length);
            for (int i = 0; i < length; i++)
            {
                sb.Append((char)('a' + random.Next(26)));
            }
            return sb.ToString();
        }

        #endregion
    }
}
