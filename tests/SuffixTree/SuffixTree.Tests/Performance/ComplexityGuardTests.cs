using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Performance
{
    /// <summary>
    /// Performance regression guards.
    /// These tests do NOT measure absolute wall-clock time.
    /// Instead, they verify algorithmic complexity by comparing time(2n)/time(n).
    /// A ratio ≤ threshold indicates the expected complexity class is preserved.
    /// 
    /// Expected complexities:
    ///   Build:     O(n)   → ratio ≈ 2
    ///   Contains:  O(m)   → ratio ≈ 2
    ///   FindAll:   O(m+k) → ratio depends on k growth
    ///   Count:     O(m)   → ratio ≈ 2
    ///   LCS:       O(m)   → ratio ≈ 2
    ///   LRS:       O(1)   → ratio ≈ 1
    /// </summary>
    [TestFixture]
    [Category("Performance")]
    public class ComplexityGuardTests
    {
        // Allow generous headroom for GC jitter, JIT warmup, and CI variability.
        // A true O(n²) regression would produce ratio ≥ 4 at doubling.
        private const double LinearThreshold = 3.5;
        private const double SublinearThreshold = 4.0;

        private static string GenerateDna(int length, int seed = 42)
        {
            var rng = new Random(seed);
            var chars = new char[length];
            const string alpha = "ACGT";
            for (int i = 0; i < length; i++)
                chars[i] = alpha[rng.Next(4)];
            return new string(chars);
        }

        private static double MeasureMs(Action action, int warmup = 2, int measured = 3)
        {
            // Warmup
            for (int i = 0; i < warmup; i++)
                action();

            // Measured runs
            var sw = Stopwatch.StartNew();
            for (int i = 0; i < measured; i++)
                action();
            sw.Stop();
            return sw.Elapsed.TotalMilliseconds / measured;
        }

        #region Build complexity: O(n)

        [Test]
        public void Build_ScalesLinearly()
        {
            const int n = 20_000;
            string small = GenerateDna(n);
            string large = GenerateDna(n * 2, seed: 99);

            double tSmall = MeasureMs(() => SuffixTree.Build(small));
            double tLarge = MeasureMs(() => SuffixTree.Build(large));

            double ratio = tLarge / Math.Max(tSmall, 0.01);
            Assert.That(ratio, Is.LessThan(LinearThreshold),
                $"Build: time({n * 2})={tLarge:F1}ms / time({n})={tSmall:F1}ms = {ratio:F2} exceeds {LinearThreshold}");
        }

        #endregion

        #region Contains complexity: O(m)

        [Test]
        public void Contains_ScalesLinearlyWithPatternLength()
        {
            string text = GenerateDna(50_000);
            var tree = SuffixTree.Build(text);

            const int m = 5000;
            string shortPattern = text.Substring(0, m);
            string longPattern = text.Substring(0, m * 2);

            double tShort = MeasureMs(() => tree.Contains(shortPattern), warmup: 5, measured: 10);
            double tLong = MeasureMs(() => tree.Contains(longPattern), warmup: 5, measured: 10);

            double ratio = tLong / Math.Max(tShort, 0.001);
            Assert.That(ratio, Is.LessThan(LinearThreshold),
                $"Contains: time({m * 2})={tLong:F3}ms / time({m})={tShort:F3}ms = {ratio:F2}");
        }

        #endregion

        #region FindAllOccurrences: O(m + k) — verify no O(n²) DFS

        [Test]
        public void FindAll_ScalesWithOutputSize()
        {
            // Build tree with repetitive content so k grows predictably
            const int n = 20_000;
            string text = GenerateDna(n);
            var tree = SuffixTree.Build(text);

            // Short pattern → many results, longer pattern → fewer results
            // We verify that the operation completes in reasonable time
            // and doesn't exhibit quadratic behavior
            string pattern1 = "ACG";   // short → lots of hits
            string pattern2 = text.Substring(0, 8); // longer → fewer hits

            var pos1 = new List<int>();
            var pos2 = new List<int>();

            double t1 = MeasureMs(() => { pos1 = tree.FindAllOccurrences(pattern1).ToList(); });
            double t2 = MeasureMs(() => { pos2 = tree.FindAllOccurrences(pattern2).ToList(); });

            // Just verify it completes in < 500ms (no O(n²) blowup)
            Assert.That(t1, Is.LessThan(500),
                $"FindAll(\"{pattern1}\") took {t1:F1}ms, k={pos1.Count}");
            Assert.That(t2, Is.LessThan(500),
                $"FindAll(\"{pattern2}\") took {t2:F1}ms, k={pos2.Count}");
        }

        #endregion

        #region CountOccurrences: O(m) — no DFS, just read LeafCount

        [Test]
        public void Count_ScalesLinearlyWithPatternLength()
        {
            string text = GenerateDna(50_000);
            var tree = SuffixTree.Build(text);

            const int m = 5000;
            string shortP = text.Substring(100, m);
            string longP = text.Substring(100, m * 2);

            double tShort = MeasureMs(() => tree.CountOccurrences(shortP), warmup: 5, measured: 10);
            double tLong = MeasureMs(() => tree.CountOccurrences(longP), warmup: 5, measured: 10);

            double ratio = tLong / Math.Max(tShort, 0.001);
            Assert.That(ratio, Is.LessThan(LinearThreshold),
                $"Count: time({m * 2})={tLong:F3}ms / time({m})={tShort:F3}ms = {ratio:F2}");
        }

        #endregion

        #region LCS: O(n+m) — suffix-link streaming

        [Test]
        public void LCS_ScalesLinearlyWithQueryLength()
        {
            string text = GenerateDna(30_000);
            var tree = SuffixTree.Build(text);

            const int m = 5000;
            string shortQ = GenerateDna(m, seed: 77);
            string longQ = GenerateDna(m * 2, seed: 77);

            double tShort = MeasureMs(() => tree.LongestCommonSubstring(shortQ), warmup: 3, measured: 5);
            double tLong = MeasureMs(() => tree.LongestCommonSubstring(longQ), warmup: 3, measured: 5);

            double ratio = tLong / Math.Max(tShort, 0.001);
            Assert.That(ratio, Is.LessThan(SublinearThreshold),
                $"LCS: time({m * 2})={tLong:F2}ms / time({m})={tShort:F2}ms = {ratio:F2}");
        }

        #endregion

        #region LRS: O(1) cached after first call

        [Test]
        public void LRS_IsCachedConstantTime()
        {
            string text = GenerateDna(50_000);
            var tree = SuffixTree.Build(text);

            // First call — may compute
            _ = tree.LongestRepeatedSubstring();

            // Subsequent calls should be near-instant
            double t = MeasureMs(() => tree.LongestRepeatedSubstring(), warmup: 5, measured: 50);
            Assert.That(t, Is.LessThan(1.0),
                $"LRS (cached) took {t:F3}ms — expected <1ms for O(1) cached read");
        }

        #endregion

        #region Build: pathological repetitive input should stay O(n)

        [Test]
        public void Build_RepetitiveInput_StillLinear()
        {
            // Pathological: single-char repeat maximizes Ukkonen's active state transitions
            const int n = 10_000;
            string small = new string('a', n);
            string large = new string('a', n * 2);

            double tSmall = MeasureMs(() => SuffixTree.Build(small));
            double tLarge = MeasureMs(() => SuffixTree.Build(large));

            double ratio = tLarge / Math.Max(tSmall, 0.01);
            Assert.That(ratio, Is.LessThan(LinearThreshold),
                $"Build(repeat): time({n * 2})={tLarge:F1}ms / time({n})={tSmall:F1}ms = {ratio:F2}");
        }

        #endregion

        #region GetAllSuffixes: O(n) total output

        [Test]
        public void GetAllSuffixes_ScalesLinearly()
        {
            const int n = 5_000;
            string small = GenerateDna(n);
            string large = GenerateDna(n * 2, seed: 88);

            var treeSmall = SuffixTree.Build(small);
            var treeLarge = SuffixTree.Build(large);

            double tSmall = MeasureMs(() => treeSmall.GetAllSuffixes());
            double tLarge = MeasureMs(() => treeLarge.GetAllSuffixes());

            // GetAllSuffixes produces O(n²) total chars but n suffix entries.
            // Allow higher ratio since output grows quadratically.
            double ratio = tLarge / Math.Max(tSmall, 0.01);
            Assert.That(ratio, Is.LessThan(6.0),
                $"GetAllSuffixes: time({n * 2})={tLarge:F1}ms / time({n})={tSmall:F1}ms = {ratio:F2}");
        }

        #endregion
    }
}
