using System;
using System.Diagnostics;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    [TestFixture]
    [Category("Algorithms")]
    [Category("Performance")]
    public class PerformanceGuardTests
    {
        // 6.0 is a safe upper bound against OS/GC noise.
        // O(N^2) regressions at N=100_000 would take multiple seconds instead of milliseconds.
        private const double ExpectedLinearScalingBound = 6.0;
        private const int BaseN = 100_000;
        
        private string GenerateRandomDna(int length, int seed = 42)
        {
            var random = new Random(seed);
            const string chars = "ACGT";
            return new string(Enumerable.Repeat(chars, length).Select(s => s[random.Next(s.Length)]).ToArray());
        }

        [Test]
        public void Build_ShouldBeLinear()
        {
            var textN = GenerateRandomDna(BaseN);
            var text2N = GenerateRandomDna(BaseN * 2);

            // Warmup
            _ = SuffixTree.Build(GenerateRandomDna(1000));

            var sw = Stopwatch.StartNew();
            _ = SuffixTree.Build(textN);
            sw.Stop();
            var timeN = sw.Elapsed.TotalMilliseconds;

            sw.Restart();
            _ = SuffixTree.Build(text2N);
            sw.Stop();
            var time2N = sw.Elapsed.TotalMilliseconds;

            // Enforce O(n) constraint manually: doubled input shouldn't exceed 2x time (+ some buffer for GC/system noise)
            Assert.That(time2N, Is.LessThanOrEqualTo(timeN * ExpectedLinearScalingBound),
                $"SuffixTree.Build took {time2N}ms for 2N, but {timeN}ms for N. Expected <= {timeN * ExpectedLinearScalingBound}ms");
        }

        [Test]
        public void Contains_ShouldBeLinearInPatternLength()
        {
            var text = GenerateRandomDna(BaseN * 2);
            var tree = SuffixTree.Build(text);
            
            int patternLen = 1000;
            var patternN = text.Substring(500, patternLen);
            var pattern2N = text.Substring(500, patternLen * 2);

            // Warmup
            _ = tree.Contains(text.Substring(0, 10));

            var sw = Stopwatch.StartNew();
            for (int i = 0; i < 10000; i++) tree.Contains(patternN);
            sw.Stop();
            var timeN = sw.Elapsed.TotalMilliseconds;

            sw.Restart();
            for (int i = 0; i < 10000; i++) tree.Contains(pattern2N);
            sw.Stop();
            var time2N = sw.Elapsed.TotalMilliseconds;

            Assert.That(time2N, Is.LessThanOrEqualTo(timeN * ExpectedLinearScalingBound),
                $"tree.Contains took {time2N}ms for 2M, but {timeN}ms for M. Expected <= {timeN * ExpectedLinearScalingBound}ms");
        }

        [Test]
        public void FindAll_ShouldBeLinearInPatternLengthPlusOccurrences()
        {
            // Create a highly repetitive text so that the occurrences count 'k' scales up
            // textN has length ~BaseN, text2N has length ~BaseN*2
            string unit = "ACGTACGTACGT";
            var textN = string.Concat(Enumerable.Repeat(unit, BaseN / unit.Length));
            var text2N = string.Concat(Enumerable.Repeat(unit, (BaseN * 2) / unit.Length));

            var treeN = SuffixTree.Build(textN);
            var tree2N = SuffixTree.Build(text2N);
            
            var pattern = "ACGT"; // finding this in repetitive text will yield Many occurrences. Occurrences double in tree2N.

            // Warmup
            _ = treeN.FindAllOccurrences("A");

            var sw = Stopwatch.StartNew();
            for (int i = 0; i < 500; i++) treeN.FindAllOccurrences(pattern);
            sw.Stop();
            var timeN = sw.Elapsed.TotalMilliseconds;

            sw.Restart();
            for (int i = 0; i < 500; i++) tree2N.FindAllOccurrences(pattern);
            sw.Stop();
            var time2N = sw.Elapsed.TotalMilliseconds;

            Assert.That(time2N, Is.LessThanOrEqualTo(timeN * ExpectedLinearScalingBound),
                $"tree.FindAllOccurrences took {time2N}ms for 2K occurrences, but {timeN}ms for K. Expected <= {timeN * ExpectedLinearScalingBound}ms");
        }

        [Test]
        public void LCS_ShouldBeLinear()
        {
            var text1_N = GenerateRandomDna(BaseN);
            var text2_N = GenerateRandomDna(BaseN);
            
            var text1_2N = GenerateRandomDna(BaseN * 2);
            var text2_2N = GenerateRandomDna(BaseN * 2);

            var treeN = SuffixTree.Build(text1_N);
            var tree2N = SuffixTree.Build(text1_2N);

            // Warmup
            _ = treeN.LongestCommonSubstring("ACGT");

            var sw = Stopwatch.StartNew();
            _ = treeN.LongestCommonSubstring(text2_N);
            sw.Stop();
            var timeN = sw.Elapsed.TotalMilliseconds;

            sw.Restart();
            _ = tree2N.LongestCommonSubstring(text2_2N);
            sw.Stop();
            var time2N = sw.Elapsed.TotalMilliseconds;

            // LCS scales with O(n + m), so doubling both inputs means T(2N) ~ 2 * T(N)
            Assert.That(time2N, Is.LessThanOrEqualTo(timeN * ExpectedLinearScalingBound),
                $"tree.LongestCommonSubstring took {time2N}ms for 2N+2M, but {timeN}ms for N+M. Expected <= {timeN * ExpectedLinearScalingBound}ms");
        }
    }
}

