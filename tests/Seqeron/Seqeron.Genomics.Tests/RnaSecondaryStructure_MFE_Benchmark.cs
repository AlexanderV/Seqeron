using System;
using System.Diagnostics;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
[Category("Benchmark")]
public class RnaSecondaryStructure_MFE_Benchmark
{
    /// <summary>
    /// Benchmark: Classic O(L³) vs Optimized MFE for various RNA sequence lengths.
    /// Optimizations: flat ArrayPool buffer, pre-indexed base positions,
    /// sliding lower-bound pointers, pair-type lookup table.
    /// </summary>
    [Test]
    [Explicit("Benchmark — run manually")]
    public void MFE_Benchmark_AllScenarios()
    {
        var scenarios = new (string Name, int Length, double GcContent)[]
        {
            ("L=100  GC=50%", 100, 0.5),
            ("L=200  GC=50%", 200, 0.5),
            ("L=300  GC=50%", 300, 0.5),
            ("L=500  GC=50%", 500, 0.5),
            ("L=500  GC=70%", 500, 0.7),
            ("L=750  GC=50%", 750, 0.5),
            ("L=1000 GC=50%", 1000, 0.5),
        };

        const int warmup = 1;
        const int runs = 3;

        Console.WriteLine("=== RNA MFE Benchmark: Classic vs Optimized ===\n");
        Console.WriteLine($"{"Scenario",-20} {"Classic ms",12} {"Optimized ms",14} {"Speedup",10} {"Match",6}");
        Console.WriteLine(new string('-', 70));

        foreach (var (name, length, gc) in scenarios)
        {
            string seq = GenerateReproducibleRna(length, gc, seed: length * 31 + (int)(gc * 100));

            // Warmup
            for (int w = 0; w < warmup; w++)
            {
                RnaSecondaryStructure.CalculateMinimumFreeEnergyClassic(seq);
                RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
            }

            // Benchmark classic
            var classicTimes = new double[runs];
            double classicResult = 0;
            for (int r = 0; r < runs; r++)
            {
                var sw = Stopwatch.StartNew();
                classicResult = RnaSecondaryStructure.CalculateMinimumFreeEnergyClassic(seq);
                sw.Stop();
                classicTimes[r] = sw.Elapsed.TotalMilliseconds;
            }

            // Benchmark optimized
            var optimizedTimes = new double[runs];
            double optimizedResult = 0;
            for (int r = 0; r < runs; r++)
            {
                var sw = Stopwatch.StartNew();
                optimizedResult = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
                sw.Stop();
                optimizedTimes[r] = sw.Elapsed.TotalMilliseconds;
            }

            double classicMedian = Median(classicTimes);
            double optimizedMedian = Median(optimizedTimes);
            double speedup = optimizedMedian > 0 ? classicMedian / optimizedMedian : 0;
            bool match = Math.Abs(classicResult - optimizedResult) < 0.01;

            string bar = new string('█', (int)Math.Min(speedup * 5, 50));
            Console.WriteLine($"{name,-20} {classicMedian,12:F2} {optimizedMedian,14:F2} {speedup,9:F2}× {(match ? "✓" : "✗"),5}");
            Console.WriteLine($"  {bar}");

            Assert.That(match, Is.True,
                $"MFE mismatch for {name}: classic={classicResult}, optimized={optimizedResult}");
        }
    }

    private static string GenerateReproducibleRna(int length, double gcContent, int seed)
    {
        var rng = new Random(seed);
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
        {
            double r = rng.NextDouble();
            if (r < gcContent / 2) sb.Append('G');
            else if (r < gcContent) sb.Append('C');
            else if (r < gcContent + (1 - gcContent) / 2) sb.Append('A');
            else sb.Append('U');
        }
        return sb.ToString();
    }

    private static double Median(double[] values)
    {
        var sorted = values.OrderBy(v => v).ToArray();
        int mid = sorted.Length / 2;
        return sorted.Length % 2 == 0 ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid];
    }
}
