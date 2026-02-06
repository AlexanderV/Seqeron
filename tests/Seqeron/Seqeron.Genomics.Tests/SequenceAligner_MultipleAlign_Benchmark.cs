using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Performance benchmark tests for Multiple Sequence Alignment.
/// Compares classic NW-based star alignment vs anchor-based suffix tree approach.
/// 
/// Run with: dotnet test --filter "Category=MSA-Benchmark" -v n
/// 
/// These tests output timing data to TestContext so you can track improvements
/// after each optimization step.
/// </summary>
[TestFixture]
[Category("MSA-Benchmark")]
[Category("Performance")]
public class SequenceAligner_MultipleAlign_Benchmark
{
    private static readonly ScoringMatrix Scoring = SequenceAligner.SimpleDna;

    #region Data Generators

    /// <summary>
    /// Generates a random DNA sequence with a fixed seed for reproducibility.
    /// </summary>
    private static string GenerateDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(4)];
        return new string(chars);
    }

    /// <summary>
    /// Generates sequences with controlled similarity to a reference.
    /// Introduces random point mutations at rate (1 - similarity).
    /// This simulates real biological sequences where closely related organisms
    /// share most of their DNA.
    /// </summary>
    private static List<DnaSequence> GenerateRelatedSequences(
        int count, int length, double similarity, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);

        // Generate reference
        string reference = GenerateDna(length, seed);
        var sequences = new List<DnaSequence> { new DnaSequence(reference) };

        for (int s = 1; s < count; s++)
        {
            var chars = reference.ToCharArray();
            int mutations = (int)(length * (1.0 - similarity));

            // Apply random point mutations
            for (int m = 0; m < mutations; m++)
            {
                int pos = rng.Next(length);
                char original = chars[pos];
                char replacement;
                do { replacement = bases[rng.Next(4)]; }
                while (replacement == original);
                chars[pos] = replacement;
            }

            // Occasionally add small indels (~1% of length)
            var result = new string(chars);
            int indelCount = Math.Max(1, length / 100);
            var sb = new System.Text.StringBuilder(result);

            for (int d = 0; d < indelCount && sb.Length > 10; d++)
            {
                if (rng.NextDouble() < 0.5)
                {
                    // Insertion
                    int pos = rng.Next(sb.Length);
                    sb.Insert(pos, bases[rng.Next(4)]);
                }
                else
                {
                    // Deletion
                    int pos = rng.Next(sb.Length);
                    sb.Remove(pos, 1);
                }
            }

            sequences.Add(new DnaSequence(sb.ToString()));
        }

        return sequences;
    }

    /// <summary>
    /// Generates completely random (unrelated) sequences.
    /// Worst case for anchor-based approach — few anchors expected.
    /// </summary>
    private static List<DnaSequence> GenerateRandomSequences(
        int count, int length, int seed)
    {
        return Enumerable.Range(0, count)
            .Select(i => new DnaSequence(GenerateDna(length, seed + i * 1000)))
            .ToList();
    }

    #endregion

    #region Warmup

    [OneTimeSetUp]
    public void WarmUp()
    {
        // JIT warmup — run a tiny alignment so first real test isn't skewed
        var tiny = new[] { new DnaSequence("ATGC"), new DnaSequence("ATGC") };
        SequenceAligner.MultipleAlign(tiny);
        SequenceAligner.MultipleAlignClassic(tiny.ToList(), Scoring);
    }

    #endregion

    #region Core Benchmark Method

    /// <summary>
    /// Runs both algorithms and reports timing comparison.
    /// Returns (classicMs, anchorMs, speedup).
    /// </summary>
    private (double ClassicMs, double AnchorMs, double Speedup) RunBenchmark(
        List<DnaSequence> sequences, string label, int warmupRuns = 1, int measuredRuns = 3)
    {
        var seqList = sequences.ToList();

        // Warmup
        for (int i = 0; i < warmupRuns; i++)
        {
            SequenceAligner.MultipleAlignClassic(seqList, Scoring);
            SequenceAligner.MultipleAlign(sequences, Scoring);
        }

        // Force GC before measurement
        GC.Collect();
        GC.WaitForPendingFinalizers();
        GC.Collect();

        // Measure Classic (NW-only)
        var classicTimes = new List<double>();
        for (int i = 0; i < measuredRuns; i++)
        {
            var sw = Stopwatch.StartNew();
            var classicResult = SequenceAligner.MultipleAlignClassic(seqList, Scoring);
            sw.Stop();
            classicTimes.Add(sw.Elapsed.TotalMilliseconds);
        }

        GC.Collect();
        GC.WaitForPendingFinalizers();
        GC.Collect();

        // Measure Anchor-based
        var anchorTimes = new List<double>();
        for (int i = 0; i < measuredRuns; i++)
        {
            var sw = Stopwatch.StartNew();
            var anchorResult = SequenceAligner.MultipleAlign(sequences, Scoring);
            sw.Stop();
            anchorTimes.Add(sw.Elapsed.TotalMilliseconds);
        }

        // Use median to reduce outlier impact
        double classicMs = Median(classicTimes);
        double anchorMs = Median(anchorTimes);
        double speedup = classicMs / Math.Max(anchorMs, 0.001);

        // Output results
        TestContext.WriteLine($"╔══════════════════════════════════════════════════════════");
        TestContext.WriteLine($"║ {label}");
        TestContext.WriteLine($"╠══════════════════════════════════════════════════════════");
        TestContext.WriteLine($"║ Sequences: {sequences.Count} × ~{sequences[0].Length} bp");
        TestContext.WriteLine($"║ Classic NW:     {classicMs,10:F2} ms  (median of {measuredRuns})");
        TestContext.WriteLine($"║ Anchor-based:   {anchorMs,10:F2} ms  (median of {measuredRuns})");
        TestContext.WriteLine($"║ Speedup:        {speedup,10:F2}×");
        TestContext.WriteLine($"╚══════════════════════════════════════════════════════════");

        return (classicMs, anchorMs, speedup);
    }

    private static double Median(List<double> values)
    {
        var sorted = values.OrderBy(v => v).ToList();
        int n = sorted.Count;
        return n % 2 == 0
            ? (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
            : sorted[n / 2];
    }

    #endregion

    #region Benchmark Scenarios — Related Sequences (realistic biological case)

    /// <summary>
    /// Small: 5 sequences × 200 bp, 90% similarity.
    /// Baseline scenario — should already show some improvement.
    /// </summary>
    [Test]
    [Order(1)]
    public void Benchmark_Small_Related_5x200_90pct()
    {
        var seqs = GenerateRelatedSequences(count: 5, length: 200, similarity: 0.90, seed: 42);
        var (classic, anchor, speedup) = RunBenchmark(seqs, "Small: 5 × 200bp, 90% similarity");

        // Correctness check: both produce valid MSA
        var classicResult = SequenceAligner.MultipleAlignClassic(seqs.ToList(), Scoring);
        var anchorResult = SequenceAligner.MultipleAlign(seqs, Scoring);

        Assert.Multiple(() =>
        {
            Assert.That(classicResult.AlignedSequences.Length, Is.EqualTo(5));
            Assert.That(anchorResult.AlignedSequences.Length, Is.EqualTo(5));

            // Both should produce equal-length aligned sequences
            Assert.That(classicResult.AlignedSequences.Select(s => s.Length).Distinct().Count(), Is.EqualTo(1));
            Assert.That(anchorResult.AlignedSequences.Select(s => s.Length).Distinct().Count(), Is.EqualTo(1));

            // Removing gaps recovers originals
            for (int i = 0; i < 5; i++)
            {
                Assert.That(anchorResult.AlignedSequences[i].Replace("-", ""),
                    Is.EqualTo(seqs[i].Sequence),
                    $"Anchor: sequence {i} gap removal should recover original");
            }
        });
    }

    /// <summary>
    /// Medium: 10 sequences × 500 bp, 85% similarity.
    /// Core benchmark — anchor approach should shine here.
    /// </summary>
    [Test]
    [Order(2)]
    public void Benchmark_Medium_Related_10x500_85pct()
    {
        var seqs = GenerateRelatedSequences(count: 10, length: 500, similarity: 0.85, seed: 123);
        var (classic, anchor, speedup) = RunBenchmark(seqs, "Medium: 10 × 500bp, 85% similarity");

        // Verify correctness
        var result = SequenceAligner.MultipleAlign(seqs, Scoring);
        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(10));
            for (int i = 0; i < 10; i++)
            {
                Assert.That(result.AlignedSequences[i].Replace("-", ""),
                    Is.EqualTo(seqs[i].Sequence),
                    $"Sequence {i} gap removal should recover original");
            }
        });
    }

    /// <summary>
    /// Large: 10 sequences × 1000 bp, 90% similarity.
    /// This is where anchor-based should show significant improvement.
    /// </summary>
    [Test]
    [Order(3)]
    public void Benchmark_Large_Related_10x1000_90pct()
    {
        var seqs = GenerateRelatedSequences(count: 10, length: 1000, similarity: 0.90, seed: 456);
        var (classic, anchor, speedup) = RunBenchmark(seqs, "Large: 10 × 1000bp, 90% similarity");

        var result = SequenceAligner.MultipleAlign(seqs, Scoring);
        Assert.That(result.AlignedSequences.Length, Is.EqualTo(10));

        for (int i = 0; i < 10; i++)
        {
            Assert.That(result.AlignedSequences[i].Replace("-", ""),
                Is.EqualTo(seqs[i].Sequence),
                $"Sequence {i} gap removal should recover original");
        }
    }

    /// <summary>
    /// XL: 20 sequences × 2000 bp, 90% similarity.
    /// Stress test — large sequences and many of them.
    /// </summary>
    [Test]
    [Order(4)]
    [Timeout(60_000)] // 60 second timeout
    public void Benchmark_XL_Related_20x2000_90pct()
    {
        var seqs = GenerateRelatedSequences(count: 20, length: 2000, similarity: 0.90, seed: 789);
        var (classic, anchor, speedup) = RunBenchmark(seqs,
            "XL: 20 × 2000bp, 90% similarity", warmupRuns: 0, measuredRuns: 1);

        var result = SequenceAligner.MultipleAlign(seqs, Scoring);
        Assert.That(result.AlignedSequences.Length, Is.EqualTo(20));
    }

    #endregion

    #region Benchmark Scenarios — Varying Similarity

    /// <summary>
    /// High similarity (95%) — lots of anchors expected, maximum speedup.
    /// </summary>
    [Test]
    [Order(5)]
    public void Benchmark_HighSimilarity_10x500_95pct()
    {
        var seqs = GenerateRelatedSequences(count: 10, length: 500, similarity: 0.95, seed: 111);
        RunBenchmark(seqs, "High similarity: 10 × 500bp, 95%");
    }

    /// <summary>
    /// Medium similarity (80%) — fewer anchors, smaller speedup expected.
    /// </summary>
    [Test]
    [Order(6)]
    public void Benchmark_MediumSimilarity_10x500_80pct()
    {
        var seqs = GenerateRelatedSequences(count: 10, length: 500, similarity: 0.80, seed: 222);
        RunBenchmark(seqs, "Medium similarity: 10 × 500bp, 80%");
    }

    /// <summary>
    /// Low similarity (60%) — few anchors, fallback to NW expected.
    /// Anchor approach may not help much here.
    /// </summary>
    [Test]
    [Order(7)]
    public void Benchmark_LowSimilarity_10x500_60pct()
    {
        var seqs = GenerateRelatedSequences(count: 10, length: 500, similarity: 0.60, seed: 333);
        RunBenchmark(seqs, "Low similarity: 10 × 500bp, 60%");
    }

    #endregion

    #region Benchmark Scenarios — Random Sequences (worst case)

    /// <summary>
    /// Random sequences — no biological similarity.
    /// Anchor approach should gracefully fall back to NW.
    /// Verifies we don't regress on worst-case input.
    /// </summary>
    [Test]
    [Order(8)]
    public void Benchmark_Random_10x500_NoSimilarity()
    {
        var seqs = GenerateRandomSequences(count: 10, length: 500, seed: 999);
        var (classic, anchor, speedup) = RunBenchmark(seqs, "Random (worst case): 10 × 500bp");

        // For random sequences, anchor approach should not be significantly slower
        // (graceful fallback). Allow up to 3× slower due to overhead.
        Assert.That(speedup, Is.GreaterThan(0.3),
            $"Anchor approach should not be >3× slower than classic on random input. " +
            $"Got speedup={speedup:F2}×");
    }

    #endregion

    #region Summary Report

    /// <summary>
    /// Runs all sizes and prints a summary comparison table.
    /// This is the main test to run after each optimization step.
    /// </summary>
    [Test]
    [Order(100)]
    [Timeout(120_000)] // 2 minute timeout
    public void Benchmark_SummaryReport()
    {
        var scenarios = new[]
        {
            ("5×100,  95%", GenerateRelatedSequences(5,   100, 0.95, 10)),
            ("5×200,  90%", GenerateRelatedSequences(5,   200, 0.90, 20)),
            ("10×500, 90%", GenerateRelatedSequences(10,  500, 0.90, 30)),
            ("10×500, 80%", GenerateRelatedSequences(10,  500, 0.80, 40)),
            ("10×1K,  90%", GenerateRelatedSequences(10, 1000, 0.90, 50)),
            ("20×1K,  90%", GenerateRelatedSequences(20, 1000, 0.90, 60)),
            ("10×500, rnd", GenerateRandomSequences(10,   500, 70)),
        };

        TestContext.WriteLine();
        TestContext.WriteLine("┌──────────────────┬────────────┬────────────┬──────────┐");
        TestContext.WriteLine("│ Scenario         │ Classic ms │ Anchor  ms │ Speedup  │");
        TestContext.WriteLine("├──────────────────┼────────────┼────────────┼──────────┤");

        foreach (var (label, seqs) in scenarios)
        {
            // Single warmup
            SequenceAligner.MultipleAlignClassic(seqs.ToList(), Scoring);
            SequenceAligner.MultipleAlign(seqs, Scoring);

            GC.Collect(); GC.WaitForPendingFinalizers(); GC.Collect();

            // Measure classic
            var sw = Stopwatch.StartNew();
            SequenceAligner.MultipleAlignClassic(seqs.ToList(), Scoring);
            sw.Stop();
            double classicMs = sw.Elapsed.TotalMilliseconds;

            GC.Collect(); GC.WaitForPendingFinalizers(); GC.Collect();

            // Measure anchor
            sw.Restart();
            SequenceAligner.MultipleAlign(seqs, Scoring);
            sw.Stop();
            double anchorMs = sw.Elapsed.TotalMilliseconds;

            double speedup = classicMs / Math.Max(anchorMs, 0.001);

            string bar = speedup >= 1.0
                ? new string('█', Math.Min((int)(speedup * 2), 20))
                : "▒";

            TestContext.WriteLine(
                $"│ {label,-16} │ {classicMs,10:F2} │ {anchorMs,10:F2} │ {speedup,6:F2}× {bar}│");
        }

        TestContext.WriteLine("└──────────────────┴────────────┴────────────┴──────────┘");
        TestContext.WriteLine();
        TestContext.WriteLine($"Timestamp: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
        TestContext.WriteLine($"Runtime: {System.Runtime.InteropServices.RuntimeInformation.FrameworkDescription}");
        TestContext.WriteLine($"OS: {System.Runtime.InteropServices.RuntimeInformation.OSDescription}");

        Assert.Pass("Benchmark summary completed — see output above.");
    }

    #endregion
}
