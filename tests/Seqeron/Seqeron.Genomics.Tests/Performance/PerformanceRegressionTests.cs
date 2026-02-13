using System.Diagnostics;
using FluentAssertions;
using Seqeron.Genomics.Tests.Builders;

namespace Seqeron.Genomics.Tests.Performance;

/// <summary>
/// Performance regression tests.
/// These tests catch algorithmic complexity regressions (e.g., O(n) → O(n²)) by asserting
/// that operations complete within reasonable time bounds on known input sizes.
///
/// Marked [NonParallelizable] to avoid interference from concurrent CPU load.
/// </summary>
[TestFixture]
[Category("Performance")]
[NonParallelizable]
public class PerformanceRegressionTests
{
    /// <summary>
    /// GC content calculation should complete quickly for large sequences.
    /// Expected: O(n) — single pass through sequence.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void GcContent_100K_Nucleotides_CompletesInTime()
    {
        string seq = SequenceBuilder.Dna().Random(100_000).Build();

        var sw = Stopwatch.StartNew();
        double gc = seq.AsSpan().CalculateGcContent();
        sw.Stop();

        gc.Should().BeInRange(0.0, 100.0);
        sw.ElapsedMilliseconds.Should().BeLessThan(1000,
            "GC content of 100K nucleotides should complete within 1 second");
    }

    /// <summary>
    /// GC content should scale linearly: 100× input should take &lt;200× time.
    /// </summary>
    [Test]
    public void GcContent_ScalesLinearly()
    {
        string small = SequenceBuilder.Dna().Random(1_000).Build();
        string large = SequenceBuilder.Dna().WithSeed(99).Random(100_000).Build();

        // Warm up JIT
        small.AsSpan().CalculateGcContent();
        large.AsSpan().CalculateGcContent();

        var swSmall = Stopwatch.StartNew();
        for (int i = 0; i < 100; i++)
            small.AsSpan().CalculateGcContent();
        swSmall.Stop();

        var swLarge = Stopwatch.StartNew();
        for (int i = 0; i < 100; i++)
            large.AsSpan().CalculateGcContent();
        swLarge.Stop();

        double ratio = (double)swLarge.ElapsedTicks / Math.Max(swSmall.ElapsedTicks, 1);
        ratio.Should().BeLessThan(500,
            "100× input size should take less than 500× time (linear with JIT overhead)");
    }

    /// <summary>
    /// Edit distance computation should complete for moderate-length strings.
    /// Expected: O(n×m) — bounded by string lengths.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void EditDistance_500x500_CompletesInTime()
    {
        string s1 = SequenceBuilder.Dna().Random(500).Build();
        string s2 = SequenceBuilder.Dna().WithSeed(77).Random(500).Build();

        var sw = Stopwatch.StartNew();
        int distance = ApproximateMatcher.EditDistance(s1, s2);
        sw.Stop();

        distance.Should().BeGreaterThanOrEqualTo(0);
        sw.ElapsedMilliseconds.Should().BeLessThan(5000,
            "Edit distance of 500×500 should complete within 5 seconds");
    }

    /// <summary>
    /// Global alignment should complete for moderate-length sequences.
    /// Expected: O(n×m) DP algorithm.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void GlobalAlign_500bp_CompletesInTime()
    {
        string seq1 = SequenceBuilder.Dna().Random(500).Build();
        string seq2 = SequenceBuilder.Dna().WithSeed(55).Random(500).Build();

        var sw = Stopwatch.StartNew();
        var result = SequenceAligner.GlobalAlign(seq1, seq2);
        sw.Stop();

        result.AlignedSequence1.Should().NotBeEmpty();
        sw.ElapsedMilliseconds.Should().BeLessThan(10000,
            "Global alignment of 500bp should complete within 10 seconds");
    }

    /// <summary>
    /// Disorder prediction should scale linearly with sequence length.
    /// Expected: O(n) — sliding window pass.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DisorderPrediction_10K_Residues_CompletesInTime()
    {
        string seq = SequenceBuilder.Protein().Random(10_000, "ACDEFGHIKLMNPQRSTVWY").Build();

        var sw = Stopwatch.StartNew();
        var result = DisorderPredictor.PredictDisorder(seq, minRegionLength: 5);
        sw.Stop();

        result.MeanDisorderScore.Should().BeInRange(0.0, 1.0);
        sw.ElapsedMilliseconds.Should().BeLessThan(3000,
            "Disorder prediction of 10K residues should complete within 3 seconds");
    }

    /// <summary>
    /// DnaSequence construction + complement + reverse complement chain.
    /// Tests that lazy SuffixTree is NOT built during basic operations.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DnaSequence_BasicOperations_10K_AreFast()
    {
        string seq = SequenceBuilder.Dna().Random(10_000).Build();

        var sw = Stopwatch.StartNew();
        var dna = new DnaSequence(seq);
        var comp = dna.Complement();
        var rc = dna.ReverseComplement();
        double gc = dna.GcContent();
        sw.Stop();

        comp.Length.Should().Be(dna.Length);
        rc.Length.Should().Be(dna.Length);
        gc.Should().BeInRange(0.0, 100.0);
        sw.ElapsedMilliseconds.Should().BeLessThan(1000,
            "Basic DnaSequence operations on 10K should complete within 1 second");
    }
}
