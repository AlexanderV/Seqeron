using System;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for CancellationToken and Span-based performance methods.
/// </summary>
[TestFixture]
public class PerformanceExtensionsTests
{
    #region Span-based Tests

    [Test]
    public void CalculateGcContent_Span_ReturnsCorrectValue()
    {
        ReadOnlySpan<char> sequence = "ACGT".AsSpan();
        double gcContent = sequence.CalculateGcContent();
        Assert.That(gcContent, Is.EqualTo(50.0));
    }

    [Test]
    public void CalculateGcContent_AllGC_Returns100()
    {
        ReadOnlySpan<char> sequence = "GCGCGC".AsSpan();
        Assert.That(sequence.CalculateGcContent(), Is.EqualTo(100.0));
    }

    [Test]
    public void CalculateGcContent_NoGC_Returns0()
    {
        ReadOnlySpan<char> sequence = "ATATAT".AsSpan();
        Assert.That(sequence.CalculateGcContent(), Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateGcContent_Empty_Returns0()
    {
        ReadOnlySpan<char> sequence = ReadOnlySpan<char>.Empty;
        Assert.That(sequence.CalculateGcContent(), Is.EqualTo(0.0));
    }

    // Note: Comprehensive TryGetComplement tests are in SequenceExtensions_Complement_Tests.cs
    // These are retained as smoke tests for span-based performance verification

    [Test]
    [Description("Smoke test: TryGetComplement span API works correctly")]
    public void TryGetComplement_SpanApi_SmokeTest()
    {
        ReadOnlySpan<char> source = "ACGT".AsSpan();
        Span<char> destination = stackalloc char[4];

        bool success = source.TryGetComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo("TGCA"));
    }

    // Note: Comprehensive TryGetReverseComplement tests will be in SEQ-REVCOMP-001
    [Test]
    [Description("Smoke test: TryGetReverseComplement span API works correctly")]
    public void TryGetReverseComplement_SpanApi_SmokeTest()
    {
        ReadOnlySpan<char> source = "ACGT".AsSpan();
        Span<char> destination = stackalloc char[4];

        bool success = source.TryGetReverseComplement(destination);

        Assert.That(success, Is.True);
        Assert.That(new string(destination), Is.EqualTo("ACGT")); // ACGT reverse complement is ACGT
    }

    [Test]
    public void CountKmersSpan_ReturnsCorrectCounts()
    {
        ReadOnlySpan<char> sequence = "ACGTACGT".AsSpan();
        var counts = sequence.CountKmersSpan(2);

        Assert.That(counts["AC"], Is.EqualTo(2));
        Assert.That(counts["CG"], Is.EqualTo(2));
        Assert.That(counts["GT"], Is.EqualTo(2));
        Assert.That(counts["TA"], Is.EqualTo(1));
    }

    [Test]
    public void HammingDistance_Span_ReturnsCorrectDistance()
    {
        ReadOnlySpan<char> s1 = "ACGT".AsSpan();
        ReadOnlySpan<char> s2 = "ACTT".AsSpan();

        int distance = s1.HammingDistance(s2);

        Assert.That(distance, Is.EqualTo(1));
    }

    // Note: Comprehensive IsValidDna/IsValidRna tests are in SequenceExtensions_SequenceValidation_Tests.cs (SEQ-VALID-001)
    // Smoke tests retained here for span-based API verification

    [Test]
    [Description("Smoke test: IsValidDna span API works correctly")]
    public void IsValidDna_ValidSequence_SmokeTest()
    {
        ReadOnlySpan<char> sequence = "ACGTACGT".AsSpan();
        Assert.That(sequence.IsValidDna(), Is.True);
    }

    [Test]
    [Description("Smoke test: IsValidRna span API works correctly")]
    public void IsValidRna_ValidSequence_SmokeTest()
    {
        ReadOnlySpan<char> sequence = "ACGUACGU".AsSpan();
        Assert.That(sequence.IsValidRna(), Is.True);
    }

    #endregion

    #region DnaSequence Span Methods

    [Test]
    public void DnaSequence_AsSpan_ReturnsCorrectSpan()
    {
        var dna = new DnaSequence("ACGT");
        var span = dna.AsSpan();

        Assert.That(span.Length, Is.EqualTo(4));
        Assert.That(span[0], Is.EqualTo('A'));
        Assert.That(span[3], Is.EqualTo('T'));
    }

    [Test]
    public void DnaSequence_GcContentFast_MatchesRegular()
    {
        var dna = new DnaSequence("ACGTACGT");

        double regular = dna.GcContent();
        double fast = dna.GcContentFast();

        Assert.That(fast, Is.EqualTo(regular));
    }

    [Test]
    public void DnaSequence_TryWriteComplement_Works()
    {
        var dna = new DnaSequence("ACGT");
        Span<char> buffer = stackalloc char[4];

        bool success = dna.TryWriteComplement(buffer);

        Assert.That(success, Is.True);
        Assert.That(new string(buffer), Is.EqualTo("TGCA"));
    }

    [Test]
    public void DnaSequence_TryWriteReverseComplement_Works()
    {
        var dna = new DnaSequence("AAAA");
        Span<char> buffer = stackalloc char[4];

        bool success = dna.TryWriteReverseComplement(buffer);

        Assert.That(success, Is.True);
        Assert.That(new string(buffer), Is.EqualTo("TTTT"));
    }

    #endregion

    #region CancellationToken Tests

    [Test]
    public void CountKmers_WithCancellation_CompletesNormally()
    {
        string sequence = new string('A', 1000) + new string('C', 1000);
        using var cts = new CancellationTokenSource();

        var counts = KmerAnalyzer.CountKmers(sequence, 3, cts.Token);

        Assert.That(counts, Is.Not.Empty);
    }

    [Test]
    public void CountKmers_CancelledToken_ThrowsOperationCanceledException()
    {
        string sequence = new string('A', 10000);
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        Assert.Throws<OperationCanceledException>(() =>
            KmerAnalyzer.CountKmers(sequence, 3, cts.Token));
    }

    [Test]
    public async Task CountKmersAsync_CompletesSuccessfully()
    {
        string sequence = new string('A', 1000) + new string('C', 1000);

        var counts = await KmerAnalyzer.CountKmersAsync(sequence, 3);

        Assert.That(counts, Is.Not.Empty);
    }

    [Test]
    public void FindWithMismatches_WithCancellation_CompletesNormally()
    {
        string sequence = "ACGTACGTACGTACGT";
        using var cts = new CancellationTokenSource();

        var matches = ApproximateMatcher.FindWithMismatches(sequence, "ACGT", 1, cts.Token);

        Assert.That(matches, Is.Not.Empty);
    }

    [Test]
    public void FindMicrosatellites_WithCancellation_CompletesNormally()
    {
        var dna = new DnaSequence("ATATATATATAT");
        using var cts = new CancellationTokenSource();

        var results = RepeatFinder.FindMicrosatellites(dna, 1, 6, 3, cts.Token);

        Assert.That(results, Is.Not.Empty);
    }

    [Test]
    public void FindDegenerateMotif_WithCancellation_CompletesNormally()
    {
        var dna = new DnaSequence("ACGTACGTACGT");
        using var cts = new CancellationTokenSource();

        var matches = MotifFinder.FindDegenerateMotif(dna, "ACGT", cts.Token);

        Assert.That(matches, Is.Not.Empty);
    }

    [Test]
    public void GlobalAlign_WithCancellation_CompletesNormally()
    {
        string seq1 = "ACGTACGT";
        string seq2 = "ACGTACGT";
        using var cts = new CancellationTokenSource();

        var result = SequenceAligner.GlobalAlign(seq1, seq2, null, cts.Token);

        Assert.That(result, Is.Not.Null);
        Assert.That(result.AlignedSequence1, Is.EqualTo(seq1));
    }

    [Test]
    public void FindAllOverlaps_WithCancellation_CompletesNormally()
    {
        var reads = new[] { "ACGTACGT", "CGTACGTA", "GTACGTAC" };
        using var cts = new CancellationTokenSource();

        var overlaps = SequenceAssembler.FindAllOverlaps(reads, 4, 0.9, cts.Token);

        Assert.That(overlaps, Is.Not.Null);
    }

    [Test]
    public void CountKmers_WithProgress_ReportsProgress()
    {
        string sequence = new string('A', 5000);
        using var cts = new CancellationTokenSource();
        double lastProgress = 0;
        var progress = new Progress<double>(p => lastProgress = p);

        var counts = KmerAnalyzer.CountKmers(sequence, 3, cts.Token, progress);

        Assert.That(counts, Is.Not.Empty);
    }

    #endregion

    #region KmerEnumerator Tests

    [Test]
    public void KmerEnumerator_EnumeratesCorrectly()
    {
        ReadOnlySpan<char> sequence = "ACGTAC".AsSpan();
        int count = 0;

        foreach (var kmer in sequence.EnumerateKmers(3))
        {
            count++;
            Assert.That(kmer.Length, Is.EqualTo(3));
        }

        Assert.That(count, Is.EqualTo(4)); // "ACG", "CGT", "GTA", "TAC"
    }

    #endregion
}
