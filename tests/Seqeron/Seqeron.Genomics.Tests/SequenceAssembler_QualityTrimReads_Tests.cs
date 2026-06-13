// ASSEMBLY-TRIM-001 — Quality Trimming (BWA / cutadapt running-sum)
// Evidence: docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-TRIM-001.md
// Source: Cutadapt algorithm docs (quality trimming); Li H. BWA bwa_trim_read;
//         Cock et al. (2010) NAR 38(6):1767-1771 (Phred+33).

using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceAssembler_QualityTrimReads_Tests
{
    #region QualityTrimReads

    // Phred+33 (Cock et al. 2010): char = Phred + 33.
    // Cutadapt worked example qualities 42,40,26,27,8,7,11,4,2,3 -> "KI;<)(,%#$".
    private const string ExampleQuality = "KI;<)(,%#$";

    // M1 — Cutadapt worked example: subtract cutoff 10, partial sums from end
    //      (70),(38),8,-8,-25,-23,-20,-21,-15,-7; minimum -25 at index 4 ->
    //      read trimmed to the first four bases. (Cutadapt algorithm docs.)
    [Test]
    public void QualityTrimReads_CutadaptWorkedExample_TrimsToFirstFourBases()
    {
        // Arrange
        var reads = new List<(string, string)> { ("ACGTACGTAC", ExampleQuality) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "single read survives min-length 1");
        Assert.That(result[0], Is.EqualTo("ACGT"),
            "cutadapt example: minimum partial sum -25 at index 4 keeps the first four bases");
    }

    // M2 — All-high-quality read: every (q - cutoff) >= 0, partial sums never go
    //      below 0, minimum is at the end -> nothing trimmed.
    [Test]
    public void QualityTrimReads_AllHighQuality_ReturnsUnchanged()
    {
        // Arrange: 'I' = Phred 40; cutoff 20.
        var reads = new List<(string, string)> { ("ACGTAC", "IIIIII") };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "high-quality read is kept");
        Assert.That(result[0], Is.EqualTo("ACGTAC"),
            "all qualities above cutoff -> no bases removed");
    }

    // M3 — All-low-quality read: every (q - cutoff) < 0, both end passes consume the
    //      read, trimmed length 0 < minLength -> dropped.
    [Test]
    public void QualityTrimReads_AllLowQuality_DropsRead()
    {
        // Arrange: '!' = Phred 0; cutoff 20.
        var reads = new List<(string, string)> { ("ACGTAC", "!!!!!!") };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 20, minLength: 1);

        // Assert
        Assert.That(result, Is.Empty,
            "all qualities below cutoff -> read trimmed to length 0 and dropped");
    }

    // M4 — Threshold below 1 disables trimming (BWA bwa_trim_read: trim_qual < 1 -> 0).
    [Test]
    public void QualityTrimReads_ThresholdZero_ReturnsUnchanged()
    {
        // Arrange: mixed qualities that WOULD trim under a positive cutoff.
        var reads = new List<(string, string)> { ("ACGTACGTAC", ExampleQuality) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 0, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "read kept");
        Assert.That(result[0], Is.EqualTo("ACGTACGTAC"),
            "cutoff < 1 disables trimming per BWA guard");
    }

    // M5 — Min-length filter: the cutadapt example trims to 4 bases; with minLength 5
    //      the survivor (length 4) is dropped.
    [Test]
    public void QualityTrimReads_TrimmedShorterThanMinLength_DropsRead()
    {
        // Arrange
        var reads = new List<(string, string)> { ("ACGTACGTAC", ExampleQuality) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 5);

        // Assert
        Assert.That(result, Is.Empty,
            "trimmed length 4 < minLength 5 -> read dropped");
    }

    // S1 — 5'-end trimming: low-quality prefix, high-quality suffix.
    //      Qualities 3,2,4,11,7,8,27,26,40,42 ("$#%,()<;IK"), cutoff 10.
    //      3' pass: partial sums from end never negative -> no 3' trim.
    //      5' pass: partial sums from start -7,-15,-21,-20,-23,-25,...; minimum -25 at
    //      index 5 -> cut after index 5 -> keep the last four bases.
    [Test]
    public void QualityTrimReads_LowQualityPrefix_TrimsFivePrimeEnd()
    {
        // Arrange
        var reads = new List<(string, string)> { ("ACGTACGTAC", "$#%,()<;IK") };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "read survives");
        Assert.That(result[0], Is.EqualTo("GTAC"),
            "5' running-sum minimum at index 5 keeps the last four (high-quality) bases");
    }

    // S2 — Min-length boundary: trimmed length exactly equals minLength -> kept.
    [Test]
    public void QualityTrimReads_TrimmedEqualsMinLength_KeepsRead()
    {
        // Arrange: cutadapt example trims to 4; minLength 4.
        var reads = new List<(string, string)> { ("ACGTACGTAC", ExampleQuality) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 4);

        // Assert
        Assert.That(result, Has.Count.EqualTo(1), "trimmed length 4 == minLength 4 -> kept");
        Assert.That(result[0], Is.EqualTo("ACGT"),
            "boundary: length-4 survivor retained at minLength 4");
    }

    // S3 — Multiple reads trimmed independently: one trimmable, one droppable, one clean.
    [Test]
    public void QualityTrimReads_MultipleReads_TrimsEachIndependently()
    {
        // Arrange
        var reads = new List<(string, string)>
        {
            ("ACGTACGTAC", ExampleQuality), // -> "ACGT"
            ("ACGTAC", "!!!!!!"),           // all-low -> dropped
            ("ACGTAC", "IIIIII"),           // all-high -> unchanged
        };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(2), "the all-low read is dropped");
            Assert.That(result[0], Is.EqualTo("ACGT"), "first read trimmed to four bases");
            Assert.That(result[1], Is.EqualTo("ACGTAC"), "high-quality read unchanged, order preserved");
        });
    }

    // C1 — Empty reads list -> empty result.
    [Test]
    public void QualityTrimReads_EmptyList_ReturnsEmpty()
    {
        // Arrange
        var reads = new List<(string, string)>();

        // Act
        IReadOnlyList<string> result = SequenceAssembler.QualityTrimReads(reads);

        // Assert
        Assert.That(result, Is.Empty, "no reads -> empty output");
    }

    // C2 — Empty sequence/quality read -> length 0 < minLength (default) -> dropped.
    [Test]
    public void QualityTrimReads_EmptySequence_DropsRead()
    {
        // Arrange
        var reads = new List<(string, string)> { (string.Empty, string.Empty) };

        // Act
        IReadOnlyList<string> result =
            SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);

        // Assert
        Assert.That(result, Is.Empty, "empty read has length 0 < minLength -> dropped");
    }

    // C3 — Null reads argument -> ArgumentNullException.
    [Test]
    public void QualityTrimReads_NullReads_Throws()
    {
        // Act / Assert
        Assert.That(
            () => SequenceAssembler.QualityTrimReads(null!),
            NUnit.Framework.Throws.ArgumentNullException,
            "null reads collection is rejected");
    }

    #endregion
}
