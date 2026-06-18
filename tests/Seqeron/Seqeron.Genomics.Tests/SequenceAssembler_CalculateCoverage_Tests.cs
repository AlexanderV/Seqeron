// ASSEMBLY-COVER-001 — Coverage (Depth) Calculation
// Evidence: docs/Evidence/ASSEMBLY-COVER-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-COVER-001.md
// Source: Metagenomics Wiki, "SAMtools: get breadth of coverage" (per-base depth = number of
//         reads covering a position); Cook D.E., "Calculate Depth and Breadth of Coverage From a
//         bam File" (average depth = Sum of Depths / genome size; breadth = covered bases / size);
//         Illumina, "Sequencing Coverage for NGS Experiments" (C = LN/G).

using System;
using System.Linq;
using NUnit.Framework;

using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceAssembler_CalculateCoverage_Tests
{
    #region CalculateCoverage

    // Worked dataset (Evidence §"Test Datasets"): reference ACGTTGCAAT with three 5-mer reads that
    // each match uniquely at positions 0, 3 and 5. Distinct substrings are used so each read's best
    // placement is unambiguous, isolating the source-defined depth-counting rule.
    private const string WorkedReference = "ACGTTGCAAT";              // length 10
    private static readonly string[] WorkedReads = { "ACGTT", "TTGCA", "GCAAT" }; // placed at 0, 3, 5

    // M1 — Per-base depth = number of reads spanning each position (Metagenomics Wiki).
    [Test]
    public void CalculateCoverage_WorkedDataset_ReturnsExpectedDepthArray()
    {
        int[] depth = SequenceAssembler.CalculateCoverage(WorkedReference, WorkedReads, minOverlap: 5);

        Assert.That(depth, Is.EqualTo(new[] { 1, 1, 1, 2, 2, 2, 2, 2, 1, 1 }),
            "Reads cover [0,5),[3,8),[5,10); per-base depth = count of reads spanning each position (Metagenomics Wiki).");
    }

    // M2 — Output array length always equals reference length (INV-01).
    [Test]
    public void CalculateCoverage_AnyInput_ArrayLengthEqualsReferenceLength()
    {
        var reference = "ACGTACGTAC"; // length 10
        var reads = new[] { "ACGTA", "CGTAC" };

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        Assert.That(depth.Length, Is.EqualTo(reference.Length),
            "One depth value is reported per reference position (INV-01; depth defined per reference base).");
    }

    // M3 — A single placed read of length L increments exactly the half-open interval [p, p+L) (INV-04).
    [Test]
    public void CalculateCoverage_SingleRead_IncrementsExactInterval()
    {
        var reference = "GGGGGAAAAAGGGGG"; // read AAAAA matches uniquely at position 5
        var reads = new[] { "AAAAA" };

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        Assert.That(depth, Is.EqualTo(new[] { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 }),
            "Read placed at position 5 covers exactly [5,10); all other positions remain 0 (INV-04).");
    }

    // M4 — A read longer than the reference cannot be placed and contributes 0 (INV-04/INV-05 boundary).
    [Test]
    public void CalculateCoverage_ReadLongerThanReference_ContributesZero()
    {
        var reference = "AAAA";          // length 4
        var reads = new[] { "AAAAAA" };  // length 6, cannot fit

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 4);

        Assert.That(depth, Is.EqualTo(new[] { 0, 0, 0, 0 }),
            "A read longer than the reference cannot be placed (best-match scan requires it to fit), so depth is all zero.");
    }

    // M5 — A read with fewer than minOverlap matching characters fails to place and adds 0 everywhere (INV-05).
    [Test]
    public void CalculateCoverage_UnmatchedRead_ContributesZero()
    {
        var reference = "AAAAAAAAAA";
        var reads = new[] { "GGGGG" }; // 0 matches vs an all-A reference < minOverlap 5

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        Assert.That(depth, Is.EqualTo(new[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }),
            "A read below the minOverlap match floor is not placed and contributes no depth (INV-05).");
    }

    // M6 — Empty reads list yields an all-zero array of reference length (no aligned reads → depth 0).
    [Test]
    public void CalculateCoverage_EmptyReads_ReturnsAllZeroArray()
    {
        var reference = "AAAAAAAAAA";
        var reads = Array.Empty<string>();

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        Assert.Multiple(() =>
        {
            Assert.That(depth.Length, Is.EqualTo(10),
                "Array length tracks the reference length even with no reads (INV-01).");
            Assert.That(depth.All(d => d == 0), Is.True,
                "With no reads, every per-base depth is 0 (Cook; Metagenomics Wiki).");
        });
    }

    // M7 — Average depth = Sum of per-base depths / reference length = 15 / 10 = 1.5 (Cook).
    [Test]
    public void CalculateCoverage_WorkedDataset_AverageDepthIsOnePointFive()
    {
        int[] depth = SequenceAssembler.CalculateCoverage(WorkedReference, WorkedReads, minOverlap: 5);
        double average = depth.Sum() / (double)depth.Length;

        Assert.Multiple(() =>
        {
            Assert.That(depth.Sum(), Is.EqualTo(15),
                "Sum of per-base depths = total bases mapped = 3 reads x 5 bases = 15 (INV-03).");
            Assert.That(average, Is.EqualTo(1.5).Within(1e-10),
                "Average depth = Sum of Depths / genome size = 15 / 10 = 1.5 (Cook).");
        });
    }

    // M8 — Breadth of coverage = (#positions with depth >= 1) / reference length = 10 / 10 = 1.0 (Metagenomics Wiki).
    [Test]
    public void CalculateCoverage_WorkedDataset_BreadthIsOne()
    {
        int[] depth = SequenceAssembler.CalculateCoverage(WorkedReference, WorkedReads, minOverlap: 5);
        double breadth = depth.Count(d => d >= 1) / (double)depth.Length;

        Assert.That(breadth, Is.EqualTo(1.0).Within(1e-10),
            "Breadth = covered bases / reference length = 10 / 10 = 1.0; every position has >=1 read (Metagenomics Wiki).");
    }

    // S1 — Partial coverage: one read covering half the reference gives breadth 0.5 and average 0.5.
    [Test]
    public void CalculateCoverage_PartialCoverage_BreadthAndAverageAreHalf()
    {
        var reference = "AAAAAAAAAA"; // length 10
        var reads = new[] { "AAAAA" };  // placed at position 0, covers [0,5)

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);
        double breadth = depth.Count(d => d >= 1) / (double)depth.Length;
        double average = depth.Sum() / (double)depth.Length;

        Assert.Multiple(() =>
        {
            Assert.That(depth, Is.EqualTo(new[] { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 }),
                "Single read covers exactly the first 5 positions (INV-04).");
            Assert.That(breadth, Is.EqualTo(0.5).Within(1e-10),
                "Breadth = 5 covered / 10 = 0.5 (Metagenomics Wiki).");
            Assert.That(average, Is.EqualTo(0.5).Within(1e-10),
                "Average = 5 / 10 = 0.5 (Cook).");
        });
    }

    // S2 — Property/invariant: every per-base depth is non-negative (INV-02), and the sum equals
    // total bases mapped (INV-03) for the worked dataset.
    [Test]
    public void CalculateCoverage_WorkedDataset_AllDepthsNonNegative()
    {
        int[] depth = SequenceAssembler.CalculateCoverage(WorkedReference, WorkedReads, minOverlap: 5);

        Assert.Multiple(() =>
        {
            Assert.That(depth.All(d => d >= 0), Is.True,
                "Depth is a count of reads, so every element is non-negative (INV-02).");
            Assert.That(depth.Sum(), Is.EqualTo(3 * 5),
                "Sum of per-base depths equals total overlap of placed reads with the reference (INV-03).");
        });
    }

    // C1 — Case-insensitive placement: a lowercase read maps against an uppercase reference.
    [Test]
    public void CalculateCoverage_LowercaseRead_MapsCaseInsensitively()
    {
        var reference = "AAAAAAAAAA";
        var reads = new[] { "aaaaa" }; // lowercase, placed at position 0

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        Assert.That(depth, Is.EqualTo(new[] { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 }),
            "Matching is case-insensitive; the lowercase read places and increments [0,5) like its uppercase form.");
    }

    // Edge — Empty reference yields an empty depth array (one value per reference position; zero
    // positions → length-0 array). Per-base depth is defined per reference base (Metagenomics Wiki);
    // with no reference bases there are no depth values. (Description §3.3.)
    [Test]
    public void CalculateCoverage_EmptyReference_ReturnsEmptyArray()
    {
        var reads = new[] { "AAAAA" };

        int[] depth = SequenceAssembler.CalculateCoverage(string.Empty, reads, minOverlap: 5);

        Assert.That(depth, Is.Empty,
            "An empty reference has no positions, so the per-base depth array is empty (length 0).");
    }

    // Edge — null reference throws ArgumentNullException (input validation).
    [Test]
    public void CalculateCoverage_NullReference_ThrowsArgumentNullException()
    {
        var reads = new[] { "AAAAA" };

        Assert.Throws<ArgumentNullException>(
            () => SequenceAssembler.CalculateCoverage(null!, reads, minOverlap: 5),
            "Null reference is rejected with ArgumentNullException (sibling validation convention).");
    }

    // Edge — null reads throws ArgumentNullException (input validation).
    [Test]
    public void CalculateCoverage_NullReads_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => SequenceAssembler.CalculateCoverage("AAAAAAAAAA", null!, minOverlap: 5),
            "Null reads list is rejected with ArgumentNullException (sibling validation convention).");
    }

    #endregion
}
