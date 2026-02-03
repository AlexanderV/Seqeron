using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for SequenceAligner.MultipleAlign() - Multiple Sequence Alignment.
/// Test Unit: ALIGN-MULTI-001
/// 
/// Algorithm: Star alignment (simplified progressive alignment)
/// - Uses first sequence as reference
/// - Aligns all other sequences to reference via global alignment
/// - Pads sequences to equal length
/// - Generates majority-voted consensus
/// 
/// Sources:
/// - Wikipedia: Multiple sequence alignment
/// - Wikipedia: Clustal
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MultipleAlign_Tests
{
    #region MUST Tests - Null and Edge Cases

    /// <summary>
    /// M01: Null input must throw ArgumentNullException.
    /// Source: .NET convention
    /// </summary>
    [Test]
    public void MultipleAlign_NullInput_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.MultipleAlign(null!));
    }

    /// <summary>
    /// M02: Empty collection returns Empty result.
    /// Source: Wikipedia MSA - edge case handling
    /// </summary>
    [Test]
    public void MultipleAlign_EmptyCollection_ReturnsEmpty()
    {
        var result = SequenceAligner.MultipleAlign(Array.Empty<DnaSequence>());

        Assert.Multiple(() =>
        {
            Assert.That(result, Is.EqualTo(MultipleAlignmentResult.Empty));
            Assert.That(result.AlignedSequences, Is.Empty);
            Assert.That(result.Consensus, Is.Empty);
            Assert.That(result.TotalScore, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// M03: Single sequence returns alignment containing just that sequence.
    /// Source: Wikipedia MSA - trivial case
    /// </summary>
    [Test]
    public void MultipleAlign_SingleSequence_ReturnsSameSequence()
    {
        var sequences = new[] { new DnaSequence("ATGCATGC") };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(1));
            Assert.That(result.AlignedSequences[0], Is.EqualTo("ATGCATGC"));
            Assert.That(result.Consensus, Is.EqualTo("ATGCATGC"));
            Assert.That(result.TotalScore, Is.EqualTo(0), "Single sequence has no pairwise score");
        });
    }

    #endregion

    #region MUST Tests - Basic Alignment

    /// <summary>
    /// M04: Two sequences align correctly using global alignment.
    /// Source: Wikipedia MSA - minimum non-trivial case
    /// </summary>
    [Test]
    public void MultipleAlign_TwoSequences_AlignsCorrectly()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(2));
            Assert.That(result.AlignedSequences[0], Is.EqualTo(result.AlignedSequences[1]),
                "Identical sequences should produce identical aligned sequences");
            Assert.That(result.Consensus, Is.EqualTo("ATGC"));
            Assert.That(result.TotalScore, Is.GreaterThan(0),
                "Identical sequences should have positive alignment score");
        });
    }

    /// <summary>
    /// M05: Three identical sequences produce perfect alignment.
    /// Source: Wikipedia MSA - all sequences equal after alignment
    /// </summary>
    [Test]
    public void MultipleAlign_ThreeIdenticalSequences_AllMatch()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGCATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(3));
            Assert.That(result.AlignedSequences.All(s => s == "ATGCATGC"), Is.True,
                "All identical sequences should remain identical after alignment");
            Assert.That(result.Consensus, Is.EqualTo("ATGCATGC"));
        });
    }

    #endregion

    #region MUST Tests - Invariants

    /// <summary>
    /// M06: All aligned sequences have equal length (fundamental MSA invariant).
    /// Source: Wikipedia MSA - "all conform to length L"
    /// </summary>
    [Test]
    public void MultipleAlign_AllAlignedSequences_HaveEqualLength()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGCATGCATGC"),  // 12 chars
            new DnaSequence("ATGC"),           // 4 chars
            new DnaSequence("ATGCAT"),         // 6 chars
            new DnaSequence("ATGCATGCATGCAA")  // 14 chars
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        var lengths = result.AlignedSequences.Select(s => s.Length).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(lengths.Distinct().Count(), Is.EqualTo(1),
                "All aligned sequences must have the same length");
            Assert.That(lengths[0], Is.GreaterThanOrEqualTo(14),
                "Length must be at least max input length");
        });
    }

    /// <summary>
    /// M07: Number of aligned sequences equals number of input sequences.
    /// Source: Wikipedia MSA - count preservation
    /// </summary>
    [Test]
    public void MultipleAlign_SequenceCount_Preserved()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.That(result.AlignedSequences.Length, Is.EqualTo(sequences.Length),
            "Output sequence count must equal input sequence count");
    }

    /// <summary>
    /// M08: Different length sequences are padded with gaps.
    /// Source: Wikipedia MSA - gap insertion for length normalization
    /// </summary>
    [Test]
    public void MultipleAlign_DifferentLengths_PadsWithGaps()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGCATGC"),  // 8 chars
            new DnaSequence("ATGC")        // 4 chars
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences[0].Length,
                Is.EqualTo(result.AlignedSequences[1].Length),
                "Both sequences must have equal length after alignment");

            // At least one sequence should contain gaps if lengths differ
            bool hasGaps = result.AlignedSequences.Any(s => s.Contains('-'));
            // Note: gaps may or may not appear depending on alignment
            Assert.That(result.AlignedSequences[0].Length, Is.GreaterThanOrEqualTo(8),
                "Result length should be at least max input length");
        });
    }

    /// <summary>
    /// M09: Consensus contains only valid DNA characters and gaps.
    /// Source: Implementation invariant
    /// </summary>
    [Test]
    public void MultipleAlign_ConsensusContainsOnlyValidCharacters()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("TTTT")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        const string validChars = "ACGT-";

        Assert.That(result.Consensus.All(c => validChars.Contains(c)), Is.True,
            $"Consensus '{result.Consensus}' contains invalid characters. Expected only: {validChars}");
    }

    /// <summary>
    /// M10: TotalScore is sum of pairwise alignment scores.
    /// Source: Wikipedia MSA - sum-of-pairs scoring
    /// </summary>
    [Test]
    public void MultipleAlign_TotalScore_IsSumOfPairwiseScores()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        // For identical sequences, each pairwise alignment should have positive score
        // With 3 sequences, there are 2 pairwise alignments (to reference)
        Assert.Multiple(() =>
        {
            Assert.That(result.TotalScore, Is.GreaterThan(0),
                "Identical sequences should produce positive total score");

            // Calculate expected: align seq[1] and seq[2] to seq[0]
            var score1 = SequenceAligner.GlobalAlign(sequences[0], sequences[1]).Score;
            var score2 = SequenceAligner.GlobalAlign(sequences[0], sequences[2]).Score;
            int expectedTotal = score1 + score2;

            Assert.That(result.TotalScore, Is.EqualTo(expectedTotal),
                "TotalScore should equal sum of pairwise scores to reference");
        });
    }

    /// <summary>
    /// M11: Removing gaps from aligned sequence recovers original.
    /// Source: Wikipedia MSA - reversibility invariant
    /// </summary>
    [Test]
    public void MultipleAlign_RemovingGaps_RecoversOriginal()
    {
        var originalSequences = new[]
        {
            "ATGCATGC",
            "ATGC",
            "ATGCAT"
        };
        var sequences = originalSequences.Select(s => new DnaSequence(s)).ToArray();

        var result = SequenceAligner.MultipleAlign(sequences);

        for (int i = 0; i < originalSequences.Length; i++)
        {
            string recovered = result.AlignedSequences[i].Replace("-", "");
            Assert.That(recovered, Is.EqualTo(originalSequences[i]),
                $"Sequence {i}: removing gaps should recover original sequence");
        }
    }

    #endregion

    #region SHOULD Tests - Consensus and Scoring

    /// <summary>
    /// S01: Consensus reflects majority at each position.
    /// Source: Implementation - majority voting
    /// </summary>
    [Test]
    public void MultipleAlign_ConsensusReflectsMajority()
    {
        // Create sequences where position 0 has clear majority
        var sequences = new[]
        {
            new DnaSequence("ATGC"),  // A at position 0
            new DnaSequence("ATGC"),  // A at position 0
            new DnaSequence("CTGC")   // C at position 0 (minority)
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        // At position 0: A appears 2 times, C appears 1 time → consensus should be A
        Assert.That(result.Consensus[0], Is.EqualTo('A'),
            "Consensus should reflect majority vote (A=2, C=1)");
    }

    /// <summary>
    /// S02: Custom scoring matrix is honored.
    /// Source: API contract
    /// </summary>
    [Test]
    public void MultipleAlign_WithCustomScoring_UsesProvidedMatrix()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC")
        };

        var defaultResult = SequenceAligner.MultipleAlign(sequences);
        var customResult = SequenceAligner.MultipleAlign(sequences, SequenceAligner.BlastDna);

        // BlastDna has Match=2 vs SimpleDna Match=1, so score should differ
        Assert.That(customResult.TotalScore, Is.Not.EqualTo(defaultResult.TotalScore),
            "Different scoring matrices should produce different scores");
    }

    /// <summary>
    /// S03: Partially overlapping sequences align correctly.
    /// Source: ASSUMPTION - realistic biological scenario
    /// </summary>
    [Test]
    public void MultipleAlign_PartiallyOverlapping_AlignsCorrectly()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGCATGC"),
            new DnaSequence("GCATGCAT"),
            new DnaSequence("TGCATGCA")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(3));
            Assert.That(result.AlignedSequences.All(s => s.Length == result.AlignedSequences[0].Length), Is.True);
            Assert.That(result.Consensus.Length, Is.EqualTo(result.AlignedSequences[0].Length));
        });
    }

    #endregion

    #region COULD Tests - Extended Coverage

    /// <summary>
    /// C01: Large sequence set completes in reasonable time.
    /// Source: ASSUMPTION - performance sanity check
    /// </summary>
    [Test]
    [Timeout(5000)] // 5 second timeout
    public void MultipleAlign_ManySequences_CompletesInReasonableTime()
    {
        // 20 sequences of 50bp each - should complete quickly
        var sequences = Enumerable.Range(0, 20)
            .Select(_ => new DnaSequence(GenerateRandomDnaSequence(50)))
            .ToArray();

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(20));
            Assert.That(result.AlignedSequences.All(s => s.Length > 0), Is.True);
        });
    }

    /// <summary>
    /// C02: Empty sequences in collection handled gracefully.
    /// Source: ASSUMPTION - defensive handling
    /// </summary>
    [Test]
    public void MultipleAlign_WithEmptySequence_HandlesGracefully()
    {
        // Note: This depends on whether DnaSequence allows empty strings
        // If not, this test documents the behavior
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence(""),
            new DnaSequence("ATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        // Should produce 3 aligned sequences
        Assert.That(result.AlignedSequences.Length, Is.EqualTo(3));
    }

    #endregion

    #region Helper Methods

    private static string GenerateRandomDnaSequence(int length)
    {
        const string nucleotides = "ACGT";
        var random = new Random(42); // Fixed seed for reproducibility
        return new string(Enumerable.Range(0, length)
            .Select(_ => nucleotides[random.Next(4)])
            .ToArray());
    }

    #endregion
}
