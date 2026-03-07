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
/// Algorithm: Star alignment (progressive alignment variant)
/// - Selects center sequence via k-mer cosine similarity (max total similarity)
/// - Aligns all other sequences to center via anchor-based global alignment
/// - Reconciles pairwise alignments into a single MSA via gap merging
/// - Generates majority-voted consensus (gaps participate in vote)
/// - Computes true sum-of-pairs (SP) score across all C(k,2) pairs
/// 
/// Sources:
/// - Wikipedia: Multiple sequence alignment (https://en.wikipedia.org/wiki/Multiple_sequence_alignment)
/// - Wikipedia: Clustal (https://en.wikipedia.org/wiki/Clustal)
/// - Wikipedia: Consensus sequence (https://en.wikipedia.org/wiki/Consensus_sequence)
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
            Assert.That(result.TotalScore, Is.EqualTo(0), "Single sequence has no pairs for SP scoring");
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
            Assert.That(result.TotalScore, Is.EqualTo(4),
                "SP: 1 pair × 4 matches × Match(1) = 4");
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
            Assert.That(result.TotalScore, Is.EqualTo(24),
                "SP: C(3,2)=3 pairs × 8 matches × Match(1) = 24");
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
    /// M08: Different length sequences are padded with gaps; removing gaps recovers originals.
    /// Source: Wikipedia MSA - gap insertion for length normalization + reversibility
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
            Assert.That(result.AlignedSequences[0].Length, Is.GreaterThanOrEqualTo(8),
                "Aligned length ≥ max input length (Wikipedia MSA: L ≥ max{nᵢ})");
            // Shorter sequence must contain gaps after alignment
            Assert.That(result.AlignedSequences.Any(s => s.Contains('-')), Is.True,
                "Different-length input must produce gaps in alignment");
            // Reversibility (Wikipedia MSA)
            Assert.That(result.AlignedSequences[0].Replace("-", ""), Is.EqualTo("ATGCATGC"),
                "Removing gaps recovers first original sequence");
            Assert.That(result.AlignedSequences[1].Replace("-", ""), Is.EqualTo("ATGC"),
                "Removing gaps recovers second original sequence");
        });
    }

    /// <summary>
    /// M09: Consensus contains only valid DNA characters and gaps.
    /// Source: Wikipedia MSA - consensus validity
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
    /// M10: TotalScore is the true sum-of-pairs (SP) score across all C(k,2) sequence pairs.
    /// Source: Wikipedia MSA - "sum of all of the pairs of characters at each position in the alignment
    /// (the so-called sum of pair score)."
    /// Column-based scoring: match/mismatch from scoring matrix, gap-nucleotide = GapExtend,
    /// gap-gap = 0 (standard bioinformatics convention, not explicitly stated in Wikipedia).
    /// </summary>
    [Test]
    public void MultipleAlign_TotalScore_IsSumOfPairwiseScores()
    {
        // 3 same-length sequences with a mismatch at position 0.
        // Same length → deterministic alignment with no gaps.
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("CTGC")  // differs at position 0
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        // Hand-computed SP score per Wikipedia MSA (SimpleDna: Match=1, Mismatch=-1):
        // Aligned: "ATGC", "ATGC", "CTGC" (same length, no gaps needed)
        // C(3,2) = 3 pairs:
        //   (0,1): col0 A=A +1, col1 T=T +1, col2 G=G +1, col3 C=C +1 → 4
        //   (0,2): col0 A≠C -1, col1 T=T +1, col2 G=G +1, col3 C=C +1 → 2
        //   (1,2): col0 A≠C -1, col1 T=T +1, col2 G=G +1, col3 C=C +1 → 2
        // SP = 4 + 2 + 2 = 8
        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.All(s => !s.Contains('-')), Is.True,
                "Same-length input: no gaps needed");
            Assert.That(result.TotalScore, Is.EqualTo(8),
                "SP score: hand-computed per Wikipedia MSA definition");
        });
    }

    /// <summary>
    /// M11: Removing gaps from aligned sequence recovers original.
    /// Source: Wikipedia MSA - reversibility invariant: "To return from S'_i to S_i, remove all gaps."
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

    /// <summary>
    /// M12: No column in the alignment consists entirely of gaps.
    /// Source: Wikipedia MSA - "no values in the sequences of S of the same column consists of only gaps."
    /// </summary>
    [Test]
    public void MultipleAlign_NoColumnIsAllGaps()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("ATGCAA"),
            new DnaSequence("ATGCATGCGG")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        int length = result.AlignedSequences[0].Length;
        for (int col = 0; col < length; col++)
        {
            bool allGaps = result.AlignedSequences.All(s => col >= s.Length || s[col] == '-');
            Assert.That(allGaps, Is.False,
                $"Column {col} consists entirely of gaps, violating MSA invariant (Wikipedia MSA)");
        }
    }

    #endregion

    #region SHOULD Tests - Consensus and Scoring

    /// <summary>
    /// S01: Consensus reflects majority at each position.
    /// Source: Wikipedia Consensus sequence - "the calculated sequence of most frequent residues,
    /// either nucleotide or amino acid, found at each position in a sequence alignment."
    /// </summary>
    [Test]
    public void MultipleAlign_ConsensusReflectsMajority()
    {
        // Position 1 has a clear majority split: A=1, C=2 → C wins
        var sequences = new[]
        {
            new DnaSequence("AATG"),  // col1: A (minority)
            new DnaSequence("ACTG"),  // col1: C (majority)
            new DnaSequence("ACTG")   // col1: C (majority)
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        // Hand-computed consensus:
        // col0: A=3 → A; col1: A=1,C=2 → C; col2: T=3 → T; col3: G=3 → G
        Assert.That(result.Consensus, Is.EqualTo("ACTG"),
            "Consensus: col0 A(3), col1 C(2)>A(1), col2 T(3), col3 G(3)");
    }

    /// <summary>
    /// S02: Custom scoring matrix is honored.
    /// Source: Wikipedia Clustal - "gap opening penalty and gap extension penalty parameters can be adjusted."
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

        Assert.Multiple(() =>
        {
            // SimpleDna (Match=1): 1 pair × 4 matches × 1 = 4
            Assert.That(defaultResult.TotalScore, Is.EqualTo(4),
                "SimpleDna SP: 1 pair × 4 matches × Match(1) = 4");
            // BlastDna (Match=2): 1 pair × 4 matches × 2 = 8
            Assert.That(customResult.TotalScore, Is.EqualTo(8),
                "BlastDna SP: 1 pair × 4 matches × Match(2) = 8");
        });
    }

    /// <summary>
    /// S03: Partially overlapping sequences align correctly.
    /// Source: Wikipedia MSA - MSA applies to any set of biological sequences regardless of overlap pattern.
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
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(3), "Count preservation");
            int len = result.AlignedSequences[0].Length;
            Assert.That(result.AlignedSequences.All(s => s.Length == len), Is.True,
                "Equal length invariant");
            Assert.That(result.Consensus.Length, Is.EqualTo(len),
                "Consensus length = aligned length");
            Assert.That(len, Is.GreaterThanOrEqualTo(8),
                "L ≥ max{nᵢ} = 8 (Wikipedia MSA)");
            // Reversibility (Wikipedia MSA: removing gaps recovers original)
            Assert.That(result.AlignedSequences[0].Replace("-", ""), Is.EqualTo("ATGCATGC"));
            Assert.That(result.AlignedSequences[1].Replace("-", ""), Is.EqualTo("GCATGCAT"));
            Assert.That(result.AlignedSequences[2].Replace("-", ""), Is.EqualTo("TGCATGCA"));
        });
    }

    /// <summary>
    /// S04: When gap and nucleotide are tied in majority voting, nucleotide is preferred.
    /// Source: Implementation design choice. Wikipedia Consensus sequence defines consensus
    /// as "most frequent residues" but does not specify tie-breaking or gap handling.
    /// </summary>
    [Test]
    public void MultipleAlign_ConsensusTieBreaking_PrefersNucleotideOverGap()
    {
        // 4 sequences: seqs 0,1 are "ATGC" (4 chars), seqs 2,3 are "AGC" (3 chars).
        // After alignment, seqs 2,3 get a gap at position 1: "A-GC".
        // Position 1: T=2, '-'=2 → tie → nucleotide 'T' preferred.
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("AGC"),
            new DnaSequence("AGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        // Verify at least one tie column exists and is resolved to nucleotide
        int length = result.AlignedSequences[0].Length;
        bool foundTie = false;
        for (int col = 0; col < length; col++)
        {
            int gapCount = result.AlignedSequences.Count(s => s[col] == '-');
            int nucCount = result.AlignedSequences.Length - gapCount;
            if (gapCount > 0 && gapCount == nucCount)
            {
                Assert.That(result.Consensus[col], Is.Not.EqualTo('-'),
                    $"Column {col}: gap-nucleotide tie ({gapCount}:{nucCount}) must resolve to nucleotide");
                foundTie = true;
            }
        }
        Assert.That(foundTie, Is.True,
            "Test input must produce at least one gap-nucleotide tie column");
    }

    #endregion

    #region COULD Tests - Extended Coverage

    /// <summary>
    /// C01: Large sequence set completes in reasonable time.
    /// Source: Wikipedia Clustal - ClustalW complexity is O(N²); star alignment is O(k² × m).
    /// 20 diverse sequences of 50bp should complete well within 5 seconds.
    /// </summary>
    [Test]
    [CancelAfter(5000)] // 5 second timeout
    public void MultipleAlign_ManySequences_CompletesInReasonableTime()
    {
        // 20 diverse sequences of 50bp each - fixed seed for reproducibility
        var random = new Random(42);
        var sequences = Enumerable.Range(0, 20)
            .Select(_ => new DnaSequence(GenerateRandomDnaSequence(50, random)))
            .ToArray();

        // Verify test setup: sequences should be diverse (not all identical)
        Assert.That(sequences.Select(s => s.Sequence).Distinct().Count(), Is.GreaterThan(1),
            "Test setup: sequences should be diverse");

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(20), "Count preservation");
            Assert.That(result.AlignedSequences.All(s => s.Length == result.AlignedSequences[0].Length), Is.True,
                "Equal length invariant");
            Assert.That(result.Consensus.Length, Is.EqualTo(result.AlignedSequences[0].Length),
                "Consensus length = aligned length");
        });
    }

    /// <summary>
    /// C02: Empty sequences in collection handled gracefully.
    /// Source: Wikipedia MSA - edge case; sequences of length 0 represent degenerate input.
    /// </summary>
    [Test]
    public void MultipleAlign_WithEmptySequence_HandlesGracefully()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence(""),
            new DnaSequence("ATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(3), "Count preservation");
            Assert.That(result.AlignedSequences.All(s => s.Length == result.AlignedSequences[0].Length), Is.True,
                "Equal length invariant");
            // Non-empty sequences recover via gap removal (Wikipedia MSA: reversibility)
            Assert.That(result.AlignedSequences[0].Replace("-", ""), Is.EqualTo("ATGC"));
            Assert.That(result.AlignedSequences[2].Replace("-", ""), Is.EqualTo("ATGC"));
            // Empty sequence produces all gaps or empty
            Assert.That(result.AlignedSequences[1].Replace("-", ""), Is.EqualTo(""));
        });
    }

    #endregion

    #region Helper Methods

    private static string GenerateRandomDnaSequence(int length, Random random)
    {
        const string nucleotides = "ACGT";
        return new string(Enumerable.Range(0, length)
            .Select(_ => nucleotides[random.Next(4)])
            .ToArray());
    }

    #endregion
}
