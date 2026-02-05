using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for global sequence alignment (Needleman-Wunsch algorithm).
/// Test Unit: ALIGN-GLOBAL-001
/// 
/// Evidence Sources:
/// - Wikipedia: Needleman-Wunsch algorithm
/// - Wikipedia: Sequence alignment
/// - Needleman & Wunsch (1970) J Mol Biol 48:443-453
/// - Rosalind Bioinformatics Problems (GLOB, GCON)
/// </summary>
[TestFixture]
public class SequenceAligner_GlobalAlign_Tests
{
    private static readonly ScoringMatrix WikipediaSimpleScoring = new(
        Match: 1,
        Mismatch: -1,
        GapOpen: -1,
        GapExtend: -1);

    #region Reference Data Tests - Wikipedia Examples

    /// <summary>
    /// Validates alignment against Wikipedia Needleman-Wunsch example.
    /// Source: https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm
    /// 
    /// The Wikipedia example uses:
    /// - Seq1: GCATGCG (7 bp)
    /// - Seq2: GATTACA (7 bp)
    /// - Match: +1, Mismatch: -1, Gap: -1
    /// 
    /// Expected optimal alignment score: 0
    /// (This can be verified by hand calculation of the DP matrix)
    /// </summary>
    [Test]
    public void GlobalAlign_WikipediaExample_CorrectScore()
    {
        var seq1 = new DnaSequence("GCATGCG");
        var seq2 = new DnaSequence("GATTACA");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaSimpleScoring);

        // Wikipedia example with match=1, mismatch=-1, gap=-1
        // Manual calculation: optimal score = 0
        // Possible alignment: GCA-TGCG
        //                     G-ATTACA
        // Matches: G, A, T, A (4 × +1 = +4)
        // Mismatches: C vs T, G vs C (2 × -1 = -2)
        // Gaps: 2 gaps (2 × -1 = -2)
        // Total: 4 - 2 - 2 = 0

        Assert.That(result.Score, Is.EqualTo(0),
            "Wikipedia Needleman-Wunsch example: GCATGCG vs GATTACA with match=1, mismatch=-1, gap=-1 should yield score 0");
    }

    /// <summary>
    /// Validates alignment properties for Wikipedia example.
    /// </summary>
    [Test]
    public void GlobalAlign_WikipediaExample_ValidAlignment()
    {
        var seq1 = new DnaSequence("GCATGCG");
        var seq2 = new DnaSequence("GATTACA");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaSimpleScoring);

        Assert.Multiple(() =>
        {
            // Alignment type
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.Global));

            // Equal aligned lengths (fundamental invariant)
            Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length),
                "Aligned sequences must have equal length");

            // Original sequences recoverable
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("GCATGCG"),
                "Original seq1 recoverable from alignment");
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo("GATTACA"),
                "Original seq2 recoverable from alignment");

            // Score matches recalculation
            Assert.That(result.Score, Is.EqualTo(
                CalculateScore(result.AlignedSequence1, result.AlignedSequence2, WikipediaSimpleScoring)),
                "Score must match recalculation from alignment");
        });
    }

    #endregion

    #region Reference Data Tests - Classic Bioinformatics Examples

    /// <summary>
    /// Validates perfect match alignment.
    /// Source: Fundamental property - identical sequences should have maximum score.
    /// </summary>
    [Test]
    public void GlobalAlign_IdenticalSequences_MaximumScore()
    {
        var seq = new DnaSequence("ACGTACGT");

        var result = SequenceAligner.GlobalAlign(seq, seq, WikipediaSimpleScoring);

        Assert.Multiple(() =>
        {
            // Perfect match: 8 matches × +1 = 8
            Assert.That(result.Score, Is.EqualTo(8),
                "Identical 8bp sequences should score 8 with match=+1");

            // No gaps in alignment
            Assert.That(result.AlignedSequence1, Does.Not.Contain("-"));
            Assert.That(result.AlignedSequence2, Does.Not.Contain("-"));
        });
    }

    /// <summary>
    /// Validates alignment with Rosalind-style scoring.
    /// Source: Rosalind GLOB problem (Global Alignment with Scoring Matrix)
    /// 
    /// Uses BLOSUM62-like scoring concepts adapted for DNA.
    /// </summary>
    [Test]
    public void GlobalAlign_RosalindStyleScoring_CorrectBehavior()
    {
        // Rosalind typically uses:
        // Match: +5, Mismatch: -1, Gap opening: -5, Gap extension: -1
        var rosalindScoring = new ScoringMatrix(
            Match: 5,
            Mismatch: -1,
            GapOpen: -5,
            GapExtend: -1);

        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("ACGT");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, rosalindScoring);

        // Perfect match: 4 × +5 = 20
        Assert.That(result.Score, Is.EqualTo(20),
            "4 perfect matches with match=+5 should score 20");
    }

    /// <summary>
    /// Validates gap penalty behavior.
    /// Source: Needleman & Wunsch (1970) - Gap penalties discourage unnecessary gaps.
    /// </summary>
    [Test]
    public void GlobalAlign_GapPenalty_CorrectlyApplied()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("AGT");  // Missing C

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaSimpleScoring);

        // Expected alignment: ACGT
        //                     A-GT
        // Matches: A, G, T (3 × +1 = +3)
        // Gap: 1 × -1 = -1
        // Total: 3 - 1 = 2
        Assert.That(result.Score, Is.EqualTo(2),
            "ACGT vs AGT should score 2 (3 matches - 1 gap)");
    }

    #endregion

    #region Original Tests

    [Test]
    public void GlobalAlign_WikipediaExample_ReconstructsInputsAndScoreMatches()
    {
        var seq1 = new DnaSequence("GCATGCG");
        var seq2 = new DnaSequence("GATTACA");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaSimpleScoring);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.Global));
            Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo(seq1.Sequence));
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo(seq2.Sequence));
            Assert.That(result.Score, Is.EqualTo(CalculateScore(result.AlignedSequence1, result.AlignedSequence2, WikipediaSimpleScoring)));
        });
    }

    [Test]
    public void GlobalAlign_StringOverload_MatchesDnaSequenceResult()
    {
        const string seq1 = "GCATGCG";
        const string seq2 = "GATTACA";

        var dnaResult = SequenceAligner.GlobalAlign(new DnaSequence(seq1), new DnaSequence(seq2), WikipediaSimpleScoring);
        var stringResult = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaSimpleScoring);

        Assert.Multiple(() =>
        {
            Assert.That(stringResult.AlignmentType, Is.EqualTo(AlignmentType.Global));
            Assert.That(stringResult.Score, Is.EqualTo(dnaResult.Score));
            Assert.That(stringResult.AlignedSequence1, Is.EqualTo(dnaResult.AlignedSequence1));
            Assert.That(stringResult.AlignedSequence2, Is.EqualTo(dnaResult.AlignedSequence2));
        });
    }

    [Test]
    public void GlobalAlign_StringOverload_EmptyInput_ReturnsEmpty()
    {
        var result = SequenceAligner.GlobalAlign(string.Empty, "GATTACA", WikipediaSimpleScoring);

        Assert.That(result, Is.EqualTo(AlignmentResult.Empty));
    }

    [Test]
    public void GlobalAlign_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.GlobalAlign((DnaSequence)null!, new DnaSequence("GATTACA"), WikipediaSimpleScoring));
    }

    #endregion

    private static string RemoveGaps(string alignedSequence)
    {
        return new string(alignedSequence.Where(c => c != '-').ToArray());
    }

    private static int CalculateScore(string aligned1, string aligned2, ScoringMatrix scoring)
    {
        int score = 0;

        for (int i = 0; i < aligned1.Length; i++)
        {
            char a = aligned1[i];
            char b = aligned2[i];

            if (a == '-' || b == '-')
            {
                score += scoring.GapExtend;
            }
            else if (a == b)
            {
                score += scoring.Match;
            }
            else
            {
                score += scoring.Mismatch;
            }
        }

        return score;
    }
}
