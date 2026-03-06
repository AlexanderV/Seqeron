using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for global sequence alignment (Needleman–Wunsch algorithm).
/// Test Unit: ALIGN-GLOBAL-001
///
/// All tests are grounded in the standard Needleman–Wunsch algorithm as documented in:
///   https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
///   https://en.wikipedia.org/wiki/Sequence_alignment
///
/// The algorithm uses a LINEAR gap penalty d. The ScoringMatrix.GapExtend field
/// represents d. ScoringMatrix.GapOpen is NOT used by the Needleman–Wunsch linear model.
///
/// Recurrence (from Wikipedia):
///   F(0,j) = d·j,   F(i,0) = d·i
///   F(i,j) = max( F(i−1,j−1)+S(Aᵢ,Bⱼ),  F(i−1,j)+d,  F(i,j−1)+d )
/// </summary>
[TestFixture]
public class SequenceAligner_GlobalAlign_Tests
{
    /// <summary>
    /// Scoring matrix matching the Wikipedia Needleman–Wunsch example.
    /// Match = +1, Mismatch = −1, Gap penalty d = −1 (linear).
    /// Source: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm#Choosing_a_scoring_system
    /// </summary>
    private static readonly ScoringMatrix WikipediaScoring = new(
        Match: 1,
        Mismatch: -1,
        GapOpen: 0,     // Not used by Needleman–Wunsch (linear gap model)
        GapExtend: -1);

    #region Wikipedia Needleman–Wunsch Example

    /// <summary>
    /// Wikipedia example: GCATGCG vs GATTACA, match=+1, mismatch/gap=−1.
    /// Expected optimal score = 0.
    ///
    /// One possible optimal alignment (shown on the Wikipedia page):
    ///   GCATG-CG        (seq1 with gap at position 6)
    ///   G-ATTACA        (seq2 with gap at position 2)
    ///   +−++−−+−  →  4×(+1) + 2×(−1) + 2×(−1) = 0
    ///
    /// Source: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm#Choosing_a_scoring_system
    /// </summary>
    [Test]
    public void WikipediaExample_OptimalScore_IsZero()
    {
        var seq1 = new DnaSequence("GCATGCG");
        var seq2 = new DnaSequence("GATTACA");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaScoring);

        Assert.That(result.Score, Is.EqualTo(0),
            "Wikipedia NW example: GCATGCG vs GATTACA with match=1, mismatch=-1, gap=-1 → score 0");
    }

    /// <summary>
    /// Verifies the three fundamental invariants of Needleman–Wunsch alignment
    /// using the Wikipedia example dataset.
    ///
    /// INV-1: Aligned sequences have equal length.
    ///        (Global alignment pads both sequences to the same length with gaps.)
    /// INV-2: Removing gaps from aligned sequences yields the original sequences.
    ///        (Traceback only inserts gaps; it never modifies characters.)
    /// INV-3: Score equals the sum of per-position match/mismatch/gap scores.
    ///        (Alignment score is defined as this sum — Wikipedia "Choosing a scoring system".)
    ///
    /// Source: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
    ///         https://en.wikipedia.org/wiki/Sequence_alignment (global alignment definition)
    /// </summary>
    [Test]
    public void WikipediaExample_AllInvariants_Hold()
    {
        var seq1 = new DnaSequence("GCATGCG");
        var seq2 = new DnaSequence("GATTACA");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.Global));

            // INV-1: Equal aligned lengths
            Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length),
                "INV-1: Aligned sequences must have equal length");

            // INV-2: Gap removal yields originals
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("GCATGCG"),
                "INV-2: Removing gaps from aligned seq1 must yield original GCATGCG");
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo("GATTACA"),
                "INV-2: Removing gaps from aligned seq2 must yield original GATTACA");

            // INV-3: Score = recalculated sum of per-position scores
            Assert.That(result.Score, Is.EqualTo(
                RecalculateLinearScore(result.AlignedSequence1, result.AlignedSequence2, WikipediaScoring)),
                "INV-3: Reported score must equal recalculated per-position sum");
        });
    }

    #endregion

    #region Linear Gap Penalty — Border Initialization Verification

    /// <summary>
    /// Verifies correct linear gap penalty for unequal-length sequences
    /// where the optimal alignment MUST use border cells.
    ///
    /// Seq1 = "T", Seq2 = "ACGT", match=+1, mismatch=−1, gap=−1.
    ///
    /// The only alignment that achieves a match is to align T with the final T in ACGT:
    ///   ---T
    ///   ACGT
    /// Score = 3×(−1) + 1×(+1) = −2
    ///
    /// Standard NW initialization: F(0,3)=−3, so F(1,4)=max(F(0,3)+1,...)=−2.
    /// A faulty initialization that adds extra penalty to the border would yield −3.
    ///
    /// Source: F(0,j) = d·j from Wikipedia pseudocode.
    /// </summary>
    [Test]
    public void UnequalLengths_ShortVsLong_CorrectLinearBorderScore()
    {
        var seq1 = new DnaSequence("T");
        var seq2 = new DnaSequence("ACGT");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(-2),
                "T vs ACGT: optimal alignment ---T/ACGT scores 3×(−1)+1×(+1)=−2");

            // INV-1, INV-2
            Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("T"));
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo("ACGT"));

            // INV-3
            Assert.That(result.Score, Is.EqualTo(
                RecalculateLinearScore(result.AlignedSequence1, result.AlignedSequence2, WikipediaScoring)));
        });
    }

    /// <summary>
    /// Verifies correct linear gap penalty when seq1 is longer than seq2.
    /// Symmetric case to the above test.
    ///
    /// Seq1 = "ACGT", Seq2 = "T", match=+1, mismatch=−1, gap=−1.
    ///
    /// Optimal alignment:
    ///   ACGT
    ///   ---T
    /// Score = 3×(−1) + 1×(+1) = −2
    ///
    /// Standard NW initialization: F(3,0)=−3, so F(4,1)=max(F(3,0)+1,...)=−2.
    ///
    /// Source: F(i,0) = d·i from Wikipedia pseudocode.
    /// </summary>
    [Test]
    public void UnequalLengths_LongVsShort_CorrectLinearBorderScore()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("T");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(-2),
                "ACGT vs T: optimal alignment ACGT/---T scores 3×(−1)+1×(+1)=−2");

            Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("ACGT"));
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo("T"));
            Assert.That(result.Score, Is.EqualTo(
                RecalculateLinearScore(result.AlignedSequence1, result.AlignedSequence2, WikipediaScoring)));
        });
    }

    /// <summary>
    /// Verifies that a purely-gap alignment (no matches possible) scores correctly.
    ///
    /// Seq1 = "AAAA", Seq2 = "TTTT", match=+1, mismatch=−1, gap=−1.
    ///
    /// All-mismatch alignment: AAAA / TTTT → 4×(−1) = −4.
    /// Any gap introduction costs the same or more. Score = −4.
    ///
    /// Source: Corner case from Wikipedia Sequence Alignment article —
    ///         "Completely different sequences" should produce a negative score
    ///         with negative mismatch/gap penalties.
    /// </summary>
    [Test]
    public void CompletelyDifferent_EqualLength_AllMismatchScore()
    {
        var seq1 = new DnaSequence("AAAA");
        var seq2 = new DnaSequence("TTTT");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            // All mismatches: 4 × (−1) = −4
            Assert.That(result.Score, Is.EqualTo(-4),
                "AAAA vs TTTT: 4 mismatches × (−1) = −4");
            Assert.That(result.AlignedSequence1, Does.Not.Contain("-"),
                "No gaps needed for equal-length all-mismatch");
            Assert.That(result.AlignedSequence2, Does.Not.Contain("-"),
                "No gaps needed for equal-length all-mismatch");
        });
    }

    #endregion

    #region Fundamental Alignment Properties

    /// <summary>
    /// Identical sequences: perfect diagonal alignment with no gaps.
    /// Score = n × Match.
    ///
    /// Source: Trivially follows from NW recurrence — diagonal path with all matches
    ///         dominates any path introducing gaps or mismatches.
    /// </summary>
    [Test]
    public void IdenticalSequences_PerfectScore_NoGaps()
    {
        var seq = new DnaSequence("ACGTACGT");

        var result = SequenceAligner.GlobalAlign(seq, seq, WikipediaScoring);

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(8),
                "8 identical positions × (+1) = 8");
            Assert.That(result.AlignedSequence1, Does.Not.Contain("-"));
            Assert.That(result.AlignedSequence2, Does.Not.Contain("-"));
            Assert.That(result.AlignedSequence1, Is.EqualTo("ACGTACGT"));
            Assert.That(result.AlignedSequence2, Is.EqualTo("ACGTACGT"));
        });
    }

    /// <summary>
    /// Single internal deletion: ACGT vs AGT.
    /// Optimal alignment: ACGT / A-GT → 3 matches + 1 gap = 3×(+1) + 1×(−1) = 2.
    ///
    /// Source: Linear gap penalty — a single gap costs d (here −1).
    /// </summary>
    [Test]
    public void SingleDeletion_CorrectGapPenalty()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("AGT");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(2),
                "ACGT vs AGT: 3 matches − 1 gap = 2");
            Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("ACGT"));
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo("AGT"));
            Assert.That(result.Score, Is.EqualTo(
                RecalculateLinearScore(result.AlignedSequence1, result.AlignedSequence2, WikipediaScoring)));
        });
    }

    /// <summary>
    /// Verifies exact alignment statistics for a single-deletion alignment.
    ///
    /// ACGT vs AGT → aligned: ACGT / A-GT (4 positions).
    /// Matches = 3 (A, G, T), Mismatches = 0, Gaps = 1 (position 2).
    /// Identity = 3/4 × 100 = 75.0%.
    /// Similarity = (3+0)/4 × 100 = 75.0% (non-gap positions).
    /// GapPercent = 1/4 × 100 = 25.0%.
    ///
    /// Source: Statistics are defined as per-position counts over the aligned sequences,
    ///         which follow deterministically from the NW traceback.
    /// </summary>
    [Test]
    public void SingleDeletion_Statistics_ExactValues()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("AGT");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaScoring);
        var stats = SequenceAligner.CalculateStatistics(result);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Matches, Is.EqualTo(3), "3 matching positions: A, G, T");
            Assert.That(stats.Mismatches, Is.EqualTo(0), "No mismatches");
            Assert.That(stats.Gaps, Is.EqualTo(1), "1 gap in seq2");
            Assert.That(stats.AlignmentLength, Is.EqualTo(4), "Aligned length = 4");
            Assert.That(stats.Identity, Is.EqualTo(75.0), "3/4 × 100 = 75.0%");
            Assert.That(stats.Similarity, Is.EqualTo(75.0), "3/4 × 100 = 75.0%");
            Assert.That(stats.GapPercent, Is.EqualTo(25.0), "1/4 × 100 = 25.0%");
        });
    }

    /// <summary>
    /// Verifies scoring with a different scoring matrix (higher match reward).
    /// ACGT vs ACGT with match=+5, gap=−1 → 4 × 5 = 20.
    ///
    /// Source: Wikipedia "Basic scoring schemes" — arbitrary match/mismatch/gap values.
    /// </summary>
    [Test]
    public void DifferentScoringMatrix_HighMatchReward_CorrectScore()
    {
        var scoring = new ScoringMatrix(Match: 5, Mismatch: -1, GapOpen: 0, GapExtend: -1);
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("ACGT");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, scoring);

        Assert.That(result.Score, Is.EqualTo(20),
            "4 matches × (+5) = 20");
    }

    /// <summary>
    /// NW score is symmetric: GlobalAlign(A,B).Score == GlobalAlign(B,A).Score.
    /// This follows from the NW recurrence — S(a,b) = S(b,a) by construction
    /// (match/mismatch depends only on equality), and gap penalties are uniform.
    ///
    /// Source: NW recurrence symmetry (Wikipedia "Advanced presentation of algorithm").
    /// </summary>
    [Test]
    public void ScoreSymmetry_ReversedInputs_SameScore()
    {
        var seq1 = new DnaSequence("GCATGCG");
        var seq2 = new DnaSequence("GATTACA");

        var resultAB = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaScoring);
        var resultBA = SequenceAligner.GlobalAlign(seq2, seq1, WikipediaScoring);

        Assert.That(resultAB.Score, Is.EqualTo(resultBA.Score),
            "NW score must be symmetric: score(A,B) == score(B,A)");
    }

    #endregion

    #region API Contract

    [Test]
    public void StringOverload_ProducesSameResultAsDnaSequenceOverload()
    {
        const string seq1 = "GCATGCG";
        const string seq2 = "GATTACA";

        var dnaResult = SequenceAligner.GlobalAlign(
            new DnaSequence(seq1), new DnaSequence(seq2), WikipediaScoring);
        var stringResult = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            Assert.That(stringResult.AlignmentType, Is.EqualTo(AlignmentType.Global));
            Assert.That(stringResult.Score, Is.EqualTo(dnaResult.Score));
            Assert.That(stringResult.AlignedSequence1, Is.EqualTo(dnaResult.AlignedSequence1));
            Assert.That(stringResult.AlignedSequence2, Is.EqualTo(dnaResult.AlignedSequence2));
        });
    }

    [Test]
    public void EmptyStringInput_ReturnsEmptyResult()
    {
        var result = SequenceAligner.GlobalAlign(string.Empty, "GATTACA", WikipediaScoring);

        Assert.That(result, Is.EqualTo(AlignmentResult.Empty));
    }

    [Test]
    public void NullDnaSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.GlobalAlign((DnaSequence)null!, new DnaSequence("GATTACA"), WikipediaScoring));
    }

    #endregion

    #region Helpers

    private static string RemoveGaps(string alignedSequence)
    {
        return new string(alignedSequence.Where(c => c != '-').ToArray());
    }

    /// <summary>
    /// Recalculates alignment score using the linear gap model:
    /// each gap position costs GapExtend (= d), each match costs Match,
    /// each mismatch costs Mismatch.
    ///
    /// This matches the NW scoring definition from Wikipedia:
    /// "The score of the whole alignment candidate is the sum of the scores of all the pairings."
    /// </summary>
    private static int RecalculateLinearScore(string aligned1, string aligned2, ScoringMatrix scoring)
    {
        int score = 0;

        for (int i = 0; i < aligned1.Length; i++)
        {
            char a = aligned1[i];
            char b = aligned2[i];

            if (a == '-' || b == '-')
                score += scoring.GapExtend;
            else if (a == b)
                score += scoring.Match;
            else
                score += scoring.Mismatch;
        }

        return score;
    }

    #endregion
}
