using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for ALIGN-LOCAL-001: Local Alignment (Smith–Waterman)
/// Evidence: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
///
/// Wikipedia example parameters:
///   Sequences: A = TGTTACGG, B = GGTTGACTA
///   Substitution: s(a,b) = +3 if a==b, −3 if a≠b
///   Linear gap penalty: W₁ = 2 (i.e., gap cost = −2 per position)
///
/// Hand-verified scoring matrix maximum: H[6,7] = 13 (row=C, col=C)
/// Traceback: (6,7)→diag(5,6)→diag(4,5)→left(4,4)→diag(3,3)→diag(2,2)→diag→H[1,1]=0 STOP
/// Result: aligned1="GTT-AC", aligned2="GTTGAC", score=13
/// </summary>
[TestFixture]
public class SequenceAligner_LocalAlign_Tests
{
    // Wikipedia example scoring: match +3, mismatch −3, linear gap penalty W₁=2
    // Source: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Example
    private static readonly ScoringMatrix WikipediaScoring = new(
        Match: 3,
        Mismatch: -3,
        GapOpen: -2,
        GapExtend: -2);

    #region M1: Wikipedia Example - Exact verification against source

    /// <summary>
    /// M1: Wikipedia example with sequences TGTTACGG and GGTTGACTA.
    /// Evidence: Smith–Waterman algorithm, Wikipedia "Example" section.
    ///   https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Example
    ///
    /// Expected result (from Wikipedia):
    ///   Alignment: "GTT-AC" / "GTTGAC", Score = 13
    ///
    /// Hand-verified scoring matrix (8×9, linear gap W₁=2):
    ///        -  G  G  T  T  G  A  C  T  A
    ///   -  [ 0  0  0  0  0  0  0  0  0  0 ]
    ///   T  [ 0  0  0  3  3  1  0  0  3  1 ]
    ///   G  [ 0  3  3  1  1  6  4  2  1  0 ]
    ///   T  [ 0  1  1  6  4  4  3  1  5  3 ]
    ///   T  [ 0  0  0  4  9  7  5  3  4  2 ]
    ///   A  [ 0  0  0  2  7  6 10  8  6  7 ]
    ///   C  [ 0  0  0  0  5  4  8 13 11  9 ]
    ///   G  [ 0  3  3  1  3  8  6 11 10  8 ]
    ///   G  [ 0  3  6  4  2  6  5  9  8  7 ]
    ///
    /// Maximum = 13 at H[6,7]. Traceback yields GTT-AC / GTTGAC.
    /// Validates INV-1 through INV-5 plus exact expected values.
    /// Traceback uses "left" branch at H[4,5]→H[4,4] (gap in seq1).
    /// </summary>
    [Test]
    public void LocalAlign_WikipediaExample_ValidatesAllInvariants()
    {
        var seq1 = new DnaSequence("TGTTACGG");
        var seq2 = new DnaSequence("GGTTGACTA");

        var result = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            // Exact expected values from Wikipedia example
            Assert.That(result.Score, Is.EqualTo(13),
                "Wikipedia example score: max H[i,j] = 13 at H[6,7]");
            Assert.That(result.AlignedSequence1, Is.EqualTo("GTT-AC"),
                "Wikipedia example aligned seq1");
            Assert.That(result.AlignedSequence2, Is.EqualTo("GTTGAC"),
                "Wikipedia example aligned seq2");

            // INV-1: Score ≥ 0 (zero floor property)
            Assert.That(result.Score, Is.GreaterThanOrEqualTo(0),
                "INV-1: Local alignment score must be non-negative (zero floor)");

            // INV-3: AlignmentType is Local
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.Local),
                "INV-3: AlignmentType must be Local");

            // INV-4: Removing gaps from aligned sequences yields substrings of originals
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("GTTAC"),
                "INV-4: Aligned seq1 with gaps removed must equal GTTAC");
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo("GTTGAC"),
                "INV-4: Aligned seq2 with gaps removed must equal GTTGAC");
            Assert.That(seq1.Sequence, Does.Contain(RemoveGaps(result.AlignedSequence1)),
                "INV-4: GTTAC must be substring of TGTTACGG");
            Assert.That(seq2.Sequence, Does.Contain(RemoveGaps(result.AlignedSequence2)),
                "INV-4: GTTGAC must be substring of GGTTGACTA");

            // INV-5: Exact positions from traceback
            Assert.That(result.StartPosition1, Is.EqualTo(1),
                "INV-5: StartPosition1 = 1 (0-indexed: seq1[1..5] = GTTAC)");
            Assert.That(result.EndPosition1, Is.EqualTo(5),
                "INV-5: EndPosition1 = 5");
            Assert.That(result.StartPosition2, Is.EqualTo(1),
                "INV-5: StartPosition2 = 1 (0-indexed: seq2[1..6] = GTTGAC)");
            Assert.That(result.EndPosition2, Is.EqualTo(6),
                "INV-5: EndPosition2 = 6");
        });
    }

    #endregion

    #region M2: Swapped Wikipedia Example - Gap in seq2 (up-branch)

    /// <summary>
    /// M2: Wikipedia sequences swapped: seq1=GGTTGACTA, seq2=TGTTACGG.
    /// Evidence: Derived from Wikipedia §Example via S-W transposition property.
    ///   H^{AB}[i,j] = H^{BA}[j,i] — transposed scoring matrix has same maximum.
    ///   Score = 13 (same as M1).
    ///
    /// Hand-verified transposed matrix H^T (rows=GGTTGACTA, cols=TGTTACGG):
    ///        -   T   G   T   T   A   C   G   G
    ///   -  [ 0   0   0   0   0   0   0   0   0 ]
    ///   G  [ 0   0   3   1   0   0   0   3   3 ]
    ///   G  [ 0   0   3   1   0   0   0   3   6 ]
    ///   T  [ 0   3   1   6   4   2   0   1   4 ]
    ///   T  [ 0   3   1   4   9   7   5   3   2 ]
    ///   G  [ 0   1   6   4   7   6   4   8   6 ]
    ///   A  [ 0   0   4   3   5  10   8   6   5 ]
    ///   C  [ 0   0   2   1   3   8  13  11   9 ]
    ///   T  [ 0   3   1   5   4   6  11  10   8 ]
    ///   A  [ 0   1   0   3   2   7   9   8   7 ]
    ///
    /// Maximum = 13 at H^T[7,6]. Traceback:
    ///   (7,6)→diag(6,5)→diag(5,4)→UP(4,4)→diag(3,3)→diag(2,2)→diag→H[1,1]=0 STOP
    ///   The "up" step at (5,4): score[5,4]=7 = score[4,4]+gap = 9−2 = 7.
    ///   This produces the gap in seq2 (char from seq1, dash in seq2).
    ///
    /// Result: aligned1="GTTGAC", aligned2="GTT-AC", score=13
    ///   start1=1, end1=6 (seq1[1..6] = GTTGAC from GGTTGACTA)
    ///   start2=1, end2=5 (seq2[1..5] = GTTAC from TGTTACGG)
    /// </summary>
    [Test]
    public void LocalAlign_SwappedWikipediaExample_GapInSeq2()
    {
        var seq1 = new DnaSequence("GGTTGACTA");
        var seq2 = new DnaSequence("TGTTACGG");

        var result = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            // Exact expected values from transposed Wikipedia matrix
            Assert.That(result.Score, Is.EqualTo(13),
                "Transposed Wikipedia example: score = 13 at H^T[7,6]");
            Assert.That(result.AlignedSequence1, Is.EqualTo("GTTGAC"),
                "Swapped: aligned seq1 has no gap");
            Assert.That(result.AlignedSequence2, Is.EqualTo("GTT-AC"),
                "Swapped: gap is in seq2 (up-branch in traceback)");

            // INV-3: AlignmentType is Local
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.Local));

            // INV-4: Gap removal yields substrings
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("GTTGAC"));
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo("GTTAC"));
            Assert.That(seq1.Sequence, Does.Contain(RemoveGaps(result.AlignedSequence1)));
            Assert.That(seq2.Sequence, Does.Contain(RemoveGaps(result.AlignedSequence2)));

            // INV-5: Exact positions from traceback
            Assert.That(result.StartPosition1, Is.EqualTo(1),
                "INV-5: seq1[1..6] = GTTGAC");
            Assert.That(result.EndPosition1, Is.EqualTo(6));
            Assert.That(result.StartPosition2, Is.EqualTo(1),
                "INV-5: seq2[1..5] = GTTAC");
            Assert.That(result.EndPosition2, Is.EqualTo(5));
        });
    }

    #endregion

    #region M3: String Overload Parity

    /// <summary>
    /// M3: String overload returns same result as DnaSequence overload.
    /// Evidence: Implementation contract — string overload delegates to same
    /// LocalAlignCore after ToUpperInvariant(), verified by code inspection.
    /// Validated here against the Wikipedia example expected values.
    /// </summary>
    [Test]
    public void LocalAlign_StringOverload_MatchesDnaSequenceResult()
    {
        const string seq1 = "TGTTACGG";
        const string seq2 = "GGTTGACTA";

        var stringResult = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        // Verify against exact Wikipedia expected values (same as M1)
        Assert.Multiple(() =>
        {
            Assert.That(stringResult.AlignmentType, Is.EqualTo(AlignmentType.Local));
            Assert.That(stringResult.Score, Is.EqualTo(13));
            Assert.That(stringResult.AlignedSequence1, Is.EqualTo("GTT-AC"));
            Assert.That(stringResult.AlignedSequence2, Is.EqualTo("GTTGAC"));
        });
    }

    #endregion

    #region M4: Empty Input Returns Empty

    /// <summary>
    /// M4: Empty string input returns AlignmentResult.Empty.
    /// Evidence: Smith–Waterman definition — when one sequence has length 0,
    /// the scoring matrix is (n+1)×1 or 1×(m+1) with only the initialized
    /// first row/column (all zeros). No cells to fill → maxScore = 0,
    /// traceback produces empty alignment.
    /// </summary>
    [Test]
    public void LocalAlign_StringOverload_EmptyInput_ReturnsEmpty()
    {
        var result = SequenceAligner.LocalAlign(string.Empty, "ACGT", WikipediaScoring);

        Assert.That(result, Is.EqualTo(AlignmentResult.Empty));
    }

    #endregion

    #region M5: Null Throws ArgumentNullException

    /// <summary>
    /// M5: Null DnaSequence throws ArgumentNullException.
    /// Evidence: .NET convention — ArgumentNullException.ThrowIfNull is the
    /// standard guard for null reference parameters in public API methods.
    /// </summary>
    [Test]
    public void LocalAlign_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.LocalAlign((DnaSequence)null!, new DnaSequence("ACGT")));
    }

    #endregion

    #region S1: Identical Sequences Full Match

    /// <summary>
    /// S1: Identical sequences produce full-length alignment.
    /// Evidence: Derivable from Smith–Waterman recurrence.
    ///   For identical sequences of length n with match score M, mismatch P, gap G
    ///   where M > 0 and M > |G|:
    ///   - H[i,i] = i × M for all i (diagonal always dominates)
    ///   - H[i,j] ≤ H[min(i,j), min(i,j)] for off-diagonal (mismatches/gaps reduce score)
    ///   - Maximum at H[n,n] = n × M
    ///   - Traceback follows full diagonal → complete alignment with no gaps
    ///
    ///   With WikipediaScoring (M=3): score = 8 × 3 = 24.
    /// </summary>
    [Test]
    public void LocalAlign_IdenticalSequences_FullLengthMatch()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("ACGTACGT");

        var result = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(24),
                "Identical 8-char sequences: score = 8 × match(3) = 24");
            Assert.That(result.AlignedSequence1, Is.EqualTo("ACGTACGT"),
                "Full diagonal alignment — no gaps");
            Assert.That(result.AlignedSequence2, Is.EqualTo("ACGTACGT"),
                "Full diagonal alignment — no gaps");
        });
    }

    #endregion

    #region S2: Dissimilar Sequences - Zero Floor

    /// <summary>
    /// S2: Completely dissimilar sequences produce score exactly 0.
    /// Evidence: Derivable from Smith–Waterman recurrence and zero floor.
    ///   Sequences "AAAA" vs "TTTT" with mismatch=−3, gap=−2:
    ///   - Every diagonal: previous + (−3) ≤ 0 − 3 = −3
    ///   - Every gap: previous + (−2) ≤ 0 − 2 = −2
    ///   - Zero floor: max(0, −3, −2, −2) = 0 for every cell
    ///   - All H[i,j] = 0 → maxScore = 0, alignment is empty.
    /// </summary>
    [Test]
    public void LocalAlign_DissimilarSequences_ScoreIsZero()
    {
        var seq1 = new DnaSequence("AAAA");
        var seq2 = new DnaSequence("TTTT");

        var result = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(0),
                "No matches possible: every H[i,j] = max(0, neg, neg, neg) = 0");
            Assert.That(result.AlignedSequence1, Is.Empty,
                "Score 0 → traceback starts at 0 → empty alignment");
            Assert.That(result.AlignedSequence2, Is.Empty,
                "Score 0 → traceback starts at 0 → empty alignment");
        });
    }

    #endregion

    #region Helper Methods

    private static string RemoveGaps(string alignedSequence)
    {
        return new string(alignedSequence.Where(c => c != '-').ToArray());
    }

    #endregion
}
