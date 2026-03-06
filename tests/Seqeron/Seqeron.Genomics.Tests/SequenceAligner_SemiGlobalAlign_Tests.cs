using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for ALIGN-SEMI-001: Semi-Global Alignment (Fitting / Query-in-Reference)
///
/// This is the query-in-reference variant of semi-global alignment, also known as
/// "fitting alignment" (Rosalind: SIMS). It aligns a short query (seq1) globally against
/// the best-matching substring of a longer reference (seq2).
///
/// Algorithm (modification of Needleman–Wunsch):
///   - First row F(0,j) = 0          → free leading gaps in reference
///   - First column F(i,0) = i * d   → query must be fully aligned
///   - Recurrence: standard NW (no zero floor, unlike Smith–Waterman)
///   - Traceback start: max_j F(m,j) → free trailing gaps in reference
///
/// Evidence:
///   - Wikipedia: Sequence alignment (semi-global / glocal definition)
///   - Wikipedia: Needleman–Wunsch algorithm (recurrence, initialization, traceback)
///   - Rosalind: Finding a Motif with Modifications / SIMS (fitting alignment definition)
///   - Rosalind: Semiglobal Alignment / SMGB (semiglobal glossary)
///   - Brudno et al. (2003), Bioinformatics 19 Suppl 1:i54-62
///
/// All expected scores are hand-computed from DP matrices using the stated scoring.
/// </summary>
[TestFixture]
public class SequenceAligner_SemiGlobalAlign_Tests
{
    // Linear gap cost: match=+1, mismatch=−1, gap=−1 (GapExtend used as linear penalty).
    // GapOpen is unused by the linear-cost NW recurrence.
    private static readonly ScoringMatrix SimpleDna = new(
        Match: 1,
        Mismatch: -1,
        GapOpen: -2,
        GapExtend: -1);

    #region M1: Short Query Embedded in Long Reference — exact score

    /// <summary>
    /// M1: Short query embedded in long reference.
    /// query="ATGC" (len 4), ref="AAAATGCAAA" (len 10).
    /// Perfect match at ref positions 4–7 (1-indexed). Score = 4 (4 matches × +1).
    ///
    /// Hand-computed DP last row (match=1, mismatch=−1, gap=−1):
    ///   F(4,j): [−4, −2, −2, −2, −2, 0, 2, 4, 3, 2, 1]
    ///   max at j=7, score=4.
    ///
    /// Evidence: Wikipedia (Sequence alignment) — semi-global is useful when one sequence
    /// is short and the other is long. Score derived from NW recurrence with first row = 0.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_ShortQueryInLongReference_FindsMatch()
    {
        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("AAAATGCAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal),
                "INV-1: AlignmentType must be SemiGlobal");

            string queryFromAlignment = RemoveGaps(result.AlignedSequence1);
            Assert.That(queryFromAlignment, Is.EqualTo("ATGC"),
                "INV-3: Removing gaps from aligned seq1 must yield original query");

            Assert.That(result.Score, Is.EqualTo(4),
                "Score = 4 matches × +1 (query perfectly embedded in reference)");
        });
    }

    #endregion

    #region M2: AlignmentType is SemiGlobal

    /// <summary>
    /// M2: AlignmentType property is set to SemiGlobal for any valid input.
    /// Evidence: Implementation contract — AlignmentType tag distinguishes alignment variant.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_ValidInput_AlignmentTypeIsSemiGlobal()
    {
        var query = new DnaSequence("ACGT");
        var reference = new DnaSequence("ACGT");

        var result = SequenceAligner.SemiGlobalAlign(query, reference);

        Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal));
    }

    #endregion

    #region M3: Aligned Sequences Have Equal Length

    /// <summary>
    /// M3: Aligned sequences have equal length (gaps inserted to equalize).
    /// Evidence: Fundamental alignment property — all pairwise alignment algorithms produce
    /// equal-length aligned strings. (Wikipedia: Sequence alignment, representation section)
    /// </summary>
    [Test]
    public void SemiGlobalAlign_ValidInput_AlignedSequencesHaveEqualLength()
    {
        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("AAAATGCAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length),
            "INV-2: Aligned sequences must have equal length");
    }

    #endregion

    #region M4: Query Fully Represented

    /// <summary>
    /// M4: Query is fully represented in aligned output (removing gaps yields original query).
    /// In fitting alignment, the query is globally aligned — every base of the query
    /// participates in the alignment.
    ///
    /// query="ATGC", ref="ATGCAAAA". Perfect match at ref[1..4], score=4.
    ///
    /// Evidence: Rosalind SIMS definition — "an alignment of a substring of s against ALL of t."
    /// </summary>
    [Test]
    public void SemiGlobalAlign_QueryFullyRepresented()
    {
        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("ATGCAAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            string queryFromAlignment = RemoveGaps(result.AlignedSequence1);
            Assert.That(queryFromAlignment, Is.EqualTo("ATGC"),
                "INV-3: Removing gaps from aligned seq1 must yield original query");

            Assert.That(result.Score, Is.EqualTo(4),
                "Score = 4 matches × +1 (perfect match at start of reference)");
        });
    }

    #endregion

    #region M5: Reference is Substring After Gap Removal

    /// <summary>
    /// M5: Removing gaps from aligned reference yields a substring of the original reference.
    /// In fitting alignment, only a substring of the reference participates in the scored region.
    ///
    /// Evidence: Rosalind SIMS — fitting alignment is "alignment of a substring of s against all of t."
    /// Wikipedia (Sequence alignment) — semi-global query-in-reference variant.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_ShortQueryInReference_ReferenceIsSubstring()
    {
        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("AAAATGCAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        string refFromAlignment = RemoveGaps(result.AlignedSequence2);
        Assert.That(reference.Sequence, Does.Contain(refFromAlignment),
            "INV-4: Removing gaps from aligned seq2 must yield substring of original reference");
    }

    #endregion

    // M6 removed: was duplicate of M4 (same inputs query="ATGC", ref="ATGCAAAA", same Score=4 assertion).

    #region M7-M8: Null Sequence Throws ArgumentNullException

    /// <summary>
    /// M7: Null seq1 throws ArgumentNullException.
    /// Evidence: .NET ArgumentNullException.ThrowIfNull API convention.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_NullSequence1_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.SemiGlobalAlign(null!, new DnaSequence("ATGC")));
    }

    /// <summary>
    /// M8: Null seq2 throws ArgumentNullException.
    /// Evidence: .NET ArgumentNullException.ThrowIfNull API convention.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_NullSequence2_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.SemiGlobalAlign(new DnaSequence("ATGC"), null!));
    }

    #endregion

    #region S1: Identical Sequences Full Alignment — exact score

    /// <summary>
    /// S1: Identical sequences produce full alignment with score = len × match.
    /// query=ref="ATGCATGC" (len 8). All diagonal matches → score = 8 × 1 = 8.
    /// maxJ = n = 8, so the entire reference is used (no trailing gaps).
    ///
    /// Evidence: NW recurrence — identical sequences yield pure diagonal traceback with
    /// score = sum of match scores. First row = 0 initialization does not affect the
    /// result when sequences are identical, because the diagonal path dominates.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_IdenticalSequences_FullAlignmentMaxScore()
    {
        var query = new DnaSequence("ATGCATGC");
        var reference = new DnaSequence("ATGCATGC");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("ATGCATGC"));
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo("ATGCATGC"));
            Assert.That(result.Score, Is.EqualTo(8),
                "Score = 8 matches × +1 for identical sequences");
        });
    }

    #endregion

    #region S2: Query at Start of Reference — exact score

    /// <summary>
    /// S2: Query matches the start of reference, with unmatched trailing reference.
    /// query="ATG" (len 3), ref="ATGCCCCC" (len 8).
    /// Perfect match at ref[1..3] → score = 3.
    /// Trailing ref positions 4–8 are appended as free end gaps.
    ///
    /// Hand-computed DP last row (match=1, mismatch=−1, gap=−1):
    ///   F(3,j): [−3, −1, 1, 3, 2, 1, 0, −1, −2]
    ///   max at j=3, score=3.
    ///
    /// Evidence: NW recurrence with fitting initialization. Trailing reference is free
    /// per query-in-reference variant definition.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_QueryAtStart_ExactScore()
    {
        var query = new DnaSequence("ATG");
        var reference = new DnaSequence("ATGCCCCC");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("ATG"),
                "Query must be fully represented");
            Assert.That(result.Score, Is.EqualTo(3),
                "Score = 3 matches × +1 (query matches start of reference)");
        });
    }

    #endregion

    #region S3: Query at End of Reference — exact score

    /// <summary>
    /// S3: Query matches the end of reference, with unmatched leading reference.
    /// query="CCC" (len 3), ref="ATGCCC" (len 6).
    /// Perfect match at ref[4..6] → score = 3.
    /// Leading ref positions 1–3 appear as gaps in aligned seq1.
    ///
    /// Hand-computed DP last row (match=1, mismatch=−1, gap=−1):
    ///   F(3,j): [−3, −3, −3, −3, −1, 1, 3]
    ///   max at j=6, score=3.
    ///
    /// Evidence: NW recurrence with fitting initialization. Leading reference gaps are
    /// handled by the traceback continuing left after query is exhausted (i=0, j>0).
    /// </summary>
    [Test]
    public void SemiGlobalAlign_QueryAtEnd_ExactScore()
    {
        var query = new DnaSequence("CCC");
        var reference = new DnaSequence("ATGCCC");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("CCC"),
                "Query must be fully represented");
            Assert.That(result.Score, Is.EqualTo(3),
                "Score = 3 matches × +1 (query matches end of reference)");
        });
    }

    #endregion

    #region S4: Custom Scoring Matrix — fitting context, exact score

    /// <summary>
    /// S4: Custom scoring matrix is applied correctly in a fitting alignment context.
    /// query="ATGC" (len 4), ref="AATGCCC" (len 7).
    /// Match=5, Mismatch=-3, GapExtend=-2. GapOpen unused (linear model).
    ///
    /// Perfect match at ref[1..4] (0-indexed) → score = 4 × 5 = 20.
    ///
    /// Hand-computed DP last row (match=5, mismatch=-3, gapExtend=-2):
    ///   F(4,j): [-8, -1, -1, 4, 11, 20, 18, 16]
    ///   max at j=5, score=20.
    ///
    /// Evidence: NW recurrence with fitting initialization and custom scoring.
    /// Wikipedia (Needleman–Wunsch): "Scoring systems" section.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_CustomScoring_ExactScore()
    {
        var customScoring = new ScoringMatrix(
            Match: 5,
            Mismatch: -3,
            GapOpen: -4,
            GapExtend: -2);

        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("AATGCCC");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, customScoring);

        Assert.Multiple(() =>
        {
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("ATGC"),
                "INV-3: Query fully represented");
            Assert.That(result.Score, Is.EqualTo(20),
                "Score = 4 × 5 = 20 (custom match score in fitting context)");
        });
    }

    #endregion

    #region INV: All Invariants Combined — exact score

    /// <summary>
    /// Combined invariant validation: INV-1 through INV-5.
    /// query="GCATGCG" (len 7), ref="AAAGCATGCGAAA" (len 13).
    /// Perfect match at ref[3..9] (0-indexed) "GCATGCG" → score = 7 × match(+1) = 7.
    ///
    /// Evidence: All invariants derived from fitting alignment definition (Rosalind SIMS)
    /// and NW recurrence (Wikipedia: Needleman–Wunsch algorithm).
    /// </summary>
    [Test]
    public void SemiGlobalAlign_ValidatesAllInvariants()
    {
        var query = new DnaSequence("GCATGCG");
        var reference = new DnaSequence("AAAGCATGCGAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            // INV-1: AlignmentType is SemiGlobal
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal),
                "INV-1: AlignmentType must be SemiGlobal");

            // INV-2: Aligned sequences have equal length
            Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length),
                "INV-2: Aligned sequences must have equal length");

            // INV-3: Query fully represented
            string queryFromAlignment = RemoveGaps(result.AlignedSequence1);
            Assert.That(queryFromAlignment, Is.EqualTo("GCATGCG"),
                "INV-3: Removing gaps from aligned seq1 must yield original query");

            // INV-4: Reference is substring
            string refFromAlignment = RemoveGaps(result.AlignedSequence2);
            Assert.That(reference.Sequence, Does.Contain(refFromAlignment),
                "INV-4: Removing gaps from aligned seq2 must yield substring of original reference");

            // INV-5: Score = max_j F(m,j) = 7 (7 matches × +1, perfect embedding)
            Assert.That(result.Score, Is.EqualTo(7),
                "INV-5: Score = 7 matches × +1 for perfectly embedded query");
        });
    }

    #endregion

    #region Score Verification: Score Can Be Negative — exact value

    /// <summary>
    /// NEG: The fitting alignment score can be negative when the query has no good match
    /// in the reference. Unlike Smith–Waterman (local alignment), there is no zero floor
    /// in the NW recurrence used for fitting alignment.
    ///
    /// query="AAAA" (len 4), ref="CCCC" (len 4). All mismatches, score = 4 × (−1) = −4.
    ///
    /// Hand-computed DP (match=1, mismatch=−1, gap=−1):
    ///   For all-mismatch equal-length sequences, F(i,j) = −i for all j ≥ 1.
    ///   F(4,j): [−4, −4, −4, −4, −4]
    ///   max at j=1 (leftmost tie), score = −4.
    ///
    /// Evidence: The NW recurrence F(i,j) = max(diag, up, left) has no max(0, ...) term.
    /// Wikipedia (Needleman–Wunsch): recurrence section. Wikipedia (Sequence alignment):
    /// "local alignment" uses zero floor (Smith–Waterman), but semi-global / fitting does not.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_AllMismatches_NegativeExactScore()
    {
        var query = new DnaSequence("AAAA");
        var reference = new DnaSequence("CCCC");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.That(result.Score, Is.EqualTo(-4),
            "Score = 4 × (−1) = −4 (all mismatches, no zero floor in NW recurrence)");
    }

    #endregion

    #region Score Verification: Score Reflects Alignment, Not Bottom-Right Cell

    /// <summary>
    /// Validates that the score is max_j F(m,j) (fitting alignment optimal), not F(m,n)
    /// (global alignment bottom-right cell). These differ when the query matches a
    /// substring that doesn't end at the reference end.
    ///
    /// query="ATG", ref="ATGCCC". match at ref[1..3] → score = 3.
    /// F(3,3) = 3 (optimal fitting), F(3,6) = 0 (degraded by trailing mismatches).
    ///
    /// Evidence: Fitting alignment traceback starts from max_j F(m,j), not F(m,n).
    /// Rosalind SIMS; Wikipedia (Sequence alignment, semi-global variant).
    /// </summary>
    [Test]
    public void SemiGlobalAlign_ScoreIsMaxOfLastRow_NotBottomRight()
    {
        var query = new DnaSequence("ATG");
        var reference = new DnaSequence("ATGCCC");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        // F(3,3)=3 is the correct fitting score (3 matches).
        // F(3,6)=0 would be incorrect (degraded by trailing C-vs-G mismatches).
        Assert.That(result.Score, Is.EqualTo(3),
            "Score must be max_j F(m,j) = 3, not F(m,n) = 0");
    }

    #endregion

    #region Score Verification: Manual DP with Mismatch

    /// <summary>
    /// Verifies score for a case with both matches and mismatches.
    /// query="ACG", ref="AACGG". Best alignment at ref[2..4]: A-C-G → 3 matches, score = 3.
    ///
    /// Hand-computed DP (match=1, mismatch=−1, gap=−1):
    ///         ''  A   A   C   G   G
    ///   ''     0  0   0   0   0   0
    ///   A     −1  1   1   0  −1  −1
    ///   C     −2  0   0   2   1   0
    ///   G     −3 −1  −1   1   3   2
    ///   max at j=4, score=3.
    ///
    /// Evidence: NW recurrence with fitting initialization.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_MatchWithOffset_ExactScore()
    {
        var query = new DnaSequence("ACG");
        var reference = new DnaSequence("AACGG");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("ACG"));
            Assert.That(result.Score, Is.EqualTo(3),
                "Score = 3 matches × +1 (ACG matches ref[2..4])");
        });
    }

    #endregion

    #region MIX: Mixed Matches and Mismatches in Optimal Alignment — exact score

    /// <summary>
    /// MIX: Verifies that the score is correct when the optimal fitting alignment
    /// contains both matches AND mismatches (not all-match or all-mismatch).
    ///
    /// query="AGT" (len 3), ref="AAACTAAA" (len 8).
    /// Best alignment at ref[2..4] (0-indexed) = "ACT":
    ///   A↔A(+1), G↔C(−1), T↔T(+1) → score = 1.
    ///
    /// Hand-computed DP (match=1, mismatch=−1, gap=−1):
    ///         ''  A   A   A   C   T   A   A   A
    ///   ''     0  0   0   0   0   0   0   0   0
    ///   A     −1  1   1   1   0  −1   1   1   1
    ///   G     −2  0   0   0   0  −1   0   0   0
    ///   T     −3 −1  −1  −1  −1   1   0  −1  −1
    ///   max at j=5, score=1.
    ///
    /// Evidence: NW recurrence produces both match (+1) and mismatch (−1) in the
    /// scored alignment. This exercises a code path not covered by pure-match tests.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_MixedMatchMismatch_ExactScore()
    {
        var query = new DnaSequence("AGT");
        var reference = new DnaSequence("AAACTAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("AGT"),
                "INV-3: Query fully represented");
            Assert.That(result.Score, Is.EqualTo(1),
                "Score = 2 matches(+1) + 1 mismatch(−1) = 1");
        });
    }

    #endregion

    #region GAP: Gap in Optimal Alignment — exact score

    /// <summary>
    /// GAP: Verifies that the score is correct when the optimal fitting alignment
    /// includes a gap (insertion/deletion). This exercises the gap-penalty branch
    /// of the NW recurrence that pure-match tests never reach.
    ///
    /// query="ACGT" (len 4), ref="AGT" (len 3). The query is longer than the
    /// reference, forcing at least one gap.
    /// Optimal alignment: ACGT / A-GT → A↔A(+1), C↔−(−1), G↔G(+1), T↔T(+1) = 2.
    ///
    /// Hand-computed DP (match=1, mismatch=−1, gap=−1):
    ///         ''  A   G   T
    ///   ''     0  0   0   0
    ///   A     −1  1   0  −1
    ///   C     −2  0   0  −1
    ///   G     −3 −1   1   0
    ///   T     −4 −2   0   2
    ///   max at j=3, score=2.
    ///
    /// Evidence: NW recurrence up-move = gap in reference (seq2) = deletion in query
    /// relative to reference. Wikipedia (Needleman–Wunsch): traceback section.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_GapInAlignment_ExactScore()
    {
        var query = new DnaSequence("ACGT");
        var reference = new DnaSequence("AGT");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo("ACGT"),
                "INV-3: Query fully represented");
            Assert.That(result.AlignedSequence2, Does.Contain("-"),
                "Optimal alignment must contain a gap in the reference");
            Assert.That(result.Score, Is.EqualTo(2),
                "Score = 3 matches(+1) + 1 gap(−1) = 2");
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
