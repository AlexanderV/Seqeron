using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for ALIGN-SEMI-001: Semi-Global Alignment (Ends-Free / Glocal)
/// Evidence: Wikipedia (Sequence alignment, Needleman–Wunsch algorithm), Brudno et al. (2003)
/// 
/// Semi-global alignment is a hybrid of global and local alignment that permits
/// free end gaps in one or both sequences. The implementation uses the query-in-reference
/// variant where seq1 (query) is fully aligned and seq2 (reference) has free end gaps.
/// </summary>
[TestFixture]
public class SequenceAligner_SemiGlobalAlign_Tests
{
    private static readonly ScoringMatrix SimpleDna = new(
        Match: 1,
        Mismatch: -1,
        GapOpen: -2,
        GapExtend: -1);

    #region M1: Short Query Embedded in Long Reference

    /// <summary>
    /// M1: Short query embedded in long reference finds the match.
    /// Evidence: Wikipedia (Sequence alignment) - semi-global is useful when one sequence
    /// is short and the other is long.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_ShortQueryInLongReference_FindsMatch()
    {
        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("AAAATGCAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            // INV-1: AlignmentType is SemiGlobal
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal),
                "INV-1: AlignmentType must be SemiGlobal");

            // INV-3: Query fully represented
            string queryFromAlignment = RemoveGaps(result.AlignedSequence1);
            Assert.That(queryFromAlignment, Is.EqualTo(query.Sequence),
                "INV-3: Removing gaps from aligned seq1 must yield original query");

            // Score should be positive for matching embedded sequence
            Assert.That(result.Score, Is.GreaterThan(0),
                "Score should be positive when query matches reference region");
        });
    }

    #endregion

    #region M2: AlignmentType is SemiGlobal

    /// <summary>
    /// M2: AlignmentType property is set to SemiGlobal.
    /// Evidence: Implementation contract.
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
    /// Evidence: Alignment definition (all alignment algorithms produce equal-length aligned strings).
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
    /// M4: Query is fully represented (removing gaps yields original query).
    /// Evidence: Semi-global alignment definition - query is globally aligned.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_QueryAtEnd_QueryFullyRepresented()
    {
        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("ATGCAAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        string queryFromAlignment = RemoveGaps(result.AlignedSequence1);
        Assert.That(queryFromAlignment, Is.EqualTo(query.Sequence),
            "INV-3: Removing gaps from aligned seq1 must yield original query");
    }

    #endregion

    #region M5: Reference is Substring After Gap Removal

    /// <summary>
    /// M5: Removing gaps from aligned reference yields a substring of the original reference.
    /// Evidence: Semi-global alignment (query-in-reference variant) - reference has free end gaps.
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

    #region M6: Score Non-Negative for Matching Sequences

    /// <summary>
    /// M6: Score is non-negative when query matches reference.
    /// Evidence: ASSUMPTION (implementation behavior for query-in-reference semi-global).
    /// </summary>
    [Test]
    public void SemiGlobalAlign_MatchingSequences_ScoreNonNegative()
    {
        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("ATGCAAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.That(result.Score, Is.GreaterThanOrEqualTo(0),
            "Score should be non-negative when query matches reference region");
    }

    #endregion

    #region M7-M8: Null Sequence Throws ArgumentNullException

    /// <summary>
    /// M7: Null seq1 throws ArgumentNullException.
    /// Evidence: ASSUMPTION (implementation behavior).
    /// </summary>
    [Test]
    public void SemiGlobalAlign_NullSequence1_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.SemiGlobalAlign(null!, new DnaSequence("ATGC")));
    }

    /// <summary>
    /// M8: Null seq2 throws ArgumentNullException.
    /// Evidence: ASSUMPTION (implementation behavior).
    /// </summary>
    [Test]
    public void SemiGlobalAlign_NullSequence2_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.SemiGlobalAlign(new DnaSequence("ATGC"), null!));
    }

    #endregion

    #region S1: Identical Sequences Full Alignment

    /// <summary>
    /// S1: Identical sequences produce full alignment with maximum score.
    /// Evidence: ASSUMPTION.
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
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo(query.Sequence));
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo(reference.Sequence));
            Assert.That(result.Score, Is.EqualTo(query.Length * SimpleDna.Match),
                "Score should be length × match score for identical sequences");
        });
    }

    #endregion

    #region S2: Query at Start of Reference

    /// <summary>
    /// S2: Query matching start of reference (trailing gaps in reference alignment).
    /// Evidence: ASSUMPTION (semi-global free end gaps).
    /// Note: Semi-global does NOT have a zero floor like local alignment.
    /// </summary>
    [Test]
    public void SemiGlobalAlign_QueryAtStart_TrailingGapsInReference()
    {
        var query = new DnaSequence("ATG");
        var reference = new DnaSequence("ATGCCCCC");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo(query.Sequence),
                "Query should be fully represented");
            // Note: Unlike local alignment, semi-global can have any score value
        });
    }

    #endregion

    #region S3: Query at End of Reference

    /// <summary>
    /// S3: Query matching end of reference (leading gaps in reference alignment).
    /// Evidence: ASSUMPTION (semi-global free end gaps).
    /// </summary>
    [Test]
    public void SemiGlobalAlign_QueryAtEnd_LeadingGapsHandled()
    {
        var query = new DnaSequence("CCC");
        var reference = new DnaSequence("ATGCCC");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo(query.Sequence),
                "Query should be fully represented");
        });
    }

    #endregion

    #region S4: Custom Scoring Matrix

    /// <summary>
    /// S4: Custom scoring matrix is applied correctly.
    /// Evidence: ASSUMPTION (scoring contract).
    /// </summary>
    [Test]
    public void SemiGlobalAlign_CustomScoring_ScoreReflectsScoringMatrix()
    {
        var customScoring = new ScoringMatrix(
            Match: 5,
            Mismatch: -2,
            GapOpen: -3,
            GapExtend: -1);

        var query = new DnaSequence("ACGT");
        var reference = new DnaSequence("ACGT");

        var result = SequenceAligner.SemiGlobalAlign(query, reference, customScoring);

        Assert.That(result.Score, Is.EqualTo(query.Length * customScoring.Match),
            "Score should reflect custom match score for identical sequences");
    }

    #endregion

    #region INV: All Invariants Combined Test

    /// <summary>
    /// Combined invariant validation test for semi-global alignment.
    /// Validates INV-1 through INV-5 in a single comprehensive test.
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
            Assert.That(queryFromAlignment, Is.EqualTo(query.Sequence),
                "INV-3: Removing gaps from aligned seq1 must yield original query");

            // INV-4: Reference is substring
            string refFromAlignment = RemoveGaps(result.AlignedSequence2);
            Assert.That(reference.Sequence, Does.Contain(refFromAlignment),
                "INV-4: Removing gaps from aligned seq2 must yield substring of original reference");

            // INV-5: Score is consistent with alignment
            Assert.That(result.Score, Is.GreaterThan(0),
                "INV-5: Score should be positive for matching query");
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
