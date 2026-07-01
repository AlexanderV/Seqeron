using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// Canonical tests for SequenceAligner.MultipleAlignProgressive() — guide-tree
/// progressive (Feng-Doolittle / Clustal-style) multiple sequence alignment.
/// Test Unit: ALIGN-MULTI-001 (deferred enhancement C4).
///
/// Algorithm under test (additive second aligner; star MultipleAlign() is unchanged):
/// 1. All pairwise Needleman-Wunsch alignments → distance matrix
///    (distance = 1 − fractional identity = 1 − identical columns / pairwise alignment length).
/// 2. UPGMA guide tree over that matrix (the classic Feng-Doolittle choice).
/// 3. Progressive profile-profile alignment from the tips using the NW recurrence over columns
///    with sum-of-pairs profile scoring and the "once a gap, always a gap" rule.
///
/// Sources (retrieved this session):
/// - Feng &amp; Doolittle (1987) "Progressive sequence alignment as a prerequisite to correct
///   phylogenetic trees", J Mol Evol 25:351-360 (PubMed 3118049) — progressive method, guide
///   tree, "Once a gap is introduced … it is preserved within all subsequent fusions".
/// - Wikipedia: Multiple sequence alignment (https://en.wikipedia.org/wiki/Multiple_sequence_alignment)
///   — guide tree built by NJ/UPGMA; sequences added sequentially along the guide tree.
/// - Wikipedia: UPGMA / Sokal &amp; Michener (1958) — cluster-averaging formula.
///
/// Expected aligned rows and SP scores below are hand-derived from the algorithm definition
/// (not echoed from the implementation); see per-test comments.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MultipleAlignProgressive_Tests
{
    #region MUST — Error / edge cases (0, 1, 2 sequences)

    /// <summary>M01: Null input throws ArgumentNullException (.NET convention, mirrors star MSA).</summary>
    [Test]
    public void MultipleAlignProgressive_NullInput_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.MultipleAlignProgressive(null!));
    }

    /// <summary>M02: Empty collection returns the Empty result.</summary>
    [Test]
    public void MultipleAlignProgressive_EmptyCollection_ReturnsEmpty()
    {
        var result = SequenceAligner.MultipleAlignProgressive(Array.Empty<DnaSequence>());

        Assert.Multiple(() =>
        {
            Assert.That(result, Is.EqualTo(MultipleAlignmentResult.Empty));
            Assert.That(result.AlignedSequences, Is.Empty);
            Assert.That(result.Consensus, Is.Empty);
            Assert.That(result.TotalScore, Is.EqualTo(0));
        });
    }

    /// <summary>M03: Single sequence is returned verbatim (no pairs, no gaps, SP = 0).</summary>
    [Test]
    public void MultipleAlignProgressive_SingleSequence_ReturnsSameSequence()
    {
        var result = SequenceAligner.MultipleAlignProgressive(new[] { new DnaSequence("ATGCATGC") });

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(1));
            Assert.That(result.AlignedSequences[0], Is.EqualTo("ATGCATGC"));
            Assert.That(result.Consensus, Is.EqualTo("ATGCATGC"));
            Assert.That(result.TotalScore, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// M04: Two identical sequences → both rows identical, no gaps, SP = 4.
    /// Hand-derived: NW("ATGC","ATGC") = "ATGC"/"ATGC"; 4 match columns × Match(1) × 1 pair = 4.
    /// </summary>
    [Test]
    public void MultipleAlignProgressive_TwoIdenticalSequences_ExactRows()
    {
        var result = SequenceAligner.MultipleAlignProgressive(new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC")
        });

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences, Is.EqualTo(new[] { "ATGC", "ATGC" }));
            Assert.That(result.Consensus, Is.EqualTo("ATGC"));
            Assert.That(result.TotalScore, Is.EqualTo(4));
        });
    }

    /// <summary>
    /// M05: Two sequences with one deletion → exact aligned rows with a single gap.
    /// Hand-derived NW("ACGT","ACT") with SimpleDna (match +1, mismatch −1, gap −1):
    ///   "ACGT" / "AC-T". SP: A=A,C=C,T=T → +3, one gap column → GapExtend(−1) ⇒ SP = 2.
    /// Removing gaps recovers "ACGT" and "ACT".
    /// </summary>
    [Test]
    public void MultipleAlignProgressive_TwoSequencesWithDeletion_ExactRows()
    {
        var result = SequenceAligner.MultipleAlignProgressive(new[]
        {
            new DnaSequence("ACGT"),
            new DnaSequence("ACT")
        });

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences, Is.EqualTo(new[] { "ACGT", "AC-T" }));
            Assert.That(result.TotalScore, Is.EqualTo(2));
            Assert.That(result.AlignedSequences[0].Replace("-", ""), Is.EqualTo("ACGT"));
            Assert.That(result.AlignedSequences[1].Replace("-", ""), Is.EqualTo("ACT"));
        });
    }

    #endregion

    #region MUST — Identical sequences invariant (a)

    /// <summary>
    /// M06: Identical sequences → all rows identical, no gaps. (Invariant (a).)
    /// SP = C(3,2)=3 pairs × 8 match columns × Match(1) = 24.
    /// </summary>
    [Test]
    public void MultipleAlignProgressive_ThreeIdenticalSequences_AllRowsIdenticalNoGaps()
    {
        var result = SequenceAligner.MultipleAlignProgressive(new[]
        {
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGCATGC")
        });

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences,
                Is.EqualTo(new[] { "ATGCATGC", "ATGCATGC", "ATGCATGC" }));
            Assert.That(result.AlignedSequences.Any(s => s.Contains('-')), Is.False,
                "Identical sequences must not introduce gaps");
            Assert.That(result.Consensus, Is.EqualTo("ATGCATGC"));
            Assert.That(result.TotalScore, Is.EqualTo(24));
        });
    }

    #endregion

    #region MUST — Hand-derivable 3-sequence progressive alignment (b)

    /// <summary>
    /// M07: 3-sequence case with a unique, hand-derivable progressive answer (invariant (b)).
    /// Input ["ACGT","ACGT","AGT"].
    /// Guide tree: pairwise identity distances — d(0,1)=0 (closest), so seq0+seq1 merge first
    /// into profile [ACGT / ACGT]; then "AGT" is aligned to that profile.
    /// Aligning column profile [A][C][G][T] to "AGT": the optimal places the gap opposite the
    /// 'C' column ⇒ "A-GT". Exact rows: ["ACGT","ACGT","A-GT"].
    /// SP (SimpleDna): col0 AAA→+3; col1 C,C,-→ (C=C)+1,(C,-)−1,(C,-)−1 = −1; col2 GGG→+3;
    /// col3 TTT→+3 ⇒ SP = 8.
    /// </summary>
    [Test]
    public void MultipleAlignProgressive_ThreeSequencesOneDeletion_ExactColumns()
    {
        var result = SequenceAligner.MultipleAlignProgressive(new[]
        {
            new DnaSequence("ACGT"),
            new DnaSequence("ACGT"),
            new DnaSequence("AGT")
        });

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences,
                Is.EqualTo(new[] { "ACGT", "ACGT", "A-GT" }));
            Assert.That(result.Consensus, Is.EqualTo("ACGT"));
            Assert.That(result.TotalScore, Is.EqualTo(8));
        });
    }

    #endregion

    #region MUST — Discriminating: progressive ≠ star (c)

    /// <summary>
    /// M08: DISCRIMINATING case — progressive gives the correct gap-free alignment while the
    /// star aligner inserts a spurious gap (invariant (c)). This test would FAIL against the
    /// star aligner, so it proves the two methods differ.
    ///
    /// Input ["AAGAA","AACAA","GGTGG","GGTGG"] — two tight clusters:
    ///   {AAGAA, AACAA} differ at one position; {GGTGG, GGTGG} are identical.
    /// Pairwise identity distances: d(2,3)=0, d(0,1)=0.2, all cross-cluster distances ≥ 0.8.
    /// UPGMA guide tree merges (2,3) then (0,1) then the two clusters. Each cluster aligns
    /// internally with no gaps (same length, no indel improves the diagonal), and the
    /// profile-profile step keeps the ungapped diagonal as optimal.
    ///
    /// EXACT progressive columns (length 5, no gaps):
    ///   AAGAA
    ///   AACAA
    ///   GGTGG
    ///   GGTGG
    /// SP (SimpleDna): every column has 1 mismatch among the 6 pairs except the all-different
    /// middle column. Hand sum over the 5 columns = −12.
    ///
    /// The star aligner on the SAME input produces a length-6 alignment with a leading gap
    /// (rows "AAG-AA","-AACAA","-GGTGG","-GGTGG", SP −13); asserting the gap-free length-5
    /// rows here therefore discriminates the two algorithms.
    /// </summary>
    [Test]
    public void MultipleAlignProgressive_TwoClusters_GivesGapFreeAlignment_DiscriminatesFromStar()
    {
        var input = new[]
        {
            new DnaSequence("AAGAA"),
            new DnaSequence("AACAA"),
            new DnaSequence("GGTGG"),
            new DnaSequence("GGTGG")
        };

        var progressive = SequenceAligner.MultipleAlignProgressive(input);
        var star = SequenceAligner.MultipleAlign(input);

        Assert.Multiple(() =>
        {
            // Exact progressive answer (gap-free, length 5).
            Assert.That(progressive.AlignedSequences,
                Is.EqualTo(new[] { "AAGAA", "AACAA", "GGTGG", "GGTGG" }));
            Assert.That(progressive.AlignedSequences.All(s => s.Length == 5), Is.True);
            Assert.That(progressive.AlignedSequences.Any(s => s.Contains('-')), Is.False,
                "Progressive alignment of this input needs no gaps");
            Assert.That(progressive.TotalScore, Is.EqualTo(-12));

            // The star aligner produces a DIFFERENT (gapped, longer) alignment — this is what
            // makes the assertions above discriminating.
            Assert.That(star.AlignedSequences, Is.Not.EqualTo(progressive.AlignedSequences),
                "Star and progressive must differ on this two-cluster input");
            Assert.That(star.AlignedSequences[0].Length, Is.EqualTo(6),
                "Star inserts a spurious gap, yielding length 6");
        });
    }

    #endregion

    #region MUST — Invariants (d)

    /// <summary>
    /// M09: All output rows equal length; removing gaps recovers each input; column count ≥
    /// max input length; row count preserved (invariant (d)).
    /// </summary>
    [Test]
    public void MultipleAlignProgressive_Invariants_Hold()
    {
        var originals = new[] { "ATGCATGCATGC", "ATGC", "ATGCAT", "ATGCATGCATGCAA" };
        var result = SequenceAligner.MultipleAlignProgressive(
            originals.Select(s => new DnaSequence(s)));

        int len = result.AlignedSequences[0].Length;
        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(originals.Length),
                "Row count preserved");
            Assert.That(result.AlignedSequences.All(s => s.Length == len), Is.True,
                "All rows equal length");
            Assert.That(len, Is.GreaterThanOrEqualTo(originals.Max(s => s.Length)),
                "Column count ≥ max input length");
            for (int i = 0; i < originals.Length; i++)
                Assert.That(result.AlignedSequences[i].Replace("-", ""), Is.EqualTo(originals[i]),
                    $"Removing gaps recovers original sequence {i}");
            Assert.That(result.Consensus.Length, Is.EqualTo(len),
                "Consensus length = aligned length");
        });
    }

    /// <summary>
    /// M10: No alignment column is entirely gaps (Wikipedia MSA invariant).
    /// </summary>
    [Test]
    public void MultipleAlignProgressive_NoColumnIsAllGaps()
    {
        var result = SequenceAligner.MultipleAlignProgressive(new[]
        {
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("ATGCAA"),
            new DnaSequence("ATGCATGCGG")
        });

        int length = result.AlignedSequences[0].Length;
        for (int col = 0; col < length; col++)
        {
            bool allGaps = result.AlignedSequences.All(s => s[col] == '-');
            Assert.That(allGaps, Is.False, $"Column {col} is entirely gaps");
        }
    }

    /// <summary>
    /// M11: Custom scoring matrix is honored (BlastDna doubles the per-match contribution).
    /// SimpleDna: 1 pair × 4 matches × 1 = 4. BlastDna: 1 pair × 4 matches × 2 = 8.
    /// </summary>
    [Test]
    public void MultipleAlignProgressive_WithCustomScoring_UsesProvidedMatrix()
    {
        var input = new[] { new DnaSequence("ATGC"), new DnaSequence("ATGC") };

        var simple = SequenceAligner.MultipleAlignProgressive(input);
        var blast = SequenceAligner.MultipleAlignProgressive(input, SequenceAligner.BlastDna);

        Assert.Multiple(() =>
        {
            Assert.That(simple.TotalScore, Is.EqualTo(4));
            Assert.That(blast.TotalScore, Is.EqualTo(8));
        });
    }

    /// <summary>
    /// M12: PROFILE-LEVEL "once a gap, always a gap". Two gapped profiles are merged and each
    /// profile's gap must be carried verbatim into the final alignment (never filled with a
    /// residue, never edited). Independently hand-derived (validator session).
    ///
    /// Input ["ACGT","ACGT","AGT","AGT"].
    /// Distances: d(0,1)=0, d(2,3)=0 (both minimal); deterministic lowest-index tie-break merges
    /// (0,1) first, then (2,3), then the two profiles. Profile P1=[ACGT/ACGT] (width 4) aligns to
    /// P2=[AGT/AGT] (width 3): the optimal places one new all-gap column in P2 opposite the 'C'
    /// column of P1. EXACT columns (length 4):
    ///   ACGT
    ///   ACGT
    ///   A-GT
    ///   A-GT
    /// SP (SimpleDna, 6 pairs/col): col0 AAAA=+6; col1 {C,C,-,-} = (C,C)+1 + four (C,-)·−1 + (-,-)0
    ///   = −3; col2 GGGG=+6; col3 TTTT=+6 ⇒ SP = 15. Gap-removal recovers every input.
    /// </summary>
    [Test]
    public void MultipleAlignProgressive_TwoGappedProfilesMerge_GapsCarriedVerbatim_ExactColumns()
    {
        var originals = new[] { "ACGT", "ACGT", "AGT", "AGT" };
        var result = SequenceAligner.MultipleAlignProgressive(
            originals.Select(s => new DnaSequence(s)));

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences,
                Is.EqualTo(new[] { "ACGT", "ACGT", "A-GT", "A-GT" }));
            Assert.That(result.TotalScore, Is.EqualTo(15));
            // Once-a-gap: removing gaps recovers each input exactly (no gap was ever residue-filled).
            for (int i = 0; i < originals.Length; i++)
                Assert.That(result.AlignedSequences[i].Replace("-", ""), Is.EqualTo(originals[i]));
        });
    }

    #endregion

    #region SHOULD — Additivity / non-breaking guarantee

    /// <summary>
    /// S01: Adding the progressive aligner does not change the star aligner. The existing star
    /// MultipleAlign() must still produce its documented length-6 gapped result on the
    /// discriminating input (byte-for-byte unchanged).
    /// </summary>
    [Test]
    public void MultipleAlign_StarBehaviour_Unchanged_OnDiscriminatingInput()
    {
        var star = SequenceAligner.MultipleAlign(new[]
        {
            new DnaSequence("AAGAA"),
            new DnaSequence("AACAA"),
            new DnaSequence("GGTGG"),
            new DnaSequence("GGTGG")
        });

        Assert.That(star.AlignedSequences,
            Is.EqualTo(new[] { "AAG-AA", "-AACAA", "-GGTGG", "-GGTGG" }),
            "Star MSA output must remain unchanged by the additive progressive aligner");
    }

    #endregion
}
