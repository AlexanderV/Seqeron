namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// The shared per-position site / PAM matching rule (B08 round 5 A5-13: three private copies - restriction scan,
/// methylation mask, PAM scan - became <c>IupacPattern</c>). Rows are the reference behaviours, not the code's own table:
/// Biopython 1.88 <c>Bio.Restriction</c> (<c>SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]</c>,
/// <c>EcoRI.search(Seq("AGAANTCA")) == []</c>, <c>HincII.search(Seq("AGTYRACA")) == []</c>) and CRISPOR <c>patMatch</c>
/// (a pattern N matches any sequence character; a sequence N does not satisfy G or R).
/// </summary>
[TestFixture]
public class IupacPattern_Tests
{
    [TestCase('A', 'N', true)]
    [TestCase('T', 'N', true)]
    [TestCase('N', 'N', true, TestName = "a_pattern_N_accepts_a_sequence_N")]
    [TestCase('R', 'N', true, TestName = "a_pattern_N_accepts_any_of_the_15_IUPAC_codes")]
    [TestCase('X', 'N', false, TestName = "a_pattern_N_rejects_a_non_IUPAC_character")]
    [TestCase('-', 'N', false)]
    [TestCase('a', 'N', false, TestName = "lower_case_is_the_callers_job_to_fold")]
    [TestCase('N', 'G', false, TestName = "a_sequence_N_does_not_satisfy_G")]
    [TestCase('N', 'R', false, TestName = "a_sequence_N_does_not_satisfy_a_degenerate_code")]
    [TestCase('R', 'R', false, TestName = "a_sequence_ambiguity_code_blocks_a_degenerate_code")]
    [TestCase('G', 'R', true)]
    [TestCase('C', 'R', false)]
    [TestCase('A', 'V', true)]
    [TestCase('T', 'V', false)]
    public void MatchesCode_FollowsBiopythonAndCrisporPatMatch(char sequenceBase, char patternCode, bool expected)
    {
        Assert.That(IupacPattern.MatchesCode(sequenceBase, patternCode), Is.EqualTo(expected));
    }

    [Test]
    public void MatchesAt_ReferenceRows()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IupacPattern.MatchesAt("GGCCNNNNNGGCC", 0, "GGCCNNNNNGGCC"), Is.True, "SfiI on its own N-spacer");
            Assert.That(IupacPattern.MatchesAt("AGAANTCA", 2, "GAATTC"), Is.False, "EcoRI blocked by a sequence N");
            Assert.That(IupacPattern.MatchesAt("AGTYRACA", 1, "GTYRAC"), Is.False, "HincII: the sequence's own Y/R do not satisfy the pattern's Y/R");
            Assert.That(IupacPattern.MatchesAt("AGTCGACA", 1, "GTYRAC"), Is.True, "HincII on a concrete GTCGAC");
            Assert.That(IupacPattern.MatchesAt("ACGTAGG", 4, "NGG"), Is.True);
            Assert.That(IupacPattern.MatchesAt("ACGTAGG", 5, "NGG"), Is.False, "window past the end does not fit");
            Assert.That(IupacPattern.MatchesAt("ACGTAGG", -1, "NGG"), Is.False, "negative start does not fit");
            Assert.That(IupacPattern.MatchesAt("", 0, "NGG"), Is.False);
        });
    }

    [Test]
    public void MatchesAtCircular_WrapsModuloTheSequenceLength()
    {
        // GAATTC read across the origin of the 8-mer "ATTCAAGA": positions 6,7,0,1,2,3 = G,A,A,T,T,C
        Assert.Multiple(() =>
        {
            Assert.That(IupacPattern.MatchesAtCircular("ATTCAAGA", 6, "GAATTC"), Is.True);
            Assert.That(IupacPattern.MatchesAtCircular("ATTCAAGA", 5, "GAATTC"), Is.False);
            Assert.That(IupacPattern.MatchesAt("ATTCAAGA", 6, "GAATTC"), Is.False, "the linear form does not wrap");
        });
    }
}
