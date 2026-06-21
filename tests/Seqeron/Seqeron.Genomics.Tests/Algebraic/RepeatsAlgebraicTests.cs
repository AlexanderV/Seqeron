using System.Linq;
using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Repeats area (inverted repeats, palindromes).
///
/// Algebraic testing pins the structural equations every detected repeat must
/// satisfy — here the defining relations between the two arms of a hairpin and
/// the self-reverse-complementarity of a DNA palindrome — together with the
/// identity behaviour on inputs that contain no such structure.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 15, 17.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Repeats")]
public class RepeatsAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: REP-INV-001 — Inverted repeats (Repeats)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 15.
    //
    // Model: an inverted repeat is a left arm followed (after a loop) by a right
    //        arm equal to the reverse complement of the left arm — the building
    //        block of a hairpin/stem-loop.
    //   — RepeatFinder.FindInvertedRepeats; docs/algorithms/Repeat_Analysis.
    //
    // Laws under test (checklist row 15):
    //   • ID   — no reverse-complement arm present → no inverted repeats. A pure
    //            homopolymer (all 'A') cannot contain its own complement ('T'),
    //            so the result set is empty.
    //   • DIST — defining relation: for every reported repeat,
    //            RightArm = reverseComplement(LeftArm).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: a homopolymer has no inverted repeats — its arm's reverse complement
    /// (all 'T') never occurs, so the detector returns the empty set.
    /// </summary>
    [Test]
    public void InvertedRepeats_Identity_HomopolymerHasNone()
    {
        var seq = new DnaSequence(new string('A', 40));
        RepeatFinder.FindInvertedRepeats(seq).Should().BeEmpty();
    }

    /// <summary>
    /// DIST: every reported inverted repeat satisfies RightArm = revcomp(LeftArm),
    /// the structural definition of the repeat. Asserted over random DNA so the
    /// invariant is checked against whatever the detector actually emits.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property InvertedRepeats_Distributive_RightArmIsRevCompOfLeft()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 14), seq =>
        {
            var repeats = RepeatFinder.FindInvertedRepeats(seq).ToList();
            bool ok = repeats.All(r =>
                r.RightArm == DnaSequence.GetReverseComplementString(r.LeftArm)
                && r.ArmLength == r.LeftArm.Length
                && r.RightArm.Length == r.ArmLength);
            return ok.Label($"a right arm was not revcomp(left arm) in \"{seq}\"");
        });
    }

    /// <summary>
    /// DIST witness: a constructed hairpin GGGG-AAA-CCCC has right arm CCCC =
    /// revcomp(GGGG), so the law is exercised non-vacuously on a real detection.
    /// </summary>
    [Test]
    public void InvertedRepeats_Distributive_WorkedHairpin()
    {
        var seq = new DnaSequence("GGGGAAACCCC");
        var repeats = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4, minLoopLength: 3).ToList();
        repeats.Should().NotBeEmpty();
        repeats.Should().OnlyContain(r =>
            r.RightArm == DnaSequence.GetReverseComplementString(r.LeftArm));
    }
}
