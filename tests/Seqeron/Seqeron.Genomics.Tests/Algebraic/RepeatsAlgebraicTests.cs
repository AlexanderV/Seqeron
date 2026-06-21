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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: REP-PALIN-001 — DNA palindromes (Repeats)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 17.
    //
    // Model: a DNA palindrome reads the same 5'→3' on both strands, i.e. it equals
    //        its own reverse complement. These are the recognition sites of many
    //        restriction enzymes (e.g. EcoRI GAATTC).
    //   — RepeatFinder.FindPalindromes; docs/algorithms/Repeat_Analysis.
    //
    // Laws under test (checklist row 17):
    //   • ID   — defining fixpoint: every palindrome equals revcomp(self).
    //   • DIST — a DNA palindrome has EVEN length: position i must pair with the
    //            base at the mirror position via complement, so an odd centre base
    //            would have to be its own complement, which no base is.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: every reported palindrome is a fixpoint of reverse complement —
    /// Sequence = revcomp(Sequence). Asserted over random DNA.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Palindromes_Identity_EqualOwnReverseComplement()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 8), seq =>
        {
            var palindromes = RepeatFinder.FindPalindromes(seq).ToList();
            bool ok = palindromes.All(p =>
                p.Sequence == DnaSequence.GetReverseComplementString(p.Sequence));
            return ok.Label($"a reported palindrome was not self-revcomp in \"{seq}\"");
        });
    }

    /// <summary>
    /// DIST: every reported DNA palindrome has even length (no self-complementary
    /// centre base exists).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Palindromes_Distributive_EvenLength()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 8), seq =>
        {
            var palindromes = RepeatFinder.FindPalindromes(seq).ToList();
            bool ok = palindromes.All(p => p.Length % 2 == 0 && p.Sequence.Length == p.Length);
            return ok.Label($"an odd-length palindrome was reported in \"{seq}\"");
        });
    }

    /// <summary>
    /// Witness: the EcoRI site GAATTC is a length-6 palindrome (revcomp = GAATTC),
    /// so both laws are exercised non-vacuously on a real detection.
    /// </summary>
    [Test]
    public void Palindromes_WorkedExample_EcoRiSite()
    {
        DnaSequence.GetReverseComplementString("GAATTC").Should().Be("GAATTC");
        var seq = new DnaSequence("AAAGAATTCAAA");
        var palindromes = RepeatFinder.FindPalindromes(seq).ToList();
        palindromes.Should().Contain(p => p.Sequence == "GAATTC" && p.Length == 6);
        palindromes.Should().OnlyContain(p =>
            p.Length % 2 == 0 && p.Sequence == DnaSequence.GetReverseComplementString(p.Sequence));
    }
}
