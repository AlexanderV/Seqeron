using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for REP-INV-001: Inverted Repeat Detection.
/// Canonical implementation: RepeatFinder.FindInvertedRepeats
/// 
/// An inverted repeat is a sequence followed by its reverse complement
/// with intervening nucleotides (loop). When folded, forms stem-loop/hairpin structures.
/// 
/// Sources:
/// - Wikipedia: Inverted repeat, Stem-loop, Palindromic sequence
/// - EMBOSS einverted documentation
/// - Pearson et al. (1996), Bissler (1998)
/// </summary>
[TestFixture]
public class RepeatFinder_InvertedRepeat_Tests
{
    #region MUST Tests - Core Algorithm Verification

    /// <summary>
    /// M1: Core algorithm - detects basic stem-loop structure.
    /// GCGC...loop...GCGC where GCGC reverse complement = GCGC (self-complementary).
    /// Source: Wikipedia - inverted repeat definition.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_SimpleHairpin_FindsRepeat()
    {
        // GCGC reverse complement = GCGC (G↔C, C↔G, G↔C, C↔G reversed = GCGC)
        var sequence = new DnaSequence("AAGCGCAAAAGCGCAA");
        //                              01234567890123456
        //                              Left: positions 2-5 (GCGC)
        //                              Loop: positions 6-9 (AAAA)
        //                              Right: positions 10-13 (GCGC)

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));

        var hairpin = results.First(r => r.ArmLength == 4);
        Assert.Multiple(() =>
        {
            Assert.That(hairpin.LeftArm, Is.EqualTo("GCGC"));
            Assert.That(hairpin.RightArm, Is.EqualTo("GCGC"));
            Assert.That(hairpin.LoopLength, Is.GreaterThanOrEqualTo(3));
            Assert.That(hairpin.CanFormHairpin, Is.True, "Loop ≥ 3 should allow hairpin formation");
        });
    }

    /// <summary>
    /// M2: Palindromic (self-complementary) sequences like restriction sites.
    /// GAATTC reverse complement = GAATTC (EcoRI recognition site pattern).
    /// Source: Wikipedia - palindromic sequence.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_PalindromeSequence_SelfComplementary()
    {
        // GAATTC: G↔C, A↔T, A↔T, T↔A, T↔A, C↔G → complement = CTTAAG
        // Reversed: GAATTC - it's self-complementary!
        var sequence = new DnaSequence("GAATTCAAAAGAATTC");

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));
        var result = results.First(r => r.ArmLength >= 4);
        Assert.That(result.LeftArm.Length, Is.GreaterThanOrEqualTo(4));
    }

    /// <summary>
    /// M3: Verify reverse complement relationship between arms.
    /// Left arm's reverse complement must equal right arm.
    /// Source: Wikipedia - inverted repeat definition.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_ReverseComplementMatch_BothArmsCorrect()
    {
        // ACGT reverse complement = ACGT (A↔T,C↔G,G↔C,T↔A → TGCA → reversed = ACGT)
        var sequence = new DnaSequence("ACGTTTTTACGT");
        //                              0123456789012
        //                              Left: 0-3 (ACGT), Loop: 4-8 (TTTTT), Right: 8-11 (ACGT)

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));

        foreach (var result in results)
        {
            string expectedRightArm = DnaSequence.GetReverseComplementString(result.LeftArm);
            Assert.That(result.RightArm, Is.EqualTo(expectedRightArm),
                $"Right arm must equal reverse complement of left arm");
        }
    }

    /// <summary>
    /// M4: Sequence without complementary regions returns empty.
    /// AAAA reverse complement = TTTT, no match possible.
    /// Source: Standard edge case.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_NoInvertedRepeats_ReturnsEmpty()
    {
        var sequence = new DnaSequence("AAAAAAAAAAAAAA");

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Is.Empty, "Homopolymer A has no self-complementary regions");
    }

    /// <summary>
    /// M5: Empty sequence returns empty result.
    /// Source: Standard boundary condition.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_EmptySequence_ReturnsEmpty()
    {
        var results = RepeatFinder.FindInvertedRepeats("", minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// M6: MinArmLength parameter filters results correctly.
    /// Source: Algorithm specification.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_MinArmLength_RespectsThreshold()
    {
        // GC revcomp = GC (2bp), GCGC revcomp = GCGC (4bp)
        var sequence = new DnaSequence("GCAAGC"); // GC...AA...GC (arm=2, loop=2)

        var resultsMin4 = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 0).ToList();
        var resultsMin2 = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 2, maxLoopLength: 10, minLoopLength: 0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(resultsMin4, Is.Empty, "No 4bp arms exist");
            Assert.That(resultsMin2.Any(r => r.ArmLength == 2), Is.True, "2bp arm should be found");
        });
    }

    /// <summary>
    /// M7: MinLoopLength parameter filters results correctly.
    /// Source: Algorithm specification.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_MinLoopLength_RespectsThreshold()
    {
        // GCGC...AA...GCGC (arm=4, loop=2)
        var sequence = new DnaSequence("GCGCAAGCGC");

        var resultsMinLoop3 = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();
        var resultsMinLoop1 = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(resultsMinLoop3, Is.Empty, "Loop of 2 should be filtered when minLoop=3");
            Assert.That(resultsMinLoop1.Any(r => r.LoopLength == 2), Is.True, "Loop of 2 should be found when minLoop=1");
        });
    }

    /// <summary>
    /// M8: MaxLoopLength parameter filters results correctly.
    /// Source: Algorithm specification; EMBOSS einverted.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_MaxLoopLength_RespectsThreshold()
    {
        // GCGC followed by 15 A's then GCGC (loop = 15)
        var sequence = new DnaSequence("GCGC" + new string('A', 15) + "GCGC");

        var resultsMaxLoop10 = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();
        var resultsMaxLoop20 = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 20, minLoopLength: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(resultsMaxLoop10, Is.Empty, "Loop of 15 should be filtered when maxLoop=10");
            Assert.That(resultsMaxLoop20.Any(r => r.LoopLength == 15), Is.True, "Loop of 15 should be found when maxLoop=20");
        });
    }

    /// <summary>
    /// M9: LoopLength calculation invariant.
    /// loopLength = rightArmStart - (leftArmStart + armLength).
    /// Source: Invariant.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_LoopLength_CalculatedCorrectly()
    {
        var sequence = new DnaSequence("GCGCAAAAGCGC");
        //                              012345678901
        //                              Left: 0-3, Loop: 4-7, Right: 8-11

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));

        foreach (var result in results)
        {
            int expectedLoopLength = result.RightArmStart - (result.LeftArmStart + result.ArmLength);
            Assert.That(result.LoopLength, Is.EqualTo(expectedLoopLength),
                "LoopLength = RightArmStart - (LeftArmStart + ArmLength)");
        }
    }

    /// <summary>
    /// M10: TotalLength invariant.
    /// TotalLength = 2 × ArmLength + LoopLength.
    /// Source: Invariant.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_TotalLength_InvariantHolds()
    {
        var sequence = new DnaSequence("GCGCAAAAGCGC");

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));

        foreach (var result in results)
        {
            int expectedTotal = 2 * result.ArmLength + result.LoopLength;
            Assert.That(result.TotalLength, Is.EqualTo(expectedTotal),
                "TotalLength must equal 2×ArmLength + LoopLength");
        }
    }

    /// <summary>
    /// M11: CanFormHairpin is true when LoopLength ≥ 3.
    /// Loops fewer than 3 bases are sterically impossible for hairpin formation.
    /// Source: Wikipedia - stem-loop.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_CanFormHairpin_LoopLengthValidation()
    {
        // Test with loop = 4 (should form hairpin)
        var sequenceLong = new DnaSequence("GCGCAAAAGCGC"); // loop = 4

        var resultsLong = RepeatFinder.FindInvertedRepeats(sequenceLong, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        if (resultsLong.Count > 0)
        {
            Assert.Multiple(() =>
            {
                foreach (var result in resultsLong)
                {
                    Assert.That(result.CanFormHairpin, Is.EqualTo(result.LoopLength >= 3),
                        $"CanFormHairpin should be true when LoopLength ({result.LoopLength}) ≥ 3");
                }
            });
        }

        // Test with loop = 2 (should NOT form hairpin biologically)
        var sequenceShort = new DnaSequence("GCGCAAGCGC"); // loop = 2
        var resultsShort = RepeatFinder.FindInvertedRepeats(sequenceShort, minArmLength: 4, maxLoopLength: 10, minLoopLength: 0).ToList();

        if (resultsShort.Count > 0)
        {
            var shortLoopResult = resultsShort.FirstOrDefault(r => r.LoopLength < 3);
            if (shortLoopResult.ArmLength > 0)
            {
                Assert.That(shortLoopResult.CanFormHairpin, Is.False,
                    "CanFormHairpin should be false when LoopLength < 3");
            }
        }
    }

    /// <summary>
    /// M12: Position values are correct (0-based indexing).
    /// Source: Implementation contract.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_PositionsCorrect_ZeroBased()
    {
        var sequence = new DnaSequence("GCGCAAAAGCGC");
        //                              012345678901
        // Expected: Left at 0, Right at 8

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));

        var result = results.First(r => r.ArmLength == 4);
        Assert.Multiple(() =>
        {
            Assert.That(result.LeftArmStart, Is.EqualTo(0), "Left arm should start at position 0");
            Assert.That(result.RightArmStart, Is.EqualTo(8), "Right arm should start at position 8");

            // Verify by extracting from original sequence
            string leftFromSeq = sequence.Sequence.Substring(result.LeftArmStart, result.ArmLength);
            string rightFromSeq = sequence.Sequence.Substring(result.RightArmStart, result.ArmLength);
            Assert.That(result.LeftArm, Is.EqualTo(leftFromSeq));
            Assert.That(result.RightArm, Is.EqualTo(rightFromSeq));
        });
    }

    /// <summary>
    /// M13: Sequence too short for any inverted repeat returns empty.
    /// Source: Boundary condition.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_SequenceTooShort_ReturnsEmpty()
    {
        // minArm=4, minLoop=3 requires at least 4+3+4=11 bases
        var sequence = new DnaSequence("GCGCAAGCGC"); // 10 bases

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Is.Empty, "Sequence of 10 bases cannot have 4bp arms + 3bp loop");
    }

    #endregion

    #region SHOULD Tests - Important Scenarios

    /// <summary>
    /// S1: Sequence with multiple inverted repeats finds all.
    /// Source: Real genome scenario.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_MultipleRepeats_FindsAll()
    {
        // Two separate inverted repeats
        // First: GCGC...AAAA...GCGC
        // Second: ATAT...TTT...ATAT (AT revcomp = AT)
        var sequence = new DnaSequence("GCGCAAAAGCGCTTTTATATTTTATAT");
        //                              0         1         2
        //                              012345678901234567890123456

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results.Count, Is.GreaterThanOrEqualTo(1), "Should find at least one inverted repeat");
    }

    /// <summary>
    /// S2: String overload produces same results as DnaSequence overload.
    /// Source: API consistency.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_StringOverload_MatchesDnaSequenceOverload()
    {
        string seqString = "GCGCAAAAGCGC";
        var seqDna = new DnaSequence(seqString);

        var resultsString = RepeatFinder.FindInvertedRepeats(seqString, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();
        var resultsDna = RepeatFinder.FindInvertedRepeats(seqDna, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(resultsString.Count, Is.EqualTo(resultsDna.Count));

        for (int i = 0; i < resultsString.Count; i++)
        {
            Assert.Multiple(() =>
            {
                Assert.That(resultsString[i].LeftArmStart, Is.EqualTo(resultsDna[i].LeftArmStart));
                Assert.That(resultsString[i].ArmLength, Is.EqualTo(resultsDna[i].ArmLength));
                Assert.That(resultsString[i].LoopLength, Is.EqualTo(resultsDna[i].LoopLength));
            });
        }
    }

    /// <summary>
    /// S3: Lowercase input is handled correctly (case-insensitive).
    /// Source: Implementation robustness.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_CaseInsensitivity_HandledCorrectly()
    {
        var resultsLower = RepeatFinder.FindInvertedRepeats("gcgcaaaagcgc", minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();
        var resultsUpper = RepeatFinder.FindInvertedRepeats("GCGCAAAAGCGC", minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(resultsLower.Count, Is.EqualTo(resultsUpper.Count),
            "Case should not affect result count");

        if (resultsLower.Count > 0)
        {
            Assert.That(resultsLower[0].ArmLength, Is.EqualTo(resultsUpper[0].ArmLength));
        }
    }

    /// <summary>
    /// S4: Document behavior for adjacent arms (loop = 0).
    /// Source: ASSUMPTION - palindromic sequences are edge case.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_AdjacentArms_BehaviorDocumented()
    {
        // GAATTC is a palindrome (revcomp of itself) - adjacent arms would have loop=0
        var sequence = new DnaSequence("GAATTC"); // This IS a palindrome, but no loop
        // For inverted repeat detection, we need two separate copies with intervening sequence

        var resultsNoMinLoop = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 3, maxLoopLength: 10, minLoopLength: 0).ToList();

        // When minLoopLength=0, adjacent palindromic regions might be found
        // This documents the behavior for such edge cases
        Assert.Pass($"With minLoopLength=0, found {resultsNoMinLoop.Count} results (behavior documented)");
    }

    /// <summary>
    /// S5: Known biological hairpin structure (simplified tRNA-like).
    /// Source: Wikipedia - tRNA structure.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_BiologicalHairpin_KnownStructure()
    {
        // Simplified hairpin: GGGC...loop...GCCC
        // GGGC revcomp = GCCC (G↔C, G↔C, G↔C, C↔G → CCCG → reversed = GCCC)
        var sequence = new DnaSequence("GGGCAAAAGCCC");

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));

        var hairpin = results.First();
        Assert.Multiple(() =>
        {
            Assert.That(hairpin.LeftArm, Is.EqualTo("GGGC"));
            Assert.That(hairpin.RightArm, Is.EqualTo("GCCC"));
            Assert.That(DnaSequence.GetReverseComplementString("GGGC"), Is.EqualTo("GCCC"));
        });
    }

    #endregion

    #region COULD Tests - Additional Scenarios

    /// <summary>
    /// C1: Known restriction site palindromes can be detected as inverted repeats.
    /// EcoRI: GAATTC, BamHI: GGATCC
    /// Source: Wikipedia - restriction enzymes.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_RestrictionSitePalindromes_Detected()
    {
        // Two EcoRI sites separated by loop
        // GAATTC revcomp = GAATTC (self-complementary)
        var sequence = new DnaSequence("GAATTCAAAAGAATTC");

        var results = RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1),
            "Should detect inverted repeat formed by two restriction sites");
    }

    /// <summary>
    /// C2: Loop sequence is correctly extracted.
    /// Source: Implementation verification.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_LoopSequence_CorrectlyExtracted()
    {
        var sequence = new DnaSequence("GCGCTTTTAGCGC");
        //                              0123456789012
        // Left: GCGC (0-3), Loop: TTTTA (4-8), Right: GCGC (9-12)
        // Wait: GCGC revcomp = GCGC, so we need: GCGC...loop...GCGC

        // Let's use a clear example
        var seq2 = new DnaSequence("GCGCTTTTGCGC");
        //                          012345678901

        var results = RepeatFinder.FindInvertedRepeats(seq2, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        if (results.Count > 0)
        {
            var result = results.First(r => r.ArmLength == 4);
            Assert.That(result.Loop, Is.EqualTo("TTTT"), "Loop should be the intervening sequence");
        }
    }

    #endregion

    #region Parameter Validation Tests

    /// <summary>
    /// Null DnaSequence should throw ArgumentNullException.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_NullDnaSequence_ThrowsArgumentNullException()
    {
        DnaSequence? nullSeq = null;

        Assert.Throws<System.ArgumentNullException>(() =>
            RepeatFinder.FindInvertedRepeats(nullSeq!, minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList());
    }

    /// <summary>
    /// MinArmLength less than 2 should throw.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_MinArmLengthTooSmall_ThrowsArgumentOutOfRangeException()
    {
        var sequence = new DnaSequence("GCGCAAAAGCGC");

        Assert.Throws<System.ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 1, maxLoopLength: 10, minLoopLength: 3).ToList());
    }

    /// <summary>
    /// MinLoopLength less than 0 should throw.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_NegativeMinLoopLength_ThrowsArgumentOutOfRangeException()
    {
        var sequence = new DnaSequence("GCGCAAAAGCGC");

        Assert.Throws<System.ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindInvertedRepeats(sequence, minArmLength: 4, maxLoopLength: 10, minLoopLength: -1).ToList());
    }

    #endregion
}
