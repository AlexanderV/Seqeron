using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for ProbeDesigner.ValidateProbe and CheckSpecificity.
/// Test Unit: PROBE-VALID-001
/// 
/// Evidence Sources:
/// - Wikipedia: Hybridization probe (cross-hybridization, stringency)
/// - Wikipedia: DNA microarray (probe specificity)
/// - Wikipedia: Off-target genome editing (mismatch tolerance 1-5 bp)
/// - Wikipedia: BLAST (approximate matching algorithms)
/// </summary>
[TestFixture]
public class ProbeDesigner_ProbeValidation_Tests
{
    #region Test Data

    // Standard probe for validation tests
    private static readonly string StandardProbe = "ACGTACGTACGTACGTACGT";

    // Self-complementary (palindromic) probe
    private static readonly string PalindromicProbe = "GCGCGCGCGCGCGCGCGCGC";

    // Unique probe that appears once in reference
    private static readonly string UniqueProbe = "ATCGATCGATCGATCGATCG";

    // Reference containing the unique probe once
    private static readonly string[] SingleMatchReference = new[]
    {
        "NNNNNATCGATCGATCGATCGATCGNNNN"
    };

    // Reference with repeated sequence (multiple matches)
    private static readonly string[] MultipleMatchReference = new[]
    {
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    };

    #endregion

    #region ValidateProbe - Boundary Conditions (Must)

    [Test]
    public void ValidateProbe_EmptyProbe_ReturnsValidationResult()
    {
        // M1: Empty probe boundary condition
        // An empty probe has no sequence to hybridize — specificity must be 0.0
        // Source: Invariant #5 — offTargetHits == 0 → specificityScore == 0.0
        var validation = ProbeDesigner.ValidateProbe("", SingleMatchReference);

        Assert.Multiple(() =>
        {
            Assert.That(validation.SpecificityScore, Is.EqualTo(0.0),
                "Empty probe cannot hybridize — specificity must be 0.0");
            Assert.That(validation.OffTargetHits, Is.EqualTo(0),
                "Empty probe should report 0 off-target hits");
            Assert.That(validation.SelfComplementarity, Is.EqualTo(0.0),
                "Empty probe should have 0.0 self-complementarity");
            Assert.That(validation.IsValid, Is.False,
                "Empty probe should be invalid");
            Assert.That(validation.Issues, Has.Count.GreaterThan(0),
                "Empty probe should report issues");
        });
    }

    [Test]
    public void ValidateProbe_EmptyReferences_ReturnsValidationWithNoOffTargetHits()
    {
        // M2: Empty references — no sequences to search means 0 off-target hits
        // Per Invariant #5: offTargetHits == 0 → specificityScore == 0.0
        // (probe hasn't been shown to hybridize to any target)
        var validation = ProbeDesigner.ValidateProbe(StandardProbe, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(0),
                "No references means no off-target hits");
            Assert.That(validation.SpecificityScore, Is.EqualTo(0.0),
                "Zero hits means zero specificity per Invariant #5");
        });
    }

    [Test]
    public void ValidateProbe_NullReferences_ThrowsArgumentNullException()
    {
        // Evidence: Methods should validate null parameters
        // Null references array is invalid input - should throw ArgumentNullException
        Assert.Throws<ArgumentNullException>(() =>
            ProbeDesigner.ValidateProbe(StandardProbe, null!),
            "Null references should throw ArgumentNullException");
    }

    #endregion

    #region ValidateProbe - Specificity Invariants (Must)

    [Test]
    public void ValidateProbe_UniqueProbe_HasSpecificityScoreOne()
    {
        // M3: Unique probe (1 hit) should have specificity = 1.0
        var validation = ProbeDesigner.ValidateProbe(UniqueProbe, SingleMatchReference);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(1),
                "Unique probe should have exactly 1 hit");
            Assert.That(validation.SpecificityScore, Is.EqualTo(1.0),
                "Single hit should give specificity of 1.0");
        });
    }

    [Test]
    public void ValidateProbe_MultipleHits_ReducesSpecificityByHitCount()
    {
        // M4: Multiple hits reduce specificity to 1.0/hitCount
        // Probe "AAAAAAAAAA" (10×A) in "AAA...A" (34×A): exact match at every position 0..24 = 25 hits
        // Specificity = 1.0/25 = 0.04
        string probe = "AAAAAAAAAA";
        var validation = ProbeDesigner.ValidateProbe(probe, MultipleMatchReference);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(25),
                "10-mer in 34-mer poly-A: 34-10+1 = 25 exact match positions");
            Assert.That(validation.SpecificityScore, Is.EqualTo(1.0 / 25).Within(0.0001),
                "Specificity = 1.0/25 = 0.04 (Invariant #6)");
        });
    }

    // M5, M6, M7 range/non-negative invariants: deleted as duplicates.
    // Range is verified by AllInvariants test (Invariant Group) and implicitly
    // by every exact-value test (M1-M4, M8-M12, S1-S4).

    #endregion

    #region ValidateProbe - Self-Complementarity (Must)

    [Test]
    public void ValidateProbe_HighSelfComplementarity_ReportsInIssues()
    {
        // M8: High self-complementarity (>30%) should be reported in issues
        // PalindromicProbe "GCGCGCGCGCGCGCGCGCGC" is its own reverse complement → selfComp = 1.0
        // 1.0 > default threshold 0.3 → issue must be generated
        var validation = ProbeDesigner.ValidateProbe(PalindromicProbe, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(validation.SelfComplementarity, Is.EqualTo(1.0),
                "GC-repeat palindrome: every position matches its reverse complement");
            Assert.That(validation.Issues, Has.Some.Contain("Self-complementarity"),
                "Issues must report self-complementarity when above threshold");
        });
    }

    #endregion

    #region ValidateProbe - Case Sensitivity (Must)

    [Test]
    public void ValidateProbe_MixedCaseProbe_HandledCaseInsensitively()
    {
        // M12: Case-insensitive probe handling
        string upperProbe = UniqueProbe.ToUpperInvariant();
        string lowerProbe = UniqueProbe.ToLowerInvariant();
        string mixedProbe = "AtCgAtCgAtCgAtCgAtCg";

        var validationUpper = ProbeDesigner.ValidateProbe(upperProbe, SingleMatchReference);
        var validationLower = ProbeDesigner.ValidateProbe(lowerProbe, SingleMatchReference);
        var validationMixed = ProbeDesigner.ValidateProbe(mixedProbe, SingleMatchReference);

        Assert.Multiple(() =>
        {
            Assert.That(validationLower.OffTargetHits, Is.EqualTo(validationUpper.OffTargetHits),
                "Case should not affect off-target hit count");
            Assert.That(validationMixed.OffTargetHits, Is.EqualTo(validationUpper.OffTargetHits),
                "Mixed case should match upper case result");
            Assert.That(validationLower.SpecificityScore, Is.EqualTo(validationUpper.SpecificityScore).Within(0.001),
                "Case should not affect specificity score");
        });
    }

    #endregion

    #region CheckSpecificity - Suffix Tree (Must)

    [Test]
    public void CheckSpecificity_UniqueSequence_ReturnsOne()
    {
        // M10: Unique match returns 1.0
        string genome = "NNNNNATCGATCGATCGATCGATCGNNNN";
        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        double specificity = ProbeDesigner.CheckSpecificity(UniqueProbe, genomeIndex);

        Assert.That(specificity, Is.EqualTo(1.0),
            "Unique probe should have specificity 1.0");
    }

    [Test]
    public void CheckSpecificity_MultipleOccurrences_ReturnsOneOverCount()
    {
        // M11: Multiple matches returns 1.0/count
        string repeatedSequence = "ACGT";
        string genome = "ACGTACGTACGTACGT"; // Contains ACGT 4 times
        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        double specificity = ProbeDesigner.CheckSpecificity(repeatedSequence, genomeIndex);
        var positions = genomeIndex.FindAllOccurrences(repeatedSequence);

        Assert.Multiple(() =>
        {
            Assert.That(positions.Count, Is.GreaterThan(1),
                "Should find multiple occurrences");
            Assert.That(specificity, Is.EqualTo(1.0 / positions.Count).Within(0.001),
                "Specificity should be 1.0 / count");
        });
    }

    [Test]
    public void CheckSpecificity_NoMatch_ReturnsZero()
    {
        // M9 variant: No match returns 0.0
        string genome = "AAAAAAAAAAAAAAAAAAAAAA";
        string probe = "GCGCGCGCGC"; // Won't match in all-A genome
        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        double specificity = ProbeDesigner.CheckSpecificity(probe, genomeIndex);

        Assert.That(specificity, Is.EqualTo(0.0),
            "Non-matching probe should have specificity 0.0");
    }

    // M9 range test: deleted as duplicate of M10 (unique→1.0), M11 (multi→1/N), NoMatch (→0.0).

    #endregion

    #region ValidateProbe - Secondary Structure (Should)

    [Test]
    public void ValidateProbe_PotentialHairpin_DetectsSecondaryStructure()
    {
        // S1: Secondary structure potential detected for hairpin sequences
        // Stem-loop: GCGC (stem, 4nt) + TTT (loop, 3nt) + GCGC (stem, 4nt) + filler
        // HasSecondaryStructurePotential checks inverted repeats with stemLen≥4, gap=3
        // revComp("GCGC") = "GCGC" → 4/4 = 100% match ≥ 80% threshold → detected
        string hairpinProbe = "GCGCTTTGCGCAAAAAAAAA"; // 20 chars

        var validation = ProbeDesigner.ValidateProbe(hairpinProbe, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(validation.HasSecondaryStructure, Is.True,
                "Hairpin stem GCGC-TTT-GCGC must be detected as secondary structure");
            Assert.That(validation.Issues, Has.Some.Contain("secondary structure"),
                "Issues must report secondary structure potential");
        });
    }

    #endregion

    #region ValidateProbe - Issues List (Should)

    [Test]
    public void ValidateProbe_ProblematicProbe_PopulatesIssuesList()
    {
        // S2: Issues list populated for problematic probes
        // 10-mer poly-A in 25-mer poly-A: 25-10+1 = 16 exact match positions → offTargetHits = 16
        // Implementation adds "{N} potential off-target sites" when offTargetHits > 1
        string probe = "AAAAAAAAAA";
        var references = new[] { "AAAAAAAAAAAAAAAAAAAAAAAAA" }; // 25 A's

        var validation = ProbeDesigner.ValidateProbe(probe, references);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(16),
                "10-mer in 25-mer poly-A: 25-10+1 = 16 positions");
            Assert.That(validation.Issues, Has.Some.Contain("16 potential off-target sites"),
                "Issues must report exact off-target count");
        });
    }

    [Test]
    public void ValidateProbe_MultipleProblems_IsValidFalse()
    {
        // S3: IsValid false when multiple issues exist
        // "GCGCGCGCGC" (10-mer) → selfComp = 1.0 (palindrome)
        // In 32-char GC-repeat, only even positions match (odd positions are shifted by 1 → 10 mismatches)
        // Even positions 0,2,4,...,22 = 12 hits
        // isValid formula: issues.Count==0 || (offTargetHits<=1 && selfComp<=0.4)
        //   → false || (12<=1 && 1.0<=0.4) → false
        string probe = "GCGCGCGCGC";
        var references = new[] { "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" }; // 32 chars

        var validation = ProbeDesigner.ValidateProbe(probe, references);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(12),
                "10-mer GC-repeat in 32-char GC-repeat: matches at 12 even positions");
            Assert.That(validation.SelfComplementarity, Is.EqualTo(1.0),
                "GC-repeat is its own reverse complement");
            Assert.That(validation.Issues.Count, Is.EqualTo(2),
                "Should report off-target + self-complementarity issues");
            Assert.That(validation.IsValid, Is.False,
                "offTargetHits>1 AND selfComp>0.4 → IsValid must be false");
        });
    }

    #endregion

    #region ValidateProbe - Approximate Matching (Should)

    [Test]
    public void ValidateProbe_ApproximateMatching_FindsNearMatches()
    {
        // S4: Approximate matching with maxMismatches works correctly
        // Reference has ACGAACGAACGAACGT at position 5 — differs from probe at positions 3,7,11 (3 mismatches)
        // With maxMismatches=0: no match (3 mismatches > 0)
        // With maxMismatches=3: 1 match → specificity = 1.0
        string probe = "ACGTACGTACGTACGT"; // 16-mer
        var references = new[] { "TTTTTACGAACGAACGAACGTTTTT" }; // near-match at pos 5

        var strict = ProbeDesigner.ValidateProbe(probe, references, maxMismatches: 0);
        var approx = ProbeDesigner.ValidateProbe(probe, references, maxMismatches: 3);

        Assert.Multiple(() =>
        {
            Assert.That(strict.OffTargetHits, Is.EqualTo(0),
                "Exact matching (0 mismatches) should find no hits for 3-mismatch variant");
            Assert.That(approx.OffTargetHits, Is.EqualTo(1),
                "Approximate matching (3 mismatches) should find the near-match");
            Assert.That(approx.SpecificityScore, Is.EqualTo(1.0),
                "Single hit → specificity = 1.0");
        });
    }

    #endregion

    #region Invariant Group Assertions

    [Test]
    public void ValidateProbe_AllInvariants_HoldForTypicalProbe()
    {
        // Combined invariant test for comprehensive coverage
        var validation = ProbeDesigner.ValidateProbe(StandardProbe, SingleMatchReference);

        Assert.Multiple(() =>
        {
            // Specificity range (M5)
            Assert.That(validation.SpecificityScore, Is.InRange(0.0, 1.0),
                "Specificity out of range");

            // Self-complementarity range (M6)
            Assert.That(validation.SelfComplementarity, Is.InRange(0.0, 1.0),
                "Self-complementarity out of range");

            // Off-target non-negative (M7)
            Assert.That(validation.OffTargetHits, Is.GreaterThanOrEqualTo(0),
                "OffTargetHits is negative");

            // Issues list not null
            Assert.That(validation.Issues, Is.Not.Null,
                "Issues list should not be null");

            // Specificity formula consistency (all three invariants)
            if (validation.OffTargetHits == 0)
            {
                Assert.That(validation.SpecificityScore, Is.EqualTo(0.0),
                    "Zero hits should give specificity 0.0 (Invariant #5)");
            }
            else if (validation.OffTargetHits == 1)
            {
                Assert.That(validation.SpecificityScore, Is.EqualTo(1.0),
                    "Single hit should give specificity 1.0 (Invariant #4)");
            }
            else if (validation.OffTargetHits > 1)
            {
                Assert.That(validation.SpecificityScore, Is.EqualTo(1.0 / validation.OffTargetHits).Within(0.001),
                    "Specificity should equal 1.0 / hitCount (Invariant #6)");
            }
        });
    }

    [Test]
    public void ValidateProbe_ZeroHits_ReturnsZeroSpecificity()
    {
        // Explicit Invariant #5: offTargetHits == 0 → specificityScore == 0.0
        // A probe that matches nothing in the references has not demonstrated
        // hybridization capability → specificity is zero.
        // Consistent with CheckSpecificity which also returns 0.0 for hitCount == 0.
        string nonExistentProbe = "TTTTTTTTTTTTTTTTTTTT"; // 20× T — unlikely in reference
        string[] references = { "ACGACGACGACGACGACGACGACGACG" };

        var validation = ProbeDesigner.ValidateProbe(nonExistentProbe, references, maxMismatches: 0);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(0),
                "Probe should not match any site in references");
            Assert.That(validation.SpecificityScore, Is.EqualTo(0.0),
                "Invariant #5: zero hits must yield zero specificity");
        });
    }

    #endregion

    #region Could - Efficiency and Completeness

    [Test]
    public void ValidateProbe_LongReference_FindsProbeCorrectly()
    {
        // C1: Long reference sequences handled correctly
        // UniqueProbe embedded at position 10_000 in 20_010-char poly-T reference
        string longRef = new string('T', 10_000) + UniqueProbe + new string('T', 10_000);

        var validation = ProbeDesigner.ValidateProbe(UniqueProbe, new[] { longRef }, maxMismatches: 0);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(1),
                "Should find exactly 1 hit in long reference");
            Assert.That(validation.SpecificityScore, Is.EqualTo(1.0),
                "Single hit → specificity = 1.0");
        });
    }

    [Test]
    public void ValidateProbe_MultipleReferences_AccumulatesHits()
    {
        // C2: Multiple references all searched — hits accumulate across references
        // UniqueProbe appears once in each of 3 separate references → 3 total hits
        var ref1 = "TTTTT" + UniqueProbe + "TTTTT";
        var ref2 = "CCCCC" + UniqueProbe + "CCCCC";
        var ref3 = "GGGGG" + UniqueProbe + "GGGGG";

        var validation = ProbeDesigner.ValidateProbe(UniqueProbe, new[] { ref1, ref2, ref3 }, maxMismatches: 0);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(3),
                "Hits should accumulate across all 3 references");
            Assert.That(validation.SpecificityScore, Is.EqualTo(1.0 / 3).Within(0.0001),
                "Specificity = 1.0/3 (Invariant #6)");
        });
    }

    #endregion

    #region Integration with Probe Design

    [Test]
    public void DesignProbes_WithGenomeIndex_UsesCheckSpecificity()
    {
        // Verify that DesignProbes with suffix tree uses CheckSpecificity internally
        string uniqueRegion = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        string repeatedRegion = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        string genome = uniqueRegion + repeatedRegion + "AAAAAAAA" + repeatedRegion;

        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        var param = ProbeDesigner.Defaults.Microarray with { MinLength = 50, MaxLength = 52 };
        var probes = ProbeDesigner.DesignProbes(uniqueRegion, genomeIndex, param, maxProbes: 3, requireUnique: true).ToList();

        // All returned probes should be unique in the genome
        foreach (var probe in probes)
        {
            var positions = genomeIndex.FindAllOccurrences(probe.Sequence);
            Assert.That(positions.Count, Is.EqualTo(1),
                $"Probe should be unique but found {positions.Count} occurrences");
        }
    }

    #endregion
}
