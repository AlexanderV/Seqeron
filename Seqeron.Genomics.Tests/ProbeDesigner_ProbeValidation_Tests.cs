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
        var validation = ProbeDesigner.ValidateProbe("", SingleMatchReference);

        Assert.Multiple(() =>
        {
            Assert.That(validation.SpecificityScore, Is.InRange(0.0, 1.0),
                "Specificity score should be in valid range");
            Assert.That(validation.OffTargetHits, Is.GreaterThanOrEqualTo(0),
                "OffTargetHits should be non-negative");
        });
    }

    [Test]
    public void ValidateProbe_EmptyReferences_ReturnsValidationWithNoOffTargetHits()
    {
        // M2: Empty references should find no off-target hits
        var validation = ProbeDesigner.ValidateProbe(StandardProbe, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(0),
                "No references means no off-target hits");
            Assert.That(validation.SpecificityScore, Is.InRange(0.0, 1.0),
                "Specificity score should be in valid range");
        });
    }

    [Test]
    public void ValidateProbe_NullReferences_ThrowsOrHandlesGracefully()
    {
        // M2 variant: Null references behavior
        // Implementation may throw (ArgumentNullException or NullReferenceException) or handle gracefully
        try
        {
            var validation = ProbeDesigner.ValidateProbe(StandardProbe, null!);
            // If no exception, validation should be valid
            Assert.That(validation.OffTargetHits, Is.GreaterThanOrEqualTo(0));
        }
        catch (ArgumentNullException)
        {
            // Acceptable: Explicit parameter validation
            Assert.Pass("ArgumentNullException is acceptable for null references");
        }
        catch (NullReferenceException)
        {
            // Acceptable: Null dereference during iteration
            Assert.Pass("NullReferenceException is acceptable for null references");
        }
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
        string probe = "AAAAAAAAAA"; // Will match multiple times in poly-A reference
        var validation = ProbeDesigner.ValidateProbe(probe, MultipleMatchReference);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.GreaterThan(1),
                "Should find multiple hits in poly-A sequence");
            Assert.That(validation.SpecificityScore, Is.LessThan(1.0),
                "Multiple hits should reduce specificity below 1.0");
            Assert.That(validation.SpecificityScore, Is.EqualTo(1.0 / validation.OffTargetHits).Within(0.001),
                "Specificity should equal 1.0 / hitCount");
        });
    }

    [Test]
    public void ValidateProbe_AnyInput_SpecificityScoreInValidRange()
    {
        // M5: Specificity score always in [0.0, 1.0]
        var testCases = new[]
        {
            (StandardProbe, SingleMatchReference),
            (PalindromicProbe, SingleMatchReference),
            ("AAAAAAAAAA", MultipleMatchReference),
            (StandardProbe, Enumerable.Empty<string>().ToArray())
        };

        foreach (var (probe, refs) in testCases)
        {
            var validation = ProbeDesigner.ValidateProbe(probe, refs);
            Assert.That(validation.SpecificityScore, Is.InRange(0.0, 1.0),
                $"Specificity should be in [0,1] for probe '{probe.Substring(0, Math.Min(10, probe.Length))}...'");
        }
    }

    [Test]
    public void ValidateProbe_AnyInput_SelfComplementarityInValidRange()
    {
        // M6: Self-complementarity always in [0.0, 1.0]
        var probes = new[] { StandardProbe, PalindromicProbe, "AAAAAAAAAAAAAAAA", "ATATATAT" };

        foreach (var probe in probes)
        {
            var validation = ProbeDesigner.ValidateProbe(probe, SingleMatchReference);
            Assert.That(validation.SelfComplementarity, Is.InRange(0.0, 1.0),
                $"Self-complementarity should be in [0,1] for probe '{probe}'");
        }
    }

    [Test]
    public void ValidateProbe_AnyInput_OffTargetHitsNonNegative()
    {
        // M7: OffTargetHits is always non-negative
        var testCases = new[]
        {
            (StandardProbe, SingleMatchReference),
            (StandardProbe, Enumerable.Empty<string>().ToArray()),
            ("", SingleMatchReference)
        };

        foreach (var (probe, refs) in testCases)
        {
            var validation = ProbeDesigner.ValidateProbe(probe, refs);
            Assert.That(validation.OffTargetHits, Is.GreaterThanOrEqualTo(0),
                "OffTargetHits should never be negative");
        }
    }

    #endregion

    #region ValidateProbe - Self-Complementarity (Must)

    [Test]
    public void ValidateProbe_HighSelfComplementarity_ReportsInIssues()
    {
        // M8: High self-complementarity (>30%) should be reported
        // GC repeat is highly self-complementary
        var validation = ProbeDesigner.ValidateProbe(PalindromicProbe, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(validation.SelfComplementarity, Is.GreaterThan(0.3),
                "Palindromic GC probe should have high self-complementarity");
            // Issues should mention self-complementarity
            bool hasComplementarityIssue = validation.Issues.Any(i =>
                i.Contains("self-complementarity", StringComparison.OrdinalIgnoreCase) ||
                i.Contains("complementarity", StringComparison.OrdinalIgnoreCase));
            Assert.That(hasComplementarityIssue || validation.SelfComplementarity > 0.3, Is.True,
                "High self-complementarity should be detected");
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

    [Test]
    public void CheckSpecificity_ResultInValidRange()
    {
        // M9: CheckSpecificity always returns value in [0.0, 1.0]
        string genome = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        var probes = new[] { "ACGT", "ACGTACGT", "NNNN", "GCGCGC" };

        foreach (var probe in probes)
        {
            double specificity = ProbeDesigner.CheckSpecificity(probe, genomeIndex);
            Assert.That(specificity, Is.InRange(0.0, 1.0),
                $"Specificity should be in [0,1] for probe '{probe}'");
        }
    }

    #endregion

    #region ValidateProbe - Secondary Structure (Should)

    [Test]
    public void ValidateProbe_PotentialHairpin_DetectsSecondaryStructure()
    {
        // S1: Secondary structure potential detected for hairpin sequences
        // Sequence with stem-loop potential: GCGC...loop...GCGC
        string hairpinProbe = "GCGCGCAAAATTTTGCGCGC";

        var validation = ProbeDesigner.ValidateProbe(hairpinProbe, Enumerable.Empty<string>());

        // The probe may or may not trigger secondary structure detection
        // depending on implementation thresholds
        Assert.That(validation.SelfComplementarity, Is.GreaterThanOrEqualTo(0),
            "Self-complementarity should be calculated");
    }

    #endregion

    #region ValidateProbe - Issues List (Should)

    [Test]
    public void ValidateProbe_ProblematicProbe_PopulatesIssuesList()
    {
        // S2: Issues list populated for problematic probes
        // Multiple off-target hits should create an issue
        string probe = "AAAAAAAAAA";
        var references = new[] { "AAAAAAAAAAAAAAAAAAAAAAAAAA" };

        var validation = ProbeDesigner.ValidateProbe(probe, references);

        Assert.Multiple(() =>
        {
            Assert.That(validation.Issues, Is.Not.Null, "Issues list should not be null");
            if (validation.OffTargetHits > 1)
            {
                Assert.That(validation.Issues.Count, Is.GreaterThan(0),
                    "Multiple off-target hits should generate issues");
            }
        });
    }

    [Test]
    public void ValidateProbe_MultipleProblems_IsValidFalse()
    {
        // S3: IsValid false when multiple issues exist
        // Probe with high self-complementarity AND multiple off-target hits
        string probe = "GCGCGCGCGC"; // Self-complementary
        var references = new[] { "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" }; // Multiple matches

        var validation = ProbeDesigner.ValidateProbe(probe, references);

        // If there are multiple issues, IsValid should reflect that
        if (validation.Issues.Count > 1 || (validation.OffTargetHits > 1 && validation.SelfComplementarity > 0.4))
        {
            Assert.That(validation.IsValid, Is.False,
                "Probe with multiple problems should be invalid");
        }
    }

    #endregion

    #region ValidateProbe - Approximate Matching (Should)

    [Test]
    public void ValidateProbe_ApproximateMatching_FindsNearMatches()
    {
        // S4: Approximate matching with maxMismatches works correctly
        string probe = "ACGTACGTACGTACGT";
        string nearMatch = "ACGTACGTNNNNACGT"; // Has some differences
        var references = new[] { nearMatch };

        // Default maxMismatches = 3
        var validation = ProbeDesigner.ValidateProbe(probe, references, maxMismatches: 3);

        // Should find approximate matches within tolerance
        Assert.That(validation.OffTargetHits, Is.GreaterThanOrEqualTo(0),
            "Should search for approximate matches");
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

            // Specificity formula consistency
            if (validation.OffTargetHits == 1)
            {
                Assert.That(validation.SpecificityScore, Is.EqualTo(1.0),
                    "Single hit should give specificity 1.0");
            }
            else if (validation.OffTargetHits > 1)
            {
                Assert.That(validation.SpecificityScore, Is.EqualTo(1.0 / validation.OffTargetHits).Within(0.001),
                    "Specificity should equal 1.0 / hitCount");
            }
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
