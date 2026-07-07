namespace Seqeron.Genomics.Tests.Unit.MolTools;

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
    private const string StandardProbe = "ACGTACGTACGTACGTACGT";

    // Self-complementary (palindromic) probe
    private const string PalindromicProbe = "GCGCGCGCGCGCGCGCGCGC";

    // Unique probe that appears once in reference
    private const string UniqueProbe = "ATCGATCGATCGATCGATCG";

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

    #region ScanOffTargetsGapped - Gapped (Smith-Waterman) off-target scan (Must)

    // PROBE-VALID-001 — gapped off-target scan + on/off-target separation
    // Evidence: docs/Evidence/PROBE-VALID-001-Evidence.md
    // Sources:
    //   Smith TF, Waterman MS (1981) J Mol Biol 147(1):195-197 (local-alignment recurrence)
    //   Altschul SF et al. (1990) J Mol Biol 215(3):403-410 (gapped local alignment finds indels the ungapped scan misses)
    //   Kane MD et al. (2000) Nucleic Acids Res 28(22):4552-4557 (>75% identity over the probe → off-target)

    // Probe used across the gapped-scan tests.
    private const string GappedProbe = "ACGTACGTACGT"; // 12 nt

    // Reference containing the EXACT on-target at start 5 and an indel off-target
    // ("ACGTACTGTACGT" = probe with a 'T' inserted after position 6) at start 27.
    private static readonly string[] IndelOffTargetReference = new[]
    {
        "NNNNN" + GappedProbe + "NNNNNNNNNN" + "ACGTACTGTACGT" + "NNNNN"
    };

    [Test]
    public void ScanOffTargetsGapped_IndelOffTarget_FoundByGappedScanMissedByHammingScan()
    {
        // MG1: An off-target reachable ONLY through a single insertion is found by the gapped
        // scan but missed by the ungapped Hamming scan (Altschul 1990: gapped finds indels).
        // The indel region "ACGTACTGTACGT" has >=6 mismatches in every fixed 12-window, so the
        // default Hamming scan (maxMismatches=3) cannot reach it; only the exact on-target at 5.
        var hamming = ProbeDesigner.ValidateProbe(GappedProbe, IndelOffTargetReference, maxMismatches: 3);
        var gapped = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, IndelOffTargetReference);

        Assert.Multiple(() =>
        {
            // Ungapped Hamming scan: only the exact on-target match (pooled into OffTargetHits).
            Assert.That(hamming.OffTargetHits, Is.EqualTo(1),
                "Ungapped Hamming scan (maxMismatches=3) finds only the exact on-target; the indel site has >=6 mismatches per window");

            // Gapped scan: 1 on-target (exact, no gap) + 1 genuine indel off-target.
            Assert.That(gapped.OnTargetHits, Has.Count.EqualTo(1),
                "Gapped scan must identify exactly one intended on-target (the perfect exact match)");
            Assert.That(gapped.OffTargetCount, Is.EqualTo(1),
                "Gapped scan must find the indel off-target that the Hamming scan misses");
            Assert.That(gapped.OffTargetHits[0].HasGaps, Is.True,
                "The off-target is reachable only via an insertion → its alignment contains a gap");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_OnTargetExactMatch_NotCountedAsOffTarget()
    {
        // MG2: The intended on-target (perfect, ungapped, full-coverage exact match) is reported
        // separately and excluded from the off-target count — fixing the OffTargetHits pooling.
        var gapped = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, IndelOffTargetReference);

        var onTarget = gapped.OnTargetHits[0];

        Assert.Multiple(() =>
        {
            Assert.That(onTarget.Start, Is.EqualTo(5),
                "On-target exact match begins at reference position 5");
            Assert.That(onTarget.End, Is.EqualTo(16),
                "On-target exact 12-mer match ends (inclusive) at reference position 16");
            Assert.That(onTarget.Identity, Is.EqualTo(1.0).Within(1e-10),
                "On-target is a 12/12 exact match → identity 1.0");
            Assert.That(onTarget.Coverage, Is.EqualTo(1.0).Within(1e-10),
                "On-target covers the full probe → coverage 1.0");
            Assert.That(onTarget.HasGaps, Is.False,
                "On-target exact match has no gaps");
            Assert.That(gapped.OffTargetHits, Has.None.Matches<ProbeDesigner.GappedProbeHit>(
                h => h.Start == 5),
                "The on-target site must NOT appear among off-target hits");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_IndelOffTarget_HasExactHandDerivedIdentity()
    {
        // MG3: Exact identity/coverage on the hand-derived indel alignment.
        // probe "ACGTAC-GTACGT" vs ref "ACGTACTGTACGT": 12/12 identical aligned columns,
        // one insertion gap → identity 1.0, coverage 1.0, HasGaps true. Site starts at 27.
        var gapped = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, IndelOffTargetReference);

        var off = gapped.OffTargetHits[0];

        Assert.Multiple(() =>
        {
            Assert.That(off.Start, Is.EqualTo(27),
                "Indel off-target begins at reference position 27");
            Assert.That(off.Identity, Is.EqualTo(1.0).Within(1e-10),
                "Indel off-target: 12 identical aligned columns / 12 probe length = 1.0");
            Assert.That(off.Coverage, Is.EqualTo(1.0).Within(1e-10),
                "Indel off-target: 12 ungapped columns / 12 probe length = 1.0");
            Assert.That(off.HasGaps, Is.True,
                "The off-target alignment contains the insertion gap");
            Assert.That(off.AlignedProbe, Is.EqualTo("ACGTAC-GTACGT"),
                "Probe side of the local alignment carries the gap at the insertion point");
            Assert.That(off.AlignedReference, Is.EqualTo("ACGTACTGTACGT"),
                "Reference side of the local alignment is the indel region");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_IndelPlusMismatch_IdentityIsHandDerivedFraction()
    {
        // MG4: An off-target with one insertion AND a trailing mismatch.
        // Region "ACGTACTGTACTT": SW (zero-floor) trims the mismatched "TT" tail →
        // probe "ACGTAC-GTAC" vs ref "ACGTACTGTAC" = 10 identical columns / 12 = 0.8333.
        var references = new[] { "NNNNNNNNNN" + "ACGTACTGTACTT" + "NNNNN" };

        var gapped = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, references, minIdentity: 0.75);

        Assert.Multiple(() =>
        {
            Assert.That(gapped.OnTargetHits, Is.Empty,
                "No perfect exact match in this reference → no on-target");
            Assert.That(gapped.OffTargetCount, Is.EqualTo(1),
                "The indel+mismatch site exceeds the 0.75 threshold → one off-target");
            Assert.That(gapped.OffTargetHits[0].Identity, Is.EqualTo(10.0 / 12.0).Within(1e-10),
                "10 identical aligned columns / 12 probe length = 0.8333... (SW trims the mismatched tail)");
            Assert.That(gapped.OffTargetHits[0].HasGaps, Is.True,
                "The hit required an insertion");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_IdentityThreshold_GatesHits()
    {
        // SG1: The minIdentity threshold (Kane et al. 2000, default 0.75) gates hits.
        // The 0.8333-identity indel+mismatch site is admitted at 0.75 but rejected at 0.90.
        var references = new[] { "NNNNNNNNNN" + "ACGTACTGTACTT" + "NNNNN" };

        var lenient = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, references, minIdentity: 0.75);
        var strict = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, references, minIdentity: 0.90);

        Assert.Multiple(() =>
        {
            Assert.That(lenient.OffTargetCount, Is.EqualTo(1),
                "0.8333 identity >= 0.75 → admitted");
            Assert.That(strict.OffTargetCount, Is.EqualTo(0),
                "0.8333 identity < 0.90 → rejected");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_SpecificProbe_IsSpecificTrue()
    {
        // SG2: A probe with exactly one on-target and no off-targets reports IsSpecific.
        var references = new[] { "GGGGG" + GappedProbe + "GGGGG" };

        var gapped = ProbeDesigner.ScanOffTargetsGapped(GappedProbe, references);

        Assert.Multiple(() =>
        {
            Assert.That(gapped.OnTargetHits, Has.Count.EqualTo(1),
                "Single exact on-target");
            Assert.That(gapped.OffTargetCount, Is.EqualTo(0),
                "No off-targets in flanking poly-G");
            Assert.That(gapped.IsSpecific, Is.True,
                "Exactly one on-target and zero off-targets → specific");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_OnTargetPlusHighIdentityAndIndelOffTargets_SeparatedAndLowIdentityExcluded()
    {
        // SG3: full Kane (2000) separation cross-check on a hand-constructed target set:
        //   ref0 = exact on-target           → 1 on-target (identity 1.0, no gap)
        //   ref1 = 17/20 = 0.85 ungapped     → off-target ABOVE the 0.75 Kane threshold (mismatched, no gap)
        //   ref2 = ~scrambled (<0.75)        → must NOT be called (best local block far below 0.75)
        //   ref3 = one inserted base         → off-target reachable ONLY via a gap (identity 1.0, HasGaps)
        // Hand-derived with an independent Smith-Waterman (BlastDna +2/-3, gap -2):
        //   ref1 alignment ACGTGGCATTACGGCATTCA / ACATGGCATAACGGCAATCA → 17 identical / 20 = 0.85
        //   ref3 alignment ACGTGGCATT-ACGGCATTCA / ACGTGGCATTAACGGCATTCA → 20 identical / 20 = 1.0, one gap
        //   ref2 best local block only 4 identical columns / 20 = 0.20 (< 0.75) → rejected
        const string probe = "ACGTGGCATTACGGCATTCA"; // 20 nt
        var references = new[]
        {
            "TTTTT" + probe + "TTTTT",                          // exact on-target
            "GGGGG" + "ACATGGCATAACGGCAATCA" + "GGGGG",         // 0.85 mismatched off-target
            "CCCCC" + "TTAATTAATTAATTAATTAA" + "CCCCC",         // low-identity, must be excluded
            "AAAAA" + "ACGTGGCATTAACGGCATTCA" + "AAAAA",        // indel-only off-target (1 insertion)
        };

        var gapped = ProbeDesigner.ScanOffTargetsGapped(probe, references, minIdentity: 0.75);

        Assert.Multiple(() =>
        {
            Assert.That(gapped.OnTargetHits, Has.Count.EqualTo(1),
                "Exactly one intended on-target (the perfect exact match in ref0)");
            Assert.That(gapped.OnTargetHits[0].ReferenceIndex, Is.EqualTo(0),
                "The on-target lives in reference 0");

            Assert.That(gapped.OffTargetCount, Is.EqualTo(2),
                "Two genuine off-targets: the 0.85 mismatched site and the indel site");

            var hi = gapped.OffTargetHits.Single(h => h.ReferenceIndex == 1);
            Assert.That(hi.Identity, Is.EqualTo(17.0 / 20.0).Within(1e-10),
                "Mismatched off-target identity = 17/20 = 0.85 (Kane 2000: >0.75 cross-hybridizes)");
            Assert.That(hi.HasGaps, Is.False, "The 0.85 off-target is ungapped (pure substitutions)");

            var indel = gapped.OffTargetHits.Single(h => h.ReferenceIndex == 3);
            Assert.That(indel.Identity, Is.EqualTo(1.0).Within(1e-10),
                "Indel off-target: 20 identical aligned columns / 20 = 1.0");
            Assert.That(indel.HasGaps, Is.True,
                "Indel off-target reachable only via a gap (insertion) — missed by the ungapped Hamming scan");

            Assert.That(gapped.OffTargetHits.Any(h => h.ReferenceIndex == 2), Is.False,
                "Low-identity reference (best block 0.20 < 0.75) must NOT be called an off-target");
            Assert.That(gapped.IsSpecific, Is.False,
                "Two off-targets present → the probe is not specific");
        });
    }

    [Test]
    public void ScanOffTargetsGapped_NullProbe_ThrowsArgumentNullException()
    {
        // Guard: null probe must throw.
        Assert.Throws<ArgumentNullException>(() =>
            ProbeDesigner.ScanOffTargetsGapped(null!, IndelOffTargetReference),
            "Null probe should throw ArgumentNullException");
    }

    [Test]
    public void ScanOffTargetsGapped_NullReferences_ThrowsArgumentNullException()
    {
        // Guard: null references must throw.
        Assert.Throws<ArgumentNullException>(() =>
            ProbeDesigner.ScanOffTargetsGapped(GappedProbe, null!),
            "Null references should throw ArgumentNullException");
    }

    [Test]
    public void ScanOffTargetsGapped_EmptyProbe_ReturnsNoHits()
    {
        // Guard: empty probe yields no on/off-target hits.
        var gapped = ProbeDesigner.ScanOffTargetsGapped("", IndelOffTargetReference);

        Assert.Multiple(() =>
        {
            Assert.That(gapped.OnTargetHits, Is.Empty, "Empty probe → no on-target hits");
            Assert.That(gapped.OffTargetHits, Is.Empty, "Empty probe → no off-target hits");
        });
    }

    #endregion

    #region Karlin–Altschul E-value / bit-score for off-target hits (Must)

    // KA1–KA7 — Karlin–Altschul statistics (Karlin & Altschul 1990, PNAS 87:2264;
    // Altschul et al. 1990, J Mol Biol 215:403).
    // Verbatim formulas (retrieved 2026-06-24):
    //   E = K·m·n·e^{−λS};  S' = (λS − ln K)/ln 2;  E = m·n·2^{−S'};
    //   λ = unique positive root of Σ p_i p_j e^{λ s_ij} = 1.
    //   Sources: NCBI "The Statistics of Sequence Similarity Scores" (Altschul),
    //   https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html ; Durand,
    //   "BLAST (Karlin–Altschul) Statistics", CMU 03-711 (cites both 1990 papers).

    // +1/−3 nucleotide scoring scheme (the scheme NCBI blastn reports λ ≈ 1.37, K ≈ 0.711 for).
    private static readonly Seqeron.Genomics.Infrastructure.ScoringMatrix MatchMismatch1_3 =
        new(Match: 1, Mismatch: -3, GapOpen: -5, GapExtend: -2);

    // KA1 — λ for the +1/−3, uniform-0.25 scheme equals the published NCBI blastn value ≈ 1.374.
    // This pins the numeric root-solver to a sourced value: a wrong solver fails this assertion.
    [Test]
    public void ComputeLambdaNucleotide_Plus1Minus3_UniformFrequencies_MatchesPublishedValue()
    {
        // 0.25·e^{λ·1} + 0.75·e^{λ·(−3)} = 1 → λ ≈ 1.3740631 (NCBI blastn +1/−3).
        double lambda = ProbeDesigner.ComputeLambdaNucleotide(match: 1, mismatch: -3);

        Assert.That(lambda, Is.EqualTo(1.3740631224599755).Within(1e-6),
            "λ for +1/−3 with uniform 0.25 base frequencies must equal the published NCBI blastn value ≈ 1.374");
    }

    // KA2 — the solved λ actually satisfies the defining equation Σ p_i p_j e^{λ s_ij} = 1.
    [Test]
    public void ComputeLambdaNucleotide_SolvedRoot_SatisfiesDefiningEquation()
    {
        double lambda = ProbeDesigner.ComputeLambdaNucleotide(match: 1, mismatch: -3);

        // p(match)=0.25, p(mismatch)=0.75 for four equiprobable bases.
        double f = 0.25 * Math.Exp(lambda * 1) + 0.75 * Math.Exp(lambda * -3);

        Assert.That(f, Is.EqualTo(1.0).Within(1e-9),
            "The returned λ must be a root of Σ p_i p_j e^{λ s_ij} = 1 (Karlin & Altschul 1990)");
    }

    // KA3 — bit score and E-value for a hand-derived (S, m, n) with the +1/−3 scheme.
    // S=30, m=20, n=1000, K=0.711, λ=1.3740631224599755:
    //   S' = (λ·30 − ln 0.711)/ln 2 = 59.9627001142850
    //   E  = 0.711·20·1000·e^{−λ·30} = 1.78015836860839e-14
    [Test]
    public void ComputeKarlinAltschul_HandDerivedExample_MatchesBitScoreAndEValue()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 30, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Lambda, Is.EqualTo(1.3740631224599755).Within(1e-6),
                "λ must be the +1/−3 published value");
            Assert.That(stats.BitScore, Is.EqualTo(59.962700114285006).Within(1e-6),
                "S' = (λS − ln K)/ln 2 (Altschul et al. 1990)");
            Assert.That(stats.EValue, Is.EqualTo(1.7801583686083893e-14).Within(1e-24),
                "E = K·m·n·e^{−λS} (Karlin & Altschul 1990)");
        });
    }

    // KA4 — the two equivalent E-value forms agree: E = K·m·n·e^{−λS} = m·n·2^{−S'}.
    [Test]
    public void ComputeKarlinAltschul_EValue_EqualsSearchSpaceTimesTwoToMinusBitScore()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 18, queryLength: 25, databaseLength: 5000, scoring: MatchMismatch1_3, k: 0.711);

        double fromBits = (double)stats.QueryLength * stats.DatabaseLength * Math.Pow(2.0, -stats.BitScore);

        Assert.That(stats.EValue, Is.EqualTo(fromBits).Within(stats.EValue * 1e-9),
            "E = m·n·2^{−S'} must equal E = K·m·n·e^{−λS} (Altschul et al. 1990)");
    }

    // KA5 — E-value strictly decreases as the raw score increases (a better hit is less likely by chance).
    [Test]
    public void ComputeKarlinAltschul_EValue_DecreasesAsScoreIncreases()
    {
        var low = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 20, queryLength: 30, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);
        var high = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 21, queryLength: 30, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.That(high.EValue, Is.LessThan(low.EValue),
            "Higher score → lower E (E = K·m·n·e^{−λS} is monotonically decreasing in S)");
    }

    // KA6 — E-value scales linearly with the search space m·n (double n → double E).
    [Test]
    public void ComputeKarlinAltschul_EValue_ScalesLinearlyWithSearchSpace()
    {
        var baseStats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);
        var doubledN = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 2000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.That(doubledN.EValue, Is.EqualTo(2.0 * baseStats.EValue).Within(baseStats.EValue * 1e-9),
            "Doubling n doubles E (E is linear in the search space m·n; Karlin & Altschul 1990)");
    }

    // KA7 (guards) — the Karlin–Altschul preconditions and argument validation.
    [Test]
    public void ComputeLambdaNucleotide_NonPositiveMatch_Throws()
    {
        // No positive score → λ undefined (Altschul et al. 1990).
        Assert.Throws<ArgumentOutOfRangeException>(
            () => ProbeDesigner.ComputeLambdaNucleotide(match: 0, mismatch: -3),
            "A scheme with no positive score has no λ");
    }

    [Test]
    public void ComputeLambdaNucleotide_NonNegativeExpectedScore_Throws()
    {
        // match=3, mismatch=−1 → expected = 0.25·3 + 0.75·(−1) = 0 → not negative → λ undefined.
        Assert.Throws<ArgumentOutOfRangeException>(
            () => ProbeDesigner.ComputeLambdaNucleotide(match: 3, mismatch: -1),
            "Expected per-pair score must be negative for λ to be defined");
    }

    [Test]
    public void ComputeKarlinAltschul_NonPositiveLength_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ProbeDesigner.ComputeKarlinAltschul(10, 0, 1000, MatchMismatch1_3),
                "Query length m must be positive");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ProbeDesigner.ComputeKarlinAltschul(10, 20, 0, MatchMismatch1_3),
                "Database length n must be positive");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ProbeDesigner.ComputeKarlinAltschul(10, 20, 1000, MatchMismatch1_3, k: 0),
                "K must be positive");
        });
    }

    // KA8 — score S = 0 boundary: E collapses to the raw search space scaled by K (E = K·m·n·e^0 = K·m·n),
    // and the bit score is the pure normalization offset S' = −ln K / ln 2 (Karlin & Altschul 1990).
    // m=20, n=1000, K=0.711 → E = 0.711·20·1000 = 14220; S' = −ln(0.711)/ln 2 = 0.4920785350426718.
    [Test]
    public void ComputeKarlinAltschul_ScoreZero_EValueEqualsKTimesSearchSpace()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 0, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.Multiple(() =>
        {
            Assert.That(stats.EValue, Is.EqualTo(14220.0).Within(1e-6),
                "At S=0, E = K·m·n·e^0 = K·m·n (Karlin & Altschul 1990)");
            Assert.That(stats.BitScore, Is.EqualTo(0.4920785350426718).Within(1e-9),
                "At S=0, S' = (0 − ln K)/ln 2 = −ln K/ln 2 (Altschul et al. 1990)");
        });
    }

    // KA9 — the K parameter: E is linear in K (E = K·m·n·e^{−λS}), and the bit score shifts by −log2(K),
    // so doubling K doubles E and lowers the bit score by exactly 1 bit (Karlin & Altschul 1990; Altschul et al. 1990).
    [Test]
    public void ComputeKarlinAltschul_DoublingK_DoublesEValueAndLowersBitByOneBit()
    {
        var baseStats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);
        var doubledK = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 1.422);

        Assert.Multiple(() =>
        {
            Assert.That(doubledK.EValue, Is.EqualTo(2.0 * baseStats.EValue).Within(baseStats.EValue * 1e-9),
                "E is linear in K: doubling K doubles E (E = K·m·n·e^{−λS})");
            Assert.That(doubledK.BitScore, Is.EqualTo(baseStats.BitScore - 1.0).Within(1e-9),
                "S' = (λS − ln K)/ln 2: doubling K lowers the bit score by exactly log2(2)=1 bit");
        });
    }

    // KA10 — E increases with the query length m as well (E linear in the search space m·n; KA6 covered n).
    // Doubling m must double E, the symmetric counterpart of doubling n (Karlin & Altschul 1990).
    [Test]
    public void ComputeKarlinAltschul_DoublingQueryLength_DoublesEValue()
    {
        var baseStats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 20, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);
        var doubledM = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 22, queryLength: 40, databaseLength: 1000, scoring: MatchMismatch1_3, k: 0.711);

        Assert.That(doubledM.EValue, Is.EqualTo(2.0 * baseStats.EValue).Within(baseStats.EValue * 1e-9),
            "Doubling m doubles E (E is linear in the search space m·n; Karlin & Altschul 1990)");
    }

    // KA11 — λ root-finder convergence is precise: f(λ) = Σ p_i p_j e^{λ s_ij} − 1 must vanish to
    // near machine precision (the bisection runs 200 iterations on [0,100], far past double resolution).
    // Independently re-solved oracle: λ(1,−3, uniform 0.25) = 1.3740631224599755.
    [Test]
    public void ComputeLambdaNucleotide_RootFinder_ConvergesToMachinePrecision()
    {
        double lambda = ProbeDesigner.ComputeLambdaNucleotide(match: 1, mismatch: -3);
        double residual = 0.25 * Math.Exp(lambda * 1) + 0.75 * Math.Exp(lambda * -3) - 1.0;

        Assert.Multiple(() =>
        {
            Assert.That(lambda, Is.EqualTo(1.3740631224599755).Within(1e-15),
                "Bisection must converge to the exact double-precision root λ = 1.3740631224599755");
            Assert.That(residual, Is.EqualTo(0.0).Within(1e-12),
                "Residual of the defining equation must vanish to near machine precision");
        });
    }

    // KA12 — the DEFAULT scheme path (no scoring arg → BlastDna +2/−3) computes λ from that matrix
    // under uniform 0.25 frequencies: independently re-solved oracle λ(2,−3) = 0.6337314430979077.
    // NOTE: this is the UNIFORM-0.25 two-term root, not NCBI's published ungapped 2/−3 value (λ≈0.55,
    // K≈0.21), which derives from the full score lattice; the documented contract is the uniform-0.25 model.
    [Test]
    public void ComputeKarlinAltschul_DefaultScheme_UsesBlastDna2_3UniformLambda()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 30, queryLength: 20, databaseLength: 1000);

        Assert.That(stats.Lambda, Is.EqualTo(0.6337314430979077).Within(1e-9),
            "Default scoring is BlastDna (+2/−3); λ is its uniform-0.25 root 0.6337314430979077");
    }

    #endregion
}
