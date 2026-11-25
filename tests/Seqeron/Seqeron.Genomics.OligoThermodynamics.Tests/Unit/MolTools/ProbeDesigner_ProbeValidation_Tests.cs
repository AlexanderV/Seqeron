// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for ThermodynamicProbeDesigner.ValidateProbe and CheckSpecificity.
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



    [Test]
    public void ValidateProbe_EmptyProbe_ReturnsValidationResult()
    {
        // M1: Empty probe boundary condition
        // An empty probe has no sequence to hybridize — specificity must be 0.0
        // Source: Invariant #5 — offTargetHits == 0 → specificityScore == 0.0
        var validation = ThermodynamicProbeDesigner.ValidateProbe("", SingleMatchReference);

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
        var validation = ThermodynamicProbeDesigner.ValidateProbe(StandardProbe, Enumerable.Empty<string>());

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
            ThermodynamicProbeDesigner.ValidateProbe(StandardProbe, null!),
            "Null references should throw ArgumentNullException");
    }



    [Test]
    public void ValidateProbe_UniqueProbe_HasSpecificityScoreOne()
    {
        // M3: Unique probe (1 hit) should have specificity = 1.0
        var validation = ThermodynamicProbeDesigner.ValidateProbe(UniqueProbe, SingleMatchReference);

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
        var validation = ThermodynamicProbeDesigner.ValidateProbe(probe, MultipleMatchReference);

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



    [Test]
    public void ValidateProbe_HighSelfComplementarity_ReportsInIssues()
    {
        // M8: High self-complementarity (>30%) should be reported in issues
        // PalindromicProbe "GCGCGCGCGCGCGCGCGCGC" is its own reverse complement → selfComp = 1.0
        // 1.0 > default threshold 0.3 → issue must be generated
        var validation = ThermodynamicProbeDesigner.ValidateProbe(PalindromicProbe, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(validation.SelfComplementarity, Is.EqualTo(1.0),
                "GC-repeat palindrome: every position matches its reverse complement");
            Assert.That(validation.Issues, Has.Some.Contain("Self-complementarity"),
                "Issues must report self-complementarity when above threshold");
        });
    }



    [Test]
    public void ValidateProbe_MixedCaseProbe_HandledCaseInsensitively()
    {
        // M12: Case-insensitive probe handling
        string upperProbe = UniqueProbe.ToUpperInvariant();
        string lowerProbe = UniqueProbe.ToLowerInvariant();
        string mixedProbe = "AtCgAtCgAtCgAtCgAtCg";

        var validationUpper = ThermodynamicProbeDesigner.ValidateProbe(upperProbe, SingleMatchReference);
        var validationLower = ThermodynamicProbeDesigner.ValidateProbe(lowerProbe, SingleMatchReference);
        var validationMixed = ThermodynamicProbeDesigner.ValidateProbe(mixedProbe, SingleMatchReference);

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

    // M9 range test: deleted as duplicate of M10 (unique→1.0), M11 (multi→1/N), NoMatch (→0.0).



    [Test]
    public void ValidateProbe_PotentialHairpin_DetectsSecondaryStructure()
    {
        // S1: Secondary structure potential detected for hairpin sequences
        // Stem-loop: GCGC (stem, 4nt) + TTT (loop, 3nt) + GCGC (stem, 4nt) + filler
        // A 20-nt A/C/G/T probe gets the default Primer3 ntthal hairpin screen (Tm > 47 °C); the fallback
        // stem-loop screen PrimerDesigner.HasHairpinPotential (exact ≥ 4-bp stem, loop ≥ 3) flags it as well.
        string hairpinProbe = "GCGCTTTGCGCAAAAAAAAA"; // 20 chars

        var validation = ThermodynamicProbeDesigner.ValidateProbe(hairpinProbe, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(validation.HasSecondaryStructure, Is.True,
                "Hairpin stem GCGC-TTT-GCGC must be detected as secondary structure");
            Assert.That(validation.Issues, Has.Some.Contain("secondary structure"),
                "Issues must report secondary structure potential");
        });
    }



    [Test]
    public void ValidateProbe_ProblematicProbe_PopulatesIssuesList()
    {
        // S2: Issues list populated for problematic probes
        // 10-mer poly-A in 25-mer poly-A: 25-10+1 = 16 exact match positions → offTargetHits = 16
        // Implementation adds "{N} potential off-target sites" when offTargetHits > 1
        string probe = "AAAAAAAAAA";
        var references = new[] { "AAAAAAAAAAAAAAAAAAAAAAAAA" }; // 25 A's

        var validation = ThermodynamicProbeDesigner.ValidateProbe(probe, references);

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
        // "GCGCGCGCGC" (10-mer) → fold-back fraction 1.0 (palindrome)
        // In 32-char GC-repeat, only even positions match (odd positions are shifted by 1 → 10 mismatches)
        // Even positions 0,2,4,...,22 = 12 hits
        // Thermodynamic screen (Primer3 probe conditions 50 nM / 50 mM / 0 Mg / 0 dNTP; primer3-py 2.3.1):
        //   calc_homodimer Tm 52.763 °C and calc_hairpin Tm 55.851 °C both exceed 47 °C → 2 structure issues.
        // IsValid = no recorded issue → false
        string probe = "GCGCGCGCGC";
        var references = new[] { "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" }; // 32 chars

        var validation = ThermodynamicProbeDesigner.ValidateProbe(probe, references);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(12),
                "10-mer GC-repeat in 32-char GC-repeat: matches at 12 even positions");
            Assert.That(validation.SelfComplementarity, Is.EqualTo(1.0),
                "GC-repeat is its own reverse complement");
            Assert.That(validation.Issues.Count, Is.EqualTo(3),
                "Should report off-target + ntthal self-dimer + ntthal hairpin issues");
            Assert.That(validation.HasSecondaryStructure, Is.True, "ntthal hairpin Tm 55.85 °C > 47 °C");
            Assert.That(validation.IsValid, Is.False,
                "Issues recorded → IsValid must be false");
        });
    }



    [Test]
    public void ValidateProbe_ApproximateMatching_FindsNearMatches()
    {
        // S4: Approximate matching with maxMismatches works correctly
        // Reference has ACGAACGAACGAACGT at position 5 — differs from probe at positions 3,7,11 (3 mismatches)
        // With maxMismatches=0: no match (3 mismatches > 0)
        // With maxMismatches=3: 1 match → specificity = 1.0
        string probe = "ACGTACGTACGTACGT"; // 16-mer
        var references = new[] { "TTTTTACGAACGAACGAACGTTTTT" }; // near-match at pos 5

        var strict = ThermodynamicProbeDesigner.ValidateProbe(probe, references, maxMismatches: 0);
        var approx = ThermodynamicProbeDesigner.ValidateProbe(probe, references, maxMismatches: 3);

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



    [Test]
    public void ValidateProbe_AllInvariants_HoldForTypicalProbe()
    {
        // Combined invariant test for comprehensive coverage
        var validation = ThermodynamicProbeDesigner.ValidateProbe(StandardProbe, SingleMatchReference);

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

        var validation = ThermodynamicProbeDesigner.ValidateProbe(nonExistentProbe, references, maxMismatches: 0);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(0),
                "Probe should not match any site in references");
            Assert.That(validation.SpecificityScore, Is.EqualTo(0.0),
                "Invariant #5: zero hits must yield zero specificity");
        });
    }



    [Test]
    public void ValidateProbe_LongReference_FindsProbeCorrectly()
    {
        // C1: Long reference sequences handled correctly
        // UniqueProbe embedded at position 10_000 in 20_010-char poly-T reference
        string longRef = new string('T', 10_000) + UniqueProbe + new string('T', 10_000);

        var validation = ThermodynamicProbeDesigner.ValidateProbe(UniqueProbe, new[] { longRef }, maxMismatches: 0);

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

        var validation = ThermodynamicProbeDesigner.ValidateProbe(UniqueProbe, new[] { ref1, ref2, ref3 }, maxMismatches: 0);

        Assert.Multiple(() =>
        {
            Assert.That(validation.OffTargetHits, Is.EqualTo(3),
                "Hits should accumulate across all 3 references");
            Assert.That(validation.SpecificityScore, Is.EqualTo(1.0 / 3).Within(0.0001),
                "Specificity = 1.0/3 (Invariant #6)");
        });
    }



    [Test]
    public void DesignProbes_WithGenomeIndex_UsesCheckSpecificity()
    {
        // Verify that DesignProbes with suffix tree uses CheckSpecificity internally
        string uniqueRegion = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        string repeatedRegion = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        string genome = uniqueRegion + repeatedRegion + "AAAAAAAA" + repeatedRegion;

        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        var param = ThermodynamicProbeDesigner.Defaults.Microarray with { MinLength = 50, MaxLength = 52 };
        var probes = ThermodynamicProbeDesigner.DesignProbes(uniqueRegion, genomeIndex, param, maxProbes: 3, requireUnique: true).ToList();

        // All returned probes should be unique in the genome
        foreach (var probe in probes)
        {
            var positions = genomeIndex.FindAllOccurrences(probe.Sequence);
            Assert.That(positions.Count, Is.EqualTo(1),
                $"Probe should be unique but found {positions.Count} occurrences");
        }
    }

    [Test]
    public void ScanOffTargetsGapped_IndelOffTarget_FoundByGappedScanMissedByHammingScan()
    {
        // MG1: An off-target reachable ONLY through a single insertion is found by the gapped
        // scan but missed by the ungapped Hamming scan (Altschul 1990: gapped finds indels).
        // The indel region "ACGTACTGTACGT" has >=6 mismatches in every fixed 12-window, so the
        // default Hamming scan (maxMismatches=3) cannot reach it; only the exact on-target at 5.
        var hamming = ThermodynamicProbeDesigner.ValidateProbe(GappedProbe, IndelOffTargetReference, maxMismatches: 3);
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
    public void ComputeBlastnStatistics_ProbeVsSubject_AlignsWithLocalAlignAffineAndMatchesBlastn()
    {
        var s = ProbeDesigner.ComputeBlastnStatistics(EValueProbe, EValueSubject);
        var kane = ThermodynamicProbeDesigner.AssessCrossHybridization(EValueProbe, new[] { EValueSubject }, bothStrands: false)[0];

        Assert.Multiple(() =>
        {
            Assert.That(EValueSubject.Length, Is.EqualTo(279));
            Assert.That(s.RawScore, Is.EqualTo(66));
            Assert.That(kane.AlignmentScore, Is.EqualTo(66), "same canonical LocalAlignAffine score");
            Assert.That(s.EffectiveSearchSpace, Is.EqualTo(8672));
            Assert.That(s.EValue, Is.EqualTo(4.327686048582086e-15).Within(1e-24));
            Assert.That(s.BitScore, Is.EqualTo(60.79747462182638).Within(1e-9));
        });
    }



    // primer3-py 2.3.1 calc_homodimer / calc_end_stability / calc_hairpin Tm at the Primer3 probe conditions
    // (mv 50 mM, dv 0, dntp 0, dna 50 nM) — the ValidateProbe default (Defaults.Microarray).

    [Test]
    public void ValidateProbe_ThermodynamicScreen_ReportsPrimer3NtthalTm()
    {
        var v = ThermodynamicProbeDesigner.ValidateProbe(PalindromicProbe, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(v.ThermodynamicScreen, Is.True);
            Assert.That(v.SelfDimerTm!.Value, Is.EqualTo(78.85652531616256).Within(1e-9));
            Assert.That(v.SelfEndDimerTm!.Value, Is.EqualTo(78.85652531616256).Within(1e-9));
            Assert.That(v.HairpinTm!.Value, Is.EqualTo(87.30265612393043).Within(1e-9));
            Assert.That(v.HasSecondaryStructure, Is.True, "hairpin Tm 87.3 °C > PRIMER_INTERNAL_MAX_HAIRPIN_TH 47 °C");
            Assert.That(v.Issues, Has.Some.StartsWith("Self-complementarity: ntthal self-dimer Tm 78.9"));
            Assert.That(v.Issues, Has.Some.StartsWith("Potential secondary structure formation: ntthal hairpin Tm 87.3"));
            Assert.That(v.IsValid, Is.False);
        });
    }

    [Test]
    public void ValidateProbe_ThermodynamicScreen_UsesStatedConditions()
    {
        // calc_homodimer / calc_end_stability / calc_hairpin at mv 100, dv 2, dntp 0.2, dna 250:
        // 69.17069845823409 / 69.17069845823409 / 74.99462150250321.
        var conditions = ThermodynamicProbeDesigner.Defaults.Microarray with
        {
            MonovalentMillimolar = 100, DivalentMillimolar = 2.0, DntpMillimolar = 0.2, DnaConcentrationNanomolar = 250,
        };
        var v = ThermodynamicProbeDesigner.ValidateProbe("ACGTACGTACGTACGTACGTACGT", Enumerable.Empty<string>(), conditions: conditions);

        Assert.Multiple(() =>
        {
            Assert.That(v.SelfDimerTm!.Value, Is.EqualTo(69.17069845823409).Within(1e-9));
            Assert.That(v.SelfEndDimerTm!.Value, Is.EqualTo(69.17069845823409).Within(1e-9));
            Assert.That(v.HairpinTm!.Value, Is.EqualTo(74.99462150250321).Within(1e-9));
        });
    }

    [Test]
    public void ValidateProbe_HighFoldBackFractionWithoutStableStructure_PassesThermodynamicScreen()
    {
        // Fold-back fraction 16/25 = 0.64 > 0.3, but primer3-py: homodimer Tm −6.43, end −99.94, hairpin 0
        // (no stable structure; Primer3 reports negative Tm as 0) → no self-structure issue.
        const string probe = "CTAGAAATGCTGTCGGGACTTCTAC";
        var thermo = ThermodynamicProbeDesigner.ValidateProbe(probe, Enumerable.Empty<string>());
        var heuristic = ThermodynamicProbeDesigner.ValidateProbe(probe, Enumerable.Empty<string>(),
            conditions: ThermodynamicProbeDesigner.Defaults.Microarray with { StructureScreen = ProbeDesigner.ProbeStructureScreen.Heuristic });

        Assert.Multiple(() =>
        {
            Assert.That(thermo.SelfComplementarity, Is.EqualTo(0.64).Within(1e-12));
            Assert.That(thermo.SelfDimerTm, Is.EqualTo(0.0));
            Assert.That(thermo.SelfEndDimerTm, Is.EqualTo(0.0));
            Assert.That(thermo.HairpinTm, Is.EqualTo(0.0));
            Assert.That(thermo.Issues, Has.None.Contain("Self-complementarity"));
            Assert.That(thermo.HasSecondaryStructure, Is.False);

            // Fallback self-dimer criterion = Primer3 alignment-mode internal-oligo oligo_compl (dpal.c, compiled):
            // self_any 9.00, self_end 7.00 ≤ PRIMER_INTERNAL_MAX_SELF_ANY/_END 12.00 → no issue (the fold-back
            // fraction 0.64 is only a reported library metric now).
            Assert.That(heuristic.ThermodynamicScreen, Is.False);
            Assert.That(heuristic.SelfDimerTm, Is.Null);
            Assert.That(heuristic.SelfComplementarity, Is.EqualTo(0.64).Within(1e-12));
            Assert.That(heuristic.SelfAny, Is.EqualTo(9.0));
            Assert.That(heuristic.SelfEnd, Is.EqualTo(7.0));
            Assert.That(thermo.SelfAny, Is.EqualTo(9.0), "alignment-mode values are reported for every probe");
            Assert.That(heuristic.Issues, Has.None.Contain("Self-complementarity"));
        });
    }

    [Test]
    public void ValidateProbe_ProbeLongerThan60nt_UsesPrimer3AlignmentSelfDimerFallback()
    {
        // thal.c THAL_MAX_ALIGN = 60: a 64-nt probe has no ntthal self-structure; the fallback self-dimer criterion is
        // Primer3 alignment-mode oligo_compl (no length limit): (ACGT)16 is its own reverse complement, dpal.c
        // self_any = self_end = 64.00 > PRIMER_INTERNAL_MAX_SELF_ANY 12.00.
        string probe = string.Concat(Enumerable.Repeat("ACGT", 16));
        var v = ThermodynamicProbeDesigner.ValidateProbe(probe, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(v.ThermodynamicScreen, Is.False);
            Assert.That(v.SelfDimerTm, Is.Null);
            Assert.That(v.HairpinTm, Is.Null);
            Assert.That(v.SelfComplementarity, Is.EqualTo(1.0));
            Assert.That(v.SelfAny, Is.EqualTo(64.0));
            Assert.That(v.SelfEnd, Is.EqualTo(64.0));
            Assert.That(v.Issues, Has.Some.EqualTo("Self-complementarity: Primer3 self_any 64.00 exceeds 12.00"));
        });
    }

    [Test]
    public void AssessCrossHybridization_MatchesBiopythonLocalAlignmentAndLcs()
    {
        var r = ThermodynamicProbeDesigner.AssessCrossHybridization(KaneProbe,
            new[] { KaneNonTargetA, KaneNonTargetB, KaneNonTargetC, KaneNonTargetD, KaneNonTargetE });

        Assert.That(r, Has.Count.EqualTo(10), "5 non-targets × 2 strands");
        // (index, reverse, score, identical, lcs) — Biopython: A fwd 53/40/6, A rc 12/–/6, B fwd 37/28/18, B rc 12/6/6,
        // C fwd 10/5/5, C rc 16/8/8, D fwd 12/6/6, D rc 100/50/50, E fwd 30/15/15, E rc 10/5/5.
        var expected = new (int Score, int? Identical, int Lcs)[]
        {
            (53, 40, 6), (12, null, 6), (37, 28, 18), (12, 6, 6), (10, 5, 5),
            (16, 8, 8), (12, 6, 6), (100, 50, 50), (30, 15, 15), (10, 5, 5),
        };
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                Assert.That(r[k].NonTargetIndex, Is.EqualTo(k / 2));
                Assert.That(r[k].ReverseComplementStrand, Is.EqualTo(k % 2 == 1));
                Assert.That(r[k].AlignmentScore, Is.EqualTo(expected[k].Score), $"score {k}");
                Assert.That(r[k].LongestContiguousMatch, Is.EqualTo(expected[k].Lcs), $"lcs {k}");
                if (expected[k].Identical is int id)
                {
                    Assert.That(r[k].IdenticalColumns, Is.EqualTo(id), $"identical {k}");
                    Assert.That(r[k].Identity, Is.EqualTo(id / 50.0).Within(1e-12));
                }
                else
                {
                    // Biopython co-optimal alignments have 9 or 6 identical columns.
                    Assert.That(r[k].IdenticalColumns, Is.AnyOf(9, 6), $"identical {k}");
                }
            }
        });
    }

    [Test]
    public void AssessCrossHybridization_AppliesKaneThresholds()
    {
        var r = ThermodynamicProbeDesigner.AssessCrossHybridization(KaneProbe,
            new[] { KaneNonTargetA, KaneNonTargetB, KaneNonTargetC, KaneNonTargetD, KaneNonTargetE });

        Assert.Multiple(() =>
        {
            // A: identity 0.80 > 0.75 (contiguous 6).
            Assert.That(r[0].ExceedsIdentityThreshold, Is.True);
            Assert.That(r[0].ExceedsContiguousThreshold, Is.False);
            // B: identity 0.56 but an 18-nt identical stretch > 15.
            Assert.That(r[2].ExceedsIdentityThreshold, Is.False);
            Assert.That(r[2].ExceedsContiguousThreshold, Is.True);
            Assert.That(r[2].CrossHybridizes, Is.True);
            // C: unrelated.
            Assert.That(r[4].CrossHybridizes || r[5].CrossHybridizes, Is.False);
            // D: probe's reverse complement → reverse strand identical.
            Assert.That(r[6].CrossHybridizes, Is.False);
            Assert.That(r[7].Identity, Is.EqualTo(1.0));
            Assert.That(r[7].CrossHybridizes, Is.True);
            // E: exactly 15 contiguous identical bases is not > 15.
            Assert.That(r[8].LongestContiguousMatch, Is.EqualTo(15));
            Assert.That(r[8].CrossHybridizes, Is.False);
        });

        // Thresholds are strict: identity 0.80 is not > 0.80; 18 nt is > 17 but not > 18.
        var strict = ThermodynamicProbeDesigner.AssessCrossHybridization(KaneProbe, new[] { KaneNonTargetA, KaneNonTargetB },
            maxIdentity: 0.80, maxContiguousMatch: 18, bothStrands: false);
        Assert.Multiple(() =>
        {
            Assert.That(strict, Has.Count.EqualTo(2));
            Assert.That(strict[0].ExceedsIdentityThreshold, Is.False);
            Assert.That(strict[1].ExceedsContiguousThreshold, Is.False);
        });
    }

    [Test]
    public void AssessCrossHybridization_SiteDuplexTm_MatchesPrimer3CalcHeterodimer()
    {
        // Non-target A: best local alignment covers strand[20..68] (49 nt; the mismatched last probe base is trimmed).
        // primer3-py 2.3.1 calc_heterodimer(probe, revcomp(site), mv 50, dv 0, dntp 0, dna 50).tm = 36.11423712379826;
        // non-target D (reverse strand = the probe itself at 15..64): 66.04038852959525.
        var r = ThermodynamicProbeDesigner.AssessCrossHybridization(KaneProbe, new[] { KaneNonTargetA, KaneNonTargetD });

        Assert.Multiple(() =>
        {
            Assert.That((r[0].SiteStart, r[0].SiteEnd), Is.EqualTo((20, 68)));
            Assert.That(r[0].DuplexTm!.Value, Is.EqualTo(36.11423712379826).Within(1e-9));
            Assert.That((r[3].SiteStart, r[3].SiteEnd), Is.EqualTo((15, 64)));
            Assert.That(r[3].DuplexTm!.Value, Is.EqualTo(66.04038852959525).Within(1e-9));
        });

        // Optional OligoArray-style duplex-Tm threshold (identity criterion relaxed to isolate it): 36.11 °C > 30 °C.
        var t30 = ThermodynamicProbeDesigner.AssessCrossHybridization(KaneProbe, new[] { KaneNonTargetA }, maxIdentity: 0.9,
            bothStrands: false, maxDuplexTm: 30);
        var t40 = ThermodynamicProbeDesigner.AssessCrossHybridization(KaneProbe, new[] { KaneNonTargetA }, maxIdentity: 0.9,
            bothStrands: false, maxDuplexTm: 40);
        var v = ThermodynamicProbeDesigner.ValidateProbe(KaneProbe, Array.Empty<string>(), nonTargetSequences: new[] { KaneNonTargetA },
            maxNonTargetIdentity: 0.9, maxDuplexTm: 30);
        Assert.Multiple(() =>
        {
            Assert.That(t30[0].ExceedsIdentityThreshold, Is.False);
            Assert.That(t30[0].ExceedsDuplexTmThreshold, Is.True);
            Assert.That(t30[0].CrossHybridizes, Is.True);
            Assert.That(t40[0].CrossHybridizes, Is.False);
            Assert.That(v.Issues, Has.Some.EqualTo(
                "Cross-hybridization risk with non-target 0: identity 80%, longest contiguous match 6 nt (Kane 2000), site duplex Tm 36.1°C > 30°C"));
        });

        // 62-nt probe, 49-nt site [20, 68]: thal.c thal_check_errors refuses only when BOTH strands exceed THAL_MAX_ALIGN
        // (60), so the duplex Tm is computed (A3-27; primer3-py calc_heterodimer(probe, revcomp(site)) at 50/0/0/50 nM:
        // 36.11423712379826); empty strand → no site.
        string longProbe = KaneProbe + "ACGTACGTACGT";
        var l = ThermodynamicProbeDesigner.AssessCrossHybridization(longProbe, new[] { KaneNonTargetA, "" }, bothStrands: false);
        Assert.Multiple(() =>
        {
            Assert.That((l[0].SiteStart, l[0].SiteEnd), Is.EqualTo((20, 68)));
            Assert.That(l[0].DuplexTm, Is.EqualTo(36.11423712379826).Within(1e-9));
            Assert.That((l[1].SiteStart, l[1].SiteEnd, l[1].DuplexTm), Is.EqualTo((-1, -1, (double?)null)));
        });
    }

    [Test]
    public void AssessCrossHybridization_LongNonTarget_ChunkedScoreEqualsCanonicalWholeStrandAlignment()
    {
        // 9 kb strand (> 4096-nt chunk) with a mutated probe copy straddling the first chunk boundary.
        var rng = new Random(11);
        char[] bases = { 'A', 'C', 'G', 'T' };
        string Random(int n) => new string(Enumerable.Range(0, n).Select(_ => bases[rng.Next(4)]).ToArray());
        char[] mutated = KaneProbe.ToCharArray();
        mutated[7] = 'G'; mutated[22] = 'A';
        string strand = Random(4070) + new string(mutated) + Random(4880);

        var r = ThermodynamicProbeDesigner.AssessCrossHybridization(KaneProbe, new[] { strand }, bothStrands: false)[0];
        var whole = Seqeron.Genomics.Alignment.SequenceAligner.LocalAlignAffine(
            KaneProbe, strand, Seqeron.Genomics.Alignment.SequenceAligner.BlastDna);

        Assert.That(r.AlignmentScore, Is.EqualTo(whole.Score));
        Assert.That(r.IdenticalColumns, Is.GreaterThanOrEqualTo(48));
    }

    [Test]
    public void AssessCrossHybridization_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => ThermodynamicProbeDesigner.AssessCrossHybridization(null!, new[] { "ACGT" }));
        Assert.Throws<ArgumentNullException>(() => ThermodynamicProbeDesigner.AssessCrossHybridization("ACGT", null!));
        Assert.Throws<ArgumentException>(() => ThermodynamicProbeDesigner.AssessCrossHybridization("", new[] { "ACGT" }));
        Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicProbeDesigner.AssessCrossHybridization("ACGT", new[] { "ACGT" }, maxIdentity: 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicProbeDesigner.AssessCrossHybridization("ACGT", new[] { "ACGT" }, maxContiguousMatch: -1));
        var empty = ThermodynamicProbeDesigner.AssessCrossHybridization("ACGT", new[] { "" });
        Assert.That(empty.All(x => x.AlignmentScore == 0 && x.LongestContiguousMatch == 0 && !x.CrossHybridizes), Is.True);
    }

    [Test]
    public void ValidateProbe_WithNonTargets_RecordsKaneCrossHybridizationIssues()
    {
        var v = ThermodynamicProbeDesigner.ValidateProbe(KaneProbe, new[] { "GGGG" + KaneProbe + "GGGG" },
            nonTargetSequences: new[] { KaneNonTargetA, KaneNonTargetB, KaneNonTargetC });
        var none = ThermodynamicProbeDesigner.ValidateProbe(KaneProbe, new[] { "GGGG" + KaneProbe + "GGGG" },
            nonTargetSequences: new[] { KaneNonTargetC });

        Assert.Multiple(() =>
        {
            Assert.That(v.CrossHybridization, Has.Count.EqualTo(6));
            Assert.That(v.Issues.Count(i => i.StartsWith("Cross-hybridization risk")), Is.EqualTo(2));
            Assert.That(v.Issues, Has.Some.EqualTo(
                "Cross-hybridization risk with non-target 0: identity 80%, longest contiguous match 6 nt (Kane 2000)"));
            Assert.That(v.Issues, Has.Some.EqualTo(
                "Cross-hybridization risk with non-target 1: identity 56%, longest contiguous match 18 nt (Kane 2000)"));
            Assert.That(v.IsValid, Is.False);
            Assert.That(none.Issues, Has.None.StartsWith("Cross-hybridization risk"));
            Assert.That(ThermodynamicProbeDesigner.ValidateProbe(KaneProbe, Array.Empty<string>()).CrossHybridization, Is.Empty);
        });
    }

    [Test]
    public void ValidateProbe_ShortProbe_ThreeMismatchSiteAtSeventyFivePercent_IsNotAnOffTargetIssue()
    {
        // 12-mer: exact site (1.0), a 3-mismatch site (9/12 = 0.75, longest run 3 — no Kane criterion) and a
        // 2-mismatch site (10/12 = 0.8333 > 0.75). Ungapped hits 3 → uniqueness score 1/3 (library convention).
        var noIssue = ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe12,
            new[] { "TTTTTGATCCGACGCTATTTTTGCTCCTACGGTATTTTT" });
        var issue = ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe12,
            new[] { "TTTTTGATCCGACGCTATTTTTGCTCCTACGGTATTTTT", "GGGGGGAACCGAGGCTAGGGGG" });

        Assert.Multiple(() =>
        {
            Assert.That(noIssue.OffTargetHits, Is.EqualTo(2));
            Assert.That(noIssue.SpecificityScore, Is.EqualTo(0.5).Within(1e-12));
            Assert.That(noIssue.CrossHybridizingHits, Is.EqualTo(1));
            Assert.That(noIssue.Issues, Has.None.Contain("off-target"));
            Assert.That(issue.OffTargetHits, Is.EqualTo(3));
            Assert.That(issue.SpecificityScore, Is.EqualTo(1.0 / 3).Within(1e-12));
            Assert.That(issue.CrossHybridizingHits, Is.EqualTo(2));
            Assert.That(issue.Issues, Has.Some.EqualTo(
                "2 potential off-target sites (Kane 2000: identity > 75% or > 15 contiguous identical nt)"));
        });
    }

    [Test]
    public void ValidateProbe_ReferenceSites_FollowKaneIdentityAndContiguityCriteria()
    {
        var k3 = ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, new[] { KaneSiteReference40 });
        var k12 = ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, new[] { KaneSiteReference40 }, maxMismatches: 12);
        var identity069 = ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, new[] { KaneSiteReference40 }, maxMismatches: 12,
            maxNonTargetIdentity: 0.69);
        var contiguous28 = ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, new[] { KaneSiteReference40 }, maxMismatches: 12,
            maxContiguousMatch: 28);

        Assert.Multiple(() =>
        {
            Assert.That((k3.OffTargetHits, k3.CrossHybridizingHits), Is.EqualTo((1, 1)), "3 mismatches: exact site only");
            // 12 mismatches: the clustered site meets the contiguity criterion (28 > 15), the spread one neither.
            Assert.That((k12.OffTargetHits, k12.CrossHybridizingHits), Is.EqualTo((3, 2)));
            Assert.That(k12.Issues, Has.Some.StartWith("2 potential off-target sites"));
            Assert.That((identity069.OffTargetHits, identity069.CrossHybridizingHits), Is.EqualTo((3, 3)), "0.70 > 0.69");
            Assert.That(identity069.Issues, Has.Some.EqualTo(
                "3 potential off-target sites (Kane 2000: identity > 69% or > 15 contiguous identical nt)"));
            Assert.That((contiguous28.OffTargetHits, contiguous28.CrossHybridizingHits), Is.EqualTo((3, 1)), "28 is not > 28");
            Assert.That(contiguous28.Issues, Has.None.Contain("off-target"));
            Assert.That(contiguous28.SpecificityScore, Is.EqualTo(1.0 / 3).Within(1e-12), "uniqueness score is reported only");
        });
    }

    [Test]
    public void ValidateProbe_InvalidKaneThresholdsOrConditions_Throw()
    {
        string[] refs = { KaneSiteReference40 };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, refs, maxNonTargetIdentity: 1.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, refs, maxNonTargetIdentity: double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, refs, maxContiguousMatch: -1));
            // Primer3 _pr_data_control: internal-oligo salt / DNA concentration > 0, divalent / dNTP >= 0.
            var p = ThermodynamicProbeDesigner.Defaults.Microarray;
            Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, refs, conditions: p with { MonovalentMillimolar = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, refs, conditions: p with { DnaConcentrationNanomolar = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, refs, conditions: p with { DivalentMillimolar = -1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicProbeDesigner.ValidateProbe(KaneSiteProbe40, refs, conditions: p with { DntpMillimolar = -1 }));
        });
    }

    [Test]
    public void ValidateProbe_BothStrands_CountsReverseComplementSitesByKaneCriteria()
    {
        string[] refs = { BothStrandsReference };
        var single3 = ThermodynamicProbeDesigner.ValidateProbe(BothStrandsProbe, refs);
        var both3 = ThermodynamicProbeDesigner.ValidateProbe(BothStrandsProbe, refs, bothStrands: true);
        var single5 = ThermodynamicProbeDesigner.ValidateProbe(BothStrandsProbe, refs, maxMismatches: 5);
        var both5 = ThermodynamicProbeDesigner.ValidateProbe(BothStrandsProbe, refs, maxMismatches: 5, bothStrands: true);

        Assert.Multiple(() =>
        {
            // Default (backward compatible): the given strand only — the exact site.
            Assert.That((single3.OffTargetHits, single3.CrossHybridizingHits), Is.EqualTo((1, 1)));
            Assert.That(single3.Issues, Has.None.Contain("off-target"));
            // Both strands: + the 2-mismatch reverse-complement site (18/20 = 0.90 > 0.75) → off-target issue.
            Assert.That((both3.OffTargetHits, both3.CrossHybridizingHits), Is.EqualTo((2, 2)));
            Assert.That(both3.SpecificityScore, Is.EqualTo(0.5).Within(1e-12));
            Assert.That(both3.Issues, Has.Some.EqualTo(
                "2 potential off-target sites (Kane 2000: identity > 75% or > 15 contiguous identical nt)"));
            Assert.That(both3.IsValid, Is.False);
            // Radius 5 reaches the 5-mismatch reverse-complement site (15/20 = 0.75, longest run 3 — no Kane criterion).
            Assert.That((single5.OffTargetHits, single5.CrossHybridizingHits), Is.EqualTo((1, 1)));
            Assert.That((both5.OffTargetHits, both5.CrossHybridizingHits), Is.EqualTo((3, 2)));
        });
    }

    [Test]
    public void ValidateProbe_BothStrands_PalindromicProbeOrSiteCountedOnce()
    {
        // Reverse-palindromic probe (rc = itself): its reverse-complement hits are the same sites — counted once,
        // as CheckSpecificity(bothStrands) counts a palindromic probe once.
        const string palindrome = "GAATTCCGGAATTC";
        Assert.That(DnaSequence.GetReverseComplementString(palindrome), Is.EqualTo(palindrome));
        // Non-palindromic probe = the reverse-palindromic site GACGTCAGCTGACGTC with position 0 substituted: the site is
        // 1 mismatch from the probe AND from its reverse complement (same position) → one site, not two.
        const string palindromicSite = "GACGTCAGCTGACGTC";
        const string nearPalindromicProbe = "TACGTCAGCTGACGTC";
        Assert.That(DnaSequence.GetReverseComplementString(palindromicSite), Is.EqualTo(palindromicSite));

        Assert.Multiple(() =>
        {
            foreach (int k in new[] { 0, 2 })
            {
                var single = ThermodynamicProbeDesigner.ValidateProbe(palindrome, new[] { "TTTTT" + palindrome + "AAAAA" }, k);
                var both = ThermodynamicProbeDesigner.ValidateProbe(palindrome, new[] { "TTTTT" + palindrome + "AAAAA" }, k, bothStrands: true);
                Assert.That((single.OffTargetHits, single.CrossHybridizingHits), Is.EqualTo((1, 1)), $"palindromic probe, k {k}");
                Assert.That((both.OffTargetHits, both.CrossHybridizingHits), Is.EqualTo((1, 1)), $"palindromic probe both strands, k {k}");
            }
            string[] siteRef = { "CCCCC" + palindromicSite + "CCCCC" };
            var siteSingle = ThermodynamicProbeDesigner.ValidateProbe(nearPalindromicProbe, siteRef, maxMismatches: 1);
            var siteBoth = ThermodynamicProbeDesigner.ValidateProbe(nearPalindromicProbe, siteRef, maxMismatches: 1, bothStrands: true);
            Assert.That((siteSingle.OffTargetHits, siteSingle.CrossHybridizingHits), Is.EqualTo((1, 1)));
            Assert.That((siteBoth.OffTargetHits, siteBoth.CrossHybridizingHits), Is.EqualTo((1, 1)), "palindromic site once");
            Assert.That(siteBoth.SpecificityScore, Is.EqualTo(1.0));
        });
    }

    [TestCaseSource(nameof(BothStrandsOracleCases))]
    public void ValidateProbe_BothStrands_MatchesPythonBruteForce(
        string probe, string[] references, int maxMismatches,
        int singleHits, int singleKane, int bothHits, int bothKane)
    {
        var single = ThermodynamicProbeDesigner.ValidateProbe(probe, references, maxMismatches);
        var both = ThermodynamicProbeDesigner.ValidateProbe(probe, references, maxMismatches, bothStrands: true);
        Assert.Multiple(() =>
        {
            Assert.That((single.OffTargetHits, single.CrossHybridizingHits), Is.EqualTo((singleHits, singleKane)), "given strand");
            Assert.That((both.OffTargetHits, both.CrossHybridizingHits), Is.EqualTo((bothHits, bothKane)), "both strands");
            Assert.That(both.SpecificityScore, Is.EqualTo(bothHits == 0 ? 0.0 : 1.0 / bothHits).Within(1e-12));
            Assert.That(both.Issues.Any(i => i.Contains("off-target")), Is.EqualTo(bothKane > 1));
        });
    }

    [Test]
    public void ValidateProbe_LongProbeFallback_FlagsStemLoopWithFourNtLoop()
    {
        var fallback = ThermodynamicProbeDesigner.ValidateProbe(StemLoop74, new[] { StemLoop74 });
        var noTargets = new ThermodynamicProbeDesigner.ProbeParameters(20, 120, -1000, 1000, 0, 1, 100, true, 0.3);
        // Opt-in ntthal (F55): thal.c (primer3-py 2.3.1 sources) compiled with -DTHAL_MAX_ALIGN=10000, calc_hairpin
        // arguments: 50 mM / 0 / 0 / 50 nM, 37 °C, max loop 30 → Tm 77.95325865166825 °C (dG −10369.44 cal/mol; the same
        // value primer3-py calc_hairpin gives for the 60-nt window ACACAACC…CACCCC); at the Microarray preset's
        // OligoArray 1 M / 1 µM → 92.272558368972682 °C.
        var thermo = ThermodynamicProbeDesigner.ValidateProbe(StemLoop74, new[] { StemLoop74 },
            conditions: noTargets with { ThermodynamicScreenMaxLength = 74 });
        var thermoMicroarray = ThermodynamicProbeDesigner.ValidateProbe(StemLoop74, new[] { StemLoop74 },
            conditions: ThermodynamicProbeDesigner.Defaults.Microarray with { ThermodynamicScreenMaxLength = 74 });
        Assert.Multiple(() =>
        {
            Assert.That(fallback.ThermodynamicScreen, Is.False, "74 nt > THAL_MAX_ALIGN 60 → fallback screens");
            Assert.That(fallback.HasSecondaryStructure, Is.True, "10-bp stem + 4-nt loop must be flagged (was False)");
            // A5-5: the fallback stem flag is a warning, not an issue.
            Assert.That(fallback.Warnings, Is.EqualTo(new[] { "Potential secondary structure formation" }));
            Assert.That(fallback.Issues, Has.None.Contain("secondary structure"));
            Assert.That(thermo.HairpinTm, Is.EqualTo(77.95325865166825).Within(1e-6));
            Assert.That(thermo.HasSecondaryStructure, Is.True);
            Assert.That(thermoMicroarray.HairpinTm, Is.EqualTo(92.272558368972682).Within(1e-6));
        });
    }

    [Test]
    public void ValidateProbe_LongProbeFallback_StemFlagIsWarning_NotValidityCriterion()
    {
        // Audit round 5, A5-5: validity rests on sourced criteria only (F52). The fallback stem-loop presence test
        // (stem ≥ 4 bp, library convention) is reported but does not make IsValid false.
        var fallback = ThermodynamicProbeDesigner.ValidateProbe(Random61, new[] { Random61 });
        var heuristic = ThermodynamicProbeDesigner.ValidateProbe(Random61, new[] { Random61 },
            conditions: ThermodynamicProbeDesigner.Defaults.Microarray with
            {
                MonovalentMillimolar = 50, DnaConcentrationNanomolar = 50,
                StructureScreen = ProbeDesigner.ProbeStructureScreen.Heuristic,
            });
        var thermo = ThermodynamicProbeDesigner.ValidateProbe(Random61, new[] { Random61 },
            conditions: ThermodynamicProbeDesigner.Defaults.Microarray with
            {
                MonovalentMillimolar = 50, DivalentMillimolar = 0, DntpMillimolar = 0, DnaConcentrationNanomolar = 50,
                ThermodynamicScreenMaxLength = 61,
            });
        Assert.Multiple(() =>
        {
            Assert.That(fallback.ThermodynamicScreen, Is.False, "61 nt > THAL_MAX_ALIGN 60");
            Assert.That(fallback.HasSecondaryStructure, Is.True, "exact 4-bp stem + ≥ 3-nt loop present");
            Assert.That(fallback.Warnings, Is.EqualTo(new[] { "Potential secondary structure formation" }));
            Assert.That(fallback.SelfAny, Is.EqualTo(7.0), "independent Smith–Waterman vs reverse complement (+1/−1, gap −2): 7");
            Assert.That(fallback.SelfEnd, Is.LessThanOrEqualTo(12.0));
            Assert.That(fallback.Issues, Is.Empty);
            Assert.That(fallback.IsValid, Is.True, "was false before A5-5 (the stem flag was an issue)");
            Assert.That(heuristic.HasSecondaryStructure, Is.True);
            Assert.That(heuristic.IsValid, Is.True);
            // Opt-in ntthal (F55): the sourced hairpin decision, 33.98 °C ≤ PRIMER_INTERNAL_MAX_HAIRPIN_TH 47 °C.
            Assert.That(thermo.ThermodynamicScreen, Is.True);
            Assert.That(thermo.HairpinTm, Is.EqualTo(33.980529935122263).Within(1e-6));
            Assert.That(thermo.HasSecondaryStructure, Is.False);
            Assert.That(thermo.Warnings, Is.Empty);
            Assert.That(thermo.IsValid, Is.True);
        });
    }

    [Test]
    public void ValidateProbe_LongProbeFallback_StemLoopRepro_StemFlagIsNotAnIssue_NtthalHairpinDecidesWithOptIn()
    {
        // F61 repro (74 nt, 10-bp stem + TTTT loop): by default the fallback still reports the structure as a warning;
        // the probe is invalid only by sourced criteria — Primer3 alignment-mode self_any 16.00 > 12.00 (independent
        // Python Smith–Waterman of the probe vs its reverse complement, +1/−1, gap −2: 16) and, with the opt-in, the
        // ntthal hairpin Tm (thal.c 77.95325865166825 °C > 47 °C).
        var fallback = ThermodynamicProbeDesigner.ValidateProbe(StemLoop74, new[] { StemLoop74 });
        var thermo = ThermodynamicProbeDesigner.ValidateProbe(StemLoop74, new[] { StemLoop74 },
            conditions: ThermodynamicProbeDesigner.Defaults.Microarray with
            {
                MonovalentMillimolar = 50, DivalentMillimolar = 0, DntpMillimolar = 0, DnaConcentrationNanomolar = 50,
                ThermodynamicScreenMaxLength = 100,
            });
        Assert.Multiple(() =>
        {
            Assert.That(fallback.HasSecondaryStructure, Is.True);
            Assert.That(fallback.Warnings, Is.EqualTo(new[] { "Potential secondary structure formation" }));
            Assert.That(fallback.SelfAny, Is.EqualTo(16.0));
            Assert.That(fallback.Issues, Is.EqualTo(new[] { "Self-complementarity: Primer3 self_any 16.00 exceeds 12.00" }),
                "the stem flag is not an issue; the sourced self_any limit is");
            Assert.That(fallback.IsValid, Is.False);
            Assert.That(thermo.ThermodynamicScreen, Is.True);
            Assert.That(thermo.HairpinTm, Is.EqualTo(77.95325865166825).Within(1e-6));
            Assert.That(thermo.HasSecondaryStructure, Is.True);
            Assert.That(thermo.Issues, Has.Some.EqualTo("Potential secondary structure formation: ntthal hairpin Tm 78.0°C exceeds 47°C"));
            Assert.That(thermo.Warnings, Is.Empty);
            Assert.That(thermo.IsValid, Is.False);
        });
    }

    [TestCase("TTT")]
    [TestCase("TTTT")]
    [TestCase("TTTTTTTT")]
    public void ValidateProbe_LongProbeFallback_AnyLoopOfAtLeastThreeIsFlagged(string loop)
    {
        string probe = StemLoop74.Replace("CAGTTTTCTG", "CAG" + loop + "CTG");
        var v = ThermodynamicProbeDesigner.ValidateProbe(probe, Enumerable.Empty<string>());
        Assert.Multiple(() =>
        {
            Assert.That(v.ThermodynamicScreen, Is.False);
            Assert.That(v.HasSecondaryStructure, Is.True);
            Assert.That(v.HasSecondaryStructure, Is.EqualTo(PrimerDesigner.HasHairpinPotential(probe)));
        });
    }

    [Test]
    public void ValidateProbe_LongProbeFallback_NoMismatchTolerance()
    {
        // The canonical screen requires an exactly complementary ≥ 4-bp stem (F49 library convention); the former
        // unsourced 80 % match tolerance is dropped.
        var v = ThermodynamicProbeDesigner.ValidateProbe(MismatchedStem73, Enumerable.Empty<string>());
        Assert.Multiple(() =>
        {
            Assert.That(v.ThermodynamicScreen, Is.False);
            Assert.That(v.HasSecondaryStructure, Is.False);
            Assert.That(v.Issues, Has.None.Contain("secondary structure"));
        });
    }

    [Test]
    public void FallbackHairpinScreen_EqualsHasHairpinPotential_ForLongNonAcgtAndHeuristicProbes()
    {
        var rng = new Random(20261008);
        var probes = new List<string> { StemLoop74, MismatchedStem73, "GCGCTTTGCGCAAAAAAAAANAAAAAAAA", "ACGTNNNNACGTAAAAAAA" };
        for (int i = 0; i < 20; i++)
            probes.Add(new string(Enumerable.Range(0, 61 + rng.Next(140)).Select(_ => "ACGT"[rng.Next(4)]).ToArray()));
        for (int i = 0; i < 10; i++)
            probes.Add(new string(Enumerable.Range(0, 61 + rng.Next(140)).Select(_ => "AC"[rng.Next(2)]).ToArray()));
        var heuristic = new ThermodynamicProbeDesigner.ProbeParameters(20, 300, -1000, 1000, 0, 1, 300, true, 0.3)
            { StructureScreen = ProbeDesigner.ProbeStructureScreen.Heuristic };
        Assert.Multiple(() =>
        {
            foreach (string p in probes)
            {
                bool expected = PrimerDesigner.HasHairpinPotential(p);
                Assert.That(ThermodynamicProbeDesigner.ValidateProbe(p, Enumerable.Empty<string>()).HasSecondaryStructure,
                    Is.EqualTo(expected), p);
                Assert.That(ThermodynamicProbeDesigner.ValidateProbe(p, Enumerable.Empty<string>(), conditions: heuristic).HasSecondaryStructure,
                    Is.EqualTo(expected), "Heuristic " + p);
                var designed = ThermodynamicProbeDesigner.DesignProbes(p, heuristic with { MinLength = p.Length, MaxLength = p.Length }, 1).ToList();
                if (designed.Count == 1)
                    Assert.That(designed[0].Warnings.Contains("Potential secondary structure"), Is.EqualTo(expected), "DesignProbes " + p);
            }
        });
    }


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

    // KA20 — end-to-end through the canonical affine aligner: probe vs a 279-nt subject carrying the probe with one
    // mismatch and a 1-nt deletion. blastn -task blastn (2/−3, 5/2) bl2seq: score 66, 60.8 bits, E 4.33e-15,
    // effective search space 8672; Biopython PairwiseAligner local (2, −3, open −7, extend −2) score 66.0.
    private const string EValueProbe = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACT";
    private const string EValueSubject =
        "CTTGTCTCCAAGTACCCATTTAGTAGACAAATCGTTCCATCACCAATTCGCTGGTTGTTGAACTATACGACCGGGGCACACTGCACTCAGTTCCCATTTAGAGGATCC" +
        "TAGCCTAGCTACGCTAAAGACAATGACATAACATACAGTCAGCACGAAACTGCGTTTGCGCATCAGGCTGTCCCATACATCAAGCGGTTCCCCTCAAATTATCCGGAC" +
        "TCGGTAAGGGCAGCGAGTAAATATTTTACAATACGTTTCTTGTCAATCTGCTGCTTTGTACGC";



    // Kane et al. (2000) NAR 28:4552: non-targets > 75 % similar over the probe, or sharing a stretch of
    // > 15 contiguous identical bases, may cross-hybridize. Reference values: Biopython 1.88
    // PairwiseAligner(mode='local', match 2, mismatch −3, open_gap_score −7, extend_gap_score −2) = BLAST+ blastn
    // scoring (existence 5, extension 2) — score and the identical-column count of the (unique) optimal alignment;
    // longest common substring by dynamic programming.
    private const string KaneProbe = "TATGCCTCCGGTACATCAACTACAGTTAGCCTTAAGAGAAAAATCCCAAA";
    // Probe with a substitution at every 5th position (40/50 identical), random flanks.
    private const string KaneNonTargetA = "CCGCACCATGAGACTGTTTCTATGGCTCCTGTACCTCAAGTACATTTAGGCTTACGAGACAAATGCCAACCACATCGGCTTCGCACGTCT";
    // Probe[10..28) (18 nt) embedded in random sequence.
    private const string KaneNonTargetB = "GGTCCCACTGATAACGTGTTACCGGCTCTAGTACATCAACTACAGTTACCTAATGCAAAAAACTGTTAACACTTTAAA";
    // Unrelated random sequence.
    private const string KaneNonTargetC = "AATGTGATAGGATGTTAAAAGCGCCGAGACGGCGGTCTGCGATGTACCCCGCAACTGGTTCTTCCCCAGCCGCGGGGGTA";
    // Reverse complement of the probe embedded in random sequence.
    private const string KaneNonTargetD = "CCCCCGGCATTGTTCTTTGGGATTTTTCTCTTAAGGCTAACTGTAGTTGATGTACCGGAGGCATACGGCGCGGAATGACG";
    // Probe[0..15) (exactly 15 nt) embedded in random sequence — at, not above, the contiguous threshold.
    private const string KaneNonTargetE = "AAATTCCCGAAGAAGCGACTTGGTAGGGGATATGCCTCCGGTACACTAGTTCGTTATCCGTCTCGTATTCTGCTC";



    // Fixtures (random.seed(20261008)); expected values from an independent Python brute-force Hamming scan
    // (identity (L - d) / L, longest run of identical positions; Kane et al. 2000: > 0.75 or > 15 nt).
    private const string KaneSiteProbe12 = "GATCCGACGCTA";
    private const string KaneSiteProbe40 = "TATGCCGTACAGTTTTAAGATAGAGCGAAAGCGCAGACAA";

    // Exact site at 10; probe with positions 0..11 substituted at 60 (identity 0.70, run 28); probe with positions
    // 1, 4, …, 34 substituted at 110 (identity 0.70, longest run 5).
    private const string KaneSiteReference40 =
        "AAAAAAAAAATATGCCGTACAGTTTTAAGATAGAGCGAAAGCGCAGACAAAAAAAAAAAAACATGGTACGCTTTTTAAGATAGAGCGAAAGCGCAGACAA"
        + "AAAAAAAAAATCTGGCGAACCGTATTCAGCTATAGGGACAGGGCCGACAAAAAAAAAAAA";



    // Expected values: independent Python brute force over both strands (oracle.py — every reference window compared
    // position by position with the probe and with its Biopython reverse complement; the union of the hit positions per
    // reference is the site set, a position cross-hybridizing when either orientation meets a Kane criterion:
    // identity (L - d) / L > 0.75 or a run of identical positions > 15). Fixture random.seed(20261008):
    // reference = A + P + B + rc(P with positions 3, 11 substituted) + C + rc(P with 1, 5, 9, 13, 17 substituted) + D.
    private const string BothStrandsProbe = "GATCCGACGCTATATGCCGT";
    private const string BothStrandsReference =
        "ACAGTTTTAAGATAGGATCCGACGCTATATGCCGTAGCGAAAGCGCAGACACGGCATAGAGCGTCGCATCAATAAATAATCCGTAACCGCAGATACCGTAGGAGCGGGAGACCTGGCACA";

    // Random cases (Python random.Random(7): probe 8–30 nt, 1–3 references of 0–120 random nt with 0–3 planted copies of
    // the probe or its reverse complement carrying 0–4 substitutions; radius 0–4) where the two modes differ. Columns:
    // probe, references, radius, (hits, Kane sites) given strand, (hits, Kane sites) both strands — from oracle.py.
    private static IEnumerable<TestCaseData> BothStrandsOracleCases()
    {
        yield return new TestCaseData("CCCCCCAATGCCCCGC", new[] { "TAGGGCGGTGCGCGCCCAATGCGCGGGGCATTGGGGGGGCCGGATTTGGTGGGTAA" }, 4, 0, 0, 1, 1);
        yield return new TestCaseData("ACTAAAGCAAGCTCCCTTGGACTA", new[] { "TTCCGTTCCCTAGCAGTCGGCGCTAGTCCAAGTGCGCTTGCTTTAGTTAACGAGAAG" }, 4, 0, 0, 1, 1);
        yield return new TestCaseData("AGCGGAGACGGTAG", new[] { "GAACGGCTATAATCTACCGTCTCCGCTAAGCCGTCGGTAAGCTTAAACTTCTTCAGGCG" }, 1, 0, 0, 1, 1);
        yield return new TestCaseData("AGGTTCTAAAGGCTATGC", new[] { "GTGAGTAACATTGCATAGCCTTTAGAACCTCGCGCCACAGGTTCAAAAGGCAATGCATGAGCACG" }, 1, 0, 0, 1, 1);
        yield return new TestCaseData("CGATGCAAATTCCTCTGTTTCTAGT", new[] { "ACTAGAAACAGAGGAATTTGCATCGACATGACTAGAAACCGAGGAATTTGCATCGTCTACTTTAACTATTCGT" }, 1, 0, 0, 2, 2);
        yield return new TestCaseData("CTGCCCACCAGTCGCGAGGCAA", new[] { "TCCACTAACAGTACAGGCACGATCTCTATTCATTCACCAACAGCAGTCCCGAAGCCTTGCCTCGCGACTGGTGGGCAG" }, 0, 0, 0, 1, 1);
        yield return new TestCaseData("ACGGCGCTTTTATTTCGGGGTC", new[] { "AATACCCCGAAATCAAAGCGCCGTGGTCGTCCAAGGAGTGCAGCTATATTCATTTGCTTCAAAAAGTAGTCATTCCGGT" }, 3, 0, 0, 1, 1);
        yield return new TestCaseData("TATGAGAAAAGTTG", new[] { "CTTATTAATCCAACTTTTCTCATATCATGTAGCCGGCCCGCAGAAGCAGCCGGTTTTTGTTAGATATTAGAAAAGTAGCG" }, 2, 1, 1, 2, 2);
        yield return new TestCaseData("CTTGGCAT", new[] { "TTCGAGGTTTATTTCTTGGCATGTGAGCAGCATCGATAAGTATGCCAAGTATGGGCCTTGGCATAAAATTACGGGGGTAACGCCACCAGTTC" }, 4, 11, 2, 20, 3);
        yield return new TestCaseData("TGTTGAGC", new[] { "CGCTATGTCTAAACGCCGCGCTTAAGGCACAAGAGTTTCTTTTGAGGAGAAGTTCTATGAGTTTGTCGAGCACGGCACTCGCAAGAGAGACTCG" }, 2, 2, 1, 3, 1);
        yield return new TestCaseData("TCGCGATATAA", new[] { "TTATATTTATATCGGGACGGGAGCTGGCGATAAACTTTC", "TAATTAGTTATATCGGGAATACTCGCGATATAAGCGCGCGCTCCATTTACAGCCCAACGCTACAAGG" }, 4, 2, 1, 6, 3);
    }



    // 74-nt probe with a 10-bp perfect stem GCGCATCCAG·CTGGATGCGC closing the 4-nt loop TTTT;
    // the A/C-only flanks cannot pair. Auditor's repro: the former private scan (loop of exactly 3 nt) missed it.
    private const string StemLoop74 =
        "ACCAACCACACAACCCAAACCCAACGCGCATCCAGTTTTCTGGATGCGCAAAAACCCAACCCACCCCACCACCC";

    // 73-nt: a 5-bp "stem" ACCAC·GTAGT (4 of 5 complementary, central C·A mismatch) closing a 3-nt loop in an A/C-only
    // background, no exact 4-bp stem anywhere: flagged only by the former unsourced "≥ 80 % matched" tolerance.
    private const string MismatchedStem73 =
        "ACCAACCACACAACCCAAACCCAACACCAAACCACAAAGTAGTACCAACCACACAACCCAAACCCAACACCAA";

    // F55's random 61-mer (PrimerDesigner_NtthalMaxAlign_Tests.Random61): thal.c -DTHAL_MAX_ALIGN=10000 at
    // 50 mM / 0 / 0 / 50 nM: hairpin 33.980529935122263 °C, ANY −0.968 °C, END1 −79.6 °C (all ≤ 47 °C: no sourced
    // structure). It contains exact 4-bp stems closing ≥ 3-nt loops (independent Python scan: True), as 98 % of random
    // 61-mers do (F61 coverage).
    private const string Random61 = "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGT";
}
