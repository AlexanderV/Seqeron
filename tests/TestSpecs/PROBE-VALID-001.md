# Test Specification: PROBE-VALID-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | PROBE-VALID-001 |
| **Area** | MolTools |
| **Title** | Probe Validation |
| **Canonical Class** | `ProbeDesigner` |
| **Canonical Methods** | `ValidateProbe`, `CheckSpecificity` |
| **Complexity** | O(n × g) |
| **Status** | ☑ Complete |
| **Last Updated** | 2026-03-04 |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Hybridization probe | Academic | Probe hybridization, stringency-dependent cross-hybridization |
| Wikipedia: DNA microarray | Academic | Probe specificity, cross-hybridization in high-density arrays |
| Wikipedia: Off-target genome editing | Academic | CRISPR/Cas9 tolerates 3-5 bp mismatches per 20nt guide |
| Wikipedia: Nucleic acid thermodynamics | Academic | Mismatch destabilization, complementarity metrics |
| Wikipedia: Off-target activity | Academic | Off-target detection methods, mismatch tolerance mechanisms |

---

## Invariants

1. **Specificity Range**: 0.0 ≤ specificityScore ≤ 1.0 (Source: Implementation)
2. **Self-Complementarity Range**: 0.0 ≤ selfComplementarity ≤ 1.0 (Source: Mathematical)
3. **Off-Target Non-Negative**: offTargetHits ≥ 0 (Source: Implementation)
4. **Unique Match Specificity**: offTargetHits == 1 → specificityScore == 1.0 (Source: Implementation)
5. **No Match Zero Specificity**: offTargetHits == 0 → specificityScore == 0.0 (Source: Implementation)
6. **Multiple Hits**: offTargetHits > 1 → specificityScore == 1.0 / offTargetHits (Source: Implementation)

---

## Test Cases

### Must (Required - Evidence-Based)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | Empty probe returns validation with zero specificity | Boundary condition | Implementation spec |
| M2 | Null/empty references returns validation without off-target hits | Boundary condition | Implementation spec |
| M3 | Unique probe (1 hit) has specificity score 1.0 | Invariant #4 | Implementation |
| M4 | Multiple hits reduce specificity to 1.0/hitCount | Invariant #6 | Implementation |
| M5 | Specificity score always in [0.0, 1.0] | Invariant #1 | Implementation |
| M6 | Self-complementarity in [0.0, 1.0] | Invariant #2 | Mathematical |
| M7 | OffTargetHits is non-negative | Invariant #3 | Implementation |
| M8 | High self-complementarity (>30%) reported in issues | Quality criterion | Implementation spec |
| M9 | CheckSpecificity with suffix tree returns correct score | API contract | Implementation |
| M10 | CheckSpecificity unique match returns 1.0 | Invariant #4 | Implementation |
| M11 | CheckSpecificity multiple matches returns 1.0/count | Invariant #6 | Implementation |
| M12 | Case-insensitive probe handling | Usability | Implementation |

### Should (Important)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | Secondary structure potential detected for hairpin sequences | Quality metric | Implementation |
| S2 | Issues list populated for problematic probes | User feedback | Implementation |
| S3 | IsValid false when multiple issues exist | Validation logic | Implementation |
| S4 | Approximate matching with maxMismatches works correctly | Off-target detection | Wikipedia (Off-target) |

### Could (Optional)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | Long reference sequences handled efficiently | Performance | Implementation |
| C2 | Multiple references all searched | Completeness | Implementation |

---

## Coverage Classification

| ID | Test Method | Classification | Notes |
|----|------------|----------------|-------|
| M1 | `ValidateProbe_EmptyProbe_ReturnsValidationResult` | ✅ Covered | Exact: specificity=0.0, offTargetHits=0, IsValid=false |
| M2 | `ValidateProbe_EmptyReferences_ReturnsValidationWithNoOffTargetHits` | ✅ Covered | Exact: specificity=0.0, offTargetHits=0 (Inv. #5) |
| — | `ValidateProbe_NullReferences_ThrowsArgumentNullException` | ✅ Covered | ArgumentNullException |
| M3 | `ValidateProbe_UniqueProbe_HasSpecificityScoreOne` | ✅ Covered | Exact: offTargetHits=1, specificity=1.0 |
| M4 | `ValidateProbe_MultipleHits_ReducesSpecificityByHitCount` | ✅ Covered | Exact: offTargetHits=25, specificity=0.04 |
| M5 | ~~`ValidateProbe_AnyInput_SpecificityScoreInValidRange`~~ | 🔁 Deleted | Subsumed by AllInvariants + every exact-value test |
| M6 | ~~`ValidateProbe_AnyInput_SelfComplementarityInValidRange`~~ | 🔁 Deleted | Subsumed by AllInvariants + M8 |
| M7 | ~~`ValidateProbe_AnyInput_OffTargetHitsNonNegative`~~ | 🔁 Deleted | Subsumed by AllInvariants + all specific tests |
| M8 | `ValidateProbe_HighSelfComplementarity_ReportsInIssues` | ✅ Covered | Exact: selfComp=1.0, issues contain "Self-complementarity" |
| M9 | ~~`CheckSpecificity_ResultInValidRange`~~ | 🔁 Deleted | Subsumed by M10 + M11 + NoMatch exact tests |
| M9 | `CheckSpecificity_NoMatch_ReturnsZero` | ✅ Covered | Exact: specificity=0.0 |
| M10 | `CheckSpecificity_UniqueSequence_ReturnsOne` | ✅ Covered | Exact: specificity=1.0 |
| M11 | `CheckSpecificity_MultipleOccurrences_ReturnsOneOverCount` | ✅ Covered | Exact: 1.0/count |
| M12 | `ValidateProbe_MixedCaseProbe_HandledCaseInsensitively` | ✅ Covered | Upper/lower/mixed → same results |
| S1 | `ValidateProbe_PotentialHairpin_DetectsSecondaryStructure` | ✅ Covered | Exact: HasSecondaryStructure=true, issues contain text |
| S2 | `ValidateProbe_ProblematicProbe_PopulatesIssuesList` | ✅ Covered | Exact: offTargetHits=16, issues contain "16 potential off-target sites" |
| S3 | `ValidateProbe_MultipleProblems_IsValidFalse` | ✅ Covered | Exact: hits=12, selfComp=1.0, issues=2, IsValid=false |
| S4 | `ValidateProbe_ApproximateMatching_FindsNearMatches` | ✅ Covered | Exact: 0 hits (strict) vs 1 hit (approx) |
| C1 | `ValidateProbe_LongReference_FindsProbeCorrectly` | ✅ Covered | 20k-char ref, exact: hits=1, specificity=1.0 |
| C2 | `ValidateProbe_MultipleReferences_AccumulatesHits` | ✅ Covered | 3 refs, exact: hits=3, specificity=1/3 |
| — | `ValidateProbe_AllInvariants_HoldForTypicalProbe` | ✅ Covered | All 6 invariants verified for typical input |
| — | `ValidateProbe_ZeroHits_ReturnsZeroSpecificity` | ✅ Covered | Explicit Invariant #5 |
| — | `DesignProbes_WithGenomeIndex_UsesCheckSpecificity` | ✅ Covered | Integration: suffix tree + DesignProbes |

---

## Evidence-Backed Parameters

All configurable parameters have external evidence justification. No assumptions remain.

| Parameter | Default | Evidence | Source |
|-----------|---------|----------|--------|
| `maxMismatches` | 3 | Lower bound of CRISPR/Cas9 off-target tolerance: 3-5 bp mismatches per 20nt guide (Hsu et al. 2013, Fu et al. 2013). Configurable per caller. | Wikipedia: Off-target genome editing |
| `selfComplementarityThreshold` | 0.3 | For random DNA (uniform base distribution), expected self-complementarity ≈ 0.25. Threshold of 0.3 (20% above baseline) detects statistically elevated palindromic character. Per-application defaults: qPCR=0.25, Microarray=0.3, NorthernBlot=0.35, FISH/SouthernBlot=0.4. | Wikipedia: Nucleic acid thermodynamics; `ProbeParameters.Defaults` |


