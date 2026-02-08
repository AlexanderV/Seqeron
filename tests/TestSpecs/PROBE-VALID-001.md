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
| **Last Updated** | 2026-01-23 |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Hybridization probe | Academic | Probe hybridization, stringency, cross-hybridization |
| Wikipedia: DNA microarray | Academic | Probe specificity, cross-hybridization issues |
| Wikipedia: Off-target genome editing | Academic | Mismatch tolerance (1-5 bp), off-target detection methods |
| Wikipedia: BLAST | Academic | Approximate matching algorithms, sequence alignment |
| Amann & Ludwig (2000) | Research | rRNA probe specificity, limitations |

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

## Audit Results

### Existing Test Coverage (ProbeDesignerTests.cs)

| Test | Status | Notes |
|------|--------|-------|
| ValidateProbe_UniqueProbe_HighSpecificity | Covered | M3 partially |
| ValidateProbe_MultipleHits_LowSpecificity | Covered | M4 partially |
| ValidateProbe_HighSelfComplementarity_ReportsIssue | Covered | M8 partially |

### Missing Tests (All Closed)

| ID | Test Case | Status |
|----|-----------|--------|
| M1 | Empty probe handling | ✅ Covered |
| M2 | Empty references handling | ✅ Covered |
| M5 | Specificity range invariant | ✅ Covered |
| M6 | Self-complementarity range invariant | ✅ Covered |
| M7 | OffTargetHits non-negativity | ✅ Covered |
| M9-M11 | CheckSpecificity with suffix tree | ✅ Covered |
| M12 | Case-insensitive handling | ✅ Covered |
| S1-S4 | Should tests | ✅ Covered |

### Weak Tests

| Test | Issue | Fix |
|------|-------|-----|
| ValidateProbe_UniqueProbe_HighSpecificity | Doesn't verify all invariants | Add Assert.Multiple with all fields |
| ValidateProbe_MultipleHits_LowSpecificity | Incomplete | Verify exact formula |

---

## Consolidation Plan

1. **Canonical File**: Create `ProbeDesigner_ProbeValidation_Tests.cs` (new file)
2. **Smoke Tests**: Keep 1-2 validation tests in `ProbeDesignerTests.cs` as smoke tests
3. ~~**Add Missing Tests**: M1-M12, S1-S4~~ ✅ Done
4. **Strengthen Tests**: Use Assert.Multiple for invariant grouping
5. **Remove**: Move validation tests from ProbeDesignerTests.cs to canonical file
6. **Naming Convention**: `Method_Scenario_ExpectedResult`

---

## Open Questions

None - behavior is well-documented in implementation.

---

## Assumptions

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | maxMismatches=3 is appropriate default | Aligns with CRISPR off-target tolerance (1-5 mismatches tolerated) |
| A2 | Self-complementarity threshold of 0.3 for warnings | Common practice for probe design |
