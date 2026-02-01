# Test Specification: PROBE-DESIGN-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | PROBE-DESIGN-001 |
| **Area** | MolTools |
| **Title** | Hybridization Probe Design |
| **Canonical Class** | `ProbeDesigner` |
| **Canonical Methods** | `DesignProbes`, `DesignTilingProbes`, `ScoreProbe` (via EvaluateProbe) |
| **Complexity** | O(n²) |
| **Status** | ☑ Complete |
| **Last Updated** | 2026-01-23 |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Nucleic acid thermodynamics | Academic | Tm calculation, nearest-neighbor method, GC content effects |
| Wikipedia: Hybridization probe | Academic | Probe design principles, applications (15-10000 nt) |
| Wikipedia: FISH | Academic | Probe size 10-25 nt for specificity, 200+ bp for FISH |
| Wikipedia: DNA microarray | Academic | Microarray probe design, 25-60 mer oligos |
| SantaLucia (1998) | Research | Unified nearest-neighbor thermodynamics for Tm |
| Breslauer et al. (1986) | Research | Predicting DNA duplex stability |

---

## Invariants

1. **Score Range**: 0.0 ≤ score ≤ 1.0 (Source: Implementation)
2. **GC Range**: 0.0 ≤ GC content ≤ 1.0 (Source: Mathematical definition)
3. **Tm Positivity**: Tm > 0 for valid probes (Source: Physical law)
4. **Coordinate Validity**: 0 ≤ Start < End < sequence.Length (Source: Implementation)
5. **Probe Substring**: probe.Sequence == input.Substring(probe.Start, probe.End - probe.Start + 1) (Source: Implementation)

---

## Test Cases

### Must (Required - Evidence-Based)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | Empty sequence returns empty result | Boundary condition | Implementation spec |
| M2 | Null sequence returns empty result | Boundary condition | Implementation spec |
| M3 | Sequence shorter than MinLength returns empty | Length constraint | Implementation spec |
| M4 | Valid sequence produces probes with score in [0,1] | Invariant #1 | Implementation |
| M5 | All probes have GC content in [0,1] | Invariant #2 | Mathematical |
| M6 | All probes have Tm > 0 | Invariant #3 | Physical law |
| M7 | Probe coordinates are valid (Start ≥ 0, End < seq.Length) | Invariant #4 | Implementation |
| M8 | Probe sequence matches substring at coordinates | Invariant #5 | Implementation |
| M9 | Tiling probes cover expected positions | Coverage guarantee | Algorithm spec |
| M10 | Tiling probes all have Type = Tiling | Type consistency | Implementation |
| M11 | Microarray defaults: length 50-70 bp | Application param | Wikipedia (DNA microarray) |
| M12 | FISH defaults: length 200-500 bp | Application param | Wikipedia (FISH) |
| M13 | High GC content (100%) results in GcContent ≈ 1.0 | Edge case | Mathematical |
| M14 | Low GC content (all A/T) results in low GcContent | Edge case | Mathematical |
| M15 | maxProbes parameter limits returned count | API contract | Implementation |

### Should (Important)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | Homopolymer runs generate warnings | Quality check | General practice |
| S2 | Case-insensitive input handling | Usability | Implementation |
| S3 | DesignAntisenseProbes returns Antisense type | Type correctness | Implementation |
| S4 | MolecularBeacon has stem sequences | Structure check | Implementation |
| S5 | Tiling probes calculate mean Tm correctly | Statistics | Implementation |
| S6 | Probes are sorted by score descending | Ranking | Implementation |

### Could (Optional)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | qPCR defaults produce 20-30 bp probes | Application param | Standard practice |
| C2 | Self-complementarity detection works correctly | Quality metric | Implementation |
| C3 | Secondary structure detection identifies hairpins | Quality metric | Implementation |

---

## Audit Results

### Existing Test Coverage (ProbeDesignerTests.cs)

| Test | Status | Notes |
|------|--------|-------|
| DesignProbes_ValidSequence_ReturnsProbes | Covered | Basic functionality |
| DesignProbes_ShortSequence_ReturnsEmpty | Covered | M3 |
| DesignProbes_EmptySequence_ReturnsEmpty | Covered | M1 |
| DesignProbes_ProbesHaveValidCoordinates | Covered | M7 partially |
| DesignProbes_ProbesHaveCalculatedProperties | Covered | M4, M5, M6 partially |
| DesignProbes_MicroarrayDefaults_ReturnsValidProbes | Covered | M11 |
| DesignProbes_FISHDefaults_AllowsLongerProbes | Weak | Only checks param, not result |
| DesignProbes_qPCRDefaults_ShorterProbes | Weak | Only checks param |
| DesignTilingProbes_CoversSequence | Covered | M9 partially |
| DesignTilingProbes_CalculatesTmStatistics | Covered | S5 |
| DesignTilingProbes_ProbesAreTilingType | Covered | M10 |
| DesignAntisenseProbes_ReturnsAntisenseType | Covered | S3 |
| DesignMolecularBeacon_CreatesBeaconWithStem | Covered | S4 |
| DesignMolecularBeacon_ShortSequence_ReturnsNull | Covered | Boundary |
| ValidateProbe_* tests | Covered | PROBE-VALID-001 scope |
| AnalyzeOligo_CalculatesAllProperties | Covered | Related utility |
| Edge cases (AllGC, AllAT, Homopolymer, CaseInsensitive) | Partially | Need strengthening |

### Missing Tests

| ID | Test Case | Priority |
|----|-----------|----------|
| M2 | Null sequence handling | Must |
| M8 | Probe sequence matches substring | Must |
| M15 | maxProbes limits result count | Must |
| - | Score range invariant assertion | Must |
| - | Assert.Multiple for invariant grouping | Refactor |

### Weak Tests

| Test | Issue | Fix |
|------|-------|-----|
| DesignProbes_FISHDefaults_AllowsLongerProbes | Only checks params | Verify actual probe lengths |
| DesignProbes_qPCRDefaults_ShorterProbes | Only checks params | Verify actual probe lengths |
| Edge case tests | Incomplete assertions | Add invariant assertions |

---

## Consolidation Plan

1. **Canonical File**: `ProbeDesignerTests.cs` → rename to `ProbeDesigner_ProbeDesign_Tests.cs`
2. **Add Missing Tests**: M2, M8, M15, invariant assertions
3. **Strengthen Weak Tests**: FISH/qPCR tests to verify actual probe lengths
4. **Use Assert.Multiple**: Group invariant assertions for clarity
5. **Remove**: None (no duplicates found)
6. **Naming Convention**: `Method_Scenario_ExpectedResult`

---

## Open Questions

None - behavior is well-documented in implementation and sources.

---

## Assumptions

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Wallace rule applies for probes ≤14 bp | Implementation uses ThermoConstants.WallaceMaxLength |
| A2 | Score penalties are implementation-specific | No standardized scoring system exists |
