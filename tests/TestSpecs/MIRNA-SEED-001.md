# MIRNA-SEED-001: Seed Sequence Analysis — TestSpec

## Test Unit Overview
| Field | Value |
|-------|-------|
| **ID** | MIRNA-SEED-001 |
| **Area** | MiRNA |
| **Class** | MiRnaAnalyzer |
| **Methods** | GetSeedSequence, CreateMiRna, CompareSeedRegions |
| **Canonical Test File** | MiRnaAnalyzer_SeedAnalysis_Tests.cs |
| **Status** | ☑ Complete |

## Consolidation Plan
- **Canonical file**: `MiRnaAnalyzer_SeedAnalysis_Tests.cs` — 21 tests (29 with parametric cases)
- **Migration**: Complete — all seed tests from `MiRnaAnalyzerTests.cs` migrated
- **Property tests**: 4 MIRNA-SEED-001 tests removed from `MiRnaProperties.cs` (3 duplicates + 1 weak)
- **Snapshot tests**: `CreateMiRna_Snapshot` in `MiRnaSnapshotTests.cs` (orthogonal layer, kept)
- **Out of scope in MiRnaAnalyzerTests.cs**: ReverseComplement, CanPair/IsWobblePair, Target/Alignment/PreMiRNA/Context/Family/Utility

## Must Tests

### M-001: GetSeedSequence — Known miRNA returns correct seed (miRBase reference)
**Evidence**: miRBase, Bartel (2009): seed = positions 2-8 of mature miRNA
**Input**: hsa-let-7a-5p (UGAGGUAGUAGGUUGUAUAGUU)
**Expected**: GAGGUAG (7 nt, positions 2-8)

### M-002: GetSeedSequence — Multiple known miRNAs validate seed extraction
**Evidence**: miRBase reference sequences
**Input**: hsa-miR-21-5p, hsa-miR-155-5p, hsa-miR-1-3p
**Expected**: Known seed sequences from miRBase

### M-003: GetSeedSequence — Seed is always 7 nucleotides
**Evidence**: Bartel (2009): family-defining region = positions 2-8 (7 nt); TargetScan site types confirm
**Invariant**: GetSeedSequence returns exactly 7 characters for any valid (≥8 nt) input

### M-004: GetSeedSequence — Empty/null input returns empty string
**Evidence**: Defensive programming requirement
**Input**: null, "", short sequences (<8 nt)
**Expected**: ""

### M-005: GetSeedSequence — Short sequence (< 8 nt) returns empty
**Evidence**: Cannot extract seed from sequence shorter than 8 nt
**Input**: "UAGCA" (5 nt)
**Expected**: ""

### M-006: GetSeedSequence — Case normalization to uppercase
**Evidence**: Implementation convention
**Input**: "ugagguaguagguuguauaguu"
**Expected**: "GAGGUAG" (uppercased)

### M-007: CreateMiRna — Factory produces correct record for known miRNA
**Evidence**: miRBase let-7a-5p
**Input**: name="let-7a", sequence="UGAGGUAGUAGGUUGUAUAGUU"
**Expected**: Name="let-7a", Sequence=uppercase RNA, SeedSequence="GAGGUAG", SeedStart=1, SeedEnd=7

### M-008: CreateMiRna — DNA input is converted to RNA (T→U)
**Evidence**: Standard conversion convention
**Input**: DNA sequence with T
**Expected**: Sequence contains U, no T; seed correctly extracted

### M-009: CompareSeedRegions — Identical seeds yield 0 mismatches and same family
**Evidence**: Family definition from TargetScan: same seed = same family
**Input**: Two let-7 family members (same seed GAGGUAG)
**Expected**: Matches=7, Mismatches=0, IsSameFamily=true

### M-010: CompareSeedRegions — Different seeds yield correct mismatch count
**Evidence**: Hamming distance between seed strings
**Input**: let-7a (seed GAGGUAG) vs miR-21 (seed AGCUUAU)
**Expected**: Computed Hamming distance, IsSameFamily=false

### M-011: CompareSeedRegions — Single mismatch correctly counted
**Evidence**: Bartel (2009): family = identical sequence at positions 2-8; any mismatch = different family
**Input**: Two miRNAs with 1 seed difference
**Expected**: Matches=6, Mismatches=1, IsSameFamily=false

### M-012: CompareSeedRegions — Empty seed handling
**Evidence**: Bartel (2009): miRNAs are ~23 nt; sequences <8 nt are not valid miRNAs — defensive guard for invalid input
**Input**: MiRna with empty seed (from short sequence)
**Expected**: Graceful result with 0 matches, IsSameFamily=false

## Should Tests

### S-001: GetSeedSequence — Exactly 8 nt input (boundary)
**Evidence**: Minimum viable miRNA for seed extraction
**Input**: "ABCDEFGH" (8 characters)
**Expected**: "BCDEFGH" (positions 2-8)

### S-002: CreateMiRna — Mixed case DNA input
**Evidence**: Robustness requirement
**Input**: "tAgCaGcAcGuAaAuAuUgGcG"
**Expected**: Correct uppercase RNA sequence, correct seed

### S-003: let-7 family members share identical seed
**Evidence**: miRBase: let-7a/b/c share seed GAGGUAG
**Input**: let-7a, let-7b, let-7c
**Expected**: All produce seed "GAGGUAG"

### S-004: CompareSeedRegions — Completely different seeds
**Evidence**: Mathematical: max mismatches = seed length
**Input**: Two miRNAs with no shared seed nucleotides
**Expected**: Matches=0, Mismatches=7

## Could Tests

### C-001: Seed extraction is a pure function (deterministic)
**Evidence**: Implementation requirement
**Input**: Same input called multiple times
**Expected**: Same output every time

## Coverage Classification

**Total: 21 tests in canonical file (29 including parametric cases)**

### Summary

| Category | Count |
|----------|-------|
| ❌ Missing → Implemented | 0 |
| ⚠ Weak → Removed | 1 (Property: `IsSubstringOfMiRna` — `Does.Contain` can't distinguish extraction position) |
| 🔁 Duplicate → Removed | 3 (Property: `HasLength7` ⊂ M-003, `PreservesNameAndSequence` ⊂ M-007, `IdenticalSeeds_AllMatches` ⊂ M-009) |
| ✅ Covered | 21 canonical tests |

### Classification Detail

#### ⚠ Weak (Removed)

| Test | File | Before | Action |
|------|------|--------|--------|
| `GetSeedSequence_IsSubstringOfMiRna` | MiRnaProperties.cs | `Does.Contain(seed)` — passes even if seed extracted from wrong position (e.g., pos 3-9 for miR-21) | Removed: canonical M-001/M-002 verify exact seed values against miRBase |

#### 🔁 Duplicate (Removed)

| Test | File | Reason |
|------|------|--------|
| `GetSeedSequence_HasLength7` | MiRnaProperties.cs | Subsumed by M-003 (5 parametric cases all assert length=7) |
| `CreateMiRna_PreservesNameAndSequence` | MiRnaProperties.cs | Subsumed by M-007 (checks all 5 fields: Name, Sequence, SeedSequence, SeedStart, SeedEnd) |
| `CompareSeedRegions_IdenticalSeeds_AllMatches` | MiRnaProperties.cs | Subsumed by M-009 (uses real let-7a/b family members; checks Matches=7, Mismatches=0, IsSameFamily) |

### Canonical File (`MiRnaAnalyzer_SeedAnalysis_Tests.cs`) — 21 tests (29 with parametric cases)

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `GetSeedSequence_Let7a_ReturnsGAGGUAG` | M-001 | ✅ |
| 2 | `GetSeedSequence_KnownMiRNAs_ReturnsExpectedSeed` (3 cases) | M-002 | ✅ |
| 3 | `GetSeedSequence_ValidInput_ReturnsExactly7Characters` (5 cases) | M-003 | ✅ |
| 4 | `GetSeedSequence_NullInput_ReturnsEmpty` | M-004 | ✅ |
| 5 | `GetSeedSequence_EmptyString_ReturnsEmpty` | M-004 | ✅ |
| 6 | `GetSeedSequence_ShortSequence_ReturnsEmpty` (4 cases) | M-005 | ✅ |
| 7 | `GetSeedSequence_LowercaseInput_ReturnsUppercaseSeed` | M-006 | ✅ |
| 8 | `GetSeedSequence_MixedCaseInput_ReturnsUppercaseSeed` | M-006 | ✅ |
| 9 | `CreateMiRna_Let7a_ProducesCorrectRecord` | M-007 | ✅ |
| 10 | `CreateMiRna_DnaInput_ConvertsToRna` | M-008 | ✅ |
| 11 | `CompareSeedRegions_IdenticalSeeds_ZeroMismatchesSameFamily` | M-009 | ✅ |
| 12 | `CompareSeedRegions_Let7a_Vs_MiR21_CorrectHammingDistance` | M-010 | ✅ |
| 13 | `CompareSeedRegions_SingleMismatch_CorrectlyReported` | M-011 | ✅ |
| 14 | `CompareSeedRegions_EmptySeed_ReturnsGracefulResult` | M-012 | ✅ |
| 15 | `CompareSeedRegions_BothEmptySeeds_ReturnsGracefulResult` | M-012 | ✅ |
| 16 | `GetSeedSequence_Exactly8Nucleotides_ReturnsSeed` | S-001 | ✅ |
| 17 | `CreateMiRna_MixedCaseDna_ConvertsCorrectly` | S-002 | ✅ |
| 18 | `GetSeedSequence_Let7Family_AllShareSameSeed` | S-003 | ✅ |
| 19 | `CompareSeedRegions_CompletelyDifferentSeeds_MaxMismatches` | S-004 | ✅ |
| 20 | `GetSeedSequence_CalledMultipleTimes_ReturnsSameResult` | C-001 | ✅ |

### Supplementary Tests (out of canonical scope)

| Test | File | Status |
|------|------|--------|
| `CreateMiRna_Snapshot` | MiRnaSnapshotTests.cs | ✅ Orthogonal (snapshot/golden-master) |

### Classification Summary

- ✅ Covered: 21 tests (20 methods, 12 parametric cases)
- ❌ Missing: 0
- ⚠ Weak: 0
- 🔁 Duplicate: 0

### Theory Verification

| Check | Result |
|-------|--------|
| Seed = positions 2–8 (Bartel 2009) | `Substring(1, 7)` = index 1–7 (0-based) = positions 2–8 (1-based) ✅ |
| miRBase reference data | All 6 sequences verified against MIMAT accessions ✅ |
| Hamming distance (M-010) | GAGGUAG vs AGCUUAU = 2 matches, 5 mismatches (hand-verified) ✅ |
| Single mismatch (M-011) | GAGGUAG vs GAGUUAG = 6 matches, 1 mismatch (hand-verified) ✅ |
| Family = identical pos 2–8 | `IsSameFamily = seed1 == seed2` ✅ |
| Wrong impl would fail tests | `Substring(0,7)` → let-7a seed "UGAGGUA" ≠ "GAGGUAG" → M-001 fails ✅ |

## Open Questions
None — all behavior well-defined by authoritative sources.

## Validation Checklist

- [x] All MUST tests have evidence source (miRBase, Bartel 2009, TargetScan)
- [x] Invariants specified and tested (seed length=7, family=identical seed)
- [x] Edge cases documented (empty, null, short, boundary 8 nt)
- [x] Reference data from miRBase (6 miRNAs with MIMAT accessions)
- [x] Hamming distances hand-verified
- [x] No assumptions — all design decisions backed by external sources
- [x] No duplicates — each test serves a distinct purpose
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
- [x] Tests passing (29/29 parametric cases)
