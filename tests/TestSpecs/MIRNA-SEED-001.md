# MIRNA-SEED-001: Seed Sequence Analysis — TestSpec

## Test Unit Overview
| Field | Value |
|-------|-------|
| **ID** | MIRNA-SEED-001 |
| **Area** | MiRNA |
| **Class** | MiRnaAnalyzer |
| **Methods** | GetSeedSequence, CreateMiRna, CompareSeedRegions |
| **Canonical Test File** | MiRnaAnalyzer_SeedAnalysis_Tests.cs |

## Consolidation Plan
- **Existing tests**: `MiRnaAnalyzerTests.cs` contains seed-related tests mixed with target/precursor/alignment tests
- **Action**: Extract seed-specific tests into dedicated `MiRnaAnalyzer_SeedAnalysis_Tests.cs`
- **Removed from existing file**: Seed Sequence Tests region, CreateMiRna tests (moved to canonical file)
- **Kept in existing file**: Target site, alignment, pre-miRNA, context, family, utility tests (belong to other test units)

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
**Evidence**: TargetScan, implementation convention
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
**Evidence**: miRNA family grouping allows 0-1 mismatches
**Input**: Two miRNAs with 1 seed difference
**Expected**: Matches=6, Mismatches=1, IsSameFamily=false

### M-012: CompareSeedRegions — Empty seed handling
**Evidence**: Defensive programming (ASSUMPTION)
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

## Audit Notes (Existing Tests)

| Existing Test | Classification | Action |
|---------------|---------------|--------|
| GetSeedSequence_ValidMiRNA_ReturnsSeed | Weak (uses non-standard miRNA name) | Replace with miRBase reference data |
| GetSeedSequence_ShortSequence_ReturnsEmpty | Covered | Migrate to canonical file |
| GetSeedSequence_EmptySequence_ReturnsEmpty | Covered | Migrate to canonical file |
| CreateMiRna_ValidSequence_CreatesMiRna | Weak (doesn't verify seed against miRBase) | Replace with reference-validated version |
| CreateMiRna_DNASequence_ConvertsToRNA | Covered | Migrate to canonical file |
| ReverseComplement tests | Not in scope for MIRNA-SEED-001 | Keep in existing file |
| CanPair/IsWobblePair tests | Not in scope | Keep in existing file |
| Target/Alignment/PreMiRNA/Context/Family/Utility tests | Not in scope | Keep in existing file |

## Open Questions
None — all behavior well-defined by authoritative sources.
