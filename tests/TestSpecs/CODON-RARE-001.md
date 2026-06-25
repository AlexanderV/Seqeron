# Test Specification: CODON-RARE-001

## Test Unit Information
- **ID**: CODON-RARE-001
- **Title**: Rare Codon Detection
- **Area**: Codon Optimization
- **Complexity**: O(n)
- **Status**: ☑ Complete

## Method Under Test

```csharp
public static IEnumerable<(int Position, string Codon, string AminoAcid, double Frequency)> 
    FindRareCodons(string codingSequence, CodonUsageTable table, double threshold = 0.15)
```

## Evidence Summary

Sources: Wikipedia (Codon usage bias), Kazusa Codon Usage Database, GenScript GenRCA, Shu et al. (2006), Sharp & Li (1987)

Key findings:
- Rare codons (AGA, AGG, CGA, CUA) have low tRNA abundance in E. coli
- Default threshold 0.15 is reasonable for most applications
- Position, codon, amino acid, and frequency are all valuable output fields

---

## Test Categories

### MUST Tests (Evidence-backed, Critical)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| M01 | Empty sequence returns empty | `FindRareCodons("", table)` → empty enumerable | Implementation contract |
| M02 | Single rare codon detection | Sequence with one AGA (0.04) detected below threshold 0.10 | Shu et al., Kazusa MG1655 |
| M03 | Multiple rare codons detection | Sequence with AGA, AGG, CGA all detected | Kazusa frequencies |
| M04 | Position is nucleotide index | AGA at codon index 1 reports position 3 | Implementation spec |
| M05 | Threshold boundary - below | Codon at freq 0.07 detected with threshold 0.10 | Math invariant |
| M06 | Threshold boundary - at | Codon at freq 0.10 NOT detected with threshold 0.10 | Math invariant (< not ≤) |
| M07 | No rare codons found | Sequence with only common codons returns empty | Logic correctness |
| M08 | Amino acid translation correct | Reported AA matches genetic code | Sharp & Li |
| M09 | Frequency value matches table | Reported frequency equals table lookup | Implementation spec |
| M10 | DNA T converted to RNA U | Input with T handled correctly | Implementation spec |

### SHOULD Tests (Important behavior)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| S01 | Default threshold is 0.15 | Single-arg overload uses 0.15 | Implementation default |
| S02 | Incomplete codons ignored | 10-nucleotide sequence processes only 3 codons | Implementation spec |
| S03 | Case insensitivity | Lowercase input handled | Robustness |
| S04 | All codons rare scenario | Sequence of all rare codons returns all positions | Stress case |
| S05 | Mixed rare/common codons | CUG (common) + CUA (rare) - only CUA reported | Kazusa frequencies |

### COULD Tests (Edge cases, nice-to-have)

| ID | Test Name | Description | Evidence |
|----|-----------|-------------|----------|
| C01 | Threshold zero | All codons with freq > 0 NOT reported as rare | Math edge |
| C02 | Threshold one | All codons reported as rare | Math edge |
| C03 | Different organism tables | Yeast table has different rare codons than E. coli | Kazusa data |
| C04 | Unknown codon handling | Non-standard codon gets freq 0 | Implementation behavior |

---

## Test Data

### E. coli K12 Rare Codons (freq < 0.10) — Kazusa MG1655 (species=316407)
- AGA: 0.04 (Arginine)
- AGG: 0.02 (Arginine)  
- CGA: 0.06 (Arginine)
- CUA: 0.04 (Leucine)
- UAG: 0.07 (Stop)

### E. coli K12 Common Codons (freq > 0.40)
- CUG: 0.50 (Leucine)
- AUG: 1.00 (Methionine)
- UGG: 1.00 (Tryptophan)

### Test Sequences

| Name | Sequence | Description |
|------|----------|-------------|
| SingleRare | AUGAGA | M + R(rare) |
| MultiRare | AUGAGAAGGCGA | M + R(rare) + R(rare) + R(rare) |
| NoRare | AUGCUGUGG | M + L(common) + W(common) |
| AllRare | AGAAGGCGA | All rare arginine codons |
| MixedLeu | CUGCUA | L(common 0.50) + L(rare 0.04) |
| Borderline | CUGAUA | L(common) + I(0.07 near threshold) |

---

## Invariants

1. All reported positions are multiples of 3
2. All reported positions are in range [0, sequence.Length - 3]
3. All reported frequencies are < threshold
4. All reported codons have length 3
5. Same input → same output (deterministic)
6. Protein sequence is preserved (amino acid translation correct)

---

## Audit Notes

### Existing Tests (CodonOptimizerTests.cs)
- `FindRareCodons_IdentifiesRare` - basic detection ✓
- `FindRareCodons_ReturnsPosition` - position check ✓
- `FindRareCodons_EmptySequence_ReturnsEmpty` - empty input ✓

### Coverage Gaps (All Closed)
- ~~No threshold boundary tests~~ ✅ Covered
- ~~No threshold=0 or threshold=1 edge cases~~ ✅ Covered
- ~~No multi-rare codon sequence test~~ ✅ Covered
- ~~No case insensitivity test~~ ✅ Covered
- ~~No DNA (T) conversion test~~ ✅ Covered
- ~~No frequency value verification~~ ✅ Covered
- ~~No amino acid translation verification~~ ✅ Covered
- ~~No incomplete codon test~~ ✅ Covered

### Consolidation Plan
- ~~Create dedicated test file: `CodonOptimizer_FindRareCodons_Tests.cs`~~ ✅ Done
- ~~Move existing 3 tests to new file with consistent naming~~ ✅ Done
- ~~Add missing MUST tests~~ ✅ Done
- ~~Add selected SHOULD/COULD tests~~ ✅ Done

---

## Open Questions / Decisions

1. **Q**: Should threshold comparison be `<` or `<=`?  
   **A**: Implementation uses `<` (strict less than). Documented.

2. **Q**: What happens with non-standard codons?  
   **A**: `GetValueOrDefault` returns 0, so they're always flagged as rare. Documented.

---

# Addendum (2026-06-24): Rare-Codon Cluster / Run Detection

Closes the documented limitation that per-codon detection does not find consecutive
rare-codon **clusters / runs / ramps**. Two opt-in methods added; per-codon `FindRareCodons`
unchanged. Sources retrieved this session (see `docs/Evidence/CODON-RARE-001-Evidence.md`
Addendum): Clarke & Clark (2008) %MinMax; Chartier et al. (2012) Sherlocc RCC.

## Methods Under Test (clusters)

```csharp
public static IReadOnlyList<MinMaxWindow> CalculateMinMaxProfile(
    string codingSequence, CodonUsageTable table, int windowSize = 18);

public static IReadOnlyList<RareCodonCluster> FindRareCodonClusters(
    string codingSequence, CodonUsageTable table,
    double rareThreshold = 0.15, int windowSize = 7, int minRareCodons = 4);
```

| Method | Type | Source / parameters |
|--------|------|---------------------|
| `CalculateMinMaxProfile` | Canonical | Clarke & Clark (2008), PLoS ONE 3:e3412; window 18 codons |
| `FindRareCodonClusters` | Canonical | Chartier et al. (2012), Bioinformatics 28(11):1438; 7-codon window, ≥4 rare |

## Evidence-derived expected values (E. coli K12 Arg family; Xavg=0.16667, Xmax=0.40, Xmin=0.02)

| Case | Input | Method / params | Expected |
|------|-------|-----------------|----------|
| MM1 | AGA×3 | %MinMax, w=3 | −86.363636…% |
| MM2 | CGC×3 | %MinMax, w=3 | +100% |
| MM3 | CUG·AGA | %MinMax, w=2 | +36.470588…% |
| C1 | AGA×7 | RCC default | 1 cluster, codons 0–6, 7 rare |
| C2 | AGA×3+CGC×4 | RCC default | no cluster (3 < 4) |
| C3 | AGA×4+CGC×3 | RCC default | 1 cluster, 4 rare (boundary ≥4) |
| C8 | AGA×7 (yeast) | RCC default | no cluster (AGA 0.48 common) |

## Cluster Test Cases (file `CodonOptimizer_RareCodonClusters_Tests.cs`)

%MinMax: MM1 all-rare %Min (exact), MM2 all-common +100, MM3 mixed %Max (exact), MM4 sliding
windows count+order, MM5 short-sequence empty, MM6 empty/null, MM7 single-codon-AA 0 (no NaN),
MM8 bound `[-100,100]`, MM9 DNA→RNA, MM10 windowSize<1 throws.

RCC: C1 7-rare cluster, C2 below-threshold (3<4) no cluster, C3 exactly-4 boundary, C4 isolated
rare codons → per-codon flags 3 but zero clusters (the closed gap), C5 long run merges to one
cluster, C6 two separated runs → two clusters, C7 default-threshold = 0.15, C8 organism-specific,
C9 short-sequence empty, C10 empty/null, C11 DNA→RNA, C12 tunable window/threshold short run,
C13 invalid params throw, C14 determinism.

## Invariants (clusters)

- INV-05: every %MinMax value ∈ `[-100, 100]`.
- INV-06: `CalculateMinMaxProfile` yields `codonCount − w + 1` windows (or none).
- INV-07: single-codon AA contributes 0 to %MinMax (no NaN).
- INV-08: every cluster contains ≥ `minRareCodons` rare codons.

## Phase 7 Work Queue (clusters) — all ✅ Done

24 cluster cases (MM1–MM10, C1–C14) implemented and executed: **Failed: 0, Passed: 24**.
Remaining: 0.

## Status note

Per the validation protocol, the ROOT `ALGORITHMS_CHECKLIST_V2.md` Status for CODON-RARE-001
is reset to ☐ pending independent re-validation of this new capability.
