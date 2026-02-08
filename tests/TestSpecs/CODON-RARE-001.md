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
| M02 | Single rare codon detection | Sequence with one AGA (0.07) detected below threshold 0.10 | Shu et al., Kazusa |
| M03 | Multiple rare codons detection | Sequence with AGA, AGG, CGA all detected | Kazusa frequencies |
| M04 | Position is nucleotide index | AGA at codon index 1 reports position 3 | Implementation spec |
| M05 | Threshold boundary - below | Codon at freq 0.09 detected with threshold 0.10 | Math invariant |
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

### E. coli K12 Rare Codons (freq < 0.10)
- AGA: 0.07 (Arginine)
- AGG: 0.04 (Arginine)  
- CGA: 0.07 (Arginine)
- CUA: 0.04 (Leucine)
- UAG: 0.09 (Stop)

### E. coli K12 Common Codons (freq > 0.40)
- CUG: 0.47 (Leucine)
- AUG: 1.00 (Methionine)
- UGG: 1.00 (Tryptophan)

### Test Sequences

| Name | Sequence | Description |
|------|----------|-------------|
| SingleRare | AUGAGA | M + R(rare) |
| MultiRare | AUGAGAAGGCGA | M + R(rare) + R(rare) + R(rare) |
| NoRare | AUGCUGUGG | M + L(common) + W(common) |
| AllRare | AGAAGGCGA | All rare arginine codons |
| MixedLeu | CUGCUA | L(common 0.47) + L(rare 0.04) |
| Borderline | CUGAUA | L(common) + I(0.11 near threshold) |

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
