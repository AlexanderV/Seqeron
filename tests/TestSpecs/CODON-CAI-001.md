# Test Specification: CODON-CAI-001

## Test Unit Information
- **ID:** CODON-CAI-001
- **Title:** Codon Adaptation Index (CAI) Calculation
- **Canonical Method:** `CodonOptimizer.CalculateCAI(string, CodonUsageTable, bool excludeSingleCodonAminoAcids = true)` → delegates to the canonical core `CodonUsageAnalyzer.CalculateCai(string, IReadOnlyDictionary<string,double>, GeneticCode)` (review 2026-09: one CAI implementation)
- **Area:** Codon Optimization
- **Complexity:** O(n)
- **Status:** ☑ Re-validated 2026-09 (review campaign B02, findings F12–F14)

## Method Under Test

```csharp
public static double CalculateCAI(string codingSequence, CodonUsageTable table, bool excludeSingleCodonAminoAcids = true)
// canonical core
public static double CodonUsageAnalyzer.CalculateCai(string sequence, IReadOnlyDictionary<string, double> referenceRscu, GeneticCode code)
```

## Algorithm Summary

CAI = geometric mean of relative adaptiveness values:
- `w_i = f_i / max(f_j)` for synonymous codons
- `CAI = exp((1/L) × Σ ln(w_i))`, L = scored codons
- Not scored: stop codons and single-codon families (Met/Trp in table 1; genetic-code dependent) — Sharp & Li 1987 (quoted by Xia 2007), CodonW `cai_out`, seqinr `cai`, Biopython `CodonAdaptationIndex`
- `w < 0.0001` → `0.01` (CodonW `cai_out`; seqinr `zero.to`; Bulmer 1988)
- Triplets with non-nucleotide symbols skipped without frame shift; DNA/RNA, any case
- Opt-in `excludeSingleCodonAminoAcids:false`: Met/Trp scored with w = 1 (EMBOSS `ajCodCalcCaiSeq`-style)

## Test Categories

### MUST Tests (Required for Completion)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| M1 | CalculateCAI_EmptySequence_ReturnsZero | Empty input returns 0 | Edge case convention |
| M2 | CalculateCAI_SingleMetCodon_ExcludedByDefault_ReturnsZero | AUG not scored → 0; opt-in inclusion → 1.0 | Sharp & Li / CodonW (`ATGTGG` → 0.000) |
| M3 | CalculateCAI_SingleTrpCodon_ExcludedByDefault_ReturnsZero | UGG not scored → 0; opt-in inclusion → 1.0 | Sharp & Li / CodonW |
| M4 | CalculateCAI_AllOptimalCodons_ReturnsOne | All optimal codons → CAI=1.0 | Sharp & Li (1987) |
| M5 | CalculateCAI_RareCodons_ReturnsLow | Rare codons → CAI < 0.5 | Sharp & Li (1987) |
| M6 | CalculateCAI_RangeIsZeroToOne | CAI always in [0, 1] | CAI definition |
| M7 | CalculateCAI_DifferentOrganisms_DifferentResults | Same sequence, different CAI per organism | Organism-specific bias |
| M8 | CalculateCAI_DnaInput_HandledCorrectly | T→U conversion works | Implementation requirement |
| M9 | CalculateCAI_LowercaseInput_Handled | Case-insensitive processing | Robustness |
| M10 | CalculateCAI_ExcludesStopCodons | Stop codons not counted in calculation | Standard practice |
| M11 | CalculateCAI_GeometricMeanProperty | Single rare codon significantly lowers CAI | Mathematical property |
| M12 | CalculateCAI_HandCalculatedValue_Matches | Verify against hand-calculated example | Validation |
| M13 | CalculateCAI_DefaultMode_ExcludesSingleCodonAminoAcids | Default == explicit-true: AUGUGG→0; explicit-false → 1.0 | Sharp & Li / Xia 2007 / CodonW |
| M14 | CalculateCAI_ExcludeMode_AllSingleCodonAA_ReturnsZero | excludeSingleCodonAminoAcids:true; AUGUGG→0 (no scored codons) | Sharp & Li (1987) / Xia (2007) |
| M15 | CalculateCAI_ExcludeMode_DropsMetFromGeometricMean | AUGCUACUA: incl=0.18566355334451112, excl=0.08 exact | Sharp & Li (1987) / Xia (2007) |
| M16 | CalculateCAI_ExcludeMode_DropsBothMetAndTrp_ScoresOnlyRemainder | AUGUGGCUA: incl=0.43088693800637673, excl=0.08 exact | Sharp & Li (1987) / Xia (2007) |
| M17 | CalculateCAI_ExcludeMode_NoSingleCodonAA_UnchangedFromDefault | CUGCUA: incl==excl=0.28284271247461906 | Exclusion only affects Met/Trp |

| M18 | CalculateCAI_AbsentCodonWithPresentSynonym_UsesCodonWZeroSubstitute | table {CUG:1}; CUACUG → 0.1; CUA → 0.01 | CodonW 1.4.4 binary (0.100 / 0.010) |
| M19 | CalculateCAI_ZeroSubstituteThreshold_IsCodonW0_0001 | w = 0.0001 kept; w = 0.00009 → 0.01 | CodonW `if (w < 0.0001) w = 0.01` |
| M20 | CalculateCAI_DelegatesToCanonicalCodonUsageAnalyzerCore | CodonOptimizer == CodonUsageAnalyzer on same values | No duplication |
| M21 | CalculateCAI_AmbiguousTriplet_SkippedWithoutFrameShift | CUGNNNCUA, CUGCURCUA → √0.08 | CodonW `ident_codon` |
| M22 | CalculateCai_SharpLiEColiIndex_MatchesBiopythonAndCodonW (CodonUsageAnalyzer) | 7 genes vs `EColiOptimalCodons` | Biopython 1.88 (full precision), CodonW 1.4.4 (3 dp) |
| M23 | CalculateCai_GeneticCode_MatchesCodonW | tables 1/2: ATAATG, TGATGG, AGAAGGCTG | CodonW `-code` |
| M24 | CalculateCai_ZeroWeightCodon_ScoredAsCodonW001 | E. coli w with CTA=0: CTGCTA 0.1 | CodonW `-cai_file` |

### SHOULD Tests (Recommended)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| S1 | CalculateCAI_LongSequence_Completes | Performance for long sequences | Usability |
| S2 | CalculateCAI_MixedOptimalRare_IntermediateValue | Mix yields intermediate CAI | Mathematical property |
| S3 | CalculateCAI_OnlyStopCodons_ReturnsZero | Sequence of only stops → 0 | Edge case |

### COULD Tests (Optional)

| # | Test Name | Description | Evidence |
|---|-----------|-------------|----------|
| C1 | CalculateCAI_AllThreeOrganismTables_Valid | E. coli, Yeast, Human all work | API coverage |

## Invariants to Verify

1. **Range Invariant:** `0 ≤ CAI ≤ 1`
2. **Single-Codon AA:** by default (`excludeSingleCodonAminoAcids:true`) Met/Trp codons are excluded from the geometric mean entirely (Sharp & Li 1987; Xia 2007; CodonW); with `false` they contribute w=1.0
3. **Monotonicity:** Replacing rare codons with optimal ones increases CAI
4. **Idempotence:** Calculating twice gives same result

## Edge Cases

| Case | Input | Expected Output |
|------|-------|-----------------|
| Empty | `""` | 0 |
| Single Met | `"AUG"` | 0 (not scored); 1.0 with `excludeSingleCodonAminoAcids:false` |
| Single Trp | `"UGG"` | 0 (not scored); 1.0 with `excludeSingleCodonAminoAcids:false` |
| DNA format | `"ATGCTG"` | Same as `"AUGCUG"` |
| Lowercase | `"augcug"` | Same as `"AUGCUG"` |
| Incomplete codon | `"CUAC"` | Based on `"CUA"` only (0.08) |
| Absent codon (w=0) | table {CUG:1}, `"CUACUG"` | 0.1 (w_CUA = 0.01) |
| Ambiguous triplet | `"CUGNNNCUA"` | √0.08 (NNN skipped, frame kept) |
| Only stop codons | `"UAAUAGUGA"` | 0 (no codons to evaluate) |

## Test Data

### Hand-Calculated Reference Values (E. coli K12, Kazusa MG1655)

**Test 1: All Optimal Codons**
```
Sequence: CUGCCGACC (Leu-Pro-Thr)
Codons:   CUG(0.50), CCG(0.53), ACC(0.44)
w values: 0.50/0.50=1.0, 0.53/0.53=1.0, 0.44/0.44=1.0
CAI = (1.0 × 1.0 × 1.0)^(1/3) = 1.0
```

**Test 2: Mixed Codons**
```
Sequence: AUGCUGACC (Met-Leu-Thr)
AUG: not scored (single-codon family)
CUG: 0.50/0.50 = 1.0
ACC: 0.44/0.44 = 1.0
CAI = (1.0 × 1.0)^(1/2) = 1.0
```

**Test 3: Rare Codons**
```
Sequence: CUAACU (Leu-Thr)
CUA: 0.04/0.50 = 0.08
ACU: 0.16/0.44 = 0.36364
CAI = (0.08 × 0.36364)^(1/2) = 0.17056
```

### Exclusion-Mode Reference Values (`excludeSingleCodonAminoAcids: true` = default, E. coli K12; "incl" = opt-in `false`)

CUA(Leu) w = 0.04/0.50 = 0.08; AUG(Met)/UGG(Trp) excluded.

```
AUGUGG    : incl=1.0                  ; excl=0 (L=0, all excluded)
AUGCUACUA : incl=0.18566355334451112  ; excl=exp((ln0.08+ln0.08)/2)=0.08
AUGUGGCUA : incl=0.43088693800637673  ; excl=exp(ln0.08/1)=0.08
CUGCUA    : incl=0.28284271247461906  ; excl=same (no Met/Trp)
```

### Codon Frequencies (E. coli K12, Kazusa MG1655)

| AA | Codon | Freq | Max | w |
|----|-------|------|-----|---|
| Leu | CUG | 0.50 | 0.50 | 1.00 |
| Leu | CUA | 0.04 | 0.50 | 0.08 |
| Pro | CCG | 0.53 | 0.53 | 1.00 |
| Pro | CCA | 0.19 | 0.53 | 0.358 |
| Thr | ACC | 0.44 | 0.44 | 1.00 |
| Thr | ACU | 0.16 | 0.44 | 0.364 |

## Audit Notes

### Current Test Coverage (30 tests in CodonOptimizer_CAI_Tests.cs)

| Category | Count | Tests |
|----------|-------|-------|
| M1 (Empty) | 2 | EmptySequence_ReturnsZero, NullSequence_ReturnsZero |
| M2/M3 (Single-codon AA) | 3 | SingleMetCodon, SingleTrpCodon, MetAndTrp |
| M4 (All optimal) | 2 | AllOptimalCodonsEColi, OptimalCodonsWithMet |
| M5 (Rare codons) | 2 | RareCodonsEColi, RareArginineCodonsEColi |
| M6 (Range) | 1 | AnyValidSequence_RangeIsZeroToOne |
| M7 (Organism diff) | 2 | SameSequence_DifferentOrganisms, YeastPreferredCodons |
| M8 (DNA input) | 1 | DnaInputWithThymine_ConvertsToUracil |
| M9 (Lowercase) | 1 | LowercaseInput_HandledCorrectly |
| M10 (Stop codons) | 3 | SequenceWithStopCodon_Excludes, OnlyStopCodons, StopCodonInMiddle |
| M11 (Geometric mean) | 2 | SingleRareCodon_SignificantlyLowers, MoreRareCodons_LowerCAI |
| M12 (Hand-calculated) | 1 | HandCalculatedRareCodons_MatchesExpected |
| S1 (Performance) | 1 | LongSequence_CompletesInReasonableTime |
| S2 (Mixed) | 1 | MixedOptimalAndRare_MatchesHandCalculated |
| S3 (Only stops) | 1 | OnlyStopCodons_ReturnsZero (shared with M10) |
| C1 (All organisms) | 1 | AllThreeOrganismTables_MatchHandCalculated |
| Edge cases | 2 | IncompleteFinalCodon, TwoIncompleteBases |
| M13–M17 (Single-codon AA exclusion) | 5 | DefaultMode_ExcludesSingleCodonAminoAcids, ExcludeMode_AllSingleCodonAA_ReturnsZero, ExcludeMode_DropsMetFromGeometricMean, ExcludeMode_DropsBothMetAndTrp_ScoresOnlyRemainder, ExcludeMode_NoSingleCodonAA_UnchangedFromDefault |

### Consolidation Status

- **Canonical file:** `CodonOptimizer_CAI_Tests.cs` — all CAI tests consolidated here
- **CodonUsageAnalyzer tests:** Separate class, kept as smoke tests in their own file

### Missing Coverage

None — all test categories fully covered.

## Deviations and Assumptions

**Resolved 2026-09 (review campaign B02):**
- **F12 — default scored Met/Trp.** Sharp & Li (1987) exclude single-codon families (quoted by Xia 2007,
  Evol. Bioinform. 3:53-58, PMC2684136 — formerly mis-cited here as "Jansen et al. 2003"); CodonW, seqinr
  and Biopython implement it. Default is now `excludeSingleCodonAminoAcids: true`; `false` is the
  EMBOSS-style opt-in.
- **F13 — zero-w handling.** The unsourced `1e-6` clamp (CodonOptimizer) and the silent drop of `w = 0`
  codons (CodonUsageAnalyzer, which *raised* CAI) are replaced by CodonW's `w < 0.0001 → 0.01`
  (Bulmer 1988; seqinr `zero.to`).
- **F14 — two CAI implementations.** CodonOptimizer.CalculateCAI now delegates to the canonical
  CodonUsageAnalyzer core (genetic-code aware; new `CalculateCai(..., GeneticCode)` overloads).

**Conventions kept (documented):**
- An amino acid with no reference data at all (family maximum 0) is not scored — no w is defined
  (CodonW refuses to build such a w table; Biopython's 0.5 pseudo-count would give w = 1).
- Reference values are rescaled by the family maximum, so RSCU or w tables are both accepted
  (a proper w table — family maximum 1, as CodonW `-cai_file` expects — is used unchanged).
- Biopython 1.88 `CodonAdaptationIndex.calculate` scores stop codons when the index contains them
  (built from sequences) — not followed (Sharp & Li / CodonW exclude stops). Biopython ≤1.79
  `cai_for_gene` divides by `L − 1` (bug) — not followed.
- CAI against `CreateCodonTableFromSequence` tables uses 0.01 for codons absent from the reference
  sequence, whereas Sharp & Li / Biopython / CodonW w-generation use a 0.5 pseudo-count; building w
  from counts belongs to the table builder (CODON-OPT-001 lead).

## Open Questions

None — algorithm well-documented in Sharp & Li (1987).

## Decisions

1. Create new canonical test file: `CodonOptimizer_CAI_Tests.cs`
2. Consolidate CAI tests from `CodonOptimizerTests.cs`
3. Keep `CodonUsageAnalyzerTests.cs` CAI tests as smoke verification (different class)
