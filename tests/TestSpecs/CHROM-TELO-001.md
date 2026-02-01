# TestSpec: CHROM-TELO-001

## Test Unit Information
- **ID:** CHROM-TELO-001
- **Title:** Telomere Analysis
- **Area:** Chromosome
- **Created:** 2026-02-01
- **Status:** Complete

---

## Methods Under Test

### Canonical Methods
| Method | Class | Type | Test Depth |
|--------|-------|------|------------|
| `AnalyzeTelomeres(chrName, seq, repeat, ...)` | ChromosomeAnalyzer | Canonical | Deep |
| `EstimateTelomereLengthFromTSRatio(tsRatio, ...)` | ChromosomeAnalyzer | Canonical | Deep |

### Constants
| Constant | Value | Description |
|----------|-------|-------------|
| `HumanTelomereRepeat` | "TTAGGG" | Default vertebrate repeat |

---

## Test Categories

### MUST Tests (Required for DoD)

#### AnalyzeTelomeres

| ID | Test Name | Description | Source |
|----|-----------|-------------|--------|
| M1 | 3PrimeTelomere_Detected | TTAGGG repeats at 3' end detected | Wikipedia, Meyne (1989) |
| M2 | 5PrimeTelomere_Detected | CCCTAA repeats at 5' end detected | Wikipedia |
| M3 | BothEnds_BothDetected | Telomeres at both ends detected independently | Wikipedia |
| M4 | EmptySequence_NoTelomere | Empty string returns no telomere, critically short | Edge case |
| M5 | NoRepeats_NoTelomere | Random sequence returns no telomere | Edge case |
| M6 | CriticallyShort_Flagged | Telomere below critical threshold flagged | Wikipedia |
| M7 | Length_MatchesRepeatCount | Measured length = repeats × 6 | Algorithm definition |
| M8 | Purity_HighForPerfectRepeats | Perfect repeats → purity > 0.95 | Algorithm definition |
| M9 | MinThreshold_Respected | Telomere below min threshold → HasTelomere = false | API contract |

#### EstimateTelomereLengthFromTSRatio

| ID | Test Name | Description | Source |
|----|-----------|-------------|--------|
| M10 | Ratio1_ReturnsReference | T/S=1.0 → returns reference length | Cawthon (2002) |
| M11 | HigherRatio_LongerTelomere | T/S > 1 → proportionally longer | Cawthon (2002) |
| M12 | LowerRatio_ShorterTelomere | T/S < 1 → proportionally shorter | Cawthon (2002) |
| M13 | ZeroRatio_ReturnsZero | T/S = 0 → length = 0 | Edge case |
| M14 | LinearProportionality | length = refLen × (ratio / refRatio) | Cawthon (2002) |

---

### SHOULD Tests (High Value)

| ID | Test Name | Description | Source |
|----|-----------|-------------|--------|
| S1 | CustomRepeat_Detected | Non-default repeat (e.g., TTTAGGG) works | Species variation |
| S2 | CaseInsensitive_Works | Lowercase sequence handled correctly | Robustness |
| S3 | DivergentRepeats_LowerPurity | Imperfect repeats reduce purity score | Algorithm behavior |
| S4 | LongTelomere_FullyMeasured | Very long telomere (>10kb) measured correctly | Boundary |
| S5 | SearchLength_Limits | Only searches within searchLength from ends | API contract |

---

### COULD Tests (Nice to Have)

| ID | Test Name | Description | Source |
|----|-----------|-------------|--------|
| C1 | ChromosomeName_PreservedInResult | Name passed through to result | API contract |
| C2 | VertebrateRepeats_AllSupported | TTAGGG works for mouse, Xenopus, etc. | Wikipedia table |
| C3 | VeryShortSequence_HandledGracefully | Sequence < repeat length → no crash | Edge case |

---

## Test Invariants

### AnalyzeTelomeres Invariants
1. `TelomereLength5Prime >= 0`
2. `TelomereLength3Prime >= 0`
3. `0 <= RepeatPurity5Prime <= 1`
4. `0 <= RepeatPurity3Prime <= 1`
5. `Has5PrimeTelomere == (TelomereLength5Prime >= minTelomereLength)`
6. `Has3PrimeTelomere == (TelomereLength3Prime >= minTelomereLength)`
7. `IsCriticallyShort == true` when empty sequence

### EstimateTelomereLengthFromTSRatio Invariants
1. `result >= 0` when `tsRatio >= 0`
2. `result = referenceLength * tsRatio / referenceRatio`

---

## Test Data

### Standard Test Sequences
```
3' Telomere: new string('A', 1000) + string.Concat(Enumerable.Repeat("TTAGGG", 200))
5' Telomere: string.Concat(Enumerable.Repeat("CCCTAA", 200)) + new string('A', 1000)
Both ends: [CCCTAA×200] + [A×1000] + [TTAGGG×200]
No telomere: new string('A', 1000)
Empty: ""
```

### T/S Ratio Test Cases
| tsRatio | refRatio | refLength | Expected |
|---------|----------|-----------|----------|
| 1.5 | 1.0 | 7000 | 10500 |
| 0.5 | 1.0 | 7000 | 3500 |
| 2.0 | 1.0 | 7000 | 14000 |
| 1.0 | 1.0 | 7000 | 7000 |
| 0.0 | 1.0 | 7000 | 0 |

---

## Audit Results

### Existing Tests (ChromosomeAnalyzerTests.cs)
| Test | Status | Action |
|------|--------|--------|
| AnalyzeTelomeres_With3PrimeTelomere_DetectsTelomere | Covered | Keep, enhance assertions |
| AnalyzeTelomeres_With5PrimeTelomere_DetectsTelomere | Covered | Keep |
| AnalyzeTelomeres_CriticallyShort_DetectsCriticalLength | Covered | Keep |
| AnalyzeTelomeres_NoTelomere_ReturnsNoTelomere | Covered | Keep |
| AnalyzeTelomeres_EmptySequence_ReturnsNoTelomere | Covered | Keep |
| EstimateTelomereLengthFromTSRatio_CalculatesCorrectly | Covered | Expand with more cases |

### Missing Tests
| Test | Priority | Added |
|------|----------|-------|
| Both ends detected | Must | Yes |
| Case insensitivity | Should | Yes |
| Length matches repeat count | Must | Yes |
| T/S ratio edge cases | Must | Yes |
| Invariants validation | Must | Yes |

### Consolidation Plan
- **Canonical file:** ChromosomeAnalyzer_Telomere_Tests.cs (new)
- **Action:** Extract telomere tests from ChromosomeAnalyzerTests.cs into dedicated file
- **Duplicates:** None found
- **Enhancements:** Add missing Must tests, improve assertions with Assert.Multiple

---

## Open Questions / Decisions

None - all behavior is well-documented.

---

## ASSUMPTIONS

None - all test rationale backed by cited sources.
