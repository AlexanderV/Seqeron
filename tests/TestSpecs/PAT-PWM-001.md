# Test Specification: PAT-PWM-001

**Test Unit ID:** PAT-PWM-001
**Area:** Pattern Matching
**Algorithm:** Position Weight Matrix (PWM)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-01-22

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Position weight matrix | https://en.wikipedia.org/wiki/Position_weight_matrix | 2026-01-22 |
| Kel et al. (2003) MATCH: A tool for searching TF binding sites | https://pmc.ncbi.nlm.nih.gov/articles/PMC169193/ | 2026-01-22 |
| Rosalind: Consensus and Profile | https://rosalind.info/problems/cons/ | 2026-01-22 |
| Nishida et al. (2008) Pseudocounts for transcription factor binding sites | Nucleic Acids Research 37(3):939-944 | Reference |
| Stormo (2000) DNA binding sites: representation and discovery | Bioinformatics Review | Reference |

### 1.2 Algorithm Description

#### Position Weight Matrix (Wikipedia)

> A position weight matrix (PWM), also known as a position-specific weight matrix (PSWM)
> or position-specific scoring matrix (PSSM), is a commonly used representation of
> motifs (patterns) in biological sequences.

**Construction Process (Wikipedia):**

1. **Position Frequency Matrix (PFM):** Count occurrences of each nucleotide at each position
2. **Position Probability Matrix (PPM):** Normalize by dividing by number of sequences
3. **Position Weight Matrix (PWM):** Convert to log-odds using background model

**Log-odds Formula (Wikipedia):**
$$M_{k,j} = \log_2\left(\frac{M_{k,j}}{b_k}\right)$$

Where:
- $M_{k,j}$ is the probability of nucleotide k at position j
- $b_k$ is the background frequency (typically 0.25 for DNA)

**Pseudocounts (Wikipedia/Nishida et al.):**
> Pseudocounts (or Laplace estimators) are often applied when calculating PPMs if
> based on a small dataset, in order to avoid matrix entries having a value of 0.

**Scoring (Wikipedia):**
> When the PWM elements are calculated using log likelihoods, the score of a sequence
> can be calculated by adding (rather than multiplying) the relevant values at each
> position in the PWM. The sequence score gives an indication of how different the
> sequence is from a random sequence.

**Score Interpretation:**
- Score = 0: Equal probability of being functional vs random site
- Score > 0: More likely to be a functional site
- Score < 0: More likely to be a random site

#### Profile Matrix (Rosalind CONS Problem)

> A profile matrix is a 4×n matrix P in which P_{1,j} represents the number of
> times that 'A' occurs in the jth position of one of the strings.

**Consensus String (Rosalind):**
> A consensus string c is a string of length n formed from our collection by taking
> the most common symbol at each position.

### 1.3 Reference Examples from Evidence

#### Wikipedia PPM Example

**Input Sequences (7 sequences of length 9):**
```
GAGGTAAAC
TCCGTAAGT
CAGGTTGGA
ACAGTCAGT
TAGGTTCAT
TTAGGTACT
GATGGTAAC
```

**Expected PPM (Wikipedia):**
| Base | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|------|---|---|---|---|---|---|---|---|---|
| A | 0.3 | 0.6 | 0.1 | 0.0 | 0.0 | 0.6 | 0.7 | 0.2 | 0.1 |
| C | 0.2 | 0.2 | 0.1 | 0.0 | 0.0 | 0.2 | 0.1 | 0.1 | 0.2 |
| G | 0.1 | 0.1 | 0.7 | 1.0 | 0.0 | 0.1 | 0.1 | 0.5 | 0.1 |
| T | 0.4 | 0.1 | 0.1 | 0.0 | 1.0 | 0.1 | 0.1 | 0.2 | 0.6 |

#### Rosalind CONS Example

**Input Sequences (7 sequences of length 8):**
```
ATCCAGCT
GGGCAACT
ATGGATCT
AAGCAACC
TTGGAACT
ATGCCATT
ATGGCACT
```

**Expected Profile Matrix:**
| Base | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|------|---|---|---|---|---|---|---|---|
| A | 5 | 1 | 0 | 0 | 5 | 5 | 0 | 0 |
| C | 0 | 0 | 1 | 4 | 2 | 0 | 6 | 1 |
| G | 1 | 1 | 6 | 3 | 0 | 1 | 0 | 0 |
| T | 1 | 5 | 0 | 0 | 0 | 1 | 1 | 6 |

**Expected Consensus:** ATGCAACT

### 1.4 Edge Cases from Evidence

| Edge Case | Expected Behavior | Source |
|-----------|-------------------|--------|
| Empty alignment | Exception | Implementation contract |
| Single sequence | Valid PWM | Mathematical definition |
| Unequal length sequences | Exception | Wikipedia (same length required) |
| All same base at position | Log-odds = log2(1/0.25) ≈ 2.0 | Wikipedia formula |
| Pseudocount = 0 | Risk of -∞ for unseen bases | Wikipedia |
| Default pseudocount 0.25 | Avoids zero probabilities | Nishida et al. |
| Threshold at boundary | Include/exclude based on >= | Implementation |
| High threshold | Few or no matches | Definition |
| Sequence shorter than PWM | No matches | Definition |
| Non-ACGT characters | Skip position (score invalid) | Implementation |
| Perfect match | MaxScore | Definition |
| Worst match | MinScore | Definition |
| Null input | ArgumentNullException | Implementation contract |

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CreatePwm(IEnumerable<string>, double)` | MotifFinder | **Canonical** | PWM construction |
| `ScanWithPwm(DnaSequence, PWM, double)` | MotifFinder | **Canonical** | Sequence scanning |

### 2.1 Supporting Types

| Type | Description |
|------|-------------|
| `PositionWeightMatrix` | PWM data structure with Matrix, Length, Consensus, MaxScore, MinScore |
| `MotifMatch` | Match result with Position, MatchedSequence, Pattern, Score |

---

## 3. Test Classification

### 3.1 Must Tests (Evidence-Based)

| ID | Test Case | Expected | Source |
|----|-----------|----------|--------|
| M1 | Single sequence creates valid PWM | Length = sequence length, Consensus = sequence | Mathematical definition |
| M2 | Multiple identical sequences create same PWM | Consensus = common sequence | Definition |
| M3 | Consensus derives from highest scoring base at each position | Follows max rule | Wikipedia/Rosalind |
| M4 | Empty sequences throws ArgumentException | Exception | Contract |
| M5 | Unequal length sequences throws ArgumentException | Exception | Wikipedia requirement |
| M6 | PWM MaxScore > MinScore for non-uniform matrix | MaxScore > MinScore | Definition |
| M7 | Log-odds scores: perfect match gets max score at position | Verified numerically | Wikipedia formula |
| M8 | ScanWithPwm finds exact trained sequence | Match found | Definition |
| M9 | Threshold filters results correctly | Only scores >= threshold | Definition |
| M10 | Null sequence in CreatePwm throws | ArgumentNullException | Contract |
| M11 | Null sequence in ScanWithPwm throws | ArgumentNullException | Contract |
| M12 | Null PWM in ScanWithPwm throws | ArgumentNullException | Contract |
| M13 | PWM length matches input sequence length | pwm.Length = sequences[0].Length | Definition |
| M14 | Rosalind CONS consensus test case | Consensus matches expected | Rosalind |
| M15 | Non-ACGT characters in scanned sequence skip position | Match marked invalid | Implementation |
| M16 | Sequence shorter than PWM returns no matches | Empty result | Definition |

### 3.2 Should Tests (Recommended)

| ID | Test Case | Expected | Source |
|----|-----------|----------|--------|
| S1 | Pseudocount prevents zero probabilities | All scores finite | Wikipedia |
| S2 | Case-insensitive input handling | Uppercase normalization | Implementation |
| S3 | Multiple matches returned in order | Sorted by position | Implementation |
| S4 | MatchedSequence property populated correctly | Substring at match position | Implementation |
| S5 | Score property reflects log-odds sum | Numeric value | Wikipedia |
| S6 | High threshold (near MaxScore) returns few matches | Filtered correctly | Definition |

### 3.3 Could Tests (Optional)

| ID | Test Case | Expected | Notes |
|----|-----------|----------|-------|
| C1 | Large PWM performance | Reasonable time | Performance |
| C2 | Genome-scale scanning | Memory efficient | Performance |

---

## 4. Invariants

| Invariant | Description | Test Method |
|-----------|-------------|-------------|
| **PWM Length** | PWM.Length equals input sequence length | Assert.That(pwm.Length, Is.EqualTo(seqLength)) |
| **Consensus Length** | Consensus length equals PWM length | Assert.That(pwm.Consensus.Length, Is.EqualTo(pwm.Length)) |
| **MaxScore >= MinScore** | Maximum always >= minimum | Assert.That(pwm.MaxScore, Is.GreaterThanOrEqualTo(pwm.MinScore)) |
| **Score Range** | All match scores between MinScore and MaxScore | Assert within bounds |
| **Matrix Dimensions** | Matrix is 4 x Length | Assert.That(pwm.Matrix.GetLength(0), Is.EqualTo(4)) |

---

## 5. Test Consolidation Plan

### 5.1 Current State

**Existing tests in MotifFinderTests.cs (PWM region):**
- `CreatePwm_SingleSequence_CreatesMatrix` - basic happy path
- `CreatePwm_MultipleSequences_CreatesAverageMatrix` - multiple input
- `CreatePwm_EmptySequences_ThrowsException` - edge case
- `CreatePwm_UnequalLengths_ThrowsException` - edge case
- `ScanWithPwm_FindsExactMatch` - basic scanning
- `ScanWithPwm_ReturnsScores` - score verification
- `Pwm_MaxMinScore_Calculated` - property verification
- `CreatePwm_NullSequences_ThrowsException` - null handling
- `ScanWithPwm_NullSequence_ThrowsException` - null handling
- `ScanWithPwm_NullPwm_ThrowsException` - null handling

### 5.2 Consolidation Actions

| Action | From | To | Rationale |
|--------|------|-----|-----------|
| Move & expand | MotifFinderTests.cs (#region PWM Tests) | MotifFinder_PWM_Tests.cs | Dedicated canonical test file |
| Retain smoke tests | - | MotifFinderTests.cs | API verification (2-3 tests) |
| Remove duplicates | - | - | None identified |
| ~~Add missing~~ | - | ✅ Done | M14 (Rosalind), threshold tests, invariant tests — all added |

### 5.3 Final Test Structure

```
MotifFinder_PWM_Tests.cs (Canonical - PAT-PWM-001)
├── CreatePwm Construction Tests
│   ├── CreatePwm_SingleSequence_CreatesValidMatrix
│   ├── CreatePwm_MultipleIdenticalSequences_CreatesSameConsensus
│   ├── CreatePwm_MixedSequences_CreatesExpectedConsensus
│   ├── CreatePwm_RosalindCONS_TestCase
│   ├── CreatePwm_EmptySequences_ThrowsArgumentException
│   ├── CreatePwm_UnequalLengths_ThrowsArgumentException
│   ├── CreatePwm_NullSequences_ThrowsArgumentNullException
│   └── CreatePwm_CaseInsensitive_NormalizesToUppercase
├── PWM Properties Tests
│   ├── Pwm_Length_MatchesInputSequenceLength
│   ├── Pwm_Consensus_HasCorrectLength
│   ├── Pwm_MaxScore_GreaterThanOrEqualToMinScore
│   ├── Pwm_Matrix_HasCorrectDimensions
│   └── Pwm_LogOdds_PerfectMatchGetsMaxPositionalScore
├── ScanWithPwm Tests
│   ├── ScanWithPwm_FindsTrainedSequence
│   ├── ScanWithPwm_ReturnsCorrectPositions
│   ├── ScanWithPwm_ReturnsMatchedSequence
│   ├── ScanWithPwm_ScoreWithinValidRange
│   ├── ScanWithPwm_ThresholdFiltersResults
│   ├── ScanWithPwm_HighThreshold_ReturnsFewerMatches
│   ├── ScanWithPwm_SequenceShorterThanPwm_ReturnsEmpty
│   ├── ScanWithPwm_NonAcgtCharacter_SkipsPosition
│   ├── ScanWithPwm_NullSequence_ThrowsArgumentNullException
│   └── ScanWithPwm_NullPwm_ThrowsArgumentNullException
└── Invariant Tests
    └── Pwm_AllInvariants_HoldForValidInput

MotifFinderTests.cs (Smoke only - retained)
└── #region PWM Tests (Smoke)
    ├── CreatePwm_Smoke_ReturnsValidMatrix
    └── ScanWithPwm_Smoke_FindsMatch
```

---

## 6. ASSUMPTIONS

| ID | Assumption | Rationale |
|----|------------|-----------|
| A1 | Default background frequency is 0.25 for each base | Standard for DNA (Wikipedia) |
| A2 | Default pseudocount is 0.25 | Implementation default, reasonable (Nishida) |
| A3 | Non-ACGT characters invalidate the position during scanning | Implementation-specific, not in sources |
| A4 | Case-insensitive input normalization | Standard bioinformatics practice |

---

## 7. Open Questions

None - all behaviors are documented in sources or marked as ASSUMPTIONS.

---

## 8. Audit Results

| Category | Count | Notes |
|----------|-------|-------|
| Covered | 16 | All tests implemented |
| Missing | 0 | ~~Rosalind test case, threshold tests, invariants, edge cases~~ ✅ All added |
| Weak | 0 | ~~Score verification~~ ✅ Strengthened |
| Duplicate | 0 | None |
| Incorrect | 0 | None |

---

## 9. Validation Criteria

- [ ] All Must tests pass
- [ ] All Should tests pass
- [ ] Zero warnings
- [ ] Tests are deterministic
- [ ] Single canonical test file for PWM
- [ ] Smoke tests retained in MotifFinderTests.cs
