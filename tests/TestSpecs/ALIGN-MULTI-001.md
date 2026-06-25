# TestSpec: ALIGN-MULTI-001 - Multiple Sequence Alignment

## Test Unit Information
- **ID:** ALIGN-MULTI-001
- **Area:** Alignment
- **Canonical Method:** `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)`
- **Complexity:** O(k² × m) where k = sequences, m = average length
- **Algorithm:** Anchor-based star alignment (progressive alignment with suffix tree anchors)
- **Status:** ☑ Complete

---

## Evidence Sources

| Source | Type | URL |
|--------|------|-----|
| Wikipedia - Multiple Sequence Alignment | Primary | https://en.wikipedia.org/wiki/Multiple_sequence_alignment |
| Wikipedia - Clustal | Primary | https://en.wikipedia.org/wiki/Clustal |
| Wikipedia - Consensus sequence | Primary | https://en.wikipedia.org/wiki/Consensus_sequence |

---

## Test Categories

### MUST Tests (Required - Evidence-backed)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| M01 | `MultipleAlign_NullInput_ThrowsArgumentNullException` | .NET null-safety convention | .NET convention |
| M02 | `MultipleAlign_EmptyCollection_ReturnsEmpty` | Edge case: empty input | Wikipedia MSA |
| M03 | `MultipleAlign_SingleSequence_ReturnsSameSequence` | Trivial case; SP=0 for 0 pairs | Wikipedia MSA |
| M04 | `MultipleAlign_TwoSequences_AlignsCorrectly` | Minimum non-trivial; exact SP=4 | Wikipedia MSA |
| M05 | `MultipleAlign_ThreeIdenticalSequences_AllMatch` | Perfect alignment; exact SP=24 | Wikipedia MSA |
| M06 | `MultipleAlign_AllAlignedSequences_HaveEqualLength` | MSA invariant: L ≥ max{nᵢ} | Wikipedia MSA: "all conform to length L" |
| M07 | `MultipleAlign_SequenceCount_Preserved` | Invariant: output count = input count | Wikipedia MSA |
| M08 | `MultipleAlign_DifferentLengths_PadsWithGaps` | Gap insertion + reversibility | Wikipedia MSA |
| M09 | `MultipleAlign_ConsensusContainsOnlyValidCharacters` | Consensus ∈ {A,C,G,T,-} | Wikipedia MSA |
| M10 | `MultipleAlign_TotalScore_IsSumOfPairwiseScores` | Hand-computed SP=8 with mismatches | Wikipedia MSA: "sum of all of the pairs of characters at each position" |
| M11 | `MultipleAlign_RemovingGaps_RecoversOriginal` | Reversibility invariant | Wikipedia MSA: "To return from S'_i to S_i, remove all gaps" |
| M12 | `MultipleAlign_NoColumnIsAllGaps` | All-gap column prohibition | Wikipedia MSA: "no column consists of only gaps" |

### SHOULD Tests (Recommended - Quality/Robustness)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| S01 | `MultipleAlign_ConsensusReflectsMajority` | Exact consensus via majority voting | Wikipedia Consensus sequence: "most frequent residues... at each position" |
| S02 | `MultipleAlign_WithCustomScoring_UsesProvidedMatrix` | Exact SP: SimpleDna=4, BlastDna=8 | Wikipedia Clustal: "gap penalties are configurable" |
| S03 | `MultipleAlign_PartiallyOverlapping_AlignsCorrectly` | Invariants + reversibility for diverse input | Wikipedia MSA: arbitrary sequence sets |
| S04 | `MultipleAlign_ConsensusTieBreaking_PrefersNucleotideOverGap` | Gap-nucleotide tie resolved to nucleotide | Implementation design choice (nucleotide > gap on tie) |

### COULD Tests (Optional - Extended coverage)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| C01 | `MultipleAlign_ManySequences_CompletesInReasonableTime` | Performance: 20 diverse seqs, 5s timeout | Wikipedia Clustal: ClustalW complexity is O(N²) |
| C02 | `MultipleAlign_WithEmptySequence_HandlesGracefully` | Invariants + reversibility for degenerate input | Wikipedia MSA: edge case |

---

## Coverage Classification

### Discovery Summary

- **Canonical:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_MultipleAlign_Tests.cs` — 18 tests
- **Property:** `tests/Seqeron/Seqeron.Genomics.Tests/Properties/AlignmentProperties.cs` — 2 MSA tests
- **Benchmark:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_MultipleAlign_Benchmark.cs` — 9 tests

### Canonical (`SequenceAligner_MultipleAlign_Tests.cs`) — 18 test methods

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `MultipleAlign_NullInput_ThrowsArgumentNullException` | M01 | ✅ |
| 2 | `MultipleAlign_EmptyCollection_ReturnsEmpty` | M02 | ✅ |
| 3 | `MultipleAlign_SingleSequence_ReturnsSameSequence` | M03 | ✅ |
| 4 | `MultipleAlign_TwoSequences_AlignsCorrectly` | M04 | ✅ |
| 5 | `MultipleAlign_ThreeIdenticalSequences_AllMatch` | M05 | ✅ |
| 6 | `MultipleAlign_AllAlignedSequences_HaveEqualLength` | M06 | ✅ |
| 7 | `MultipleAlign_SequenceCount_Preserved` | M07 | ✅ |
| 8 | `MultipleAlign_DifferentLengths_PadsWithGaps` | M08 | ✅ |
| 9 | `MultipleAlign_ConsensusContainsOnlyValidCharacters` | M09 | ✅ |
| 10 | `MultipleAlign_TotalScore_IsSumOfPairwiseScores` | M10 | ✅ |
| 11 | `MultipleAlign_RemovingGaps_RecoversOriginal` | M11 | ✅ |
| 12 | `MultipleAlign_NoColumnIsAllGaps` | M12 | ✅ |
| 13 | `MultipleAlign_ConsensusReflectsMajority` | S01 | ✅ |
| 14 | `MultipleAlign_WithCustomScoring_UsesProvidedMatrix` | S02 | ✅ |
| 15 | `MultipleAlign_PartiallyOverlapping_AlignsCorrectly` | S03 | ✅ |
| 16 | `MultipleAlign_ConsensusTieBreaking_PrefersNucleotideOverGap` | S04 | ✅ |
| 17 | `MultipleAlign_ManySequences_CompletesInReasonableTime` | C01 | ✅ |
| 18 | `MultipleAlign_WithEmptySequence_HandlesGracefully` | C02 | ✅ |

#### Property (`AlignmentProperties.cs`) — 2 MSA tests

| # | Test Method | Status | Notes |
|---|-------------|--------|-------|
| 1 | `MultipleAlign_AllSequences_HaveEqualLength` | ✅ | Structural invariant, different inputs from M06 |
| 2 | `MultipleAlign_ConsensusLength_EqualsAlignedLength` | ✅ | Consensus length = aligned length |

#### Classification Summary

- ✅ Covered: 18 canonical + 2 property = 20 total
- ❌ Missing: 0
- ⚠ Weak: 0
- 🔁 Duplicate: 0

### Work Queue

| # | Test | Prior Status | Action Taken | Final Status |
|---|------|-------------|--------------|--------------|
| 1 | M04 | ⚠ Weak | `TotalScore > 0` → exact `== 4` (hand-computed) | ✅ Done |
| 2 | M05 | ⚠ Weak | Added exact `TotalScore == 24` (hand-computed) | ✅ Done |
| 3 | M08 | ⚠🔁 Weak+Dup | Added gap-existence + reversibility (distinct from M06) | ✅ Done |
| 4 | M10 | ⚠ Weak | Rewritten: non-trivial input with mismatches, hand-computed SP=8 | ✅ Done |
| 5 | S01 | ⚠ Weak | New input, exact full consensus "ACTG" (hand-computed) | ✅ Done |
| 6 | S02 | ⚠ Weak | `Not.EqualTo` → exact SP values: SimpleDna=4, BlastDna=8 | ✅ Done |
| 7 | S03 | ⚠ Weak | Added reversibility + L ≥ max{nᵢ} assertions | ✅ Done |
| 8 | S04 | ❌ Missing | NEW: tie-breaking test (gap vs nucleotide → nucleotide) | ✅ Done |
| 9 | C01 | ⚠ Weak | Fixed Random seed (was producing 20 identical seqs); added diversity check + invariants | ✅ Done |
| 10 | C02 | ⚠ Weak | Added equal-length + reversibility assertions | ✅ Done |
| 11 | Helper | Dead code | Removed `ComputeExpectedSPScore` (mirrored implementation) | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **Remaining:** 0

---

## Invariants Under Test

1. **Equal Length:** `result.AlignedSequences.All(s => s.Length == result.AlignedSequences[0].Length)` — M06, M08, S03, C01, C02
2. **Count Preservation:** `result.AlignedSequences.Length == inputSequences.Count()` — M07, S03, C01, C02
3. **Consensus Validity:** `result.Consensus.All(c => "ACGT-".Contains(c))` — M09
4. **Reversibility:** Removing gaps from aligned sequence recovers original — M08, M11, S03, C02
5. **No all-gap columns:** No column has gaps in all sequences — M12
6. **Sum-of-pairs score:** TotalScore = Σ column-based SP across all C(k,2) pairs — M04, M05, M10, S02
7. **Consensus majority:** Consensus character is majority-voted per column — S01, S04

---

## Test Data

### Canonical Test Sequences

```csharp
// Identical sequences — SP = C(3,2) × 8 × Match(1) = 24
["ATGCATGC", "ATGCATGC", "ATGCATGC"] → Consensus: "ATGCATGC", SP: 24

// Same-length with mismatch — hand-computed SP = 8
["ATGC", "ATGC", "CTGC"] → Consensus: "ATGC", SP: 8

// Majority voting — col1: A=1, C=2 → C
["AATG", "ACTG", "ACTG"] → Consensus: "ACTG"

// Tie-breaking — col1: T=2, '-'=2 → nucleotide T preferred
["ATGC", "ATGC", "AGC", "AGC"] → Consensus contains no gaps
```

---

## Design Decisions

| Decision | Rationale | Source |
|----------|-----------|--------|
| Hand-computed SP values (not helper) | Avoids mirroring implementation logic in test code | Testing best practice |
| Exact assertions over range checks | `Is.EqualTo(8)` instead of `Is.GreaterThan(0)` — catches regressions precisely | Testing best practice |
| Dynamic tie-column detection (S04) | Robust to alignment changes while verifying tie-breaking behavior | Wikipedia MSA consensus definition |
| Diverse random seqs (C01) | Single shared `Random(42)` instance produces 20 unique sequences | Bug fix: was creating new `Random(42)` per call |

---

## Definition of Done

- [x] Evidence document created
- [x] Algorithm documentation created
- [x] Tests implemented in canonical file
- [x] All MUST tests passing (12 tests)
- [x] All SHOULD tests passing (4 tests)
- [x] All COULD tests passing (2 tests)
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
- [x] All SP scores hand-computed from theory, not taken from implementation output
- [x] Zero warnings
- [x] Zero assumptions
- [x] Checklist updated

---

# Addendum: Iterative Refinement (`MultipleAlignIterative`) — 2026-06-23

This addendum covers the iterative-refinement aligner added to remove the progressive method's
single-pass "once a gap, always a gap" limitation.

- **Canonical Method:** `SequenceAligner.MultipleAlignIterative(IEnumerable<DnaSequence>, ScoringMatrix?, int maxIterations)`
- **Type:** Canonical (deep evidence-based testing)
- **Status:** ☐ Pending re-validation (checklists reset)
- **Canonical test file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_MultipleAlignIterative_Tests.cs`

## Evidence Sources (addendum)

| Source | Type | URL |
|--------|------|-----|
| Edgar (2004) MUSCLE, NAR 32(5):1792-1797 — Stage 3 steps 3.1-3.4 | Primary (rank 1) | https://academic.oup.com/nar/article/32/5/1792/2380623 |
| Barton & Sternberg (1987) J Mol Biol 198:327-337 | Primary (rank 1) | https://pubmed.ncbi.nlm.nih.gov/3430611/ |
| Wallace, O'Sullivan & Higgins (2005) Bioinformatics 21(8):1408 | Primary (rank 1) | https://academic.oup.com/bioinformatics/article/21/8/1408/249176 |
| Wikipedia — Multiple sequence alignment (SP score) | rank 4 | https://en.wikipedia.org/wiki/Multiple_sequence_alignment |

## Invariants (addendum)

| ID | Invariant | Evidence |
|----|-----------|----------|
| INV-06 | Iterative SP score ≥ progressive seed SP score (monotonic non-decreasing). | Edgar 2004 step 3.4 ("kept ... if the SP score is improved") |
| INV-07 | Deterministic; converges within the cap; degaps to inputs; all rows equal length; no all-gap column. | Edgar 2004 convergence clause; Wikipedia MSA invariants |

## Test Cases (addendum)

### MUST

| ID | Test Name | Description | Expected | Evidence |
|----|-----------|-------------|----------|----------|
| M01 | `MultipleAlignIterative_NullInput_ThrowsArgumentNullException` | null guard | throws | .NET convention |
| M02 | `MultipleAlignIterative_NonPositiveMaxIterations_Throws` | cap must be positive | throws ArgOutOfRange | impl contract |
| M03 | `MultipleAlignIterative_EmptyCollection_ReturnsEmpty` | empty input | `Empty`, SP=0 | Wikipedia MSA edge |
| M04 | `MultipleAlignIterative_SingleSequence_ReturnsSameSequence` | single seq | verbatim, SP=0 | Wikipedia MSA |
| M05 | `MultipleAlignIterative_TwoSequences_EqualsProgressiveSeed` | no internal edge → equals seed | `ACGT`/`AC-T`, SP=2 | hand-derived |
| M06 | `MultipleAlignIterative_CorrectsEarlyGapPlacement_RaisesSpScore` | constructed gap-relocation correction | seed SP=−8 → refined SP=−6 (+2) | Edgar 2004 + hand-derived SP |
| M07 | `MultipleAlignIterative_SpScore_NeverBelowProgressiveSeed` | monotonicity over a panel | iter SP ≥ prog SP each | Edgar 2004 step 3.4 |
| M08 | `MultipleAlignIterative_Invariants_Hold` | equal length / degap / count / no all-gap col | invariants hold | Wikipedia MSA |
| M09 | `MultipleAlignIterative_AlreadyOptimal_ChangesNothing` | idempotence on identical seqs | gap-free, SP=24 | Edgar 2004 (strict-improve filter) |
| M10 | `MultipleAlignIterative_IsFixedPoint` | refining the refined result | same SP | Edgar 2004 convergence |
| M11 | `MultipleAlignIterative_ConvergesWithinCap` | maxIter 1 == maxIter 100 | SP=−6 both | Edgar 2004 convergence/cap |
| M12 | `MultipleAlignIterative_IsDeterministic` | repeated runs identical | byte-identical | impl contract (no RNG) |

### SHOULD

| ID | Test Name | Description | Expected | Notes |
|----|-----------|-------------|----------|-------|
| S01 | `MultipleAlignIterative_WithCustomScoring_UsesProvidedMatrix` | custom matrix honored | SimpleDna=4, BlastDna=8 | exact |
| S02 | `SiblingAligners_Unchanged_OnDiscriminatingInput` | star & progressive byte-for-byte unchanged | documented rows | additivity guarantee |
| S03 | `MultipleAlignIterative_RandomInputs_NeverRegressAndStayValid` | 500 random inputs: never regress, always valid | iter ≥ prog; degaps; equal len | property test (O(k·L²)+) |

## Coverage (addendum)

- ✅ Covered: 15 (12 MUST + 3 SHOULD) — all in the canonical iterative file.
- ❌ Missing: 0 | ⚠ Weak: 0 | 🔁 Duplicate: 0 | Remaining: 0

## Assumptions (addendum)

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Edge-partition refinement (MUSCLE) chosen over remove-one-sequence (Barton-Sternberg); identical accept semantics. | implementation structure |
| 2 | SP gap-gap = 0 / residue-gap = GapExtend (inherited from `ComputeSumOfPairsScore`). | SP scoring |

---

# Addendum: Consistency-based MSA (`MultipleAlignConsistency`, T-Coffee) — 2026-06-23

Adds the consistency-based aligner that optimises the **T-Coffee consistency objective** (a distinct
objective class from the SP score), via primary library + library extension + progressive alignment
on the extended library.

- **Canonical Method:** `SequenceAligner.MultipleAlignConsistency(IEnumerable<DnaSequence>, ScoringMatrix?)`
- **Type:** Canonical (deep evidence-based testing)
- **Status:** ☐ In Progress (checklists reset)
- **Canonical test file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_MultipleAlignConsistency_Tests.cs`

## Evidence Sources (addendum)

| Source | Type | URL |
|--------|------|-----|
| Notredame, Higgins & Heringa (2000) T-Coffee, J Mol Biol 302:205–217 | Primary (rank 1) | https://doi.org/10.1006/jmbi.2000.4042 |
| T-Coffee Technical Documentation (min() rule) | Project docs (rank 3) | https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html |
| Gotoh (1982) DP, J Mol Biol 162:705–708 | Primary (rank 1) | https://doi.org/10.1016/0022-2836(82)90398-9 |

## Invariants Under Test (addendum)

| ID | Invariant | Evidence |
|----|-----------|----------|
| INV-08 | Extended-library weight = primary weight + Σ min(W1,W2) over informative triplets; ≥ primary (extension never lowers). | T-Coffee p.209 (GARFIELD: 88 → 165) |
| INV-09 | A consistency-supported residue pair has strictly greater extended weight than an unsupported pair. | T-Coffee p.209 |
| INV-10 | Result is a valid MSA: rows equal length; degap recovers each input; no all-gap column; count preserved. | MSA definition |
| INV-11 | On a case engineered to mislead a single pairwise alignment, the consistency objective of the consistency alignment ≥ that of the plain progressive seed. (Not a universal guarantee: T-Coffee is a progressive heuristic; on arbitrary inputs a different UPGMA merge order can lower the total objective.) | T-Coffee central claim |
| INV-12 | Deterministic (no RNG); identical inputs → trivial gapless exact alignment. | T-Coffee identity corner case |

## Test Cases (addendum)

### MUST

| ID | Test Name | Description | Expected | Evidence |
|----|-----------|-------------|----------|----------|
| TM01 | `MultipleAlignConsistency_NullInput_ThrowsArgumentNullException` | null guard | throws | .NET convention |
| TM02 | `MultipleAlignConsistency_EmptyCollection_ReturnsEmpty` | empty input | `Empty`, SP=0 | MSA edge |
| TM03 | `MultipleAlignConsistency_SingleSequence_ReturnsSameSequence` | single seq | verbatim, SP=0 | MSA |
| TM04 | `ExtendedWeight_ConsistencyPair_ExceedsPrimaryWeight` | GARFIELD-style integer relation on the primary→extended weight of a pair supported through an intermediate | extended == primary + min(W1,W2) and > primary | T-Coffee p.209 (88→165) |
| TM05 | `ExtendedWeight_SupportedPair_GreaterThanUnsupportedPair` | a pair backed by an intermediate vs an inconsistent pair | supported > unsupported | T-Coffee p.209 |
| TM06 | `MultipleAlignConsistency_IdenticalSequences_TrivialExactAlignment` | 3 identical seqs | rows == input, gap-free | identity corner case |
| TM07 | `MultipleAlignConsistency_ValidMsa_Invariants` | equal length / degap / count / no all-gap col | invariants hold | MSA definition |
| TM08 | `MultipleAlignConsistency_ConsistencyObjective_NotBelowProgressiveSeed` | objective(consistency) ≥ objective(progressive seed) on a misleading case | ≥ holds | T-Coffee central claim |
| TM09 | `MultipleAlignConsistency_TwoSequences_EqualsPairwiseGlobal` | k=2: no triplets → pairwise global | rows == global align rows | extended=primary |
| TM10 | `MultipleAlignConsistency_IsDeterministic` | repeat runs identical | byte-identical | no RNG |

### SHOULD

| ID | Test Name | Description | Expected | Notes |
|----|-----------|-------------|----------|-------|
| TS01 | `SiblingAligners_Unchanged_WhenConsistencyAdded` | star / progressive / iterative byte-for-byte unchanged | documented rows | additivity guarantee |
| TS02 | `MultipleAlignConsistency_RandomInputs_AlwaysValidMsa` | random panel (100 trials): always a structurally valid MSA | equal length / degap / count / no all-gap col | property test (O(N³L²) extension path); objective dominance asserted only on the engineered TM08 case |

## Coverage (addendum)

- ✅ Covered: 12 (10 MUST + 2 SHOULD) — all in the canonical consistency test file.
- ❌ Missing: 0 | ⚠ Weak: 0 | 🔁 Duplicate: 0 | Remaining: 0

## Assumptions (addendum)

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Primary weight = integer percent identity of the pairwise global alignment; local entries from Smith–Waterman added by signal addition. | primary library |
| 2 | Guide tree = existing UPGMA (NJ in paper); affects merge order only. | progressive step |
| 3 | Progressive DP gap score = 0 (T-Coffee p.210); group columns use average extended-library score. | progressive DP |
