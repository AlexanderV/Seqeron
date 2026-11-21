# TestSpec: ALIGN-MULTI-001 - Multiple Sequence Alignment

## Test Unit Information
- **ID:** ALIGN-MULTI-001
- **Area:** Alignment
- **Canonical Method:** `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)`
- **Complexity:** O(k² × m²) (exact center-star pair scores)
- **Algorithm:** Center-star alignment (Gusfield 1993) with exact pairwise Needleman-Wunsch alignments; variants: `MultipleAlignAnchored` (suffix-tree anchors, approximate), `MultipleAlignProgressive`, `MultipleAlignIterative`, `MultipleAlignConsistency`
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
- **Benchmark:** the former `[Explicit]` fixture `SequenceAligner_MultipleAlign_Benchmark` is gone (A6-8d: no tier ran it); its scenarios are the deterministic tests `SequenceAligner_MultipleAlignAnchored_Scenarios_Tests` (8) and the timing table is `apps/Seqeron.Alignment.Benchmarks`

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
- **Status:** ☑ Validated (the 2026-06 checklist reset is superseded by the audit rounds of this review; see the audit sections of this file)
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
| 1 | Primary weight = ⌊100·identical / aligned (non-gap) columns⌋ of the pairwise alignment (T-Coffee `get_seq_sim` sim1); local entries from Smith–Waterman added by signal addition. | primary library |
| 2 | Guide tree = existing UPGMA (NJ in paper); affects merge order only. | progressive step |
| 3 | Progressive DP gap score = 0 (T-Coffee p.210); group columns use average extended-library score. | progressive DP |

---

## B09 re-validation addendum (2026-11-04)

Tests added (all lock oracle values: Biopython 1.88 pair optima and independent Python references, see Evidence addendum):

| ID | Test | Expected | Source |
|---|---|---|---|
| B01-B03 | `MultipleAlign_TwoSequences_IsOptimalPairwiseAlignment` ×3 | SP = Biopython optimum 24 / 27 / −3 (previous anchored stage: 9 / 22 / −35) | Biopython 1.88; Gusfield 1993 |
| B04-B06 | `MultipleAlign_{Three,Four,Five}Sequences_…MatchesIndependentCenterStar` | exact rows; SP −3 / −5 / 10 | independent Python center-star |
| B07 | `MultipleAlign_CenterIsTheSequenceWithBestTotalPairScore_NotTheFirst` | totals −2/4/2 → center S1 | Gusfield center |
| B08-B10 | determinism, random validity (150 sets), empty member inside set | invariants INV-01..05, INV-11 | MSA definition |
| B11-B14 | `MultipleAlignAnchored_*` | valid rows, ≤ optimum, identical → gap-free (SP 30) | approximate variant contract |
| C01-C07 | `AnchorBasedAligner_Tests` (`ChainAnchors`) | overlap trimming (0,0,10)+(6,20,10) → (0,0,10),(10,24,6); coverage over count; crossing anchors; random invariants | Abouelhoda & Ohlebusch 2005; Delcher 1999 |
| C08-C12 | `AlignWithAnchors_*` | identical; empty side = all-gap alignment; linear gap scoring; related pairs ≤ `GlobalAlign`; repeated blocks | MUMmer scheme |
| C13-C14 | `GlobalAlignLinearSpace_*` | Hirschberg score = NW optimum (24/27 Biopython; 200 random pairs × 4 schemes, forced recursion) | Hirschberg 1975 |
| T01 | `LibraryWeights_PercentIdentityExcludesGapColumns_AndExtensionAddsMinOfLegs` | primary 200, extended 385 | T-Coffee `get_seq_sim`, p.209 |
| T02 | `LibraryWeights_PairWithoutDirectSupport_GetsWeightThroughIntermediate` | primary 0, extended 133 | T-Coffee `residue_pair_extended_list` |
| T03 | `MultipleAlignConsistency_MatchesIndependentTCoffeeReference` | `-CTT`, `AC-T`, `-CAA` | independent Python T-Coffee |
| I01 | `MultipleAlignIterative_CorrectsEarlyGapPlacement_RaisesSpScore` (re-engineered for the ClustalW distance) | −7 → −3 | hand SP, Python refinement |
| M01-M08 | MCP: `multiple_align*` bindings, `global/local/semi_global_align_affine` | Biopython/parasail values (−17, 12 @ (7,12,0,5), 7 `--GA-TTACA--`) | SequenceAligner unit tests |


---

# Addendum: ClustalW guide-tree distance (audit A1-15, B09 F51)

**Class:** `SequenceAligner_ClustalWDistance_Tests` (17 tests) + `Seqeron.Mcp.Alignment.Tests.MultipleAlignTests` (2).
**Subject:** `SequenceAligner.ClustalWPairwiseAlign` (port of ClustalW 2.1 `FullPairwiseAlign`) and `MsaDistanceMode.ClustalWAffine` for `MultipleAlignProgressive/Iterative/Consistency` (additive `MsaOptions`, default `Linear` = unchanged).

| ID | Test | Source / oracle |
|----|------|-----------------|
| CW1 | printed Score, identity count and DP score equal the oracle on 3 600 pairs (ACGT, IUB codes, ≤ 320 nt) | installed `clustalw` 2.1 `-TYPE=DNA -ALIGN` ("Sequences (i:j) Aligned. Score: N" = ⌊100·identities/MIN(len1,len2)⌋) + line-by-line Python port asserted equal to the binary |
| CW2 | hand-run pairs: Score 91 / 40 / 26 (identities 32 / 12 / 8) | clustalw 2.1 binary |
| CW3 | rows spell the local region `s[Start..End]`; identical columns = `Identities`; empty result has −1 coordinates | invariant |
| CW4 | rows re-score to `Score` with the IUB×100 matrix (+1000/−900) and gap 3000 + 666·k (open charged once per gap run) | FullPairwiseAlign parameters |
| CW5 | empty input → distance 1, percent identity 0; lower case = upper case; `DnaSequence` overload = string overload | `FullPairwiseAlign.cpp` (`n == 0 \|\| m == 0 → 1.0`) |
| CW6 | gap penalties are parameters (free gaps bridge an insertion); invalid symbol / negative / NaN penalty / null throw | contract |
| CW7 | progressive and iterative rows with `ClustalWAffine` equal the independent Python on 140 sets × 4 schemes; ≥ 5 sets differ from `Linear` | independent Python (NW, profile DP, exact-rational UPGMA, tree-partition refinement) |
| CW8 | default options = historical behaviour (progressive, iterative, consistency); consistency with ClustalW distances spells its inputs; undefined mode throws; two sequences identical in both modes | regression |
| CW9 (MCP) | `distance = clustalw` gives the discriminating set's ClustalW rows, `linear` the old rows, bad value throws | Python oracle |


---

# Addendum: Affine-gap MSA (audit A1-10, B09 F52)

**Class:** `SequenceAligner_MultipleAlignAffine_Tests` (16 tests) + `Seqeron.Mcp.Alignment.Tests.MultipleAlignTests` (2 `gapModel` tests).
**Subject:** `MsaOptions.GapModel = MsaGapModel.Affine` on `MultipleAlign`, `MultipleAlignAnchored`, `MultipleAlignProgressive`, `MultipleAlignIterative`, `MultipleAlignConsistency` and the objective `SequenceAligner.SumOfPairsScore(rows, scoring, gapModel)` (natural affine sum of pairs, Altschul 1989). Default (`Linear`) unchanged.

| ID | Test | Source / oracle |
|----|------|-----------------|
| AF1 | `SumOfPairsScore` (linear and affine) equals the independent Python re-score on 453 random gapped MSAs; affine SP never exceeds the sum of the Biopython 1.88 pairwise optima | Altschul 1989 definition; Biopython `PairwiseAligner(mode='global')` |
| AF2 | hand-derived example and run semantics (gap-gap columns skipped, a gap in row i followed by a gap in row j = two runs, run continues across gap-gap columns) | definition |
| AF3 | SP of the optimal `GlobalAlignAffine` rows = `Score` (k = 2 objective = Gotoh score) | F1 |
| AF4 | rows of the exact optimum of 260 tiny sets score to the optimum under the C# objective | exact DP over (i, j, l, run state per pair), validated by exhaustive enumeration (lengths ≤ 4 in the generator) |
| AF5 | no method exceeds the exact optimum (5 methods x 260 sets); k = 2: star, progressive, iterative equal the optimum; 3-sequence optimum rate recorded (170 / 189 / 189 of 227) | same |
| AF6 | affine star / progressive / iterative / consistency rows equal the independent Python implementations on 150 sets x 7 schemes | independent Python (Gotoh, SWG, Gusfield merge, UPGMA, affine profile DP, edge refinement, T-Coffee extension) |
| AF7 | two sequences: star/progressive/iterative score and progressive rows equal `GlobalAlignAffine` (120 random pairs x 4 schemes) | F1 / Biopython |
| AF8 | a 10-nt indel costs one opening (affine 122 vs linear 138 on the hand-built trio) | hand-derived |
| AF9 | default options are linear and every method is unchanged; rows spell inputs, equal length, `TotalScore` = `SumOfPairsScore`; iterative ≥ progressive; GapOpen 0 makes both objectives equal; invalid model throws; empty/single trivial | contract |
| AF10 (MCP) | `gapModel=affine` honours `gapOpen` (66 + 2·(46 − 18) vs linear 66 + 2·(46 − 10)); invalid value throws on all five tools | hand-derived |


---

# Addendum: guide trees (audit A1-12, B09 F53) and anchor chaining (audit A1-16, B09 F54)

**Classes:** `SequenceAligner_GuideTree_Tests` (18 tests), `AnchorBasedAligner_ChainEquivalence_Tests` (5), `Seqeron.Mcp.Alignment.Tests.MultipleAlignTests` (2 `guideTree` tests).
**Subject:** `MsaOptions.GuideTree` (`Upgma` default / `NeighborJoining`), the canonical `PhylogeneticAnalyzer.BuildTreeFromMatrix` + `TreeBuildOptions` (`TieTolerance`, `TieOrder`), `SequenceAligner.BuildProgressiveGuideTree` / `BuildClustalWNeighborJoiningGuideTree` (internal); `AnchorBasedAligner.ChainAnchors`.

| ID | Test | Source / oracle |
|----|------|-----------------|
| GT1 | UPGMA guide tree equals the removed private implementation on 900 random / integer-tie / tenth-tie / ClustalW-like matrices (k = 2-14) | former code kept in the test as the oracle |
| GT2 | hand-run UPGMA ties: d(0,3) = d(1,2) → ((0,3),(1,2)); rounding-noise ties go to the lowest pair | hand-derived |
| GT3 | canonical NJ = Biopython 1.88 `DistanceTreeConstructor().nj` (edge by edge, ≤ 1e-6) and canonical UPGMA = scipy `linkage('average')` on 400 random matrices each | Biopython / scipy (`generate_guidetree_biopython_oracle.py`) |
| GT4 | ClustalW NJ: the unrooted tree (|b| < 1e-4 → 0, 5 decimals) equals the `.dnd` of `clustalw -ALIGN -MAXDIV=0 -NEWTREE` and the merge sizes of the rooted tree equal "Group n: Sequences: s" on 420 DNA sets | clustalw 2.1 (`generate_clustalw_guidetree_oracle.py`) |
| GT5 | the oracle covers 2-11 sequences, > 20 negative NJ branches, > 30 sets with duplicate sequences (tied distances) | coverage guard |
| GT6 | two sequences → the single step (0,1); fewer than two → ArgumentOutOfRange | contract |
| GT7 | `TreeBuildOptions`: default = parameterless overload (100 random matrices × UPGMA/NJ); `SecondIndexMajor` changes a tied NJ topology (5-taxon integer matrix found by exhaustive search; splits {0,4},{2,3} vs {1,3},{2,4}); tolerance turns 1e-12 noise into a tie | hand-run Saitou-Nei |
| GT8 | `GuideTree` default = Upgma, explicit Upgma identical rows (40 sets); NJ gives structurally valid MSAs for progressive / iterative / consistency, both distance modes (80 sets), iterative ≥ progressive; non-ultrametric 5-taxon example ((((4,3),2),(0,1)), steps 2,3,2,5); undefined enum throws | invariants |
| GT9 | `Phylogenetics` assembly does not reference `Alignment`, `Alignment` references `Phylogenetics` (plus the solution builds with 0 errors and the Architecture tests stay green) | layering |
| GT10 (MCP) | `guideTree=nj` returns the NJ rows (hand-set where UPGMA and NJ differ: SP 32 vs 28), `upgma` the old rows, aliases accepted, bad value throws on all three tools | library call |
| CH1 | `ChainAnchors` = the O(a²) scan (kept as oracle) on 1 500 random sets (coordinate ranges 6-600, lengths 4-80: dense overlaps and ties to sparse) | former code |
| CH2 | same on 400 MEM-like sets (backbone + repeat copies on shifted diagonals + noise) and 2 × 3 000 anchors | former code |
| CH3 | same on 12 adversarial shapes (identical anchors, same start, nested in one long anchor, anti-diagonal, staircase, shifted diagonals, one-sided overlaps, duplicates, empty, single, invalid entries) | former code |
| CH4 | 10⁵ anchors: chain strictly colinear and non-overlapping, completes in seconds (the O(a²) scan needs ~5·10⁹ tests) | timing / invariant |

---

# Addendum: ClustalW fast pairwise stage and ClustalW multiple alignment (audit A1-15b / A1-10b, B09 F55 / F56)

**Classes:** `SequenceAligner_ClustalWQuick_Tests` (10 tests), `SequenceAligner_ClustalWMultiple_Tests` (9), `SequenceAligner_GuideTree_Tests` (unchanged, now run on ClustalW's NJ arithmetic), `Seqeron.Mcp.Alignment.Tests.MultipleAlignTests` (2 new) and `ToolRegistrationContractTests` (`multiple_align_clustalw`).
**Subject:** `SequenceAligner.ClustalWQuickPairwiseAlign` / `MsaDistanceMode.ClustalWQuick` / `MsaOptions.QuickPairwise`, `SequenceAligner.MultipleAlignClustalW(sequences, ClustalWOptions?, ScoringMatrix?)`, `SequenceAligner.ClustalWNeighborJoining` (internal).

| ID | Test | Source / oracle |
|----|------|-----------------|
| CQ1 | printed `Score:` (or "Not Aligned") of `clustalw -QUICKTREE -ALIGN` equals `ClustalWQuickPairwiseAlign` on 8 257 pairs of 605 sets (k-tuple 1-4, -WINDOW/-TOPDIAGS/-PAIRGAP/-SCORE variants, IUB-N sets, "(Partial alignment)" sets) | installed `clustalw` 2.1 (`generate_clustalw_quick_oracle.py`) |
| CQ2 | the oracle covers k = 1-4, absolute scores, "Not Aligned", explicit pair gaps and reaches the partial-alignment exit | coverage guard |
| CQ3 | rows run with ClustalW's defaults (2 / 4 / 4 / 5) equal a parameterless call; the pair gap of `-KTUPLE=2` alone is 6 (27.5862 vs 25.8621 on a hand-run pair) | clustalw binary |
| CQ4 | hand run: 92.5926 / 27.7778 / 27.7778 | clustalw binary |
| CQ5 | score / distance consistency, defaults, contracts (null, invalid symbol, k-tuple range, window / topdiags / pair gap ≥ 1, maximum length < sequence length, empty, shorter than k) | contract |
| CQ6 | `MsaDistanceMode.ClustalWQuick` gives valid MSAs on 60 sets and changes some guide trees; options reach the guide stage; default distance unchanged | invariants |
| CM1 | rows of `MultipleAlignClustalW` equal `clustalw -TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` on 1 000 sets with option variants | installed `clustalw` 2.1 (`generate_clustalw_msa_full_oracle.py`) |
| CM2 | coverage guards: every option (-GAPOPEN, -GAPEXT, -TRANSWEIGHT, -MAXDIV, -NOWEIGHTS, -NEGATIVE, -QUICKTREE, -KTUPLE, -WINDOW, -TOPDIAGS, -PAIRGAP, -PWGAPOPEN, -PWGAPEXT) in > 20 sets, two-sequence sets, ≥ 14-sequence sets, IUB codes, duplicates | coverage guard |
| CM3 | each mechanism changes the rows on many oracle sets when switched off (weights, delayed divergent sequences, transition weight, negative matrix, gap penalties, fast pairwise stage) — equality is not an accident of easy inputs | discrimination |
| CM4 | hand-run three-sequence example; rows spell the inputs, equal length, input order | clustalw binary / invariant |
| CM5 | contracts (null, empty, single, invalid symbol, empty sequence, negative / NaN penalties, undefined mode, k-tuple range); `DnaSequence` overload = string overload = lower case; defaults = ClustalW DNA defaults; GAPEXT stored as float | contract |
| CM6 (MCP) | `multiple_align_clustalw` delegates (defaults and every option), bad distance / sequence / negative penalty throws; registered in the contract list | library call |
| TC1 | rows of `MultipleAlignTCoffee` equal `t_coffee -seq -lib -usetree -relax_lib {1,0} -output fasta_aln -n_core 1` on 800 committed sets (400 random DNA sets x relax 1/0; k = 3-6; libraries of this repository ("rand") and of the program ("binary")) | installed `t_coffee` 13.41 (`generate_tcoffee_oracle.py`) |
| TC2 | oracle coverage guard (rows >= 780, both relax modes, both library kinds, k = 3..6); rows spell the inputs with equal length; contracts (null, library size mismatch, empty sequence, tree not a permutation, pair outside sequence, self pair) | coverage / contract |
| TC3 | `MsaOptions.ConsistencyEngine` default = `PaperExtension` and unchanged; `TCoffeeLinked` yields a valid deterministic MSA; MCP `consistencyEngine` (`paper` default, `tcoffee`, bad value throws) | invariants |

---

## 2026-11 audit round 1 (A1-11, B09 F60): exact sum-of-pairs optimum

`SequenceAligner.MultipleAlignExact` / `MultipleAlignExactWithStatistics` (+ `MsaOptions.ExactMaxStateCells`, `ExactCarrilloLipman`, record `ExactMsaResult`; MCP `multiple_align_exact`). Class `SequenceAligner_MultipleAlignExact_Tests`, oracle `TestData/Alignment/msa_exact_oracle.tsv` (804 sets, k = 2..5, both objectives; two independent Python exact programs = exhaustive enumeration = Biopython for k = 2; `scripts/alignment/generate_msa_exact_oracle.py`).

| ID | Test | Expected | Source |
|----|------|----------|--------|
| EX1 | `Oracle_OptimalRows_ScoreToTheOptimum_UnderTheObjective` | the oracle's rows score (`SumOfPairsScore`) to its optimum | independent Python |
| EX2 | `ExactOptimum_EqualsTheIndependentOracle_WithAndWithoutCarrilloLipman` | optimum of all 804 sets with and without the bound; rows spell the inputs; visited <= lattice | oracle |
| EX3 | `TwoSequences_AreGlobalAlign_AndGlobalAlignAffine_RowForRow` | k = 2 rows and score = `GlobalAlign` / `GlobalAlignAffine` | canonical |
| EX4 | `Heuristics_NeverExceedTheExactOptimum_AndTheirOptimalShareIsRecorded` | star, progressive, iterative, consistency, anchored <= optimum; optimal share recorded (Evidence) | exact optimum |
| EX5 | `ExactOptimum_EqualsExhaustiveEnumeration_OfTinySets_LinearAndAffine` | 120 tiny k = 3, 4 sets vs an in-test enumeration of every column sequence | definition |
| EX6 | `CarrilloLipman_PrunesSimilarSequences_AndKeepsTheOptimum`, `CarrilloLipman_AgreesWithTheFullLattice_OnMediumSets` | similar k = 4 set: visited < 1 % of > 2 M cells; pruned = full lattice on 12 medium sets | Carrillo-Lipman 1988 |
| EX7 | `IdenticalSequences_ScoreIsPairCountTimesLengthTimesMatch`, `Result_IsDeterministic_AndExactMsaResultCarriesTheBound` | 5 identical sequences: 10 x len x match, one diagonal corridor; k = 4 has 75 states | hand |
| EX8 | `ComplexityGuard_ThrowsClearException_InsteadOfExhaustingMemory`, `TooManySequences_AndNullAndTrivialInputs` | clear `ArgumentException` naming `ExactMaxStateCells`; > 8 / > 6 sequences; null / empty / single / empty member | contract |
| EX9 (MCP) | `GlobalCoOptimalAndExactMsaTests.MultipleAlignExact_*` | >= heuristics, lattice / visited reported, k = 2 = `global_align_affine`, guard, validation | library call |

---

## 2026-11 audit round 1 (A1-14, B09 F61): protein / substitution-matrix multiple alignment

`MultipleAlign*(IEnumerable<ProteinSequence> | IEnumerable<string>, MatrixScoring, MsaOptions?)` for the star, anchored, progressive, iterative, consistency (both engines) and exact aligners, `MatrixSumOfPairsScore`, `MultipleAlignmentResult.ExactTotalScore`, `ExactMsaResult.ScoreScale`; MCP `type = dna | protein`, `matrix` on `multiple_align*`. Class `SequenceAligner_MultipleAlignMatrix_Tests`; oracles `TestData/Alignment/msa_matrix_sp_oracle.tsv`, `msa_matrix_exact_oracle.tsv`, `msa_matrix_oracle.tsv` (`scripts/alignment/generate_msa_matrix_oracle.py`, independent Python, exact rationals; Biopython 1.88 pairwise optima) and `tcoffee_protein_oracle.tsv` (`t_coffee`).

| ID | Test | Expected | Source |
|----|------|----------|--------|
| MX1 | `MatrixSumOfPairsScore_MatchesIndependentRescore_AndRespectsPairwiseOptimumBound`, `_HandDerivedExample`, `_Validation` | 420 random gapped protein / IUPAC-DNA MSAs: linear and affine SP equal the independent re-score and stay below the sum of the Biopython pairwise optima; BLOSUM62 hand examples incl. EMBOSS 10 / 0.5 (-10, -10.5); null / length / alphabet / model errors | independent Python, Biopython |
| MX2 | `MultipleAlignExact_EqualsTheOracleOptimum_AndItsRowsAchieveIt`, `_Rows_...`, `_PruningDoesNotChangeTheOptimum_AndStatisticsAreScaled` | optimum of 308 sets (k = 2: Biopython; k = 3, 4: lattice DP; tiny sets: enumeration) = `ExactTotalScore`, rows re-score to it, with and without Carrillo-Lipman; k = 2 rows = `GlobalAlign(Affine)`; `ScoreScale` 2 for gaps in halves | oracle |
| MX3 | `Aligners_EveryRowEqualsTheIndependentPythonImplementation` | star / anchored / progressive / iterative / consistency rows of 238 sets (BLOSUM62, PAM250, EDNAFULL; linear / affine; fractional gaps) equal the Python implementation (>= 1 150 row checks); `ExactTotalScore` = `MatrixSumOfPairsScore` of the rows | independent Python |
| MX4 | `TwoSequences_EveryMethodIsAnOptimalPairwiseAlignment`, `ConsistencyEngineTCoffeeLinked_...` | k = 2: star / progressive / iterative score = pairwise optimum, star rows = the pairwise alignment; TCoffeeLinked valid on protein | canonical / Biopython |
| MX5 | `MatrixPath_WithMatchMismatchMatrix_ReproducesTheScalarAligners` | a +2 / -3 matrix over ACGT with integer gaps gives the same rows, score and consensus as the scalar aligners for all six methods (40 sets x linear / affine) | no-second-implementation guard |
| MX6 | `Overloads_...`, `ConsensusIsTheColumnMajority_...`, `TrivialInputs_...`, `InvalidInput_...`, `ScoreIsExactForFractionalGaps_...` | protein = string overload = lower-case input; consensus ties; empty / single input; null, '-', digit, positive / NaN / too fine gap, bad iterations / gap model, > 8 sequences; 8.5 -> 9 rounding, cheap gaps beat a mismatch | contract / hand |
| TC4 | `MultipleAlignTCoffee_Protein_EveryRowEqualsTheProgram` | 320 rows of `t_coffee` 13.41 on 160 random protein sets (explicit libraries incl. the program's own BLOSUM library, random trees, relax 1 / 0) reproduced | installed `t_coffee` |
| MX7 (MCP) | `MultipleAlignTests.MultipleAlign_Protein_...`, `MultipleAlignExact_Protein_...`, `MultipleAlign_TypeAndMatrix_AreValidated` | every tool delegates to the matrix overload (default BLOSUM62, PAM250, EDNAFULL IUPAC DNA), `exactTotalScore` / `scoreScale`, type / matrix combinations validated, fractional gaps need a matrix | library call |

---

## 2026-11 audit round 1 (A1-15b part b, B09 F62): ClustalW PROTEIN pairwise stages

`ClustalWPairwiseAlign(a, b, ClustalWProteinPairwiseOptions)` (matrix GONNET / BLOSUM / PAM / ID, `GapOpen` 10, `GapExtend` 0.1), `ClustalWQuickPairwiseAlign(a, b, ClustalWSequenceType.Protein, ClustalWQuickOptions?, int?)`, `SubstitutionMatrix.ClustalWGonnet250 / ClustalWBlosum30 / ClustalWPam350 / ClustalWIdentity`, `MsaOptions.ClustalWProtein` (protein `MatrixScoring` + `MsaDistanceMode.ClustalWAffine | ClustalWQuick`); MCP `pwMatrix`, `pwGapOpen`, `pwGapExtend` (+ protein k-tuple defaults) on `multiple_align_progressive | iterative | consistency`. Class `SequenceAligner_ClustalWProtein_Tests`; oracles `TestData/Alignment/clustalw_protein_pairwise_oracle.tsv`, `clustalw_protein_quick_oracle.tsv`, `msa_matrix_clustalw_oracle.tsv`, `clustalw_protein_matrices_reference.tsv` (`scripts/alignment/generate_clustalw_protein_pairwise_oracle.py`, installed `clustalw` 2.1).

| ID | Test | Expected | Source |
|----|------|----------|--------|
| CP1 | `FullStage_PrintedScoreIdentitiesAndDpScore_MatchClustalW21Oracle` | printed `Score:`, identity count and DP score of 4 204 pairs (4 matrices x 10 gap variants, B/Z/X, low-complexity / repeat pairs) | installed `clustalw -TYPE=PROTEIN -ALIGN` + Python port with `matrices.h` |
| CP2 | `FullStage_AlignedRows_RescoreBelowOrToTheDpScore`, `FullStage_OracleCoversMatricesGapsAndResidues` | the Myers-Miller rows never re-score above the DP optimum, and equal it for gap-free / one-gap paths; oracle covers every matrix, >= 9 gap variants, ambiguity residues, the identity matrix (avg 0) branch and discriminates the matrix | invariant / coverage |
| CP3 | `FullStage_HandRunExample_PrintsSameScores`, `FullStage_Contracts` | 40 / 38 / 31 on three globin fragments; null, J, '*', digit, negative / NaN gap, undefined matrix, empty (distance 1), case, `ProteinSequence` overload, O/U accepted, protein defaults | clustalw binary / contract |
| CP4 | `MatrixTables_EqualTheirCanonicalReferences`, `MatrixTables_AreSymmetricWithTheClustalWAlphabet` | BLOSUM30 = EMBOSS EBLOSUM30 (23 x 23), GONNET250 = round(10 x Biopython GONNET1992) (20 x 20) = 929 cells; all four tables symmetric over `ABCDEFGHIKLMNPQRSTVWXYZ`; anchor values | EMBOSS / Biopython |
| CP5 | `QuickStage_PrintedScore_MatchesClustalW21Oracle`, `QuickStage_DefaultOptions_...`, `QuickStage_OracleCoversTheParameterSpace`, `QuickStage_HandRunExample_...`, `QuickStage_Contracts_AndDnaDelegation` | `Score:` of 3 774 `-QUICKTREE` pairs (k-tuple 1-2, window / topdiags / pairgap / absolute score, 357 "Not Aligned"); parameterless call = defaults 1 / 5 / 5 / 3; 34.0426 / 38.2979 / 26.2295 hand run; protein k-tuple <= 2, DNA type delegates | installed `clustalw -TYPE=PROTEIN -QUICKTREE -ALIGN` |
| CP6 | `ClustalWAffineDistance_ProteinMsaRows_EqualTheIndependentPythonImplementation`, `ClustalWAffineDistance_ChangesTheRowsOfManySets`, `QuickAndNeighbourJoining_ProteinMsa_AreValid_...`, `ProteinDistance_WithResiduesClustalWDoesNotRead_...` | progressive / iterative / consistency rows with the ClustalW protein distance (UPGMA) equal the independent Python implementation on 60 sets x 3 methods; the distance changes rows; quick / NJ variants valid; '*' rejected; EDNAFULL uses the DNA stage | independent Python |
| CP7 (MCP) | `MultipleAlignTests.MultipleAlign_Protein_ClustalWDistance_DelegatesWithTheProteinPairwiseOptions` | `pwMatrix` / `pwGapOpen` / `pwGapExtend` / protein k-tuple defaults reach the library; bad `pwMatrix` throws | library call |

## 2026-11 audit round 1 (A1-10c part b, B09 F63): ClustalW PROTEIN profile stage (`MultipleAlignClustalW` for amino acids)

`MultipleAlignClustalW(IEnumerable<ProteinSequence> | IEnumerable<string>, ClustalWProteinOptions, MatrixScoring?)` (matrix series `-MATRIX` GONNET / BLOSUM / PAM / ID chosen per profile pair by percent identity, `-GAPOPEN 10 -GAPEXT 0.2`, `-NOPGAP`, `-NOHGAP`, `-NOVGAP`, `-GAPDIST`, `-HGAPRESIDUES`, `-NEGATIVE`, `-MAXDIV`, `-NOWEIGHTS`, pairwise stage `-PWMATRIX/-PWGAPOPEN/-PWGAPEXT/-QUICKTREE/-KTUPLE/-WINDOW/-TOPDIAGS/-PAIRGAP`), `SubstitutionMatrix.ClustalWBlosum40 / 45 / 80`, `ClustalWPam20 / 60 / 120`, `ClustalWGonnet40 / 80 / 120 / 160 / 350`, MCP `multiple_align_clustalw type=protein`. Fixture `SequenceAligner_ClustalWProteinMultiple_Tests` (13 tests) + `MultipleAlignTests` (MCP, +1).

| ID | Test | Expected | Source |
|----|------|----------|--------|
| PM1 | `MultipleAlignClustalW_Protein_Rows_EqualClustalW21_OnAllOracleSets` | every output row of `clustalw -TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` on 1 000 protein sets (2-25 sequences, <= 365 residues, 20 option variants; two `-NOVGAP` sets with the forced gdist of the uninitialised C++ member) | installed `clustalw` 2.1 |
| PM2 | `..._UninitialisedGdist_OnlyTheDocumentedSetsNeedAnOverride`, `..._UninitialisedGapDistances_AreIgnoredWithoutVariableRegionPenalties` | without any override exactly the sets 368 and 915 differ from the binary, both `-NOVGAP`; the override list is ignored without `-NOVGAP` | `-DEBUG=5` profile dumps of the binary |
| PM3 | `..._OracleCoversEveryOptionAndTheSpecialPaths`, `..._EachMechanism_ChangesRowsOnManyOracleSets` | > 20 oracle sets per option and per `-MATRIX`, 2-sequence / >= 14-sequence / >= 100-residue / B Z X O U / duplicate sets; for every mechanism (residue-specific, hydrophilic, variable-region penalties, gap distance, BLOSUM / PAM / ID series, negative, weights, delayed sequences, fast stage, gap penalties, pairwise matrix) many sets change rows when it is switched | coverage |
| PM4 | `..._HydrophilicResidues_AreTruncatedToNineLetters` | `MAXHYDRESIDUES` = 9: 12 letters = their first 9; empty / lower-case list keeps / equals the default | `UserParameters::setHydResidues` + clustalw.h |
| PM5 | `..._Invariants_...`, `..._ProteinSequenceOverload_EqualsStringOverload_AndLowerCase`, `..._Contracts`, `..._DefaultOptions_AreClustalWProteinDefaults`, `MultipleAlignClustalW_DnaStringOverload_StillSelectsTheDnaStage` | rows spell the inputs, equal length, input order; `ProteinSequence` = string = lower-case; null / J / `*` / empty / negative or NaN penalty / undefined enum / k-tuple 3 throw; defaults GONNET 10 / 0.2, pairwise 10 / 0.1, maxdiv 30, gapdist 4, `GPSNDQEKR`; the DNA overload is unchanged | contract |
| PM6 | `SeriesMatrices_EqualEmbossDataFiles_OnEveryCell`, `SeriesMatrices_GonnetTables_AreSymmetricAndOnTheClustalWAlphabet` | BLOSUM40 / 45 / 80, PAM20 / 60 / 120 / 350 and BLOSUM62 (= `blosum62mt2` / 2) equal EMBOSS EBLOSUM40/45/80/62, EPAM20/60/120/350 on all 8 x 529 cells; GONNET tables symmetric over `ABCDEFGHIKLMNPQRSTVWXYZ`, anchor values of `matrices.h` | EMBOSS data files |
| PM7 (MCP) | `MultipleAlignTests.MultipleAlignClustalW_Protein_DelegatesToTheProteinClustalWPort_AndParsesEveryOption` | default `type=protein` rows = clustalw (set 20); every option reaches the library; bad `type` / `matrix` / `pwMatrix` / `scoreMatrix` / residue J / negative gap throw; DNA defaults unchanged by the protein parameters | library call |


## 2026-11 audit round 1 (A1-10d file part, B09 F64): ClustalW user matrix files

`ClustalWUserMatrix.Parse` / `ClustalWMatrixFile.Parse` (`-MATRIX=file`, `-PWMATRIX=file`, `CLUSTAL_SERIES` files) feeding `ClustalWProteinOptions.UserMatrix / UserMatrixSeries` and `ClustalWProteinPairwiseOptions.UserMatrix`; test class `SequenceAligner_ClustalWUserMatrix_Tests` (15 tests), oracle `TestData/Alignment/clustalw_user_matrix_oracle.tsv` (240 sets, matrix files embedded).

| ID | Test | Expected | Source |
|----|------|----------|--------|
| UM1 | `MultipleAlignClustalW_UserMatrixFiles_Rows_EqualClustalW21_OnAllOracleSets` | every output row of `clustalw -TYPE=PROTEIN -ALIGN -MATRIX=file / -PWMATRIX=file` on 240 sets (single, series, pairwise, both; integer / float / unlabelled / tab / subset / shuffled / junk-triangle / quirk matrices; option variants) | clustalw 2.1 binary |
| UM2 | `ClustalWPairwiseAlign_UserMatrixFile_Scores_EqualClustalW21_OnAllOraclePairs` | the printed `Sequences (i:j) Aligned. Score:` of > 2 000 pairs (full stage) | clustalw 2.1 binary |
| UM3 | `UserMatrixOracle_CoversEveryFormatAndPath`, `UserMatrixFiles_ChangeRowsAgainstBuiltInMatrices_OnManyOracleSets`, `SeriesAndSingleMatrix_AreDistinctPaths_OnTheOracle` | > 60 series / > 30 single / > 60 `-PWMATRIX` sets, float / tab / comment layouts, quick and negative variants; user matrices change the rows of > 20 % of the sets; series and first-entry-only differ | oracle design |
| UM4 | `Parse_NcbiLayout_GivesTheCanonicalMatrix`, `Parse_LabelledAndUnlabelledRows_AreTheSameMatrix`, `Parse_FloatingPointRows_AreScaledByTenAndTruncated_PerRow`, `Parse_OnlyTheLowerTriangleIsRead`, `Parse_CommentsBlankLinesAndTrailer_AreSkipped` | BLOSUM62 file = `SubstitutionMatrix.Blosum62` on 23 x 23 cells; per-row x10 and truncation (25, -5, 30, 1, 0, 4); upper triangle ignored; `#` / `!` handling | `SubMatrix::readUserMatrix` |
| UM5 | `Parse_MalformedFiles_AreRejectedLikeClustalW`, `Parse_UnknownHeaderLetters_AreAcceptedWithAShiftedTable_AsInClustalW`, `MatrixFile_Series_ParsesEntriesAndFollowsClustalWsRules`, `MatrixFile_Series_ReadsAtMostTenMatrices` | no `*` row / header only / series where one matrix is expected / short row -> `FormatException` (the messages of the installed binary); J header = `Cross -1`; limits, ignored lines, 10-entry cap | binary + source |
| UM6 | `FromSubstitutionMatrix_EqualsTheParsedFile_AndDrivesTheAlignmentIdentically` | a canonical `SubstitutionMatrix` = the equivalent file | API |
| UM7 (MCP) | `MultipleAlignTests.MultipleAlignClustalW_Protein_UserMatrixFiles_EqualClustalW_AndAreValidated` | `matrixFile` + `matrixFiles` (series) = clustalw rows of oracle set 11; `matrixFile` / `pwMatrixFile` single; missing series file / no header / series as `pwMatrixFile` throw | oracle + library call |

## 2026-11 audit round 1 (A1-10d part 2, B09 F67): ClustalW end-gap options of the profile stage + the program's abort

`ClustalWOptions` / `ClustalWProteinOptions` `UseEndGaps` (`UserParameters::useEndGaps`) and `EndGapPenalties` (`endGapPenalties`), and `ClustalWAbortedException` (first-pass profile score < 0 -> ClustalW writes no output); test class `SequenceAligner_ClustalWEndGaps_Tests` (9 tests), oracle `TestData/Alignment/clustalw_end_gaps_oracle.tsv` (340 sets x 4 modes = 1 360 rows; 68 of them `ABORT`), produced by the user-run Docker package `docs/Validation/review-2026-09/handoff/A1-10d` (`clustalw` 2.1+lgpl-7 with one edit: the `-ENDGAPS` handler reads `CW_USEENDGAPS` / `CW_ENDGAPPEN`; control runs = installed binary on 340 / 340 sets).

| ID | Test | Expected | Source |
|----|------|----------|--------|
| EG1 | `MultipleAlignClustalW_EndGapModes_RowsEqualClustalW_OnAll340Sets` (control, useendgaps) | every output row of the 340 sets (120 DNA + 220 protein, option variants incl. `-QUICKTREE`) | patched clustalw 2.1 |
| EG2 | `MultipleAlignClustalW_EndGapPenaltyModes_RowsEqualClustalW_OnAll306SetsAndAbortOn34` (endgappen, both) | 306 sets: every row; the 34 sets (31 DNA + 3 protein) without an output file raise `ClustalWAbortedException` | patched clustalw 2.1 (MSA.cpp 150-155, Clustal.cpp 205-208) |
| EG3 | `Abort_IsRaisedExactlyForTheProgramsNoOutputSets_AndNeverWithoutEndGapPenalties`, `Abort_ScoresOfAbortingSetsAreNegativeIntegers` | `Score` <= -1 and `Group` >= 1 in the exception and its message; the same sets align without end-gap penalties | source + oracle |
| EG4 | `OracleHasEveryModeForEverySet`, `OracleDiscriminatesTheOptions`, `ControlMode_EqualsTheDefaultOptions` | 340 x 4 rows, 68 aborts; `useendgaps` changes > 100 sets, `endgappen` / `both` > 150; defaults = control | oracle design |
| EG5 (MCP) | `MultipleAlignTests.MultipleAlignClustalW_EndGapOptions_DelegateAndAbortIsExplicit` | `useEndGaps` / `endGapPenalties` = library rows of oracle sets; the aborting set fails with the 'abandons the alignment' error | oracle + library call |

## 2026-11 audit round 1 (A1-10d part 3, B09 F68): ClustalW nucleotide matrices (`-DNAMATRIX`, `-PWDNAMATRIX`, `-TRANSWEIGHT`)

`ClustalWOptions.DnaMatrix / PairwiseDnaMatrix / UserDnaMatrix / PairwiseUserDnaMatrix / UserDnaMatrixUse`, `ClustalWDnaPairwiseOptions`, `SequenceAligner.ClustalWPairwiseAlign(.., ClustalWDnaPairwiseOptions)`; test class `SequenceAligner_ClustalWDnaMatrix_Tests` (9 tests), oracle `TestData/Alignment/clustalw_dna_matrix_oracle.tsv` (360 sets: 125 command-line runs, 235 interactive-menu runs, matrix files embedded).

| ID | Test | Expected | Source |
|----|------|----------|--------|
| DM1 | `MultipleAlignClustalW_DnaMatrices_Rows_EqualClustalW21_OnAllOracleSets` | every output row of 345 sets (IUB / CLUSTALW built-ins, command-line files = ignored, menu files = applied, 13-17 letter headers); the short-header (undefined behaviour) sets: all but <= 2 equal | clustalw 2.1 binary |
| DM2 | `MultipleAlignClustalW_DnaMatrices_GroupScores_EqualClustalW21_OnAllOracleSets` | the 1 393 printed first-pass `Group g: Sequences: n Score: s` / `Delayed` lines (pin the profile matrix and the 0.66 gap scale) | clustalw 2.1 binary |
| DM3 | `ClustalWPairwiseAlign_DnaMatrices_Scores_EqualClustalW21_OnAllOraclePairs` | 3 485 printed pairwise `Score:` (full stage; IUB / CLUSTALW / user pairwise matrices, `-TRANSWEIGHT`, gap variants) | clustalw 2.1 binary |
| DM4 | `MalformedMatrixFiles_AreRejectedLikeClustalW` | the 15 files that end clustalw (no `*` row, short row, series file, empty, too few rows) -> `FormatException` | clustalw 2.1 binary + F64 reader |
| DM5 | `Oracle_CoversEveryModeMatrixAndLayout`, `Oracle_DiscriminatesTheMatrixOptions`, `CommandLineUse_IgnoresTheFile_AppliedUses_It` | >= 330 sets / > 100 command-line (125) / > 150 menu (235); every option > 20 sets; tab / float / comment / short-header layouts; user matrices change the rows of > 1/3 and the group scores of > 2/3 of the menu sets, command-line files change nothing (all 75), CLUSTALW changes > 1/4 | oracle design |
| DM6 | `UserMatrixEqualToIub_BehavesLikeTheBuiltInMatrix`, `PairwiseOptions_Defaults_AndContracts` | a file holding `swgapdnamt` in packing order = the built-in IUB (profile and pairwise); defaults 15 / 6.66 / 0.5 / IUB / Applied; null / negative / NaN / unknown-enum contracts | API |
| DM7 (MCP) | `MultipleAlignTests.MultipleAlignClustalW_DnaMatrices_DelegateAndFollowClustalWsCommandLineAndMenu` | `dnaMatrixFile` = clustalw menu rows of oracle set 244 (applied), `commandline` = rows without the file, `dnaMatrix` / `pwDnaMatrix` = `clustalw`; bad names / series file / malformed file throw `ArgumentException` | oracle + library call |

## 2026-11 audit round 1 (A1-10c part a, B09 F65): MUSCLE 3.8.1551 progressive alignment

`MultipleAlignMuscle(...)` (`-maxiters 1 | 2`) with `MuscleOptions`; test class `SequenceAligner_Muscle_Tests` (20 tests). **The row-for-row oracle (the MUSCLE 3.8.1551 binary, user-run Docker package `handoff/A1-10c-muscle`) was returned and is equal on 420 / 420 sets: `SequenceAligner_MuscleOracle_Tests` (below).** Locked besides it:

| ID | Test | Expected | Source |
|----|------|----------|--------|
| MU1 | `LogExpectationScore_OfSingleResidueColumns_IsLnOfTheVtml200Entry_MinusCenter`, `LogExpectationScore_IsWeightedByOccupancy` | `ln(VTML_LA[a][b]) + 0.52` (A-A 2.25080, C-C 15.79469, A-C 1.31180); occupancy 0.5 halves the score | Edgar 2004 + `scorepp.cpp` / `params.cpp` |
| MU2 | `NucleotideScore_IsTheBlastzMatrixPlusCenter`, `WildcardCounts_FollowMsa2Cpp_IncludingTheNucleotideAxRQuirk` | A-A 91 + 60 ... ; B -> D + N, Z -> E + Q, X -> 1/20; DNA X -> G + A, N / R -> 0.05 x 4 | `nucmx.cpp`, `msa2.cpp` |
| MU3 | `KimuraDistance_FollowsTheFormulaAndTheDayhoffTable`, `Upgma_OnTheSourceTestMatrix_BuildsTheExpectedTree`, `Upgmb_UsesTheBiasedLinkage`, `ClustalWWeights_OnTheSourceExample`, `KmerDistance_IsZeroForIdenticalAndThreeWhenNoSextetIsShared` | `-ln(1 - p - p²/5)`, table entries 195 / 250 / 988, `10` above 93 %; `((0,1),(2,3)),4)` with edge lengths 1 1 2 2 6 5 3 10; 0.1 avg + 0.9 min; weights 0.2 0.2 1.0 1.1 / 2.5; k-mer distance 0 / 3 | source comments (`msadistkimura.cpp`, `upgma2.cpp`, `clwwt.cpp`) |
| MU4 | `NwSmall_FindsAnOptimalPath_ComparedWithExhaustiveEnumeration`, `NwSmall_AlignsIdenticalProfilesWithoutGaps_AndPlacesASingleDeletionAsOneGapRun` | the DP path is a best alignment under the position-specific open / close / extension model on 480 random profile pairs (enumeration of all paths) | `nwsmall.cpp` semantics |
| MU5 | `MultipleAlignMuscle_Rows_SpellTheInputs_HaveEqualLength_AndAreDeterministic`, `..._ProteinAndDnaOverloads_EqualTheStringOverload`, `..._TwoSequences_Stop...`, `..._MethionineHack_LetterRules_AndTreeOrder`, `..._Iteration2_RealignsOnlyChangedSubtrees_ButKeepsEveryRowValid`, `..._Options_AreValidated`, `..._ScoreModels_AndGapParameters_ProduceValidAlignments` | invariants, overloads, two sequences, M hack, letter rules, iteration 2 changes some drafts and `MaxTrees 0` = draft, option validation (`MaxIterations 3` -> `NotSupportedException`), LE / SP / SV / RNA sets | API contract |
| MU6 (MCP) | `MultipleAlignTests.MultipleAlignMuscle_DelegatesToThePort_ParsesEveryOption_AndValidates` | tool = library for every option, `auto` type, errors | library call |

### Addendum (audit A1-10c a + a2, B09 F69 / F70)

| ID | Test | Expected | Source |
|----|------|----------|--------|
| MO1 | `SequenceAligner_MuscleOracle_Tests.OracleHas420Sets_AllRun`, `Rows_EqualMuscle3` x420, `OutputOrder_EqualsMuscle3` x420 | rows (input order) and the order MUSCLE writes equal the MUSCLE v3.8.1551 binary on 420 sets (270 protein / 110 DNA / 40 RNA; `-maxiters 1|2`, `-le/-sp/-sv`, gap / center parameters, `-distance1`, `-cluster1/-cluster2`, `-sueff`, `-maxtrees 0`) | user-run Docker oracle `muscle3_oracle.tsv` |
| MO2 | `SequenceAligner_MuscleOptions_Tests`: `HenikoffWeights_...`, `HenikoffPbWeights_...`, `NoneWeights_...`, `GscWeights_...`, `ThreeWayWeights_...` | weights 0.375 / 0.25 / 0.375 (Henikoff, Henikoff-PB), uniform, GSC 0.25 / 0.25 / 0.5, three-way 0.2 / 0.2 / 0.6 (hand traces of `henikoffweight*.cpp`, `setgscweights.cpp`, `threewaywt.cpp`) | source; the program's weights verified on the round-2 oracle (MO10) |
| MO3 | `ObjScoreSp_...`, `ObjScoreSp_TerminalGaps_...`, `ObjScore_XpSpfAndSp_AreConsistent`, `FindAnchorCols_...` | SP objective 3.25 / 0.75, `-termgaps half / full / ext`, XP / SPF = SP, anchors 4 / 12 / 20 / 28 / 36 | `objscore2.cpp`, `spfast.cpp`, `anchors.cpp` |
| MO4 | `DistKmer20_4_...`, `DistKmer20_3_...`, `DistKbit20_3_CountsTheUnionOfTheBitVectorsPlusTheSecondVector`, `ScoreDist_...` | 0.25 / 0.99, 0.8 / 1, 8/9 / 1 (union of the bit vectors + the bits of the second vector again, as `CommonBitCount` does: verified by the round-2 oracle), 0 for identical | `fastdistkmer.cpp`, `fastdistjones.cpp`, `fastdistkbit.cpp`, `scoredist.cpp` |
| MO5 | `NeighborJoining_ReproducesAnAdditiveMatrix`, `FixRoot_PreservesPathLengths_AndRootsAsDocumented`, `DiffTrees_FindsTheMarriedSubtreesOfTheSourceExample` | path lengths = input matrix; root at the midpoint of the longest span; the program's min-average relation; diffs k, C, D | `clust.cpp`, `phy3.cpp`, `difftrees.cpp` |
| MO6 | `FormatG_...`, `MuscleRun_WritesTheTreeAndDistanceFilesInTheProgramFormat`, `MuscleRun_ComputeWeights_...`, `ParseUserTree_...`, `MuscleRun_UserTree_...`, `ReadMatrixFile_...`, `MuscleRun_UserMatrix_...` | `%g`, file formats, node numbering of the Newick reader, `-usetree` skips iteration 2, `-matrix` incl. the lost last column | `phytofile.cpp`, `phyfromfile.cpp`, `readmx.cpp` |
| MO7 | `Diags_Hydro_DoubleAffine_AndResourceOptions_HaveNoEffect`, `Iterations3AndAbove_*`, `ClusteringDistanceAndRootOptions_*`, `DistanceOfTheWrongKind_*`, `Brenner_*` | rows unchanged by the compiled-out options; valid alignments for 72 objective x weighting x anchors and 17 clustering / distance / root combinations | `muscle.h` (`VER_3_52 = 0`), `refinehoriz.cpp` |
| MO8 | `ProfileAlign_*`, `Refine_*`, `ScoreAlignment_*`, `NucleotideOnlyModes_*`, `LocalAlign_*`, `RefineWindows_*` | the MSA command lines; protein `-ppscore/-sw/-refinew` rejected (the program crashes) | `profile.cpp`, `refine.cpp`, `dosp.cpp`, `ppscore.cpp`, `local.cpp`, `refinew.cpp` |
| MO9 (MCP) | `MultipleAlignTests.MultipleAlignMuscle_NewOptions_AreParsedAndDelegated`, `MuscleMsaCommand_RunsEveryMode_...` | tool = library, every option parsed, errors | library call |

| MO10 | `SequenceAligner_MuscleOracleRound2_Tests.OracleHas630Sets_AllKindsPresent`, `Align_RowsWrittenOrderAndFiles_EqualMuscle3` x441, `CommandLine_EqualsMuscle3` x189 | rows in MUSCLE's written order, names, `-tree1/-tree2/-distmx1/-distmx2/-computeweights` files, `SP=%.4g`, `Score=%.4g`, the `-maketree` file and the `-profile/-refine/-refinew/-sw` rows equal the MUSCLE v3.8.1551 binary on 630 recorded sets (every option of MO2-MO8, `-maxiters` 3-16 and the default 8); the 23 runs the program refuses or crashes on are refused by the port with the matching exception (crash / `-stable` -> `NotSupportedException`, `Invalid DISTANCE_n` -> `ArgumentException`, `Clust::SetLeafCount(1)` / `g_MuscleWeights = 0` -> `InvalidOperationException`) | user-run Docker oracle `muscle3_r2_oracle.tsv` (`handoff/A1-10c-muscle-r2`, 638 sets, 8 of them recorded in round 3) |
| MO11 | `SequenceAligner_MuscleOracleRound3_Tests` (`Oracle_Has60Sets_TheIdsOfThePackage`, `TheEightRoundTwoGapIdsAreNowInTheRound3Oracle`, `Align_RowsWrittenOrderAndFiles_EqualMuscle3` x37, `CommandLine_EqualsMuscle3` x23) + `SequenceAligner_MuscleOptions_Tests.MuscleRun_UserMatrix_OnNucleotides_...` | the 8 `-matrix blosum62nostar` sets, `-matrix` on DNA / RNA (nucleotide alphabet then amino-acid SP run, MO11 R638-R645), one- / two-sequence profiles / trees / scores / refinements / alignments (the 10 crashes / errors are refused), `-le -sp -sv -spn` on both alphabets: equal on 60 / 60 sets of the muscle3 3.8.1551 binary (a crash of the program on a nucleotide `-matrix` run must come with `MuscleRunResult.ReferenceBehaviourUndefined`) | local run of the Debian `muscle3` package: `muscle3_r3_oracle.tsv` (`handoff/A1-10c-muscle-r3`, `scripts/alignment/generate_muscle3_r3_oracle.py`); fresh sample (uncommitted comparison) `scripts/alignment/muscle3_fresh_sample.py`: 1 800 / 1 800 |
| MO12 | `SequenceAligner_ConsistencyOptions_Tests` (5), MCP `multiple_align_consistency` | `MsaOptions.PrimaryLocalAlignments` (top-N Lalign alignments at their own identity, signal-added; 1 = the old library exactly) and `MsaOptions.SequenceWeights = TCoffee` (extended weight x min(w_i, w_j); mean weight 1; both engines) | independent re-derivation from `LocalAlignTopN` (EMBOSS matcher oracle, F47) and the program path's weight function (t_coffee oracle, F66); the paper engine as a whole has no t_coffee equivalent (audit A2-16, F93) |
| MO13 | `SequenceAligner_ClustalWTypeGuess_Tests` (500 oracle sets + hand cases), `SequenceAligner_ClustalWUpgma_Tests` (5), MCP `ClustalWAutoTypeUpgmaTests` (2) | `GuessClustalWSequenceType` (second sequence decides, >= 85 % A/C/G/T/U/N) equals the unit printed by the installed clustalw 2.1 on 500 sets; `-CLUSTERING=UPGMA` source port locked to hand-derived joins / heights / float weights (NO row oracle: the installed clustalw segfaults on every UPGMA run; Docker package `handoff/A2-14-clustalw-upgma`) | clustalw 2.1 `clustalw_type_oracle.tsv`; source of `tree/UPGMA` (audit A2-14, F92) |

**The options of MO2-MO8 are verified against the binary (MO10: 630 / 630 recorded sets equal after the fixes F76).** Round 3 (MO11, audit A2-1 / F89): 60 / 60 equal, `-matrix` for nucleotides implemented; a fresh sample of 1 800 sets is 1 800 / 1 800 equal.

## Addendum (audit A1-13b, B09 F66): T-Coffee 13.41 program pipeline

Class `SequenceAligner_TCoffeeProgram_Tests` (13 tests), oracle = installed `t_coffee` 13.41 (`scripts/alignment/generate_tcoffee_program_oracle.py`; TestData `tcoffee_program_library_oracle.tsv` 820 sets, `tcoffee_program_msa_oracle.tsv` 1 480 sets, `tcoffee_program_big_oracle.tsv` 2 sets of 100 / 101 sequences).

| ID | Test | Scenario | Expected | Source |
|----|------|----------|----------|--------|
| TP01 | `OracleFiles_AreComplete` | counts per method / type / sequence number | >= 60 libraries per (type, method), >= 200 alignment sets per headline profile | oracle |
| TP02 | `PrimaryLibrary_EveryPairEqualsTheProgram` | `-lib_only -out_lib` of slow_pair, lalign_id_pair, proba_pair and merges, DNA + protein | every pair, weight, count equal | `t_coffee` 13.41 |
| TP03 | `PrimaryLibrary_TwoMethods_IsNotTheUnion_WorkedExample` | GGATCACAG / CTACACTG, `-lalign_n_top 5` | the program's 9 / 13 / 10 entries (not the union of 14) | `t_coffee` 13.41 |
| TP04 | `SlowPair_MyersMiller_PlacesTheGapInside_WorkedExample` | TTTCAACTC vs TTTTCAACTC | `T-TTCAACTC` (library pairs 1-1, 2-3, ...) | `t_coffee` 13.41 |
| TP05 | `LalignTopAlignments_AreNonIntersecting_AndInScoreOrder` | TTTTCGAACTTC vs TTTTCAACTC, 4 alignments | scores 790 / 500 / 400 / 370, first alignment rows, no pair reported twice | `t_coffee` 13.41, Waterman-Eggert definition |
| TP06 | `GuideTree_EqualsTheProgram_OnEveryOracleSet` | trees of all sets of >= 3 sequences | equal | `t_coffee` 13.41 `Group` log |
| TP07 | `MultipleAlignTCoffeeProgram_EveryRowEqualsTheProgram` | 1 480 sets | every row equal | `t_coffee` 13.41 |
| TP08 | `..._HundredSequences_FastNeighborJoining_EqualsTheProgram` | 100 / 101 sequences | tree and rows equal | `t_coffee` 13.41 |
| TP09 | `..._RowsSpellTheInputs_AndHaveEqualLength` | DNA example, plain + slow,lalign | rows spell the inputs | invariant |
| TP10 | `..._Example_EqualsTheProgram` | 18-nt example, `slow_pair,lalign_id_pair` | the program's four rows | `t_coffee` 13.41 |
| TP11 | `..._DegenerateInputs`, `..._Contract_Errors` | none / one / two sequences, empty / gapped / RNA input, bad `lalign_n_top` | defined results / exceptions | contract |
| TP12 | `TCoffeeKtupSimilarity_IdenticalSequencesAre100_DisjointTuplesAre70` | `ktup_dist_mat` formula | 100 / 70 | source |


## Addendum (audit A1-13c / A1-13d, B09 F71-F75): command-line options of the T-Coffee 13.41 program

Oracle: `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/tcoffee_program_options_oracle.tsv` (`scripts/alignment/generate_tcoffee_options_oracle.py`, the installed `t_coffee` 13.41, `-n_core 1 -no_warning`): `feature`, `type`, the flag string, the sequences, the rows (input order; the output-file order for `-outorder aligned`) and the sequence names. `TCoffeeProgramOptions.FromCommandLine(flags)` (+ `Names`) is run on every row of the file.

| ID | Test (`SequenceAligner_TCoffeeOptions_Tests`) | Sets | Expected | Oracle |
|----|-----------------------------------------------|------|----------|--------|
| OP01 | `Rna_EveryRowEqualsTheProgram` | 430 RNA sets (U), 5 method mixes | every row equal (A1-13c: RNA pair-HMM, `get_tot_prob`) | `t_coffee` 13.41 |
| OP02 | `RelaxAndFilterLib_EveryRowEqualsTheProgram` | 500 | `-relax_lib 0 / 2 / 3`, `-filter_lib 20 / 100` | `t_coffee` 13.41 |
| OP03 | `GapPenaltiesAndNoMatch_EveryRowEqualsTheProgram` | 560 | `-gapopen`, `-gapext`, `-cosmetic_penalty`, `-nomatch`, `-do_normalise` | `t_coffee` 13.41 |
| OP04 | `TreeModeAndDistanceMatrixMode_EveryRowEqualsTheProgram` | 730 | `-tree_mode upgma`, `-distance_matrix_mode fast / idscore / no`, `-distance_matrix_sim_mode` | `t_coffee` 13.41 |
| OP05 | `SequenceWeights_EveryRowEqualsTheProgram` | 470 | `-seq_weight t_coffee`, weight files, `-weight no` | `t_coffee` 13.41 |
| OP06 | `NameOrder_EveryRowEqualsTheProgram` | 560 | `-inorder aligned / input`, `-outorder aligned`, shuffled / prefix / s0..s12 names | `t_coffee` 13.41 |
| OP07 | `Matrix_EveryRowEqualsTheProgram` | 675 | `-matrix` for DNA / RNA / protein, `proba_pair` and `slow_pair` | `t_coffee` 13.41 |
| OP08 | `MethodParameterOverrides_EveryRowEqualsTheProgram` | 390 | `slow_pair@EP@GOP@..`, `lalign_id_pair@EP@MATRIX@..`; the 3 sets on which `t_coffee` core-dumps (SIM edit-script overflow, `lalign_id_pair@EP@GOP@-5`) raise `TCoffeeAbortedException` | `t_coffee` 13.41 |
| OP09 | `ClustalwPair_EveryRowEqualsTheProgram` | 400 | `clustalw_pair` through the ClustalW 2.1 ports | `t_coffee` 13.41 + `clustalw` 2.1 |
| OP10 | `FromCommandLine_ParsesMethodsAndRejectsUnknownFlags`, `..._MatrixAndMethodParameters` | - | every supported flag parsed; unknown flags / values rejected by name | contract |
| OP11 | `ClustalwPair_MixedAlphabets_EveryRowEqualsTheProgram` | 230 | `clustalw_pair` on sets from different alphabets: the external clustalw guesses the type of every pair file itself (the second sequence), a protein letter in DNA mode scores 0 (A2-7, F107) | `t_coffee` 13.41 + `clustalw` 2.1 |
| OP12 | `NoMatchOnLargeSets_FloatSumOrder_EveryRowEqualsTheProgram` | 260 | 6-10 protein sequences with `-nomatch 40 / 52 -gapext -1`, `-seq_weight t_coffee`: the single-precision pair list sums in ascending sequence-index order (A2-7, F107) | `t_coffee` 13.41 |
| OP13 | `BiphasicPair_EveryRowEqualsTheProgram` | 390 | `biphasic_pair` (DNA / protein / RNA, with `slow_pair`, `@EP@`, `proba_pair` + upgma; A2-3, F108) | `t_coffee` 13.41 |
| OP14 | `LalignRsSPair_EveryRowEqualsTheProgram` | 390 | `lalign_rs_s_pair` (SIM on every pair incl. each sequence against itself; A2-4, F108) | `t_coffee` 13.41 |
| OP15 | `PromoAndBlastrPair_EveryRowEqualsTheProgram` | 405 | `promo_pair`, `blastr_pair` (dinucleotide alphabet, lossy shared sequence array, weight `extend_sim`; A2-4, F108) | `t_coffee` 13.41 |
| OP16 | `RefusedMethods_RaiseTCoffeeAbortedException_AsTheProgramFails` | 160 | `clean_slow_pair`, `hash_pair`, `subop1_pair`, `subop2_pair`, `cfast_pair`: every set is a program failure (exit status 1) and raises `TCoffeeAbortedException` (`Reason`, `ExitStatus`); one sequence: no refusal (A2-3 / A2-4, F108) | `t_coffee` 13.41 |
| OP17 | `ExtendAndUnextendSeq_FollowReformatC` | hand-derived | `extend_seq` / `unextend_seq` of `reformat.c`: ac -> e, last residue kept, n / x, u, gaps, protein letters -> j | source |
| OP18 | `DistanceSimModeIdscore_GapsInTheFirstRowSegfault_OtherwiseEveryRowEqualsTheProgram` | 260 + 280 | `-distance_matrix_sim_mode idscore*` with `-distance_matrix_mode fast`: a gap in the FIRST aligned row = segmentation fault (exit 11; heap dependent after the first 12 pair alignments, 5 sets survived), a gap in the second row reads 0 (23 of 24 sets; 1 heap-dependent difference); all other rows equal (A2-6, F109) | `t_coffee` 13.41 |
| OP19 | `DistanceSimModeGroups_EveryRowEqualsTheProgram` | 265 | the residue groups of `make_group_aa` (simple, mafft, clustalw, polarity, vasiliky, make_all, unknown name) selected by the sim mode (A2-6, F109) | `t_coffee` 13.41 |
| OP20 | `TreeAndDistanceModesTheProgramRefuses_RaiseTCoffeeAbortedException` | 360 | `-tree_mode cwph / kmeans / unknown` (3+ sequences: exit 11 / 0 without alignment / 11; two sequences ignore the mode), `-distance_matrix_mode aln` (exit 11), `ktup2`, unknown (exit 1): every failure row raises `TCoffeeAbortedException`; `random` is not reproducible (`NotSupportedException`) (A2-6, F109) | `t_coffee` 13.41 |
| OP21 | `FastPair_EveryRowEqualsTheProgram` | 510 | `fast_pair` = `fasta_gotoh_pair_wise` (k-tuple diagonals, z-score threshold, banded Gotoh DP); `-nomatch` / `-ndiag` do not reach library methods (A2-2, F110) | `t_coffee` 13.41 |
| OP22 | `DistanceModesOfTheRelaxedLibrary_EveryRowEqualsTheProgram` | 430 + 356 + 340 + 470 | `-distance_matrix_mode slow / cscore / clustalw / very_fast` and `-quicktree`, computed between relaxation and re-weighting (A2-6, F111) | `t_coffee` 13.41 |
| OP23 | `GeometricDistance_EveryRowEqualsTheProgram` | 310 | `-distance_matrix_mode geometric` (fasta_pair_wise, pam250mt, identity vectors against the first three sequences; A2-6, F111) | `t_coffee` 13.41 |
| OP24 | `ExonPairs_EveryRowEqualsTheProgram` | 420 + 220 | `exon_pair` / `exon2_pair` / `exon3_pair` with `_G_` templates (exon2mt, subset letters, process-wide residue groups fixed by `fix_seq_seq`; 87 empty-library failures raise; A2-4, F112) + the library of `-lib_only -out_lib` of a hand-made pair | `t_coffee` 13.41 |
| OP25 | `GenericDpModes_EveryRowEqualsTheProgram` | 800 + 480 + 790 | `-dp_mode gotoh_pair_wise / ssec / gotoh_pair_wise_test / myers_miller_pair_wise / default / fasta / cfasta / very_fast_pair_wise` with `-extend_mode` very_fast / slow / triplet / fast_triplet, `-max_n_pair`, `-tg_mode`, `-ndiag`, gap penalties (A2-5, F113); 2 sets are the program's heap-overflow aborts (exit 6) of gotoh_pair_wise, checked against the necessary condition alignment length = 8 (mod 16) | `t_coffee` 13.41 |
| OP25b | `BiphasicAndGlocalDp_EveryRowEqualsTheProgram` | 550 | `-dp_mode biphasic_pair_wise` / `glocal_pair_wise` on the column scores (A2-5, F113); 6 sets are the row-writer overflow of the program | `t_coffee` 13.41 |
| OP25c | `GotohLgpDp_EveryRowEqualsTheProgram` | 620 | `-dp_mode gotoh_pair_wise_lgp` (position specific gap opening / closing penalties of `aln2local_penalties`) with `-gapopen` / `-gapext`, `-nomatch`, `-extend_mode`, weights, UPGMA (A2-5, F115); 2 sets are the row-writer overflow | `t_coffee` 13.41 |
| OP25d | `KtupAndPrecomputedDp_EveryRowEqualsTheProgram` | 250 + 220 + 120 | `ktup_pair_wise` (rows unchanged), `precomputed_pair_wise` (equal-length sets: rows unchanged; 120 unequal sets stop with exit status 1) (A2-5, F115) | `t_coffee` 13.41 |
| OP25e | `LocalAndTestDp_EveryRowEqualsTheProgram` | 500 + 420 | `gotoh_pair_wise_sw` (best local columns only; 50 undefined behaviour) and `test_pair_wise` (98 null allocations raise, 50 equal, 272 row overflows) (A2-5, F115) | `t_coffee` 13.41 |
| OP25f | `ProcoffeeDp_EveryRowEqualsTheProgram` | 460 | `-dp_mode procoffee_pair_wise` (sparse linked DP, `cl2diag_cap390` over the list with the borders, groups swapped) (A2-5, F115) | `t_coffee` 13.41 |
| OP25g | `ViterbiDp_EveryRowEqualsTheProgram` | 420 | `-dp_mode viterbi_pair_wise` (Viterbi path of the protein pair-HMM on the frequency profiles); 1 set = row-writer overflow (A2-5, F115) | `t_coffee` 13.41 |
| OP26b | `PcMixtAndGCoffeeExtendModes_EveryRowEqualsTheProgram` | 760 | `-extend_mode pc / mixt / g_coffee / g_coffee_quadruplets` on gotoh / myers_miller / fasta dp modes; 6 sets = row-writer overflow (A2-5, F115) | `t_coffee` 13.41 |
| OP26d | `QuadrupletExtendMode_EveryRowEqualsTheProgram` | 360 | `-extend_mode quadruplet` (`get_dp_cost_quadruplet`, `fast_get_dp_cost_quadruplet`, nomatch twice); 1 set = row-writer overflow (A2-5, F115) | `t_coffee` 13.41 |
| OP26c | `TestAndRnaExtendModes_EveryRowEqualsTheProgram` | 470 | `-extend_mode test` (raw sum x 10 / (ns1 x ns2)) and `rna` / `rna1-4` without a contact library (= rna0); 3 sets = row-writer overflow (A2-5, F115) | `t_coffee` 13.41 |
| OP29 | `MusclePair_EveryRowEqualsTheProgram` | 460 | `-method muscle_pair` (MUSCLE 3.8.1551 on every pair, `-maxiters 8`) with `slow_pair` / `proba_pair`, UPGMA, no relaxation, mixed alphabets (A2-11, F114) | `t_coffee` 13.41 + `muscle3` |
| OP30 | `MafftPair_EveryRowEqualsTheProgram` (was `MafftPair_OracleCommitted_PortWaitsForMafft`) | 460 | `-method mafft_pair` / `mafft_mixed` rows equal the program (the pair stage of MAFFT is ported, A2-19 / F119) | oracle |
| OP31 | `TMalignPair_EveryRowEqualsTheProgram` / `MustangPair_EveryRowEqualsTheProgram` | 160 + 160 | `TMalign_pair` / `mustang_pair` on synthetic CA-only PDB files: `fill_P_template` acceptance, the glue of `tc_generic_method.pl` (TMalign `tail -4`, mustang exit status 0), `proba_pair` fallback for pairs without two structures (A2-10, F116) | `t_coffee` 13.41 + `TMalign` + `mustang` |
| OP32 | `SapPair_IsTMalignWhenSapIsNotInstalled`, `DaliAndFuguePair_TheProgramStopsBecauseTheAlignersAreNotInstalled` | 150 + 8 + 60 + 60 | `sap_pair` = TM-align (SAP absent), `-mode 3dcoffee` 7 / 8 = `proba_pair,sap_pair` (the mode raises), `dali_pair` stops at the first pair of two structures, `fugue_pair` stops (A2-12, F116) | `t_coffee` 13.41 |
| OP33 | `TMalignAlignment_IsTheFirstAndThirdOfTheLastFourLines_EachWithAFinalGap`, `MustangAlignment_IsTheFirstTwoFastaRecords` | - | the glue parsers (contract) | script `tc_generic_method.pl` |
| OP34b | `TblastxMsa_EveryRowAndLibraryEqualsTheProgram` | 106 | `tblastx_msa`: the library from the HSPs of the local all-against-all `tblastx` XML (every nucleotide of the query range, reverse frames mirrored, hits of later sequences only); 5 sets without a hit stop the program (A2-9, F117) | `t_coffee` 13.41 + `blastall` |
| OP34 | `BlastpMsa_EveryRowEqualsTheProgram`, `BlastpMsa_LibraryEqualsTheProgram` | 170 + 85 | `blastp_msa`: the library from the HSPs of the local all-against-all BLAST XML (weights = truncated HSP identity, self hits, CONS 0); rows and `-lib_only -out_lib` libraries (A2-9, F117) | `t_coffee` 13.41 + `blastall` |
| OP26 | `ExtendModes_EveryRowEqualsTheProgram_RelativeTripletIsUndefinedBehaviour` | 340 | `-extend_mode rna0` (= slow_triplet), `pc` / `mixt` / `relative_triplet` on the linked DP (no effect); `relative_triplet` with a generic dp mode reads past a library list (UB: `NotSupportedException`, or SIGFPE exit 8 = `TCoffeeAbortedException`; equal on every set without such a read) | `t_coffee` 13.41 |
| OP27 | `DpAndExtendModesTheProgramRefuses_RaiseTCoffeeAbortedException` | 760 | the dp modes / extend modes the program refuses or dies on (exit 1 / 6 / 11), with the one-sequence exceptions | `t_coffee` 13.41 |
| OP28 | `FromCommandLine_FastaOptionsQuicktreeAndTemplates` | - | `-quicktree` anywhere on the line, `-template_file @t0,t1,..` (empty = none), `-ndiag`, `-ktuple`, `-tg_mode`; a template count that does not match the sequences is rejected | contract |

All committed sets of the file (20 011 since WP30: A2-5 added 4 600 dp / extend-mode sets (lgp 620, sw 500, ktup 250, precomputed 340, test 420, procoffee 460, viterbi 420, extend modes pc / mixt / g_coffee* 760, quadruplet 360 and test / rna* 470) + 40 refused, A2-11 920 (muscle_pair 460, mafft_pair oracle 460); the new file `tcoffee_external_oracle.tsv` holds 874; 14 451 since WP29: A2-2 / A2-4 / A2-5 / A2-6 added 6 736; 7 715 before) are compared exactly (A2-7: `KnownDifferences` and the `crashed <= 5` tolerance of the first version were removed; the three former differences were the ascending-index order of the pair-list sums (2 sets) and the sequence type of the clustalw pair file (1 set), F107). The only tolerated deviations are the documented heap-dependent undefined behaviours of the idscore sim mode (<= 6 of 540 sets, `TCoffeeAbortedException.HeapDependent` or a second-row gap) - an abort of the program is asserted as a `TCoffeeAbortedException` on every committed `FAIL` row.

---

## 2026-11 audit round 2 (A2-38, A2-39, A2-40, A2-41, A2-42; B09 F81-F86): input contract, size guards, parallel exceptions, cancellation

| ID | Test | Input | Expected | Source |
|----|------|-------|----------|--------|
| IC-1 | `SequenceAligner_InputContract_Tests.NullElement_*` | null element on every ClustalW (4), MUSCLE (3), T-Coffee program and Seqeron aligner overload; null collection | `ArgumentNullException` naming the index (was `NullReferenceException` on the ClustalW / MUSCLE DNA overloads and `DnaRows`) | contract |
| IC-2 | `MalformedUserTreeBranchLength_*`, `MalformedTCoffeeMethodAndWeightValues_*`, `ClustalWUserMatrixParse_StaysFormatException_ByDocumentedDesign` | `UserTree = "(0:abc,1:1,2:1);"`, `-method slow_pair@EP@GOP@abc`, `-seq_weight @1,x,3` | `ArgumentException` (was `FormatException`); `ClustalWUserMatrix.Parse` keeps `FormatException` (documented) | contract |
| IC-3 | `EmptyMember_ThreeProgramPorts_ShareOneMessage`, `EmptyMember_SeqeronAligners_ReturnAnAllGapRow` | one empty member | ClustalW / MUSCLE / T-Coffee: one message naming the program and the aligners that accept empty members; star / progressive / iterative / consistency: an all-gap row | the programs refuse empty input (ClustalW cannot read it, MUSCLE ignores empty FASTA records) |
| SG-1 | `TCoffeeProba_RefusesAnOversizedPair_BeforeAllocating`, `MuscleProfileDp_RefusesAnOversizedProfilePair_WithAClearArgumentException` | 9 000 x 9 000 protein (405 M floats > 2^28); 40 000 x 40 000 profile columns (> 1e9 cells) | `ArgumentException`; `TCoffeeProbaMaxMatrixFloats` = 2^28, `MuscleProfileMatrixMaxCells` = 1e9; 30 000 x 30 000 still accepted | t_coffee `util_dp_suboptimal_nw.c` (32-bit `l = (len1 + 1)(len2 + 1) NumMatrixTypes`), MUSCLE 3.8.1551 `nwsmall.cpp` (one byte per cell) |
| PX-1 | `ExceptionInTheParallelPairStage_KeepsItsType_ForEveryK` | k = 2, 3, 5, linear + affine; the range guard fires inside a worker | `ArgumentException` (never `AggregateException`) | `ExceptionDispatchInfo`, lowest failing index |
| PX-2 | `MaxDegreeOfParallelism_ChangesNothing_AndIsValidated` | 1, 2, 8, -1; 0 and -2 | identical rows / score; `ArgumentOutOfRangeException` | contract |
| CN-1 | `SequenceAligner_Cancellation_Tests.PreCancelled_*` | pre-cancelled token on star / anchored / progressive / iterative / consistency / exact, ClustalW (DNA, protein), MUSCLE, T-Coffee program, `WithCancellation` over top-N / co-optimal / affine / general-gap, the four `GenomicAnalyzer` scans | `OperationCanceledException`, no result | contract |
| CN-2 | `CancelledWhileRunning_*` | 100 ms token on calls that run for seconds-minutes (5 x 6 kb star affine, 2 x 30 kb ClustalW, 2 x 20 kb MUSCLE, 2 x 2.5 kb protein T-Coffee proba, 4 kb top-N / affine, 1.2 kb general gap, 3 Mb ORFs) | `OperationCanceledException` within 15 s | row-boundary polls |
| CN-3 | `UncancelledToken_GivesTheSameResultAsNoToken_AndTheScopeDoesNotLeak` | live token; then a cancelled call; then the plain call | identical rows / score; the scope is restored | contract |

## Round 2 additions (audit A2-36, F100)

MCP `PortStageToolsTests` (Seqeron.Mcp.Alignment.Tests): `ClustalWPairwiseAlign_Protein_ReproducesTheClustalWScores`, `ClustalWPairwiseAlign_Dna_DelegatesToTheLibrary`, `ClustalWQuickPairwiseAlign_ReproducesTheHandRun`, `TCoffeeSlowPairAlign_PlacesTheGapInside_WorkedExample`, `TCoffeeLalignTopAlignments_ScoresInOrder`, `TCoffeeKtupSimilarity_And_GuideTree`, `TCoffeePrimaryLibrary_DelegatesAndOrdersEntries`, `SumOfPairsScore_LinearAffineAndMatrix`.

## Audit A2-14 (B09 F101 / F102 / F103): ClustalW UPGMA oracle, iteration, profile commands

Tests (`Seqeron.Genomics.Tests`): `SequenceAligner_ClustalWUpgmaOracle_Tests` x3 (the rows of the DNA and PROTEIN `MultipleAlignClustalW` with `Clustering = Upgma` and the `-NEWTREE` text equal the unmodified Debian source built at `-O0` on 600 / 600 user-run sets), `SequenceAligner_ClustalWUpgma_Tests` (hand-derived joins / heights / float weights), `SequenceAligner_ClustalWIteration_Tests` x6 (`-ITERATION=tree | alignment -NUMITER` rows and "FINAL score" lines equal the installed clustalw 2.1 on 700 / 700 sets, > 300 of them differing from the non-iterated alignment; defaults 3 iterations / none; validation `NumIterations <= 0`, undefined mode; two sequences never iterated; UPGMA + iteration integrity), `SequenceAligner_ClustalWProfile_Tests` x4 (`AlignClustalWProfiles` / `AlignClustalWSequencesToProfile` rows equal clustalw `-PROFILE` / `-SEQUENCES` on 1 200 / 1 200 sets; inputs kept in order; empty profile, invalid residue and `Iteration = Tree` refused). MCP (`Seqeron.Mcp.Alignment.Tests`): `ClustalWAutoTypeUpgmaTests` (`type = auto`, `clustering = upgma`, `iteration` / `numIterations` equal the library and validated), `ClustalWProfileToolsTests` x3 (`clustalw_profile_align` / `clustalw_sequences_to_profile` equal clustalw on the documented example; protein / options; bad input). Oracles: `clustalw_upgma_oracle.tsv`, `clustalw_iteration_oracle.tsv`, `clustalw_profile_oracle.tsv` (`scripts/alignment/generate_clustalw_*_oracle.py`). Test status: the `ClustalW|MultipleAlign` targeted tier 1 396 passed, MCP.Alignment 174 passed.



## Audit A2-14b / c / e (B09 F104 / F105 / F106): guide-tree files and row order, output formats and command line, secondary-structure masks

Tests (`Seqeron.Genomics.Tests`, all against the installed `clustalw` 2.1): `SequenceAligner_ClustalWCommands_Tests` x19 - the oracle `clustalw_commands_oracle.tsv` (3 300 sets: `-ALIGN -NEWTREE`, `-NEWTREE`, `-USETREE` with valid / odd / malformed trees, `-PROFILE`, `-SEQUENCES` with `-USETREE1/2`, `-NEWTREE1/2`; rows, `-OUTORDER` order, guide-tree text and the refusals as `ClustalWTreeFileException`), the five-sequence worked example, the two-sequence tree line, case-insensitive name matching, seven typed refusals, input validation, `ClustalWFixed` ties (`0.015625` -> `0.01562`). `SequenceAligner_ClustalWOutput_Tests` x14 - `clustalw_output_oracle.tsv` (2 611 files, 373 per format (330 aligned, 43 `-CONVERT`), byte for byte), oracle coverage, MSF and CLUSTAL worked examples, the other formats, `-RANGE` / `-SEQNO_RANGE`, refusals. `SequenceAligner_ClustalWStructure_Tests` x9 - `clustalw_structure_oracle.tsv` (1 200 profile pairs, CLUSTAL and GDE profile files), a hand-run helix / strand example incl. `-HELIXGAP=9`, the readers' rules (first structure line wins, missing line in a later block), `-NOSECSTR1`. `ClustalWCommandLine_Tests` x8 - option spelling, `findMatch`, ranges, verbs, refusals, `ClustalWProgram.Run`. MCP (`Seqeron.Mcp.Alignment.Tests`): `ClustalWCommandToolTests` x3, contract test +`clustalw_command`. Undefined behaviour of the program (nested trichotomy, text ending inside a name, node / leaf overflow, single-child nodes) is classified, not claimed: 21 hang / crash rows of the oracle accept any refusal.

## Audit A2-20 (B09 F118): Kalign 3.3.5

Tests (`Seqeron.Genomics.Tests`, all against the installed Debian `kalign` 1:3.3.5-1 AVX2 build): `SequenceAligner_KalignOracle_Tests` x527 (`kalign_oracle.tsv`: 526 runs - DNA / RNA / protein families with lower-case runs, ambiguity codes, gap characters and sequences without residues; every `--type`, every `--gpo / --gpe / --tgpe` combination incl. the extremes and the quirk cells; `-n 1 / 2 / 4`; sets of more than 500 columns, 255 - 520 residues, 30 - 90 sequences, 99 / 100 / 101 / 110 / 120 / 150 / 200 / 1 000 sequences (UPGMA / bisecting k-means); the 5 refusals: rows equal, or `ArgumentException`), `SequenceAligner_KalignProgram_Tests` x88 (`kalign_command_oracle.tsv`: 83 command lines - type words, `atof` / `atoi` prefixes, `getopt_long_only` forms and abbreviations, two inputs, missing file, mismatch, one sequence, 60-column output; exit status and text equal), `SequenceAligner_Kalign_Tests` x25 (Myers bit-parallel distance against an independent semi-global edit-distance DP for pattern lengths 1 ... 300 incl. the 255 limit and 2 / 4 / 13 letters, the `--tgpe` quirk, `DnaInternal` = DNA with terminal extension 8, dropped sequences, alphabet detection, case / gap handling, contracts, the one-thread hand-run). MCP `KalignTests` x4 (`multiple_align_kalign`, `kalign_command`).


## Audit A2-20b (B09 F222): Kalign 3.3.5 formats, inputs and non-quiet output

Tests (`Seqeron.Genomics.Tests`, against the installed Debian `kalign` under the time shim `scripts/alignment/faketime.c`): `SequenceAligner_KalignFormatsOracle_Tests` x2 - `kalign_formats_oracle.tsv` (251 command lines from `scripts/alignment/generate_kalign_formats_oracle.py`: `--format msf|clu|fasta|fa` to stdout / `-o` / `-o dir/x`, names up to 300 characters, widths around multiples of 60, non-quiet runs with `-n 1 / 2`, aligned / unaligned input warnings, empty sequences, MSF / Clustal inputs written by the program and by other tools, a Clustal line without white space, MSF without `Name:` lines (the program segfaults), several inputs of two formats and alphabets, a file of unknown format, `-h`, `--showw`, `-n 0`): exit status, standard output (the text after `CPU Time:` normalised) and the `-o` file equal on 251 / 251; `SequenceAligner_KalignProgram_Tests.MsfAndClustalOutput_AreWritten_AsTheProgramWritesThem`; MCP `KalignTests` unchanged.

## Audit A2-19 (B09 F119): MAFFT 7.505, stage 1 (two sequences)

Tests (`Seqeron.Genomics.Tests`, against the installed Debian `mafft` 7.505-1): `SequenceAligner_MafftPair_Tests` x886 (`mafft_pair_oracle.tsv`: 884 runs - 220 protein and 200 nucleotide pairs, 160 with IUPAC codes / unusual symbols (digits, `*`, `#`, `.`, B Z J O U) under `--anysymbol`, mixed alphabets, `-` in the input, `--op`, `--ep`, `--bl` 30 ... 0, `--nuc` / `--amino`, `--retree`, `--fft`, `--nofft`, `--reorder`, 60 plain-command pairs, 47 pairs of 400 - 1 500 residues: rows equal; the plain command returns nucleotides in lower case), `SequenceAligner_Mafft_Tests` x6 (the scoring matrices and integer penalties against the matrices `disttbfast -d` prints: BLOSUM62 default, BLOSUM80 with `--op 2 --ep 0.1`, the Kimura matrix; hand runs of the program; free end gaps; contracts), T-Coffee `MafftPair_EveryRowEqualsTheProgram` (the 460 `mafft_pair` / `mafft_mixed` sets of the program, replaces the OP30 "port waits" test). MCP `MafftTests` x2.

## Audit A2-27 (B09 F120): EMBOSS `-aformat` writers

Tests (`Seqeron.Genomics.Tests`, against the installed EMBOSS 6.6.0): `SequenceAligner_EmbossFormats_Tests` - `emboss_aformat_oracle.tsv` (needle / water / stretcher x 29 formats x 114 alignment sets, 3 303 files byte for byte, faked run dates), `emboss_seqout_oracle.tsv` (`seqret`, 11 sequence-set formats, 1 309 files), >= 100 cases per format, format-name parsing, a hand-written `pair` layout, argument contracts. MCP `EmbossFormatTests` (`write_alignment`, `stretcher_align`).

## Audit A2-28 (B09 F121 / F122): stretcher, needle options, water

`SequenceAligner_Stretcher_Tests` (689 oracle pairs: rows and score, optimum bound, defaults, contract), `SequenceAligner_Needle_Tests` (887 pairs incl. `-endweight` / `-endopen` / `-endextend`, `ajRoundFloat`, bound against the affine optima), `SequenceAligner_Water_Tests` (860 pairs: rows, starts, scores).

## Audit A2-27b (B09 F126): Biopython Bio.Align text formats
`AlignmentFormats_Tests` (`biopython_msa_oracle.tsv`, Biopython 1.88): writers byte for byte (fasta, a2m, clustal, phylip, nexus, stockholm; 210 alignment sets each, multi-alignment files, quoted / duplicate ids, descriptions, Stockholm annotation, nexus interleave) and readers equal to `Bio.Align.parse` on every file Biopython can read (>= 100 per format); MCP `FormatParameterTests`.

## Audit A2-23 (B09 F127): ProbCons 1.12

Tests (`Seqeron.Genomics.Tests`, against the installed Debian `probcons` 1.12-14): `SequenceAligner_ProbConsOracle_Tests` x408 (`probcons_oracle.tsv`: 407 command lines run by `scripts/alignment/generate_probcons_oracle.py` - 50 default runs, `-c 0 ... 5`, `-ir 0 ... 1000`, `-clustalw`, `-a`, `-annot`, `-pairs`, `-viterbi`, `-co`, `-go` / `-ge` on equal-size merges, `-pre 1 ... 4`, `-t` / `-e` on one to three files, `-p` trained parameter files, `-v`, FASTA oddities (CRLF, no trailing newline, wrapped lines, spaces / `.` / `-` inside, empty sequences, one residue, unknown character), MSF / Clustal inputs, several files, 200 - 400 residues, tiny / identical sets, 21 error exits: exit code, stdout, stderr (heading, parameter set at `%.10g`, tree, `-v` lines) and every written file byte for byte), `SequenceAligner_ProbCons_Tests` x13 (the program's default rows and tree, `-annot` scores, `-a` order, Viterbi pairs, glibc `rand()` values 1804289383 ..., `CFormat.General`, default parameters = `Defaults.h` and the parameter-file round trip, training file layout, option validation, one sequence, the Clustal writer, the Debian patch message, `strtol` base 0), MCP `ProbConsTests` x3 (`multiple_align_probcons`, `probcons_command`). The T-Coffee `proba_pair` tests (`TCoffee`) pin the shared pair-HMM engine after its refactoring.

## Audit A2-26 (B09 F128): POA 2.0

Tests (`Seqeron.Genomics.Tests`, against the installed Debian `poa` 2.0+20060928): `SequenceAligner_PoaOracle_Tests` x360 (`poa_oracle.tsv`: 359 command lines run by `scripts/alignment/generate_poa_oracle.py` - protein families with BLOSUM80 and the truncated / decaying gap model of `blosum80_trunc`, nucleotide families with four lower-case matrices (`GAP-PENALTIES`, `GAP-PENALTIES-X`, `GAP-TRUNCATION-LENGTH`, `GAP-DECAY-LENGTH`), default (iterative local), `-do_global`, `-do_progressive` (local / global), `-read_pairscores`, `-fuse_all`, `-hb` / `-hbmin` / `-best`, `-preserve_seqorder`, `-tolower` / `-toupper`, unknown letters, odd FASTA (CRLF, comments, `*` lines, long names, empty records, one residue), `-v`, `-printmatrix`, 150 - 300 residue sets, error exits: exit code, standard output and the PO / PIR / CLUSTAL files byte for byte), `SequenceAligner_Poa_Tests` x9 (rows / PIR / CLUSTAL / PO head of the program's default run, heaviest bundle, progressive + permutation + `-v`, matrix parser and printer, letters outside the matrix, validation, usage, FASTA reader), MCP `PoaTests` x3 (`multiple_align_poa`, `poa_command`).



## Audit A2-26b (B09 F223): POA MSA / PO file input

Tests (`Seqeron.Genomics.Tests`, against the installed Debian `poa` 2.0+20060928): `SequenceAligner_PoaMsaOracle_Tests` x2 - `poa_msa_oracle.tsv` (687 command lines from `scripts/alignment/generate_poa_msa_oracle.py`: single / two MSAs, MSA + FASTA, `-read_msa_list`, `-subset` / `-remove` for both files, PO / PIR / CLUSTAL inputs the program wrote and hand-made ones (comments, digits, CRLF, wrapped rows, blank lines, 510-character `fgets` pieces, empty / garbage files, PO files with unusual fields), every option, error exits): exit status (non-zero = non-zero; rows where the program crashes or loops stop here), standard output, `...` lines of standard error and the three written files equal on 687 / 687; `SequenceAligner_Poa_Tests.Program_Usage_AndMissingMsaFileStopsTheRun`.

## Audit round 3 (B09 F129 - F135, WP35): input contract, resource limits, penalties, cancellation, MUSCLE default, T-Coffee `random`

| ID | Test | Locks |
|----|------|-------|
| A3-1a | `SequenceAligner_ProbCons_Tests.GapCharactersInTheInput_AreDroppedLikeTheProgramsReader`, `AnAllGapSequence_IsRefused...`; `probcons_oracle.tsv` rows `gapin` (7) | `-` / `.` / blanks dropped, `*` kept, all-gap record refused by the API (dropped by the program): rows = installed probcons 1.12 |
| A3-1b | `SequenceAligner_Poa_Tests.NonAsciiCharacters_AreReadAsUtf8BytesLikeTheProgram`; `poa_oracle.tsv` rows `utf8` (5) | byte view: U+017F / U+0131 = two residues each = `A`, NBSP / NEL not white space, ASCII-only case folding |
| A3-2 | `SequenceAligner_ResourceLimits_Tests` (13) | `NeedleMaxCells` 2^27, `WaterMaxCells` 2^28, `SupermatcherMaxPathCells`, `PoaMaxMatrixCells` 2^27, `ProbConsProfileMaxCells` 2^26 + pair-HMM guard with ProbCons wording (T-Coffee wording kept), `MafftMaxPairCells`, `Ksw2Aligner.MaxTracebackBytes`: `ArgumentException` of one style raised with < 64 MB allocated; wide band = no band; PHYLIP / NEXUS / Stockholm header counts |
| A3-3 | same class | NaN / +-Infinity of every EMBOSS penalty and `MatrixScoring.Emboss` -> `ArgumentOutOfRangeException`; stretcher path score < 2^28; dual-affine `int.MinValue`; ksw2 `q` 0..127, `e` 1..127 |
| A3-4 | `SequenceAligner_InputContract_Tests.EmptyList_SixProgramPorts_ShareOneMessage_MafftReturnsEmpty`, `SequenceAligner_Mafft_Tests.SingleSequence_IsWrittenBackLikeTheProgram` (7), MCP `KalignTests` / `MultipleAlignTests` | empty list: six refusals, one message; MAFFT empty result (exit 0 of the program) and single sequence written back |
| A3-5 | `SequenceAligner_Cancellation_Tests` (+10: pre-cancelled for 12 entry points, mid-run for needle / water, stretcher, supermatcher, needleall, ksw2 x4 + dual affine, T-Coffee Myers-Miller, `TcGotohProfile`, MUSCLE refinement, anchors, distmat; every mid-run test asserts the poll frame in the stack trace), MCP `CancellationTokenTests` (2) | the ambient token is polled at row / iteration boundaries of every ported program and the MCP wrappers pass it |
| A3-6 | `SequenceAligner_TCoffeeOptions_Tests.RandomDistance_EveryRowEqualsTheProgram_WithDebugSrand` (310 sets of `tcoffee_program_options_oracle.tsv`, feature `dist_random`), `RandomDistance_SeedIsAnOption_...` | `-distance_matrix_mode random` = glibc `rand() % 100` per pair, seed 10 (`DEBUG_SRAND`), `TCoffeeProgramOptions.RandomSeed` |
| A3-7 | `SequenceAligner_MuscleDefaultIterations_Tests` (`muscle3_default_iters_oracle.tsv`) | the program default of `-maxiters` is 8 (not the stale `-help` 16, not the options default 2): `MuscleOptions.ProgramDefaults` |

## Audit A3-9a (B09 F150): nucmer clustering

| ID | Test | Oracle / what | Verdict |
|---|---|---|---|
| A3-9a | `NucmerClusters_Tests.Cluster_EqualsMgaps_AllRows`, `WorkedValues_RunOnTheBinary`, `Validation` | 758 mgaps 3.23 runs (`mgaps_oracle.tsv`) | as specified |
| A3-9a | `GenomicAnalyzer_NucmerClusters_Tests.FindNucmerClusters_EqualsMummerPipedIntoMgaps_AllRows`, `ForwardOnly_AndDefaults` | 700 `mummer | mgaps` pipelines (`nucmer_clusters_oracle.tsv`), three match modes, both strands | as specified |

## Audit A3-9c (B09 F151): minimap2 chaining

| ID | Test | Oracle / what | Verdict |
|---|---|---|---|
| A3-9c | `Chain_EqualsMinimap2_AllCases` | 450 minimap2 2.26 runs (`minimap2_chain_oracle.jsonl`): chains from minimap2's own anchors, PAF `cm` / `s1` | as specified |
| A3-9c | `CollectAnchors_EqualsTheSeedArrayOfMinimap2_AllCases`, `MapChains_EqualsTheChainsOfMinimap2_FromTheSequences` | the same runs from the sequences: anchor arrays (order included) and chains | as specified |
| A3-9c | `SmallChain_ScoreIsTheSumOfAnchorSpansMinusGapPenalties`, `CollectAnchors_FrequentMinimizers_AreNotPorted_AndValidation`, `Validation_AndSorting` | hand values, the documented limitation, guards | as specified |

## Audit A4-4 (B09 F156): option boundaries

`OptionBoundaries_Oracle_Tests` (`option_boundary_oracle.json`, `generate_option_boundary_oracle.py`): 30 minimap2 2.26 option values (`MinimapChainOptions`: chain score + anchors equal to the program's PAF `s1` / `cm`, NaN / infinite / huge scales through the C float-to-int conversion), `KmerSize` outside 1 - 28 refused (`mm_check_opt`, `mm_sketch` assertion), `q_occ_frac` NaN / -1 / inf = filter off, 17 mgaps 3.23 option values (`NucmerClusterOptions`: negative / zero / NaN / infinite / huge values accepted, the C `(long)` conversion = LONG_MIN reproduced).

## Audit A4-1 (B09 F157): cancellation

`SequenceAligner_Cancellation_Tests` (`PreCancelled_A4_1_...`, one `CancelledWhileRunning_...` test per class) and the MCP `CancellationTokenTests` / `EmbossProgramToolsTests`: `Minimap2Chaining.Chain` / `CollectAnchors`, `NucmerClusters.Cluster` poll the ambient token (a mid-run cancel stops with `OperationCanceledException` inside the named poll) and the MCP tool takes the request token.

## Audit A4-2 (B09 F158): resource guards

`SequenceAligner_ResourceLimits_Tests`: `Minimap2Chaining.CollectAnchors` / `Chain` (`MinimapMaxAnchors`), `NucmerClusters.Cluster` (`NucmerMaxClusterMatches`) - documented limit constants (`SequenceAligner.ResourceLimits.cs`), the exception of the `SequenceAligner.ResourceLimit` style raised before the allocation in 64-bit arithmetic, and a low-allocation assertion per limit.

## Audit A4-7 (B09 F165 / F166): minimap2 seeding options, frequent minimizers, primary / secondary hits, long join, RMQ

| ID | Test | Locks |
|---|---|---|
| A4-7 | `Minimap2Seeds_Tests.CollectSeeds_EqualsMinimap2_AllGroups` | 1 980 `--print-seeds` runs (`minimap2_seed_oracle.jsonl.gz`, six groups of 330: default frequent minimizers, `-f` / `-U` / `-e` / `--q-occ-frac`, `-H`, `-T`, `-X`, `--for-only` / `--rev-only`): anchors, `RS` repeat length and `mid_occ` equal |
| A4-7 | `Minimap2Seeds_Tests.FrequentMinimizers_AreHandledLikeTheProgram_NotRefused` | the old `NotSupportedException` is gone: poly-A key dropped, `RepeatLength > 0`, `CollectAnchors` = `CollectSeeds` anchors |
| A4-7 | `Minimap2Map_Tests.Map_SecondaryChainMarking_EqualsMinimap2` / `Map_ChainingOptions_EqualMinimap2` / `Map_LongJoinRechaining_EqualsMinimap2` / `Map_RmqChaining_EqualsMinimap2` | 1 320 PAF runs (`minimap2_map_oracle.jsonl.gz`, four groups of 330): columns 3 - 12 and `tp cm s1 s2 dv rl` of every line (secondary marking, `-N` / `-p` / mask options, `-r b,B` long join, `--rmq`, `-g -n -m` and the chain scales) |
| A4-7 | MCP `Minimap2ChainsTests.Minimap2Chains_PafMode_IsThePafLineOfMinimap2`, `..._FrequentMinimizers_...` | the `paf` mode and the new options of `minimap2_chains` |

## Audit A5-1 / A5-2 (B09 F167 / F168): minimap2 query-name hash, T-Coffee pair-DP limits and cancellation

| ID | Test | Locks |
|---|---|---|
| A5-1 | `Minimap2Map_Tests.Map_NonAsciiAndNulQueryNames_EqualMinimap2` | group `qname` of `minimap2_map_oracle.jsonl.gz` (330 PAF runs, names with a first byte >= 0x80 (189), 3- / 4-byte UTF-8 characters, U+0000 (99); references with 2 - 4 identical copies of the query): the hash of the name (`__ac_X31_hash_string`: signed `char`, C string) orders equal-score chains like `minimap2 2.26` |
| A5-2 | `SequenceAligner_TCoffeeLimits_Tests.BiphasicPair_Above2GiB_...` / `FastPair_FlaggedDiagonalBand_IsRefused` / `GotohPairWise_TracebackMatrix_...` / `ViterbiPairWise_IsRefused_...` | `TCoffeePairDpMaxBytes` (2 GiB): one `ArgumentException` of the `ResourceLimit` style per pair DP (`biphasic_pair_wise` 48 bytes per cell, `glocal_pair_wise` 40, `gotoh_pair_wise_lgp` 32, `gotoh_pair_wise` / `_sw` traceback 4, `viterbi_pair_wise` 44, FASTA band 20 per row x flagged diagonal), 64-bit arithmetic (the largest int lengths do not overflow), raised before the allocation (< 64 MB allocated) |
| A5-2 | `LalignTopAlignments_IntMaxValue_MeansEveryAlignment`, `FracMinHashSketch_KAboveTheSequenceLength_IsAnEmptySketch_NotAnOutOfMemory` | the lazily grown vertex list of `lalign_id_pair`; `KmerAnalyzer.CreateFracMinHashSketch` with `k = int.MaxValue` |
| A5-2 | `BiphasicPair_Cancels_...`, `GlocalGotohSwAndGotohProfile_Cancel_...`, `FastPair_Cancels_...`, `LgpAndViterbiDpModes_Cancel_...`, `LalignAndSlowPair_Cancel_...` | a pre-cancelled token and a cancellation at the middle poll INSIDE each pair method (`PollCancellation`: the call is run once to record its polls, then cancelled at the poll of the named method: no sleep) |

## Audit A5-3 / A5-4 / A5-5 (B09 F169 / F171 / F170): null contract, deterministic cancellation, MCP documents

| ID | Test | Locks |
|---|---|---|
| A5-3 | `SequenceAligner_NullContract_Tests` (MUSCLE stage methods, Kalign / POA / ProbCons runners, `Minimap2Chaining.Map`, T-Coffee k-tuple / primary library / pair stages, `Ksw2Aligner.SimpleMatrix`, `AlignmentFormats` a2m, POA name / matrix) + MCP `NullContractTests` | a null list / element, a character the matrix cannot index, an out-of-range size -> `ArgumentException` family with the parameter name (before: NullReference / IndexOutOfRange / KeyNotFound / Overflow); the programs' own argument handling is untouched (usage screen exit 255, unknown option exit 1) |
| A5-4 | `SequenceAligner_Cancellation_Tests.CancelledWhileRunning_*` (every test names its poll site) | `PollCancellation`: the token is cancelled at the 3rd poll of the named method, the call stops with `OperationCanceledException` (no sleep, thread or time bound) |
| A5-4 | `Minimap2Seeds_Tests.FrequentMinimizers_...`, MCP `Minimap2ChainsTests`, `SequenceAligner_ResourceLimits_Tests` | exact values of the installed programs instead of `> 0` / `Is.Not.Empty`: `mid_occ = 10`, `RS 24`, 35 / 1 866 anchors (`minimap2 --print-seeds`), stretcher score 20, seqmatchall 282 lines, dotmatcher plot (397 segments, first segment) |
| A5-5 | `ToolDocsContractTests` (MCP.Alignment.Tests: 76 tools; MCP.Analysis.Tests: the B09 tools) over `ToolDocsContractTestsBase` | description, parameter names / order / type / default / description, `required`, the `.md` Input Schema table equal the reflection metadata of the tool; `python3 scripts/docs/sync_mcp_tool_docs.py [--check]` regenerates / checks the same fields |

## WP42 additions (A5-9 / A5-10)
- `Minimap2Presets_Tests` (11 presets x 220 oracle cases, `-F` / `-G`, the options.c tables, `FormatPaf`), `MashCommands_Oracle_Tests` (bounds / dist / triangle: 660 oracle rows) and the MCP tests `Minimap2PresetToolsTests`, `MashToolsTests`; sources and numbers in docs/Evidence/ALIGN-MULTI-001-Evidence.md.

## Audit round 6 / WP43 (B09 F177 - F183)

| ID | Test | Locks |
|---|---|---|
| A6-1 | `KmerAnalyzer_SketchStreaming_Tests` | the streaming MinHash / FracMinHash builders equal the previous string-per-k-mer implementation on 1 200 + 1 200 random inputs (k 1 .. 60, canonical / forward, case / N / IUPAC / U, null and short records) and on all 65 536 UTF-16 code units; a 4 Mbp genome allocates < 16 MB; the token is observed between records and inside a record |
| A6-2 | `EmbossAcdCommandLine_Oracle_Tests` (16), `EmbossSequenceQualifiersText_Oracle_Tests` (16), `EmbossSequenceQualifiers_CommandLine_Tests` (9), MCP `EmbossSequenceQualifiersToolsTests` | the qualifier text resolves as the ACD library of the program resolves it: 2 080 ACD-log rows (`emboss_acd_args_oracle.tsv`) and 640 report rows (`emboss_seqopts_text_oracle.tsv`) equal; the forms of the audit (`-sbegin1 5` on one input, `-sbegin 3` after the files, `=`, abbreviations, ambiguity) with the program's messages |
| A6-3 | `SequenceAligner_MatcherGap_Tests` | `matcher -gapopen` below `-gapextend`, 0 / 1 / 100, reset warnings and refusals: 872 rows of `matcher_gap_oracle.tsv` |
| A6-4 | `EmbossNullContract_Tests` (5), `MashCommands_Oracle_Tests.Dist_And_Triangle_RefuseNullInputs_...` | typed exceptions with the parameter name and the index for null set members / `Name` / `Residues`, hits without `EmbossStart`, an undefined report format, null report fields, null Mash inputs |
| A6-5 | `SequenceAligner_SupermatcherWordLength_Tests` | `-wordlen` around / above the targets: stderr `Error: no kmers found`, error file, alignments equal on 260 rows; `int.MaxValue` returns at once; sibling programs |
| A6-6 | `MashCommands_Oracle_Tests.FilesWithoutARecordOfKSymbols_AreRefused_...` | 80 `error` rows of `mash_commands_oracle.jsonl.gz` (status 1 messages, the SIGFPE of `dist -i`) |
| A6-7 | `Minimap2Presets_Tests` (11 presets x 260 cases), MCP `Minimap2PresetToolsTests` | `-x PRESET` followed by `-k -w -g -r -n -m -p -N -e -f -U --rmq --secondary -H --no-long-join` (half with the default value): 2 860 PAF texts; whole PAF lines of `-x asm5 -k 15`, `-w 10 -m 40 -n 3`, `-x ava-ont -m 40` |
| A6-8 | MCP wiring tests, `GenomicAnalyzer_ProteinPatterns_Oracle_Tests`, `SequenceAligner_MultipleAlignAnchored_Scenarios_Tests` | the program's whole result instead of `Any(...)` / `Is.Not.Empty`; exact counters |
| A6-10 | MCP `EmbossReportSingleSearchTests` | the report mode runs the search once (the same number of polls as without `reportFormat`) |

## Audit round 7 (WP45): minimap2 2.26 base-level alignment, paired reads, program options (A3-9c2, A5-9b, A5-9c; F187 - F190)

| ID | Test | Locks |
|---|---|---|
| A3-9c2 | `Minimap2Align_Tests` (12 groups + 1) | the whole stdout of `minimap2 [options] r.fa q.fa` as text on 2 700 command lines (`minimap2_align_oracle.jsonl.gz`): `-c` with `cg` / `cs` / `MD` / `--eqx` and the tags `NM ms AS nn ts de zd` for the presets `map-ont` / `map-hifi` / `map-pb` (homopolymer-compressed) / `asm5` / `asm10` / `asm20` (150 each), random alignment options (300), SAM (300: flags 0 / 16 / 256 / 2048, MAPQ, `SA`, clipping, header with `@PG CL`, `-Y -L --secondary-seq --sam-hit-only -R -y -Q`), the short-read mode (300), spliced alignment (300: `-x splice`, `-u -J -C --junc-bonus --junc-bed --splice-flank -G`, `ts`, `N` / `~` in `cg` / `cs`), repeat families with secondary hits (300), Z-drop and inversion rescue (300: `tp:A:I`, `zd`); the oracle contains the features it claims (`OracleContainsTheFeaturesItClaims`) |
| A3-9c2 | `Ksw2SpliceLocal_Tests` (6) | the C `ksw_ll_i16` (score, query end and target end including the zero-scored padding lanes), `ksw_exts2_sse`, `ksw_extz2_sse`, `ksw_extd2_sse` of minimap2 v2.26 (`gcc -msse4.1`, `scripts/alignment/ksw2/`) on 4 x 500 random inputs, every field of `ksw_extz_t` and the CIGAR (`ksw2_c_oracle.tsv`); an intron as `30M100N30M`; the local score of identical sequences |
| A3-9c2 | `SequenceAligner_ResourceLimits_Tests.Ksw2_GapPenalties_FollowTheHeaderOfKsw2` | `ExtendD2` accepts `e2 = 0` (`-E 2,0`), refuses `e2 < 0` |
| A5-9b / A5-9c | `Minimap2Pair_Tests` (7 groups + 1) | the whole stdout of the paired mode on 2 100 command lines (`minimap2_pair_oracle.jsonl.gz`): two query files and an interleaved file, PAF without / with `-c` (`mm_pair`), SAM (flags 0x1 0x2 0x8 0x20 0x40 0x80, mate position, TLEN, MAPQ), `--frag=yes` without a preset (forward-forward), options `-F -N -p -k -w -A -B -O -E -s -z --end-bonus --no-pairing -r -g -n -m -M --secondary=yes --heap-sort=no -t -K -I -2` |
| A5-9d | `Minimap2Options_Tests` (9 groups) | 2 401 command lines (`minimap2_options_oracle.jsonl.gz`): index parts (`-I`) and `--split-prefix`, chaining / seeding options with `-c`, output options and `-o FILE`, `--alt` / `--alt-drop`, `--qstrand`, all-vs-all with `-c`, unusual input (lower case, IUPAC, U, empty and very short reads, CRLF, blank lines, multi-line FASTA / FASTQ, comments), a query of 610 kb with more than 65 535 operations (`CG:B:I`), option errors / help / version (stdout, stderr at `-v 2`, exit status) |
| A3-9c2 .. A5-9d | MCP `Minimap2MapToolTests`, `ToolRegistrationContractTests`, `ToolDocsContractTests` | `minimap2_map`: `-c` PAF with `cg` and the tags (whole line prefix), SAM header with the command line, a proper pair (flags 99 / 147, TLEN +-450), exit status and error text, `-o` file, `-I` warning, argument validation; `ksw2_extend` variants `exts2` (`30M100N30M`) and `ll` (score 20 at target 9); tool documents equal the metadata of the code |


## Audit A3-9b (B09 F199): nucmer `postnuc` and the whole pipeline

| Item | Tests | Oracle | Status |
|---|---|---|---|
| A3-9b | `MummerNucmer_Oracle_Tests.Nucmer_FilesEqualTheProgram_AllRows` | 434 `nucmer` 3.23 runs (`nucmer_oracle.jsonl.gz`, `generate_nucmer_oracle.py`): `out.delta` / `out.cluster` text equal, 5 abort rows | as specified |
| A3-9b | `MummerNucmer_Oracle_Tests.LibstdcxxSort_...`, `SwAlign_MatchScore...`, `Postnuc_ReadsTheMgapsTextOfTheProgram` | libstdc++ `std::sort` order (stable below 17 elements), the +3 / -7 score, the program's output for one mgaps text | as specified |
| A3-9b | MCP `EmbossProgramToolsTests.Nucmer_WritesTheDeltaFileOfTheProgram`, `ToolRegistrationContractTests` | `nucmer -l 15` of the F150 example | as specified |

## Audit A3-9b (B09 F200 - F202): the delta programs

| Item | Tests | Oracle | Status |
|---|---|---|---|
| A3-9b | `MummerShow_Oracle_Tests.DeltaPrograms_OutputEqualsTheProgram_AllRuns` | 2 442 runs of `delta-filter`, `show-coords`, `show-aligns`, `show-snps`, `show-diff`, `show-tiling` (+ files), `dnadiff` (`mummer_show_oracle.jsonl.gz`, `generate_mummer_show_oracle.py`): every output text / file equal (dnadiff's hash-ordered lines as sets) | as specified |
| A3-9b | MCP `EmbossProgramToolsTests.MummerDeltaPrograms_...`, `Dnadiff_...`, `ToolRegistrationContractTests`, `ToolDocsContractTests` | the installed programs on the three-substitution and two-contig examples | as specified |

## Audit A3-9d / A3-9e / A3-9f (B09 F204 - F209): the other programs of the MUMmer package

| Item | Tests | Oracle | Status |
|---|---|---|---|
| A3-9e | `MummerCli_Oracle_Tests` (2) | 700 `mummer` runs (`mummer_cli_oracle.jsonl.gz`, `generate_mummer_cli_oracle.py`): stdout and the deterministic stderr lines equal, 47 error runs, option conflicts refused | as specified |
| A3-9d | `MummerPromer_Oracle_Tests.Promer_AllStagesEqualTheProgram_AllRows`, `Translate_...`, `Prepro_...`, `Postpro_...` | 320 `promer` runs (`promer_oracle.jsonl.gz`): delta / cluster / coords / stderr and the temporary files equal; the `prepro` examples (`A*AAAAA*A`, `-m 5`, `-m -3`), the `compdna [-1]` rule | as specified |
| A3-9d | `MummerPromer_Oracle_Tests.DeltaPrograms_OnPromerDeltaFiles_...` | 105 PROMER delta cases (`promer_show_oracle.jsonl.gz`): `delta-filter` 210, `show-coords` 315, `show-aligns` 210, `show-snps` 315, `show-diff` 210, `show-tiling` 234 runs equal | as specified |
| A3-9e | `MummerGaps_Oracle_Tests` (4) | `gaps` 259 + the crash, `delta2blocks` 226 + 20 crashes, `mummer-annotate` 230 (+ 14 undefined, valgrind), `run-mummer1` 55 + 60 failing (`mummer_gaps_oracle.jsonl.gz`) | as specified |
| A3-9e | `MummerDelta2Maf_Oracle_Tests` (2) | 375 `delta2maf` runs (`delta2maf_oracle.jsonl.gz`) | as specified |
| A3-9e / A3-9f | `MummerPerl_Oracle_Tests` (2) | 260 `nucmer2xfig` runs; 170 `nucmer -p` runs: files and stderr equal on 123 | as specified |
| A3-9d .. A3-9f | MCP `MummerProgramToolsTests` (9), `ToolRegistrationContractTests`, `ToolDocsContractTests` | the installed programs on the small examples (`mummer`, `promer`, `show_aligns` / `show_snps` on PROMER, `mummer_gaps`, `delta2blocks`, `mummer_annotate`, `run_mummer1`, `delta2maf`, `nucmer2xfig`, `nucmer` prefix) | as specified |

## Audit A3-9e2 / A3-9e3 / A3-9e4 (B09 F210 - F213): `combineMUMs`, `run-mummer3`, `repeat-match`, `exact-tandems`, `mapview`, `mummerplot`

| Item | Tests | Oracle | Status |
|---|---|---|---|
| A3-9e2 | `MummerCombineMums_Oracle_Tests.CombineMums_AndRunMummer3_EqualTheInstalledPrograms_AllRows` | `combinemums_oracle.jsonl.gz` (`generate_combinemums_oracle.py`): `combineMUMs` >= 200 and `run-mummer3` >= 200 runs equal (files, stdout, stderr), the failing runs fail, the undefined runs are skipped | as specified |
| A3-9e3 | `MummerRepeatMatch_Oracle_Tests.RepeatMatch_AndExactTandems_EqualTheInstalledPrograms_AllRows`, `..._ListsTheMatchesInTheOrderOfTheSuffixTreeOfTheProgram`, `ExactTandems_TheInstalledScriptIsBroken_...` | `repeatmatch_oracle.jsonl.gz` (`generate_repeatmatch_oracle.py`): `repeat-match` >= 200, `exact-tandems` installed >= 200 and intended >= 200 runs equal | as specified |
| A3-9e4 | `MummerMapview_Oracle_Tests` (2) | `mapview_oracle.jsonl.gz` (`generate_mapview_oracle.py`): >= 200 runs equal incl. >= 30 `-f pdf / ps` | as specified |
| A3-9e4 | `MummerPlot_Oracle_Tests` (3) | `mummerplot_oracle.jsonl.gz` (`generate_mummerplot_oracle.py`): >= 200 runs equal incl. the dies of the script | as specified |
| A3-9e2 .. A3-9e4 | MCP `MummerProgramToolsTests` (+ 8), `ToolRegistrationContractTests`, `ToolDocsContractTests` | the installed programs on the small examples (`combine_mums`, `run_mummer3`, `repeat_match`, `exact_tandems`, `mapview`, `mummerplot`) | as specified |


## Audit A2-27c (B09 F224): Biopython coordinate formats bed / psl / chain

Tests (`Seqeron.Genomics.Tests`, against Biopython 1.88): `AlignmentCoordinateFormats_Tests` x3 - `biopython_coord_oracle.tsv` (957 rows from `scripts/alignment/generate_biopython_coord_oracle.py`: writer rows with every `bedN`, PSL `header` / `mask` / `wildcard`, strands of query and target, protein / translated DNA, known and unknown residues, two alignments per file; reader rows over the writers' own text, PSLX, comments, blank lines between chain blocks, 13-column chains, 10 / 11-column BED and other texts Biopython refuses): writers byte for byte (390 / 390), readers equal to `Bio.Align.parse` (567 / 567); MCP `CoordinateAlignmentTests` x3.
