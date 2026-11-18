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
