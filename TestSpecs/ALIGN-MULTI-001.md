# TestSpec: ALIGN-MULTI-001 - Multiple Sequence Alignment

## Test Unit Information
- **ID:** ALIGN-MULTI-001
- **Area:** Alignment
- **Canonical Method:** `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)`
- **Complexity:** O(n² × m)
- **Algorithm:** Star alignment (simplified progressive alignment)

---

## Evidence Sources

| Source | Type | URL |
|--------|------|-----|
| Wikipedia - Multiple Sequence Alignment | Primary | https://en.wikipedia.org/wiki/Multiple_sequence_alignment |
| Wikipedia - Clustal | Primary | https://en.wikipedia.org/wiki/Clustal |

---

## Test Categories

### MUST Tests (Required - Evidence-backed)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| M01 | `MultipleAlign_NullInput_ThrowsArgumentNullException` | .NET null-safety convention | .NET convention |
| M02 | `MultipleAlign_EmptyCollection_ReturnsEmpty` | Edge case: empty input | Wikipedia MSA |
| M03 | `MultipleAlign_SingleSequence_ReturnsSameSequence` | Trivial case handling | Wikipedia MSA |
| M04 | `MultipleAlign_TwoSequences_AlignsCorrectly` | Minimum non-trivial case | Wikipedia MSA |
| M05 | `MultipleAlign_ThreeIdenticalSequences_AllMatch` | Perfect alignment scenario | Wikipedia MSA |
| M06 | `MultipleAlign_AllAlignedSequences_HaveEqualLength` | MSA invariant: equal length after alignment | Wikipedia MSA |
| M07 | `MultipleAlign_SequenceCount_Preserved` | Invariant: output count = input count | Wikipedia MSA |
| M08 | `MultipleAlign_DifferentLengths_PadsWithGaps` | Gap insertion for length normalization | Wikipedia MSA |
| M09 | `MultipleAlign_ConsensusContainsOnlyValidCharacters` | Consensus validity | Implementation |
| M10 | `MultipleAlign_TotalScore_IsSumOfPairwiseScores` | Sum-of-pairs scoring | Wikipedia MSA |
| M11 | `MultipleAlign_RemovingGaps_RecoversOriginal` | Reversibility invariant | Wikipedia MSA |

### SHOULD Tests (Recommended - Quality/Robustness)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| S01 | `MultipleAlign_ConsensusReflectsMajority` | Verify majority voting logic | Implementation |
| S02 | `MultipleAlign_WithCustomScoring_UsesProvidedMatrix` | Scoring matrix parameter honored | API contract |
| S03 | `MultipleAlign_PartiallyOverlapping_AlignsCorrectly` | Realistic biological scenario | ASSUMPTION |

### COULD Tests (Optional - Extended coverage)

| ID | Test Name | Rationale | Source |
|----|-----------|-----------|--------|
| C01 | `MultipleAlign_LargeSequenceSet_CompletesInReasonableTime` | Performance sanity | ASSUMPTION |
| C02 | `MultipleAlign_MixedCase_NormalizesInput` | Case-insensitivity | ASSUMPTION |

---

## Invariants Under Test

1. **Equal Length:** `result.AlignedSequences.All(s => s.Length == result.AlignedSequences[0].Length)`
2. **Count Preservation:** `result.AlignedSequences.Length == inputSequences.Count()`
3. **Consensus Validity:** `result.Consensus.All(c => "ACGT-".Contains(c))`
4. **Reversibility:** Removing gaps from aligned sequence recovers original
5. **Non-empty consensus for non-empty input:** `inputSequences.Any() => !string.IsNullOrEmpty(result.Consensus)`

---

## Existing Test Coverage Audit

### Location: `SequenceAlignerTests.cs`

| Test Name | Classification | Action |
|-----------|---------------|--------|
| `MultipleAlign_TwoSequences_Aligns` | Weak | Migrate and strengthen |
| `MultipleAlign_ThreeSequences_CreatesConsensus` | Weak | Migrate and strengthen |
| `MultipleAlign_DifferentLengths_PadsWithGaps` | Covered | Migrate |
| `MultipleAlign_SingleSequence_ReturnsSame` | Covered | Migrate |
| `MultipleAlign_Empty_ReturnsEmpty` | Covered | Migrate |
| `MultipleAlign_ReturnsTotalScore` | Weak | Migrate and strengthen |
| `MultipleAlign_NullSequences_ThrowsException` | Covered | Migrate |

### Consolidation Plan

1. **Create** `SequenceAligner_MultipleAlign_Tests.cs` as canonical test file
2. **Migrate** all Multiple Alignment tests from `SequenceAlignerTests.cs`
3. **Remove** migrated tests from `SequenceAlignerTests.cs`
4. **Add** missing MUST tests (M06, M07, M09, M10, M11)
5. **Strengthen** weak tests with proper assertions

---

## Open Questions / Decisions

1. **Q:** Should consensus use gap character when all positions are gaps?  
   **A:** Per Wikipedia, no column should be all-gaps, so this shouldn't occur in valid output.

2. **Q:** Is the choice of first sequence as reference documented?  
   **A:** Implementation-specific; documented in algorithm doc as limitation.

---

## Test Data

### Canonical Test Sequences

```csharp
// Identical sequences
["ATGC", "ATGC", "ATGC"] → Consensus: "ATGC"

// Different lengths  
["ATGCATGC", "ATGC", "ATGCAA"] → All aligned to max length

// Majority voting test
["ATGC", "ATGC", "ATCC"] → Consensus position 3: G (2 vs 1)
```

---

## Definition of Done

- [x] Evidence document created
- [x] Algorithm documentation created
- [x] Tests implemented in canonical file
- [x] All MUST tests passing
- [x] Existing tests migrated and consolidated
- [x] Zero warnings
- [x] Checklist updated
