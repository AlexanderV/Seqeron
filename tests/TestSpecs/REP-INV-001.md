# Test Specification: REP-INV-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | REP-INV-001 |
| **Area** | Repeats |
| **Title** | Inverted Repeat Detection |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-01-22 |

---

## Methods Under Test

| Method | Class | Type | Test Priority |
|--------|-------|------|---------------|
| `FindInvertedRepeats(DnaSequence, minArmLength, maxLoopLength, minLoopLength)` | RepeatFinder | Canonical | Deep testing |
| `FindInvertedRepeats(string, minArmLength, maxLoopLength, minLoopLength)` | RepeatFinder | Overload | Smoke testing |
| `FindInvertedRepeats(string, minLength, minSpacing, maxSpacing)` | RnaSecondaryStructure | Alternative (RNA) | Smoke testing |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia - Inverted repeat](https://en.wikipedia.org/wiki/Inverted_repeat) | Definition | Sequence followed by its reverse complement with intervening nucleotides |
| [Wikipedia - Stem-loop (Hairpin)](https://en.wikipedia.org/wiki/Stem-loop) | Structure | Stem forms via intramolecular base pairing; optimal loop 4-8 bases |
| [Wikipedia - Palindromic sequence](https://en.wikipedia.org/wiki/Palindromic_sequence) | Edge case | Inverted repeat with zero intervening nucleotides |
| [EMBOSS einverted](https://emboss.sourceforge.net/apps/cvs/emboss/apps/einverted.html) | Algorithm | Dynamic programming approach; parameters: minArmLength, maxLoopLength |
| Pearson et al. (1996) | Review | Significance for DNA replication; cruciform formation |
| Bissler (1998) | Review | DNA inverted repeats and human disease |

---

## Test Categories

### MUST Tests (Required for DoD)

All MUST tests are justified by evidence or explicitly marked.

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| M1 | SimpleHairpin_FindsRepeat | Core algorithm - detects stem-loop structure | Wikipedia - inverted repeat |
| M2 | PalindromeSequence_SelfComplementary | Self-complementary sequence (revcomp = self) | Wikipedia - palindromic sequence |
| M3 | ReverseComplementMatch_BothArmsCorrect | Left arm revcomp must equal right arm | Wikipedia - inverted repeat definition |
| M4 | NoInvertedRepeats_ReturnsEmpty | Sequence without complementary regions | Standard edge case |
| M5 | EmptySequence_ReturnsEmpty | Boundary - empty input | Standard boundary |
| M6 | MinArmLength_RespectsThreshold | Filter by minimum arm length | Algorithm specification |
| M7 | MinLoopLength_RespectsThreshold | Filter by minimum loop length (spacing) | Algorithm specification |
| M8 | MaxLoopLength_RespectsThreshold | Filter by maximum loop length | Algorithm specification; EMBOSS einverted |
| M9 | LoopLength_CalculatedCorrectly | loopLength = rightArmStart - leftArmEnd | Invariant |
| M10 | TotalLength_InvariantHolds | TotalLength = 2×ArmLength + LoopLength | Invariant |
| M11 | CanFormHairpin_LoopLengthValidation | CanFormHairpin true when loopLength ≥ 3 | Wikipedia - Stem-loop: "loops fewer than three bases long are sterically impossible" |
| M12 | LeftRightArmPositions_Correct | Verify position accuracy (0-based indexing) | Implementation contract |
| M13 | SequenceTooShort_ReturnsEmpty | Sequence shorter than minimum structure | Boundary condition |

### SHOULD Tests (Important but not blocking)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| S1 | MultipleInvertedRepeats_FindsAll | Sequence with multiple hairpin structures | Real genome scenario |
| S2 | StringOverload_MatchesDnaSequenceOverload | API consistency | Implementation contract |
| S3 | CaseInsensitivity_HandledCorrectly | Lowercase input processed | Implementation robustness |
| S4 | AdjacentArms_MinLoopZeroAllowsDetection | Palindromic (loop=0) arms allowed when minLoopLength=0 | Wikipedia - Inverted repeat ("intervening sequence...can be any length including zero") |
| S5 | BiologicalHairpin_KnownStructure | Test with known biological stem-loop | Wikipedia - tRNA |

### COULD Tests (Nice to have)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| C1 | RestrictionSitePalindromes_Detected | EcoRI, BamHI recognition sites | Wikipedia - restriction enzymes |
| C2 | LoopSequence_CorrectlyExtracted | Loop sequence matches intervening nucleotides | Implementation verification |
| C3 | OverlappingRepeats_AllReported | All matches reported including overlapping | EMBOSS einverted (behavioral difference) |

---

## Open Questions / Decisions

| Question | Decision | Justification |
|----------|----------|---------------|
| Should loop=0 be valid? | Allowed when minLoopLength=0; filtered by default (minLoopLength=3) | Wikipedia: "intervening sequence...can be any length including zero." Palindromes are inverted repeats with loop=0 |
| Report overlapping structures? | Yes (all matches) | Our algorithm uses exact matching with HashSet dedup, not DP scoring like EMBOSS einverted. All valid (leftStart, rightStart, armLength) tuples are reported, including those sharing arm positions |
| Case sensitivity? | Case-insensitive | Implementation uses ToUpperInvariant() |

---

## Definition of Done Checklist

- [x] All MUST tests implemented
- [x] Tests pass deterministically
- [x] No duplicate tests across files
- [x] Edge cases covered (empty, boundary, invalid)
- [x] Invariants verified with Assert.Multiple
- [x] Clean Code principles applied
