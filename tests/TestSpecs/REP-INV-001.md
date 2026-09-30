# Test Specification: REP-INV-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | REP-INV-001 |
| **Area** | Repeats |
| **Title** | Inverted Repeat Detection |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-09-30 |

---

## Methods Under Test

| Method | Class | Type | Test Priority |
|--------|-------|------|---------------|
| `FindInvertedRepeats(DnaSequence, minArmLength, maxLoopLength, minLoopLength, maxMismatches = 0, maxArmLength = int.MaxValue, allowWobble = false)` | RepeatFinder | Canonical | Deep testing |
| `FindInvertedRepeats(string, …same options…)` | RepeatFinder | Overload | Deep testing (options) |
| `FindInvertedRepeatsScored(string / DnaSequence, gapPenalty = 12, threshold = 50, matchScore = 3, mismatchScore = −4, maxRepeatLength = 2000)` | RepeatFinder | Variant (EMBOSS einverted) | Deep testing |
| `FindInvertedRepeats(string, minLength, minSpacing, maxSpacing)` | RnaSecondaryStructure | Alternative (RNA) | Smoke testing |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia - Inverted repeat](https://en.wikipedia.org/wiki/Inverted_repeat) | Definition | Sequence followed by its reverse complement with intervening nucleotides |
| [Wikipedia - Stem-loop (Hairpin)](https://en.wikipedia.org/wiki/Stem-loop) | Structure | Stem forms via intramolecular base pairing; optimal loop 4-8 bases |
| [Wikipedia - Palindromic sequence](https://en.wikipedia.org/wiki/Palindromic_sequence) | Edge case | Inverted repeat with zero intervening nucleotides |
| [EMBOSS einverted](https://emboss.sourceforge.net/apps/cvs/emboss/apps/einverted.html) (`einverted.c`, einverted.acd, EMBOSS 6.6.0 binary) | Reference tool (scored variant) | DP with match 3 / mismatch −4 / gap 12 / threshold 50 / maxrepeat 2000; only a/c/g/t score as match; implemented as `FindInvertedRepeatsScored` |
| [EMBOSS palindrome](https://emboss.sourceforge.net/apps/cvs/emboss/apps/palindrome.html) (`palindrome.c`, EMBOSS 6.6.0 binary) | Reference tool | Maximal stems (`-overlap Y`), `-nummismatches` (inward walk, trailing mismatches trimmed), `-maxpallen` (start-span bound + print filter); expected values produced with the binary |
| Crick 1966 (J Mol Biol 19:548); Varani & McClain 2000 (EMBO Rep 1:18) | Base pairing | G·U wobble pairs in RNA helices (`allowWobble`, via `RnaSecondaryStructure.CanPair`) |
| Pearson et al. (1996) | Review | Significance for DNA replication; cruciform formation |
| Bissler (1998) | Review | DNA inverted repeats and human disease |

---

## Test Categories

### MUST Tests (Required for DoD)

All MUST tests are justified by evidence or explicitly marked.

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| M1 | SimpleHairpin_FindsRepeat | Core algorithm - detects stem-loop structure | Wikipedia - inverted repeat |
| M2 | PalindromeSequence_SelfComplementary | Self-complementary arms; exactly one (maximal) stem (0,10,6) | Wikipedia - palindromic sequence; EMBOSS palindrome |
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
| C3 | OverlappingRepeats_AllReported | Overlapping but non-nested stems all reported | EMBOSS palindrome |
| C4 | SubStemsOfPerfectStem_NotReported | Sub-stems of a perfect stem not reported | EMBOSS palindrome |
| C5 | SlippedPairingInsideStem_NotReported | Slipped re-pairings inside a stem not reported | EMBOSS palindrome |
| C6 | InwardExtension_StopsAtMinLoop | Inward extension limited by minLoopLength | EMBOSS palindrome (minLoop 0) + definition |
| C7 | EmbossPalindromeWorkedExample | Multi-stem EMBOSS output locked | EMBOSS palindrome |
| C8 | AmbiguousAndNonDnaBases_NeverPair | N/IUPAC/U never pair | EMBOSS einverted (a/c/g/t only) |
| C9 | HugeMaxLoop_Terminates / InvalidParameters_ThrowEagerly | Degenerate params: no hang; eager validation incl. maxLoop < minLoop | Contract |

### Options and scored variant (added 2026-09-30, `RepeatFinder_InvertedRepeatOptions_Tests`)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| O1 | MaxMismatches_EmbossPalindromeWorkedExample | k = 0/1/2 outputs locked | EMBOSS palindrome binary |
| O2 | MaxMismatches_TrailingMismatchTrimmed_InteriorCounted | trim rule; `Mismatches` count | palindrome.c `mismatchAtEnd` |
| O3 | MaxMismatches_EndsPair_MismatchCountBounded / _Zero_IdenticalToDefault | invariants; default preserved | palindrome.c |
| O4 | MaxArmLength_LongerStemNotSplit / _EmbossWorkedExample | `-maxpallen` rule locked | EMBOSS palindrome binary |
| O5 | AllowWobble_GuPairCompletesStem / _GtOnlyStem / _MatchesAllStemsBruteForce | G·U pairing under the maximal-stem rule | Crick 1966; Varani & McClain 2000 |
| O6 | AllOptions_MatchLiteralPalindromeOracle | every option vs literal palindrome.c transcription (400 cases) | palindrome.c |
| O7 | Scored_EinvertedWorkedExample_Defaults / _CustomGap / _EinvertedPerfectRepeat | coordinates, score, counts, alignment rows locked | EMBOSS einverted binary |
| O8 | Scored_DefaultsAreEinverted / _NeverPairsN_CaseInsensitive / _ResultInvariants / _InvalidParameters_ThrowEagerly | defaults, alphabet, invariants, validation | einverted.acd, einverted.c |

---

## Open Questions / Decisions

| Question | Decision | Justification |
|----------|----------|---------------|
| Should loop=0 be valid? | Allowed when minLoopLength=0; filtered by default (minLoopLength=3) | Wikipedia: "intervening sequence...can be any length including zero." Palindromes are inverted repeats with loop=0 |
| Report overlapping structures? | Only maximal stems (revised 2026-09-29) | EMBOSS palindrome `-overlap Y`: a stem contained in both arms of another stem is dropped; overlapping non-nested stems are all reported. The earlier "all tuples" rule reported every sub-stem of a stem (e.g. 6 hits for GAATTC·AAAA·GAATTC vs EMBOSS's 1) |
| N / ambiguity codes? | Never pair | EMBOSS einverted scores a match only for a/c/g/t |
| Mismatch/gap-tolerant search? | Implemented (2026-09-30): `maxMismatches` (= palindrome `-nummismatches`), `maxArmLength` (= `-maxpallen`), `allowWobble`, and `FindInvertedRepeatsScored` (= einverted) | Cross-checked against the EMBOSS 6.6.0 binaries with 0 mismatches (Evidence) |
| Case sensitivity? | Case-insensitive | Implementation uses ToUpperInvariant() |

---

## Definition of Done Checklist

- [x] All MUST tests implemented
- [x] Tests pass deterministically
- [x] No duplicate tests across files
- [x] Edge cases covered (empty, boundary, invalid)
- [x] Invariants verified with Assert.Multiple
- [x] Clean Code principles applied
