# Test Specification: REP-DIRECT-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | REP-DIRECT-001 |
| **Area** | Repeats |
| **Title** | Direct Repeat Detection |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-09-30 |

---

## Methods Under Test

| Method | Class | Type | Test Priority |
|--------|-------|------|---------------|
| `FindDirectRepeats(DnaSequence, minLength, maxLength, minSpacing)` | RepeatFinder | Canonical | Deep testing |
| `FindDirectRepeats(string, minLength, maxLength, minSpacing)` | RepeatFinder | Overload | Smoke testing |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| MUMmer `repeat-match.cc` (mummer4, raw GitHub; compiled and run) | Reference tool | Reports maximal exact repeats: left-maximal `Data[i-1] != Data[j-1]`, right-maximal (full common prefix); `-f` forward only |
| MUMmer `mummer.cpp` `-n` | Reference tool | "match only the characters a, c, g, or t" |
| Gusfield 1997 §7.12; Abouelhoda et al. 2004 (WebSearch snippets) | Definition | Maximal pair = left- and right-maximal repeated pair |
| [Wikipedia - Direct repeat](https://en.wikipedia.org/wiki/Direct_repeat) | Definition | Sequence repeated with same directionality; may have intervening nucleotides |
| [Wikipedia - Repeated sequence (DNA)](https://en.wikipedia.org/wiki/Repeated_sequence_(DNA)) | Context | Direct vs inverted repeats; types: tandem, interspersed, flanking |
| Ussery et al. (2009) | Technical | Computing for Comparative Microbial Genomics, Springer, Chapter 8 |
| Richard (2021) PMC8145212 | Clinical | Trinucleotide repeat expansions and mismatch repair |

---

## Test Categories

### MUST Tests (Required for DoD)

All MUST tests are justified by evidence or explicitly marked.

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| M1 | SimpleDirectRepeat_FindsRepeat | Core algorithm - detects identical sequences at two positions | Wikipedia - Direct repeat |
| M2 | AdjacentRepeats_WithZeroSpacing_Found | Adjacent repeats (minSpacing=0) should be detected | Wikipedia - tandem direct repeats |
| M3 | NoRepeats_ReturnsEmpty | Sequence without repeated patterns | Standard edge case |
| M4 | EmptySequence_ReturnsEmpty | Boundary - empty input | Standard boundary |
| M5 | NullSequence_ThrowsArgumentNullException | Parameter validation | Implementation contract |
| M6 | MinLengthTooSmall_ThrowsException | minLength < 2 is invalid | Implementation contract |
| M7 | MaxLengthLessThanMinLength_ThrowsException | Invalid parameter combination | Implementation contract |
| M8 | SpacingCalculation_Correct | Spacing = SecondPosition - FirstPosition - Length | Invariant |
| M9 | FirstPosition_LessThanSecondPosition | FirstPosition < SecondPosition always | Invariant |
| M10 | RepeatSequence_MatchesActualSequence | RepeatSequence equals substring at FirstPosition | Invariant |
| M11 | MinLength_RespectsThreshold | `ACGTATTACGTAGGACGTC`: min 5 → (0,7,5); min 4 → (0,7,5),(0,14,4),(7,14,4) | repeat-match -f |
| M12 | MaxLength_RespectsThreshold | ACGTACGTAC+TTTT+ACGTACGTAC: max 8 → (0,18,6),(3,13,7); max 50 adds (0,14,10). A maximal repeat longer than maxLength is not reported (not truncated) | brute force = repeat-match -f |
| M13 | MinSpacing_RespectsThreshold | Only repeats with spacing ≥ minSpacing returned | Algorithm specification |
| M14 | SequenceTooShort_ReturnsEmpty | Sequence shorter than 2×minLength | Boundary condition |

### SHOULD Tests (Important but not blocking)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| S1 | MultipleDirectRepeats_FindsAllPairs | Three copies with distinct flanks (`ACGTATTACGTAGGACGTACC`) → (0,7),(0,14),(7,14) | Gusfield maximal pairs |
| S1b | PeriodicCopies_ReportsMaximalPairsOnly | `ACGTATTACGTATTACGTA`: spacing ≥ 1 → (0,14,5) only; all pairs → (0,7,12,−5),(0,14,5,9) | repeat-match -f |
| S2 | StringOverload_MatchesDnaSequenceOverload | API consistency | Implementation contract |
| S3 | CaseInsensitivity_HandledCorrectly | Lowercase input processed | Implementation robustness |
| S4 | LongSpacing_Detected | Repeats with large intervening region | Wikipedia - interspersed repeats |
| S5 | BiologicalRepeat_TrinucleotideCAG | Test with disease-relevant repeat | Richard (2021) |

### COULD Tests (Nice to have)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| C1 | OverlappingPatterns_OnlyMaximalPairs | `AAAAAATTTTAAAAAA` 4–6 → (0,10,6),(0,11,5),(0,12,4),(1,10,5),(2,10,4) (was 14 nested windows) | repeat-match -f |
| C2 | LargeSequence_Performance | 1000 bp completes in bounded time | DoD requirement |
| C3 | AllMaximalPairs_MatchRepeatMatch | 4 sequences, minSpacing = int.MinValue → exact repeat-match -f lists | repeat-match -f |
| C4 | NegativeMinSpacing_NoSelfPairs | `ACGTACGTACGT` 4,4,−4 → (0,8,4) only; minSpacing int.MaxValue → empty | definition (i < j) |
| C5 | NonAcgt_NeverMatches | N-run → empty; `acgtannnnnacgta` → (0,10,ACGTA); N splits copies | MUMmer `-n` |
| C6 | Output_SortedAndUniquePerPositionPair | sorted by (First, Second), unique pairs, left/right-maximal | definition |

---

## Test Audit

**Canonical file:** `RepeatFinder_DirectRepeat_Tests.cs` — 14 MUST, 6 SHOULD, 6 COULD (C3 is a 4-case TestCase).
Also: `RepeatsDifferentialTests` (brute-force maximal-pair oracle, 5 fixed + 300 random cases), `RepeatsCombinatorialTests`
(27-cell grid, maximality asserted), `RepeatsFuzzTests` (homopolymer exact list), `RepeatFinderProperties`
(INV-6 left/right-maximal + unique), snapshot (4 maximal pairs), MCP `FindDirectRepeatsTests`.

---

## Deviations and Assumptions

- **Reporting convention (2026-09-30, B04 F11):** maximal repeated pairs (Gusfield 1997 §7.12), identical to MUMmer
  `repeat-match -f` output, then filtered by `minLength ≤ L ≤ maxLength` and `Spacing ≥ minSpacing`. The previous
  "every (i, j, len) window" behaviour reported O(L²) nested sub-repeats per repeat and is removed.
- `maxLength` filters: a maximal repeat longer than `maxLength` is not reported (no truncation into sub-windows).
- `minSpacing` may be negative (overlap admitted); pairs always have `i < j` (no self-pairs).
- Only A/C/G/T match (MUMmer `mummer -n`); case-insensitive. `null`/empty string → empty; `null` DnaSequence throws.
- Evidence: `docs/Evidence/REP-DIRECT-001-Evidence.md` (repeat-match cross-check, 0 mismatches).
