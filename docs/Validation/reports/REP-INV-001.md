# Validation Report: REP-INV-001 — Inverted Repeat Detection

- **Validated:** 2026-06-12   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindInvertedRepeats(DnaSequence, minArmLength=4, maxLoopLength=50, minLoopLength=3)`; overload `FindInvertedRepeats(string, …)`; alternative RNA variant `RnaSecondaryStructure.FindInvertedRepeats(string, minLength, minSpacing, maxSpacing)` (smoke).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Inverted repeat** (https://en.wikipedia.org/wiki/Inverted_repeat): "An inverted repeat (or IR) is a single stranded sequence of nucleotides followed downstream by its reverse complement." "The intervening sequence of nucleotides … can be any length including zero." "When the intervening length is zero, the composite sequence is a palindromic sequence." Confirms: IR = sequence + (optional spacer) + reverse complement; spacer ≥ 0; spacer 0 ⇒ palindrome.
- **Wikipedia — Stem-loop** (https://en.wikipedia.org/wiki/Stem-loop): two complementary regions of the same strand base-pair into a stem ending in a loop; "loops that are fewer than three bases long are sterically impossible and thus do not form"; optimal loop "about 4–8 bases long." Confirms the `CanFormHairpin = (loop ≥ 3)` rule.
- **Wikipedia — Palindromic sequence** (https://en.wikipedia.org/wiki/Palindromic_sequence): "a (single-stranded) nucleotide sequence is said to be a palindrome if it is equal to its reverse complement." EcoRI `GAATTC` example given. Confirms self-complementary semantics used in M2/C1.
- **EMBOSS einverted** (emboss.sourceforge.net / bioinformatics.nl manuals): finds stem-loops as local alignments of the sequence vs. its reverse complement exceeding a threshold score; parameters are gap penalty (12), match (3), mismatch (−4), threshold (50); alignments may include mismatches/gaps (bulges). Confirms the *concept* (revcomp matching, arm/loop) and confirms the spec's documented divergence: this implementation uses **exact** revcomp matching with HashSet dedup, not DP scoring, so it reports all exact arm pairs (including overlaps) and does not score imperfect stems.
- **Pearson (1996), Bissler (1998)**: cited for biological significance (cruciform/replication, human disease). Background only; not numeric.

### Definitions & conventions
- Coordinate base: 0-based, half-open arm `[start, start+armLen)`. `LoopLength = RightArmStart − (LeftArmStart + ArmLength)`. `TotalLength = 2·ArmLength + LoopLength`. All standard and internally consistent.
- Parameters: `minArmLength` (min stem length), `minLoopLength`/`maxLoopLength` (spacer/loop bounds). Perfect (exact) revcomp only — no mismatch parameter (legitimate divergence from einverted, documented in spec Open Questions).

### Edge-case semantics (sourced)
- No IR ⇒ empty (homopolymer A: revcomp T, no match). Empty input ⇒ empty.
- Loop = 0 ⇒ palindrome, valid only when `minLoopLength = 0` (Wikipedia: spacer "any length including zero").
- Loop < 3 ⇒ `CanFormHairpin = false` (Stem-loop steric rule). Default `minLoopLength = 3` filters non-hairpins.

### Independent cross-check (hand computation)
- **GAATTC** complement `CTTAAG`, reversed `GAATTC` ⇒ self-complementary (matches Wikipedia EcoRI). 
- **Worked example (M1):** `AAGCGCAAAAGCGCAA`. revcomp(`GCGC`)=`GCGC`. Left arm at 2 (`GCGC`), loop `AAAA` (4 nt), right arm at 10 (`GCGC`). Expected: LeftArmStart 2, RightArmStart 10, ArmLength 4, LoopLength 4, TotalLength 12, CanFormHairpin true. ✓
- **Worked palindrome (S4):** `GGGCGCCC`, revcomp(`GGGC`)=`GCCC`. With `minLoopLength=0`: left at 0, right at 4, loop 0, CanFormHairpin false; with `minLoopLength=1`: empty. ✓

### Findings / divergences
PASS. Exact-match (perfect-stem) policy vs einverted's scored/imperfect alignment is an explicit, sourced design choice (spec Open Questions). No biological/mathematical error in the description.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:275-352` (public overloads + `FindInvertedRepeatsCore`); result type `InvertedRepeatResult` at `:576-590` (`TotalLength => 2*ArmLength + LoopLength`).

### Formula realised correctly?
- revcomp matching via `DnaSequence.GetReverseComplementString(leftArm)` compared `==` to `rightArm` — exact, perfect-stem. ✓
- `LoopLength = j − (i + armLen)`, `loop` substring extracted from `[i+armLen, j)`. ✓
- `CanFormHairpin: loopLength >= 3` — matches Stem-loop steric minimum. ✓
- Spacer/loop window: `minJ = i + armLen + minLoopLength`, `maxJ = min(i + armLen + maxLoopLength, seq.Length − armLen)` — enforces both min and max loop and keeps right arm in bounds. ✓
- Outer bound `i <= seq.Length − 2·minArmLength − minLoopLength` is a correct *necessary* lower bound: any valid match needs `i ≤ seq.Length − 2·armLen − minLoopLength ≤ seq.Length − 2·minArmLength − minLoopLength`, so no valid left start is skipped (verified algebraically). No off-by-one.
- Validation: `minArmLength < 2`, `minLoopLength < 0`, and null `DnaSequence` throw; empty string ⇒ empty. ✓
- Case handling: string overload `ToUpperInvariant()`; `DnaSequence` already normalised. ✓ (S3)
- Dedup: `HashSet<(i,j,armLen)>` — all distinct exact pairs reported including overlaps (C3). ✓

### Cross-verification table recomputed vs code (via tests)
| Case | Input | Expected | Result |
|------|-------|----------|--------|
| M1 hairpin | `AAGCGCAAAAGCGCAA` | 1 IR, L2/R10/arm4/loop4/total12 | PASS |
| M2/C1 EcoRI | `GAATTCAAAAGAATTC` | 6 IRs incl. arm6 L0/R10/loop4 | PASS |
| M4 no-IR | `AAAAAAAAAAAAAA` | empty | PASS |
| M7 minLoop | `GCGCAAGCGC` minLoop3→∅, minLoop1→loop2 | PASS |
| M8 maxLoop | `GCGC`+15A+`GCGC` maxLoop10→∅, 20→loop15 | PASS |
| M13 too-short | `GCGCAAGCGC` (10 nt) | empty | PASS |
| S4 palindrome loop0 | `GGGCGCCC` minLoop0→1, minLoop1→∅ | PASS |
| C3 overlapping | `GCGCAAAGCGCAAAGCGC` | 3 pairs (0-7,0-14,7-14) | PASS |

### Variant/delegate consistency
String overload mirrors `DnaSequence` overload (S2 verifies equality). RNA variant (`RnaSecondaryStructure.FindInvertedRepeats`, different signature) is a separate smoke-level alternative; not under deep test here.

### Test quality audit
43 inverted-repeat tests assert exact sourced positions/lengths/loop/total/hairpin flags and parameter thresholds (not "no-throw" tautologies); deterministic; cover empty, too-short, no-IR, loop=0 palindrome, min/max loop boundaries, min arm, overlapping, case, null/range validation. Edge cases from Stage A are all covered.

### Findings / defects
None. Code faithfully realises the validated definition.

## Verdict & follow-ups
- Stage A: PASS. Stage B: PASS. State: CLEAN — no defects; 43/43 inverted tests pass, full suite 4461/4461.
- Note (not a defect): perfect/exact-stem only; imperfect stems with mismatches/bulges (einverted-style scoring) are out of scope by design and documented in the TestSpec Open Questions.
