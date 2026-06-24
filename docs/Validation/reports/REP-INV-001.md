# Validation Report: REP-INV-001 — Inverted Repeat Detection

- **Validated:** 2026-06-24   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindInvertedRepeats(DnaSequence, minArmLength=4, maxLoopLength=50, minLoopLength=3)`; string overload `FindInvertedRepeats(string, …)`; alternative RNA variant `RnaSecondaryStructure.FindInvertedRepeats(string, minLength, minSpacing, maxSpacing)` (smoke).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Inverted repeat** (fetched live): "a single stranded sequence of nucleotides followed downstream by its reverse complement"; the intervening nucleotides "can be any length, including zero"; example `5'-TTACGnnnnnnCGTAA-3'`; "When the intervening sequence has zero length, an inverted repeat becomes a palindromic sequence." Confirms: IR = arm + (optional spacer/loop ≥ 0) + reverse complement of arm; spacer 0 ⇒ palindrome. This is the exact definition the unit must implement, and it is a **reverse-complement** match, not a direct (forward) repeat.
- **Wikipedia — Stem-loop** (fetched live): "Loops that are fewer than three bases long are sterically impossible and thus do not form"; optimal loop ≈ 4–8 bases (UUCG tetraloop example). Confirms the `CanFormHairpin = (loop ≥ 3)` rule and the sensible default `minLoopLength = 3`.
- **Wikipedia — Palindromic sequence** (per TestSpec): palindrome = sequence equal to its own reverse complement (EcoRI `GAATTC`); = IR with spacer 0. Confirms M2/C1 semantics and the REP-PALIN relationship (palindrome is the spacer=0 special case of an inverted repeat).
- **EMBOSS einverted** (per TestSpec): DP local-alignment of sequence vs. its reverse complement above a score threshold; parameters include arm/loop bounds and tolerate mismatches/bulges. Confirms the *concept* and frames the documented divergence: this implementation does **exact (perfect-stem)** revcomp matching with HashSet dedup, not scored/imperfect alignment.
- **Pearson (1996), Bissler (1998)**: biological-significance background only (cruciform/replication, human disease); non-numeric.

### Conventions confirmed
- **Stem (arm) length:** `minArmLength` (default 4); arms grow upward; both arms equal length (perfect stem).
- **Loop/spacer:** `minLoopLength` (default 3) ≤ loop ≤ `maxLoopLength` (default 50); loop=0 permitted only when `minLoopLength=0` (= palindrome).
- **Mismatches:** none allowed — exact reverse-complement only (legitimate, sourced divergence from einverted; TestSpec Open Questions).
- **Coordinates:** 0-based; arm span half-open `[start, start+armLen)`. `LoopLength = RightArmStart − (LeftArmStart + ArmLength)`; `TotalLength = 2·ArmLength + LoopLength`.
- **Difference from REP-PALIN:** palindrome = inverted repeat with spacer = 0; this unit detects the general spacer ≥ minLoopLength case.

### Edge-case semantics (sourced)
- No IR ⇒ empty (homopolymer A → revcomp T, no match). Empty/too-short input ⇒ empty.
- Loop = 0 ⇒ palindrome, valid only with `minLoopLength = 0` (Wikipedia "any length including zero").
- Loop < 3 ⇒ `CanFormHairpin = false` (steric minimum). Default `minLoopLength = 3` filters non-hairpins.

### Independent cross-check (hand computation, this session)
- **Designed spacer case (reverse-complement, NOT direct repeat):** `AACCGGTTTCCGGTT` (15 nt). rc(`AACCGG`) = `CCGGTT`. Left arm `AACCGG` at idx 0; loop `TTT` (3 nt); right arm at idx 9 = `CCGGTT` = rc(left). Note left arm ≠ right arm, so this is a genuine reverse-complement (inverted) match, not a direct repeat. Trace of code: i=0, armLen=6, minJ=9, maxJ=min(56,9)=9, j=9, rightArm `CCGGTT` == leftArmRevComp ✓ ⇒ LeftArmStart=0, RightArmStart=9, ArmLength=6, LoopLength=3, Loop=`TTT`, TotalLength=2·6+3=15, CanFormHairpin=true. All fields match the by-hand result.
- **EcoRI palindrome:** rc(`GAATTC`)=`GAATTC` (self-complementary) — matches Wikipedia EcoRI example; confirms spacer=0 palindrome semantics.

### Findings / divergences
PASS. Definition, formulas, conventions and edge cases all match authoritative sources. Exact-stem (no mismatch) policy vs. einverted's scored/imperfect alignment is an explicit, sourced design choice — not an error.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:275-355` (public overloads + `FindInvertedRepeatsCore`); result struct `InvertedRepeatResult` at `:646-660` (`TotalLength => 2*ArmLength + LoopLength`).

### Formula realised correctly?
- Reverse-complement matching: `DnaSequence.GetReverseComplementString(leftArm)` compared `==` to `rightArm` — exact perfect stem (`:320`,`:332`). ✓
- `LoopLength = j − (i + armLen)`; loop substring `[i+armLen, j)` (`:334-335`). ✓
- `CanFormHairpin: loopLength >= 3` (`:349`) — matches steric minimum. ✓
- Loop window: `minJ = i + armLen + minLoopLength`; `maxJ = min(i + armLen + maxLoopLength, seq.Length − armLen)` (`:323-324`) — enforces min and max loop and keeps right arm in bounds; `j + armLen > seq.Length` break is redundant-safe (`:328`). ✓
- Outer bound `i <= seq.Length − 2·minArmLength − minLoopLength` (`:315`) is a correct *necessary* lower bound (a valid match needs `i ≤ seq.Length − 2·armLen − minLoopLength ≤ seq.Length − 2·minArmLength − minLoopLength`) — no valid left start skipped; no off-by-one.
- Validation: `minArmLength < 2`, `minLoopLength < 0` throw; null `DnaSequence` throws; empty string ⇒ empty (`:281-301`). ✓
- Case handling: string overload `ToUpperInvariant()`; `DnaSequence` already normalised (`:303`). ✓ (S3)
- Dedup: `HashSet<(i,j,armLen)>` reports all distinct exact pairs incl. overlaps (C3) (`:313`,`:337-340`). ✓
- `TotalLength = 2·ArmLength + LoopLength` (`:659`). ✓

### Cross-verification (recomputed vs code)
| Case | Input | Expected | Result |
|------|-------|----------|--------|
| Designed spacer (this session) | `AACCGGTTTCCGGTT` | L0/R9/arm6/loop3/total15, hairpin=true | matches hand trace |
| M1 hairpin | `AAGCGCAAAAGCGCAA` | L2/R10/arm4/loop4/total12 | covered by tests |
| M2/C1 EcoRI | `GAATTCAAAAGAATTC` | arm6 L0/R10/loop4 + overlaps | covered |
| No-IR | homopolymer A | empty | covered |
| loop=0 palindrome | `GGGCGCCC` minLoop0→1, minLoop1→∅ | per Wikipedia | covered |

### Variant/delegate consistency
String overload mirrors the `DnaSequence` overload (S2 asserts equality). RNA variant (`RnaSecondaryStructure.FindInvertedRepeats`, different signature/parameters) is a separate smoke-level alternative, not under deep test here.

### Test quality audit
`RepeatFinder_InvertedRepeat_Tests` = 24 canonical tests; all repeat-suite tests = 700; both **0 failed**. Assertions check exact sourced positions/lengths/loop/total/hairpin flags and parameter thresholds (not no-throw tautologies); deterministic; cover empty, too-short, no-IR, loop=0 palindrome, min/max loop boundaries, min arm, overlapping, case, null/range validation. All Stage-A edge cases covered.

### Findings / defects
None. Code faithfully realises the validated definition.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.** No defects; no code changed. Canonical class 24/24 green, full repeat suite 700/700 green.
- Note (not a defect): perfect/exact-stem only; imperfect stems with mismatches/bulges (einverted-style scoring) are out of scope by design and documented in the TestSpec Open Questions.
