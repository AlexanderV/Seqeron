# Validation Report: RNA-STEMLOOP-001 — RNA Stem-Loop / Hairpin Detection

- **Validated:** 2026-06-24   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.FindStemLoops(rnaSequence, minStemLength=3, minLoopSize=3, maxLoopSize=10, allowWobble=true)`; supporting `CanPair` / `GetBasePairType` / `BuildPairLookup`, `TryBuildStemLoop`, `DetectPseudoknots(basePairs)`, `GenerateDotBracket`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (fetched this session, not from memory)
- **Wikipedia "Stem-loop"** (fetched 2026-06-24): a stem-loop / hairpin = "two regions of the same nucleic acid strand, usually complementary in nucleotide sequence, base-pair to form a double helix that ends in a loop of unpaired nucleotides." Confirms verbatim: "loops that are fewer than three bases long are sterically impossible and thus do not form"; "Optimal loop length tends to be about 4–8 bases long"; G-C pairs (3 H-bonds) more stable than A-U (2 H-bonds).
- **Wikipedia "Wobble base pair"** (fetched 2026-06-24): G–U is one of "the four main wobble base pairs"; quoted as "a fundamental building block of RNA structure crucial to RNA function" (Varani & McClain). Confirms G-U is a legitimate, naturally occurring RNA pair — correctly gated by `allowWobble`.
- **Tetraloop / Pseudoknot** definitions carried from the Evidence doc and prior report and re-checked: GNRA = G-N-R(A/G)-A; pseudoknot = crossing pairs (i,j),(k,l) with i<k<j<l. Both consistent with cited sources.

### Formula / definitions check
- **Stem pairing rules:** A-U, G-C (Watson-Crick) and G-U (wobble) — matches Wikipedia/IUPAC. Wobble gated by `allowWobble`. No T (RNA-only) — correct.
- **Min loop size:** ≥3 nt (steric constraint) — matches default `minLoopSize=3`; values below 3 are clamped up.
- **Tetraloop:** loop size exactly 4 nt; GNRA = G-N-R-A — matches spec/test.
- **Pseudoknot:** crossing-pairs definition i<k<j<l; correctly modelled as *detection from a base-pair list* (`DetectPseudoknots`), not prediction from sequence (the latter is NP-complete and out of scope here).
- **Coordinates:** 0-based, inclusive `Start`/`End`; stem reported as 5'/3' arm coords plus length; loop reported as start/end/size/sequence — all explicit and standard.

### Edge-case semantics check
Empty/null → empty (no throw); too short (`len < minStem*2 + minLoop`) → empty; no complement → empty; lowercase normalised to uppercase; loop <3 clamped to 3. All sourced (steric constraint) or conventional.

### Independent cross-check (hand computation)
- **`GGGAAAACCC`**: loop AAAA at idx 3–6; extend G(2)-C(7), G(1)-C(8), G(0)-C(9) → **stem 3, loop 4 (AAAA)**, dot-bracket `(((....)))`. Matches code/test.
- **`GGGCGAAAGCCC`** (tetraloop): loop GAAA at idx 4–7; extend C(3)-G(8), G(2)-C(9), G(1)-C(10), G(0)-C(11) → **stem 4, GNRA loop GAAA**, `((((....))))`. Matches.
- **`GGGGCUUUUGCCCC`** (prompt's designed hairpin, idx 0–13): loop UUUU at idx 5–8; extend C(4)-G(9) WC, G(3)-C(10), G(2)-C(11), G(1)-C(12), G(0)-C(13) → **stem 5, loop UUUU (4 nt)**, all Watson-Crick. Confirms the extend-outward model and that an internal C-G flanked by G-C pairs is handled correctly.
- **Wobble `GGUAAAGCC`**: G(0)-C(8), G(1)-C(7), U(2)-G(6) wobble → 3 bp stem incl. one wobble. Matches.
- **Pseudoknot**: pairs (0,6),(3,9): 0<3<6<9 → crossing → knot; (0,9),(3,6): nested → none. Matches.

### Findings / divergences
None. Description and Evidence are biologically and mathematically correct against re-fetched primary/encyclopedic sources.

## Stage B — Implementation

### Code path reviewed (`src/.../RnaSecondaryStructure.cs`)
- `:449-456` `BuildPairLookup`: A-U/U-A/G-C/C-G = code 1 (WC); G-U/U-G = code 2 (Wobble); no T — RNA-only, correct.
- `:465-486` `CanPair` / `GetBasePairType`: case-insensitive, bounds-guarded (`(b1|b2) < 128`), maps codes 1/2 to WatsonCrick/Wobble.
- `:500-533` `FindStemLoops`: clamps `minLoopSize` to ≥3 (`:511`); empty/too-short guard `len < minStem*2 + minLoop` (`:513`); upper-cases input (`:516`, handles lowercase EC-003); double loop over loop-start × loop-size calling `TryBuildStemLoop`.
- `:535-601` `TryBuildStemLoop`: extends stem outward from the loop, stops at first non-pair, rejects wobble when `!allowWobble`, returns null if `stemLength < minStemLength`; loop typed `LoopType.Hairpin`; reverses base pairs to 5'→3' order.
- `:1976-2014` `DetectPseudoknots`: normalises each pair to (open<close) and orders the two pairs by opening position, then applies exact crossing test `i<k && k<j && j<l`. This is **more robust** than (and equivalent to) the simpler test described in the prior report — it is direction- and storage-order-independent.
- `:2023-2038` `GenerateDotBracket`: standard `(` / `)` / `.` notation.

### Formula realised correctly?
Yes. Pairing rules, min-stem, loop bounds, wobble gating, tetraloop extraction, and pseudoknot crossing all match the Stage-A-validated description. All hand traces above reproduce exactly against the code and the test assertions.

### Cross-verification recomputed vs code
| Case | Expected | Code/test result | Status |
|------|----------|------------------|--------|
| `GGGAAAACCC` | stem 3, loop AAAA, `(((....)))` | identical (exact-value test) | PASS |
| `GGGCGAAAGCCC` | stem 4, loop GAAA, `((((....))))` | identical | PASS |
| `GGGGCUUUUGCCCC` (hand) | stem 5, loop UUUU, all WC | matches model | PASS |
| `AAAAAAAAAAAAAAA` | empty | empty | PASS |
| `GCAUC` (too short) | empty | empty | PASS |
| minLoopSize=1 | clamped to 3 → empty for 2-nt loop seq | clamped | PASS |
| wobble on/off (`GGUAAAGCC` / `GCGAAAACGU`) | G-U incl./excl. | matches | PASS |
| pairs (0,6),(3,9) | pseudoknot | detected | PASS |

### Variant/delegate consistency
Spec lists `FindHairpins` / `FindPseudoknots(sequence)` as Should variants; the spec's own audit (§4, §8) records the design decision that hairpin ≡ stem-loop (no separate method) and that pseudoknot detection takes a base-pair list (`DetectPseudoknots`). Implementation matches that decision — no inconsistency.

### Test quality audit
`RnaSecondaryStructureTests.cs` — stem-loop region tests use exact-value assertions (positions, stem arm coords, loop sequence/size, dot-bracket strings, specific wobble U-G pair identity, min-loop clamp with positive control), covering all Must/Should/edge cases. Values are derived from RNA biology, not fitted to the implementation. Filtered class run: 128 passed / 0 failed.

## Verdict & follow-ups
Stage A PASS, Stage B PASS. No defects found. No code changed. Full `Seqeron.Genomics.Tests` suite: **18213 passed / 0 failed**. **STATE: CLEAN.**

Note: the prior report described `DetectPseudoknots` via a simpler crossing test; the current code normalises pair orientation and orders by opening position. This is a strict improvement (equivalent logic, direction-independent) — not a regression — and is correctly exercised by the PK-001/PK-002 tests.
