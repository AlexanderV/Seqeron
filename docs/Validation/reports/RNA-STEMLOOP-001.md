# Validation Report: RNA-STEMLOOP-001 — RNA Stem-Loop / Hairpin Detection

- **Validated:** 2026-06-12   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.FindStemLoops(rnaSequence, minStemLength=3, minLoopSize=3, maxLoopSize=10, allowWobble=true)`; supporting `CanPair`/`GetBasePairType`, `DetectPseudoknots(basePairs)`, `GenerateDotBracket`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "Stem-loop"** (https://en.wikipedia.org/wiki/Stem-loop): a stem-loop/hairpin = two complementary regions of the same strand base-pairing into a double helix ending in a loop of unpaired nucleotides. Confirms: "loops fewer than three bases long are sterically impossible and thus do not form"; optimal loop length 4–8 bases; G-C pairs (3 H-bonds) more stable than A-U (2 H-bonds, RNA).
- **Wikipedia "Tetraloop"** (https://en.wikipedia.org/wiki/Tetraloop): tetraloop = 4-base hairpin loop. **GNRA** = G (5') + any N + purine R (G/A) + A (3'). **UNCG** = U + any N + C + G. "The UUCG tetraloop is the most stable tetraloop." UUCG and GNRA together make up 70% of all tetraloops in 16S rRNA (Woese 1990).
- **Wikipedia "Pseudoknot"** (https://en.wikipedia.org/wiki/Pseudoknot): pseudoknot = base pairs that "overlap" in sequence position (crossing pattern, pairs (i,j),(k,l) with i<k<j<l). Confirms standard secondary-structure prediction (Nussinov/DP) does **not** predict pseudoknots — general problem is NP-complete.
- Heus & Pardi (1991, Science 253:191-194) and Woese (1990, PNAS 87:8467) corroborate GNRA stability and prevalence as cited in the Evidence doc.

### Formula / definitions check
- **Stem pairing rules:** A-U, G-C (Watson-Crick) and G-U (wobble) — matches Wikipedia/IUPAC. Wobble gated by `allowWobble`.
- **Min loop size:** ≥3 nt (steric constraint) — matches spec default `minLoopSize=3`, clamped up if a smaller value is requested.
- **Tetraloop:** loop size exactly 4 nt; GNRA = G-N-R-A. Matches spec.
- **Pseudoknot:** crossing-pairs definition i<k<j<l — matches spec; spec correctly treats it as detection from a base-pair list (`DetectPseudoknots`), not prediction from sequence.

### Independent cross-check (hand computation)
- Worked example `GGGAAAACCC` (SL-001/002/003, DB): loop AAAA at indices 3–6; extend G(2)-C(7), G(1)-C(8), G(0)-C(9) → **stem length 3, loop 4 (AAAA)**, dot-bracket `(((....)))`. Matches.
- Tetraloop `GGGCGAAAGCCC` (TL-001/002): loop GAAA at indices 4–7; extend C(3)-G(8), G(2)-C(9), G(1)-C(10), G(0)-C(11) → **stem 4, loop GAAA (size 4)**, GNRA (G,N=A,R=A,A). Matches.
- Pseudoknot pairs (0,6),(3,9): 0<3<6<9 → crossing → pseudoknot. (0,9),(3,6): nested → none. Matches.

### Findings / divergences
None. The "GGGG-AAAA-CCCC" style worked example yields stem 4 + tetraloop 4 exactly as the description predicts.

## Stage B — Implementation

### Code path reviewed
- `RnaSecondaryStructure.cs:412-419` — `BuildPairLookup`: A-U/U-A/G-C/C-G = code 1 (WC), G-U/U-G = code 2 (Wobble). No T (RNA-only) — correct.
- `:463-496` `FindStemLoops` — clamps `minLoopSize` to ≥3 (`:474`); empty/too-short guard `len < minStem*2 + minLoop` (`:476`); upper-cases input (`:479`, handles EC-003 lowercase); scans loop positions × loop sizes, calls `TryBuildStemLoop`.
- `:498-565` `TryBuildStemLoop` — extends stem outward from loop, stops at first non-pair, rejects wobble when `!allowWobble`, returns null if `stemLength < minStemLength`. Loop typed `LoopType.Hairpin`.
- `:1537-1566` `DetectPseudoknots` — exact crossing test `i1<i2 && i2<j1 && j1<j2`.
- `:1575-1590` `GenerateDotBracket` — `(`/`)`/`.` standard notation.

### Formula realised correctly?
Yes. Pairing rules, min-stem, loop bounds, wobble gating, tetraloop loop extraction, and pseudoknot crossing all match the Stage-A-validated description. Hand traces above reproduce exactly against the code.

### Cross-verification recomputed vs code
| Case | Expected | Code result | Status |
|------|----------|-------------|--------|
| `GGGAAAACCC` | stem 3, loop AAAA, `(((....)))` | identical (test asserts exact) | ✅ |
| `GGGCGAAAGCCC` | stem 4, loop GAAA | identical | ✅ |
| `AAAAAAAAA` | empty | empty | ✅ |
| `GCAUC` (too short) | empty | empty | ✅ |
| minLoopSize=1 | clamped to 3 | clamped | ✅ |
| wobble on/off | G-U included/excluded | matches | ✅ |
| pairs (0,6),(3,9) | pseudoknot | detected | ✅ |

### Variant/delegate consistency
Spec lists `FindHairpins`/`FindPseudoknots(sequence)` as Should variants; the spec's own audit (§4, §8) records the decision that hairpin ≡ stem-loop (no separate method) and that pseudoknot detection takes a base-pair list (`DetectPseudoknots`). Implementation matches that decision — no inconsistency.

### Test quality audit
`RnaSecondaryStructureTests.cs` — 137 tests, all exact-value assertions (positions, stem coords, loop sequence, dot-bracket, wobble pair identity), covering all Must/Should/edge cases (empty, null, lowercase, min-loop clamp, min-stem, loop range, tetraloop, wobble in/out, pseudoknot crossing/nested). Not fitted to implementation; derived from RNA biology.

## Verdict & follow-ups
Stage A PASS, Stage B PASS. No defects. Tests: 137 RnaSecondaryStructure tests pass; full suite 4486 passed / 0 failed (baseline). No code changed. **STATE: CLEAN.**
