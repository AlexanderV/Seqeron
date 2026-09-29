# Evidence Artifact: REP-INV-001

**Test Unit ID:** REP-INV-001
**Algorithm:** Exact (perfect-stem) inverted-repeat detection — maximal stems (EMBOSS `palindrome` reporting rule, 0 mismatches) with a minimum/maximum loop window
**Date Collected:** 2026-09-29 (campaign 2026-09, batch B04)

---

## Sources opened

### EMBOSS palindrome — source code (`emboss/palindrome.c`, M. Faller)

**URL:** https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/palindrome.c (fetched, 482 lines)
**Authority rank:** 2 (original tool source code)

Key points (read from the code):

1. "Brute force inverted repeat finder. Allows mismatches but not gaps." (file header).
2. For every left position `current` and right position `rev` (descending) whose bases complement, the stem
   is extended **inward** (`ic++`, `ir--`) while `mismatches <= maxmismatches && ic < ir`; trailing mismatches
   are removed (`count -= mismatchAtEnd`); `gap = rev - current - 2*count + 1`.
3. A stem is kept when `count >= minLen && gap <= maxGap && !alln` (an all-`n` stem is rejected).
4. With `-overlap Y` (default) a new stem is dropped when `palindrome_AInB(new, old)` — "Palindrome A is within
   Palindrome B (is a subset of B) in both halves of the stem" — for any earlier stem. This is the reporting rule
   implemented by `RepeatFinder.FindInvertedRepeats`.
5. Complementarity uses `ajBaseAlphacharComp` (IUPAC complement; `n` complements `n`).

### EMBOSS palindrome — application documentation

**URL:** https://emboss.sourceforge.net/apps/cvs/emboss/apps/palindrome.html (fetched)

"palindrome finds inverted repeats (stem loops) in nucleotide sequences … It finds all possible inverted matches
satisfying the specified conditions of minimum and maximum length of palindrome, maximum gap between repeated
regions and number of mismatches allowed." Qualifiers `-minpallen`, `-maxpallen`, `-gaplimit`, `-nummismatches`
(default 0), `-[no]overlap` (default Y, "Report overlapping matches").

### EMBOSS einverted — source code (`emboss/einverted.c`, Durbin & Thierry-Mieg 1993)

**URL:** https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/einverted.c (fetched, 707 lines)

1. "Inverted repeats by dynamic programming" — local alignment of the sequence against its reverse complement with
   `match`/`mismatch`/`gap` scores and a `threshold`; mismatches and gaps (bulges) are allowed.
2. The match table (`revmatch`) is filled only for the four bases a/c/g/t (`base[] = "acgt-"`): any other code,
   including `n`, always scores `mismatch`. Used as the source for "N never pairs".
3. Conclusion: einverted is a *scored, imperfect* inverted-repeat search. The unit spec does not promise einverted
   scoring; the unit is an exact-stem finder and is declared as such (EMBOSS `palindrome` with 0 mismatches is the
   matching reference tool).

### Wikipedia (Inverted repeat, Stem-loop, Palindromic sequence)

Not re-opened in this session (blocked); definitions recorded in the prior report
`docs/Validation/reports/REP-INV-001.md` (IR = sequence followed downstream by its reverse complement, spacer of any
length including zero; loops < 3 nt sterically impossible).

---

## Reference implementations run

- **EMBOSS 6.6.0 binaries** (Ubuntu package `emboss 6.6.0+dfsg-12ubuntu2`): `palindrome -nummismatches 0 -overlap Y
  -minpallen <minArm> -maxpallen <len> -gaplimit <maxLoop>`.
  Note: EMBOSS ACD clamps `-minpallen`/`-maxpallen` to `len/2` ("integer value out of range … reset to"); comparison
  cases use `minArm ≤ len/2`.
- **Independent Python brute force** (campaign scratch `ref.py`): enumerate every exact stem `(i, j, A)` with
  `s[i+k]` Watson–Crick-complementary to `s[j+A−1−k]` (A/C/G/T only), `A ≥ minArm`, `minLoop ≤ j−(i+A) ≤ maxLoop`;
  keep stems not contained in both arms of another stem.

| Cross-check | Cases | Mismatches |
|---|---|---|
| Python brute force vs EMBOSS palindrome (minLoop 0; random ACGT/AT/GC, n 8–60, minArm 2–5, gap 0–12; seeds 1, 7) | 800 | 0 |
| C# `RepeatFinder.FindInvertedRepeats` vs EMBOSS palindrome (minLoop 0; n 8–80, minArm 2–6, gap 0–20) | 1000 | 0 |
| C# vs Python brute force (n 0–45, alphabets incl. N, lowercase, U; minArm 2–5, minLoop 0–6, maxLoop minLoop..+15) | 3000 (1138 non-empty) | 0 |

### Locked worked values (EMBOSS palindrome output → 0-based `(LeftArmStart, RightArmStart, ArmLength)`)

| Sequence | minArm / maxLoop | EMBOSS (minLoop 0) | minLoop 3 (Python ref = C#) |
|---|---|---|---|
| GGGGGAAACCCCC | 4 / 50 | (0,8,5) | (0,8,5) |
| GGGGGGAAACCCCCC | 4 / 50 | (0,9,6) | (0,9,6) |
| GAATTCAAAAGAATTC | 4 / 10 | (0,10,6) | (0,10,6) |
| ACGTAAATTTTTACGT | 4 / 50 | (0,9,7) loop 2 | (0,10,6) loop 4 |
| ACGTACGTAAAACGTACGT | 4 / 50 | (0,4,4) (0,11,8) (4,14,5) (11,15,4) | (0,11,8) (4,14,5) |
| GCGCAAAGCGCAAAGCGC | 4 / 10 | (0,7,4) (0,14,4) (7,14,4) | same |
| TTACGTAAACCCACGTTT | 4 / 50 | (0,4,4) (2,12,4) (6,14,4) | (2,12,4) (6,14,4) |
| TGTC×3·ACGTGC·TTAA·GCACGT·TGTC×3 (40 nt) | 3 / 5 | (11,14,3) (12,20,8) | — |

## Performance (Release harness)

| Input | Params (minArm, maxLoop, minLoop) | Results | Time |
|---|---|---|---|
| random 200 kb | 4, 50, 3 | 28 520 | 0.2 s |
| (AT)×5000 | 2, int.MaxValue, 0 | 9 997 | 0.5 s |
| G×5000·C×5000 | 2, int.MaxValue, 0 | 1 | 0.5 s |
| N×20000 | 2, 50, 0 | 0 | 0.1 s |
| random 5 kb | 2, int.MaxValue, 0 | 566 236 | 13 s (output-bound) |
