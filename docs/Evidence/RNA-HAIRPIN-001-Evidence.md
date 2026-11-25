# Evidence Artifact: RNA-HAIRPIN-001

**Test Unit ID:** RNA-HAIRPIN-001
**Algorithm:** Hairpin Loop and Stem Free-Energy Calculation (Turner 2004 nearest-neighbor model)
**Date Collected:** 2026-06-14

---

## Online Sources

<!-- The live NNDB webserver (rna.urmc.rochester.edu/NNDB) was down for maintenance on the
     access date; all NNDB pages below were retrieved from the Internet Archive Wayback
     Machine snapshots of the same canonical NNDB URLs, via `curl`. -->

### NNDB — Turner 2004 Hairpin Loop Parameters (rules and formula)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/hairpin.html
**Retrieved via:** `curl` of Wayback snapshot `http://web.archive.org/web/20240709061712/https://rna.urmc.rochester.edu/NNDB/turner04/hairpin.html` (HTTP 200, 34247 bytes)
**Accessed:** 2026-06-14
**Authority rank:** 2 (official parameter-set specification)

**Key Extracted Points:**

1. **Hairpin free-energy formula (loop > 3 nt):** Verbatim — "ΔG°37 hairpin (>3 nucleotides in loop) = ΔG°37 initiation (n) + ΔG°37 (terminal mismatch) + ΔG°37 (UU or GA first mismatch) + ΔG°37 (GG first mismatch) + ΔG°37 (special GU closure) + ΔG°37 penalty (all C loops)".
2. **3-nt hairpin formula:** "ΔG°37 hairpin (3 unpaired nucleotides) = ΔG°37 initiation (3) + ΔG°37 penalty (all C loops)" — three-nt loops do NOT receive a sequence-dependent first-mismatch term; an all-C 3-nt loop receives a stability penalty.
3. **special GU closure scope:** "applied only to hairpins in which a GU closing pair (not UG) is preceded by two Gs".
4. **all-C penalty (>3 nt):** linear "ΔG°37 penalty (all C loops; > 3 unpaired nucleotides) = An + B".
5. **Short loops prohibited:** "The nearest neighbor rules prohibit hairpin loops with fewer than 3 nucleotides."
6. **Length extrapolation (n>9):** "ΔG°37 initiation (n>9) = ΔG°37 initiation (9) + 1.75 RT ln(n/9)".
7. **Special hairpin loops:** "hairpin loop sequences of 3, 4, and 6 nucleotides that have stabilities poorly fit by the model ... are assigned stabilities based on experimental data."

### NNDB — Turner 2004 length-dependent loop initiation (loop.txt)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/loop.txt
**Retrieved via:** `curl` of Wayback snapshot `.../web/20240709061712/.../turner04/loop.txt` (HTTP 200)
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **HAIRPIN initiation column (kcal/mol):** size 3 → 5.4; 4 → 5.6; 5 → 5.7; 6 → 5.4; 7 → 6.0; 8 → 5.5; 9 → 6.4; 10 → 6.5 (extrapolated thereafter).

### NNDB — Turner 2004 Hairpin First-Mismatch bonus/penalty (hairpin-mismatch-parameters.html)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-mismatch-parameters.html
**Retrieved via:** `curl` of Wayback snapshot `.../web/20240709061712/.../turner04/hairpin-mismatch-parameters.html` (HTTP 200)
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **UU or GA first mismatch:** ΔG°37 = −0.9 kcal/mol.
2. **GG first mismatch:** ΔG°37 = −0.8 kcal/mol.
3. **special GU closure:** ΔG°37 = −2.2 kcal/mol.
4. **C3 loop (all-C, 3 nt):** ΔG°37 = +1.5 kcal/mol.
5. **All-C loop linear params:** A = +0.3 kcal/mol/nt, B = +1.6 kcal/mol.

### NNDB — Turner 2004 terminal-mismatch stacking table (tstack.txt)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/tstack.txt
**Retrieved via:** `curl` of Wayback snapshot `.../web/20240709061712/.../turner04/tstack.txt` (HTTP 200)
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Closing A-U, first mismatch X=A / last mismatch Y=A** (block AX/AY): −0.8 kcal/mol.
2. **Closing A-U, X=G / Y=G** (block AX/GY, row G, col G): −0.8 kcal/mol.

### NNDB — Turner 2004 Watson-Crick stacking + AU-end penalty (wc-parameters.html, stack.txt)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html
**Retrieved via:** `curl` of Wayback snapshot `.../web/20240709061712/.../turner04/wc-parameters.html` (HTTP 200)
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Stacking ΔG°37 (kcal/mol):** 5'AA3'/3'UU5' = −0.93; 5'AU3'/3'UA5' = −1.10; 5'UA3'/3'AU5' = −1.33; 5'CU3'/3'GA5' = −2.08; 5'CA3'/3'GU5' = −2.11; 5'GU3'/3'CA5' = −2.24; 5'GA3'/3'CU5' = −2.35; 5'CG3'/3'GC5' = −2.36; 5'GG3'/3'CC5' = −3.26; 5'GC3'/3'CG5' = −3.42.
2. **Per AU end:** ΔG°37 = +0.45 kcal/mol (applied once per AU/UA pair at each end of a helix).
3. **Helix stacking count:** "For helices of P uninterrupted basepairs, there are P-1 stacks of pairs."

### NNDB — Turner 2004 special hairpin loops (triloop.txt, tloop.txt, hexaloop.txt)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/triloop.txt , /tloop.txt , /hexaloop.txt
**Retrieved via:** `curl` of Wayback snapshots `.../web/20240709061712/.../turner04/{triloop,tloop,hexaloop}.txt` (HTTP 200)
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Triloops (closing+3nt+closing):** CAACG = 6.8; GUUAC = 6.9 kcal/mol (total energy, replaces model).
2. **Tetraloops:** e.g. CCUCGG = 2.5; CUACGG = 2.8; CAACGG = 5.5 kcal/mol.
3. **Hexaloops:** ACAGUGUU = 1.8; ACAGUACU = 2.8; ACAGUGCU = 2.9; ACAGUGAU = 3.6 kcal/mol.

### NNDB — Turner 2004 worked Hairpin Example 1 and Example 2

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-1.html , /hairpin-example-2.html
**Retrieved via:** `curl` of Wayback snapshots `.../web/20240709061712/.../turner04/hairpin-example-{1,2}.html` (HTTP 200)
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Example 1 (6-nt loop):** ΔG°37 = (CG followed by AU −2.11) + (AU followed by CG −2.24) + (CG followed by AU −2.11) + (AU end penalty +0.45) + (terminal mismatch AU followed by AA −0.8) + (initiation(6) +5.4) = **−1.4 kcal/mol**. The loop (initiation + terminal mismatch) component = +5.4 − 0.8 = **+4.6**; the helix (3 stacks + AU end) component = −2.11 −2.24 −2.11 +0.45 = **−6.01**.
2. **Example 2 (5-nt loop, GG first mismatch):** ΔG°37 = (−2.11) + (−2.24) + (−2.11) + (AU end +0.45) + (terminal mismatch AU followed by GG −0.8) + (GG first mismatch −0.8) + (initiation(5) +5.7) = **−1.9 kcal/mol**. Loop component = +5.7 − 0.8 − 0.8 = **+4.1**; helix component = **−6.01**.

### NNDB — Hairpin references (primary literature)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-references.html
**Retrieved via:** `curl` of Wayback snapshot `.../web/20240709061712/.../turner04/hairpin-references.html` (HTTP 200)
**Accessed:** 2026-06-14
**Authority rank:** 1 (primary peer-reviewed source for the parameters)

**Key Extracted Points:**

1. **Primary parameter source:** "Mathews, D.H., Disney, M.D., Childs, J.L., Schroeder, S.J., Zuker, M. and Turner, D.H. (2004) Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. Proc. Natl. Acad. Sci. USA, 101, 7287-7292."

---

## Documented Corner Cases and Failure Modes

### From NNDB Turner 2004 hairpin.html

1. **Loops < 3 nt prohibited:** the nearest-neighbor rules forbid hairpin loops with fewer than 3 unpaired nucleotides; no thermodynamic value is defined.
2. **3-nt loops have no first-mismatch term:** only initiation(3) plus (if applicable) the all-C penalty.
3. **Special loops override the model:** experimentally measured tri/tetra/hexaloop totals (keyed by closing pair + loop) replace the additive model calculation.
4. **special GU closure asymmetry:** the −2.2 bonus applies to a `G-U` closing pair (5'G / 3'U) preceded by two Gs, NOT to a `U-G` closing pair.

### From NNDB Turner 2004 wc-parameters.html

5. **Empty / single base pair stem:** a helix of P base pairs contributes P−1 stacks; a stem of 0 pairs contributes no stacking energy.

---

## Test Datasets

### Dataset: NNDB Hairpin Example 1 (6-nt loop)

**Source:** NNDB Turner 2004 hairpin-example-1.html (Mathews et al. 2004)

| Parameter | Value |
|-----------|-------|
| Closing base pair (5'/3') | A / U |
| Loop size (n) | 6 |
| First / last loop base | A / A |
| Initiation(6) | +5.4 kcal/mol |
| Terminal mismatch (A·A on A-U) | −0.8 kcal/mol |
| **Hairpin loop ΔG°37** | **+4.6 kcal/mol** |
| Helix: 3 stacks (CA/GU, AC/UG, CA/GU) + AU end | −2.11 −2.24 −2.11 +0.45 = **−6.01** |
| Total stem-loop ΔG°37 | **−1.4 kcal/mol** |

### Dataset: NNDB Hairpin Example 2 (5-nt loop, GG first mismatch)

**Source:** NNDB Turner 2004 hairpin-example-2.html (Mathews et al. 2004)

| Parameter | Value |
|-----------|-------|
| Closing base pair (5'/3') | A / U |
| Loop size (n) | 5 |
| First / last loop base | G / G |
| Initiation(5) | +5.7 kcal/mol |
| Terminal mismatch (G·G on A-U) | −0.8 kcal/mol |
| GG first-mismatch bonus | −0.8 kcal/mol |
| **Hairpin loop ΔG°37** | **+4.1 kcal/mol** |

### Dataset: Special / penalty hairpins

**Source:** NNDB triloop.txt, tloop.txt, hexaloop.txt, hairpin-mismatch-parameters.html

| Loop (closing+loop+closing) | ΔG°37 | Note |
|-----------------------------|-------|------|
| C **AAC** G (triloop) | 6.8 | special total replaces model |
| C **GAAA** G... → CCUCGG tetraloop | 2.5 | special total |
| all-C 3-nt loop, e.g. closing G-C, loop CCC | init(3) 5.4 + C3 penalty 1.5 = **6.9** | C3 penalty |
| all-C 4-nt loop, closing G-C, loop CCCC | init(4) 5.6 + tm(GCCC −0.7) + A·4+B (0.3·4+1.6=2.8) = **7.7** | linear all-C penalty |

---

## B12 systematic cross-check (2026-11-25): NNDB Turner 2004 vs ViennaRNA 2.7.2

### Sources opened
- NNDB Turner 2004 (live server, fetched 2026-11-25 with a neutral User-Agent): `rna_2004/rna_2004_hairpin_loops.html`
  (formulae §5.1, examples §5.3, extrapolation §5.4), `parameter_tables/turner_2004_rnastructure/turner_2004_loop_dg.txt`,
  `turner_2004_terminal_mismatch_dg.txt`, `turner_2004_triloop_dg.txt`, `turner_2004_tetraloop_dg.txt`, `turner_2004_hexaloop_dg.txt`,
  hairpin first-mismatch bonus page (UU/GA −0.9, GG −0.8, special GU −2.2, C3 +1.5, all-C A = 0.3, B = 1.6).
- Mathews et al. 2004 PNAS (special GU closure: "a 5′ closing G preceded by two G residues").
- ViennaRNA 2.7.2 sdist: `src/ViennaRNA/eval/eval_hairpin.c` (`vrna_E_hairpin`), `params/default.c` (`hairpin37`, `mismatchH37`),
  `misc/rna_turner2004.par` (Triloops/Tetraloops/Hexaloops).

### Rule-by-rule comparison
| Rule | NNDB Turner 2004 | ViennaRNA 2.7.2 default | Implemented |
|---|---|---|---|
| Special loops | 2 triloops, 16 tetraloops, 4 hexaloops (txt files) | same 22 keys, same values (×100) | NNDB = Vienna (22/22) |
| Special loop value | loop total | returned directly (no TerminalAU, no mismatch) | total; engines add no AU/GU penalty |
| Initiation 3..30 | loop.txt HAIRPIN column | `hairpin37` identical 28/28 | loop.txt |
| Initiation > 30 | init(9) + 1.75·RT·ln(n/9) | hairpin[30] + ⌊lxc37·ln(n/30)⌋, lxc37 = 107.856 | NNDB formula (Δ ≤ 0.01) |
| loop.txt 10..30 vs formula | tabulated values; 6 entries (n = 15, 20, 22, 24, 27, 29; formula 6.951, 7.261, 7.364, 7.458, 7.585, 7.662) are 0.1 below round(formula) | same table | table as published |
| n = 3 | init(3) + C3 penalty; no first mismatch | hairpin[3] + TerminalAU (type > 2) | init(3) [+1.5]; engines add +0.45 AU/GU |
| n ≥ 4 first mismatch | tstack + UU/GA −0.9 (not AG) + GG −0.8 | `mismatchH37` (TerminalAU folded in) | tstack + bonuses; engines add +0.45 |
| mismatchH entries | — | 149/150 = NNDB tstack + bonuses + 0.50 (N rows = max); UA·(A,G) = −0.5 vs NNDB −0.6 | NNDB (−1.1 tstack) |
| AU/GU closure penalty | 0.45 (helix end) | 0.50 | 0.45 |
| Special GU closure | −2.2 if G(5′)-U(3′) closing and the two 5′ neighbours of that G are G; n > 3 formula only | not implemented | NNDB (sequence context i−1, i−2) |
| All-C loops | n = 3: +1.5; n > 3: 0.3n + 1.6 | not implemented | NNDB |
| Loops < 3 nt | prohibited | hairpin[0..2] = INF | sentinel 100.0 (scalar); evaluator rejects (fix below) |
| N / unknown bases | — | first mismatch: `mismatchH[·][0][·]` = max; loop body irrelevant | max over A/C/G/U (F9) |

### Oracle method and numbers
Python model built only from the NNDB files above (independent of the C# tables) and ViennaRNA 2.7.2 `fc.eval_hp_loop(i,j)` on
`pre(2) + closing + loop + closing + post(2)` contexts. **Gotcha:** `RNA.params_load()` returns 1 but a later `RNA.fold_compound(seq)`
created *without* an explicit `md` after a default fold compound already exists silently reuses the cached default parameters; create
fold compounds with `RNA.md()` (or call `fc.params_subst()`) and verify the file is in effect (e.g. `GGGAAAACCC (((....)))` −2.10 default
vs −2.02 NNDB-precision file).
- 160,500 hairpins: all loops of 3–7 nt × 6 closing pairs (G-U cases half with a GG 5′ context), 300 random loops per size 8–60, 75,
  100, 150, 200 (⅓ with N, every 50th all-C), 12,000 N-containing 3–6-nt loops, all-C loops 3–40 × 6 pairs × GG/AA context.
- `CalculateHairpinLoopEnergy` vs the independent NNDB model: **0 mismatches / 160,500**.
- Engine hairpin term (`HairpinTerm`, used by MFE fill/traceback, PF and `CalculateFreeEnergy`) vs `eval_hp_loop` with the NNDB-precision
  file (RNA-STRUCT-001 Evidence §12.2) + NNDB-only terms (special GU −2.2, all-C): **0 mismatches** (n ≤ 30 exact; n > 30 ≤ 0.01, ViennaRNA
  truncation). With the default file: 106,913 differ, all explained by TerminalAU 0.50 vs 0.45 (−0.05) and UA·A/G (−0.10).
- `RNA.eval_structure` (NNDB-precision file) vs `CalculateFreeEnergy` on 17 hairpin structures (sizes 3–102, all-C, special GU and its
  non-GG context, N, A-U hexaloop): equal after NNDB-only terms (e.g. `CGGGGAAAGUCCG` −2.60 − 2.2 = −4.80 = NNDB §5.3 example).

### Finding (Stage B)
- **F-H1 `CalculateFreeEnergy` scored a hairpin of < 3 nt as a finite energy.** `CalculateFreeEnergy("GGGAACCC","(((..)))")` = 93.48
  (the 100.0 scalar sentinel plus stacks); `"((()))"` and `"(((.)))"` likewise. NNDB: "The nearest neighbor rules prohibit hairpin loops with
  fewer than 3 nucleotides"; ViennaRNA `eval_structure` = 100000 (INF). Fix: the evaluator rejects such a structure with `ArgumentException`
  (consistent with its rejection of non-canonical pairs); MCP `rna_structure_energy` inherits it. The MFE/PF engines never build such loops.

### Stage A corrections
- Algorithm doc §2.2 claimed init(n > 9) is the formula; the implementation (correctly) uses the NNDB loop.txt table up to 30. §3.3 claimed
  T is not converted (wrong since F10). Both corrected; AU/GU-penalty convention, special-GU context and NNDB-vs-ViennaRNA choices documented.

---

## Assumptions

<!-- All correctness-affecting parameters in scope are source-backed; no assumptions remain. -->

1. **ASSUMPTION: rounding** — NNDB tabulates final ΔG°37 to one decimal; the implementation rounds intermediate sums to two decimals (`Math.Round(.., 2)`). This is a display-precision choice that does not change any source-defined parameter; tests assert with `.Within(1e-9)` against the exact arithmetic sum of the cited parameters.

---

## Recommendations for Test Coverage

1. **MUST Test:** `CalculateHairpinLoopEnergy` on NNDB Example 1 (closing A-U, 6-nt loop A…A) returns +4.6 — Evidence: hairpin-example-1.html, loop.txt, tstack.txt.
2. **MUST Test:** `CalculateHairpinLoopEnergy` on NNDB Example 2 (closing A-U, 5-nt loop G…G) returns +4.1 (terminal mismatch + GG bonus) — Evidence: hairpin-example-2.html, hairpin-mismatch-parameters.html.
3. **MUST Test:** special triloop CAACG → 6.8 and special tetraloop CCUCGG → 2.5 (special table overrides model) — Evidence: triloop.txt, tloop.txt.
4. **MUST Test:** 3-nt loop has no first-mismatch term: closing G-C, loop "AAA" → init(3)=5.4 — Evidence: hairpin.html (3-nt formula), loop.txt.
5. **MUST Test:** all-C penalty: 3-nt all-C loop adds +1.5; >3-nt all-C loop adds 0.3·n+1.6 — Evidence: hairpin-mismatch-parameters.html.
6. **MUST Test:** loops < 3 nt are prohibited (prohibitive energy) — Evidence: hairpin.html.
7. **MUST Test:** `CalculateStemEnergy` on Example 1 helix (pairs C-G, A-U, C-G, A-U) returns −6.01 (3 stacks + one AU end penalty) — Evidence: wc-parameters.html, hairpin-example-1.html.
8. **MUST Test:** `CalculateStemEnergy` empty base-pair list returns 0 — Evidence: wc-parameters.html (P−1 stacks ⇒ 0 stacks for P≤1).
9. **SHOULD Test:** special GU closure bonus −2.2 applied only for G-U closing (not U-G) — Rationale: documented asymmetry, correctness-affecting.
10. **SHOULD Test:** UU/GA first-mismatch bonus −0.9 — Rationale: distinct additive term.
11. **COULD Test:** initiation extrapolation for n>30 monotonic increase — Rationale: Jacobson-Stockmayer term.

---

## References

1. Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH (2004). Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. Proc. Natl. Acad. Sci. USA 101:7287-7292. https://doi.org/10.1073/pnas.0401799101
2. Turner DH, Mathews DH (2010). NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure. Nucleic Acids Res. 38:D280-D282. https://doi.org/10.1093/nar/gkp892
3. NNDB Turner 2004 Hairpin Loop Parameters. https://rna.urmc.rochester.edu/NNDB/turner04/hairpin.html (accessed 2026-06-14 via Wayback snapshot 20240709061712)

---

## Change History

- **2026-06-14**: Initial documentation (RNA-HAIRPIN-001).
- **2026-11-25**: B12 systematic NNDB / ViennaRNA 2.7.2 cross-check; evaluator rejects hairpins < 3 nt.
