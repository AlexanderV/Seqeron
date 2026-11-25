# Evidence Artifact: RNA-MFE-001

**Test Unit ID:** RNA-MFE-001
**Algorithm:** Minimum Free Energy (MFE) RNA secondary structure prediction (Zuker–Stiegler dynamic programming with Turner 2004 nearest-neighbor parameters)
**Date Collected:** 2026-06-14

---

## Online Sources

<!-- Live NNDB pages were retrieved via `curl` of the Internet Archive Wayback Machine
     snapshot 20240709061712 of the canonical NNDB URLs (the live rna.urmc.rochester.edu
     server blocked the fetch tool with HTTP 404; the same approach was used for the
     sibling unit RNA-HAIRPIN-001). The canonical URL is given first; the exact Wayback
     URL fetched is recorded under "Retrieved via". -->

### Zuker & Stiegler (1981) — original MFE dynamic-programming algorithm

**URL:** https://doi.org/10.1093/nar/9.1.133
**Retrieved via:** `WebFetch` of PubMed Central full-text page `https://pmc.ncbi.nlm.nih.gov/articles/PMC326673/` (the open-access archive of the article)
**Accessed:** 2026-06-14
**Authority rank:** 1 (primary peer-reviewed paper)

**Key Extracted Points:**

1. **Citation (verbatim):** "Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information", M Zuker, P Stiegler, *Nucleic Acids Research*, 9(1):133–148, January 10 1981, DOI 10.1093/nar/9.1.133.
2. **Method:** a dynamic-programming approach for determining the RNA secondary structure of **minimum free energy** using published stacking and destabilizing thermodynamic energies.
3. **Loop decomposition:** RNA structures are decomposed into distinct loop types — **hairpin loops**, **stacking regions** (stacked base pairs), **bulge and interior loops**, and **multibranched (junction) loops**. The total free energy is the sum of these loop contributions.
4. **Scale:** demonstrated by folding a 459-nucleotide immunoglobulin mRNA fragment and *E. coli* 16S rRNA fragments — i.e. the algorithm is polynomial, not exponential.

### Lorenz et al. (2011) — ViennaRNA Package 2.0 (algorithm class & complexity)

**URL:** https://doi.org/10.1186/1748-7188-6-26
**Retrieved via:** `WebFetch` of PubMed Central full-text page `https://pmc.ncbi.nlm.nih.gov/articles/PMC3319429/`
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference-implementation paper) / cites rank-1 primary

**Key Extracted Points:**

1. **Foundational algorithm:** the RNA folding problem is solved by dynamic programming, attributed to **Zuker & Stiegler (1981)**; ViennaRNA's MFE folding is "derived from the decomposition scheme as described by Zuker and Stiegler [1981]".
2. **Time complexity:** standard thermodynamic MFE folding runs in **O(n³)** time (McCaskill's partition-function algorithm achieves the same asymptotic class).
3. **Full citation:** Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology* 6:26.

### Ward, Datta, Wise, Mathews (2017) — multi-loop MFE recurrences (matrices & complexity)

**URL:** https://doi.org/10.1093/nar/gkx512
**Retrieved via:** `WebFetch` of PubMed Central full-text page `https://pmc.ncbi.nlm.nih.gov/articles/PMC5737859/`
**Accessed:** 2026-06-14
**Authority rank:** 1 (primary peer-reviewed paper)

**Key Extracted Points:**

1. **Matrices:** the algorithm computes, for all subsequences, **C(i,j)** = the minimum free energy of the substructure *closed* (enclosed) by base pair (i,j), plus multiloop fragment matrices (M / M1) and an exterior-loop matrix F.
2. **C(i,j) decomposition (verbatim sense):** "a structure enclosed in a base pair is either a hairpin loop, delimited by an interior loop, or branches in a multiloop" — i.e. C(i,j) = min(hairpin, interior/bulge over an inner pair, multiloop).
3. **Multiloop structure:** the multiloop fragment is composed of a part with one or more components (M) and a part with exactly one component (M1).
4. **Complexity:** the **standard (linear/affine multiloop) Turner model runs in O(n³) time and O(n²) space**; the logarithmic multiloop model raises this to O(n⁴)/O(n³). Seqeron implements the standard linear multiloop model a + b·unpaired + c·branches with the Turner 2004 DP values a = 9.3, b = 0, c = −0.9 (ViennaRNA 2.7.2 `params/default.c`: `ML_closing37 = 930`, `ML_BASE37 = 0`, `ML_intern37 = −90`; Mathews et al. 2004) → O(n³)/O(n²).

### NNDB — Turner 2004 Hairpin Loop worked Example 1 (full stem-loop ΔG)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-1.html
**Retrieved via:** `curl` of Wayback snapshot `http://web.archive.org/web/20240709061712/https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-1.html` (HTTP 200, 28452 bytes)
**Accessed:** 2026-06-14
**Authority rank:** 2 (official parameter specification with worked example)

**Key Extracted Points (verbatim term list):**

1. ΔG°37 = ΔG°37(Watson-Crick Pairs) + ΔG°37(terminal mismatch) + ΔG°37 Hairpin initiation(6).
2. Terms: (CG followed by AU) **−2.11** + (AU followed by CG) **−2.24** + (CG followed by AU) **−2.11** + **AU end penalty +0.45** + (terminal mismatch AU followed by AA) **−0.8** + Hairpin initiation(6) **+5.4**.
3. **Total ΔG°37 = −1.4 kcal/mol** (arithmetic sum = −1.41).
4. Verbatim note: "for unimolecular secondary structures, the helical intermolecular initiation does not appear" — confirms no helix-initiation constant is added for an intramolecular hairpin.
5. Helix component (3 stacks + one AU end) = −2.11 −2.24 −2.11 +0.45 = **−6.01**; loop component (init6 + terminal mismatch) = +5.4 −0.8 = **+4.6**.

### NNDB — Turner 2004 Hairpin Loop worked Example 2 (GG first mismatch)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-2.html
**Retrieved via:** `curl` of Wayback snapshot `http://web.archive.org/web/20240709061712/.../hairpin-example-2.html` (HTTP 200, 29119 bytes)
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. Closing pair A-U, 5-nt loop with first/last mismatch G/G. Terms: 3 stacks (−2.11, −2.24, −2.11) + AU end (+0.45) + terminal mismatch (AU followed by GG) −0.8 + **GG first mismatch −0.8** + initiation(5) +5.7.
2. **Total ΔG°37 = −1.9 kcal/mol** (arithmetic sum = −1.91); hairpin-loop component = +5.7 −0.8 −0.8 = **+4.1**.

### NNDB — Turner 2004 Watson-Crick Helix worked Example 2 (non-self-complementary)

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/wc-nsc-example.html
**Retrieved via:** `curl` of Wayback snapshot `http://web.archive.org/web/20240709061712/.../wc-nsc-example.html` (HTTP 200, 17317 bytes)
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. Helix duplex 5'-GCACG-3' / 3'-CGUGC-5'; stacks (GC followed by CG) **−3.42** + (CG followed by AU) **−2.11** + (AU followed by CG) **−2.24** + (CG followed by GC) **−2.36**.
2. **Pure nearest-neighbor stacking sum = −10.13 kcal/mol** (both helix ends are G-C, so no AU end penalty). Their tabulated −6.04 kcal/mol additionally includes the **intermolecular initiation +4.09** that applies only to bimolecular duplexes, NOT to an intramolecular stem.

### NNDB — Turner 2004 stacking, loop-initiation, and end-penalty parameter tables

**URL:** https://rna.urmc.rochester.edu/NNDB/turner04/wc-parameters.html , /loop.txt , /hairpin-mismatch-parameters.html , /tstack.txt , /triloop.txt , /tloop.txt , /hexaloop.txt
**Retrieved via:** `curl` of Wayback snapshots `http://web.archive.org/web/20240709061712/.../turner04/<file>` (HTTP 200) — these tables are the source for the constants already embedded in `RnaSecondaryStructure.cs` and independently verified for the sibling unit RNA-HAIRPIN-001 (see `docs/Evidence/RNA-HAIRPIN-001-Evidence.md`).
**Accessed:** 2026-06-14
**Authority rank:** 2

**Key Extracted Points:**

1. **Stacking ΔG°37 (kcal/mol):** AA/UU −0.93; AU/UA −1.10; UA/AU −1.33; CU/GA −2.08; CA/GU −2.11; GU/CA −2.24; GA/CU −2.35; CG/GC −2.36; GG/CC −3.26; GC/CG −3.42.
2. **Per AU end penalty:** +0.45 kcal/mol applied at each helix end terminating in A-U/U-A (and per GU end).
3. **Hairpin initiation (kcal/mol):** size 3→5.4, 4→5.6, 5→5.7, 6→5.4, 7→6.0, 8→5.5, 9→6.4; 10–30 from NNDB `loop.txt` (six entries 0.1 below the rounded formula, RNA-HAIRPIN-001); n>30 by ΔG(9)+1.75·R·T·ln(n/9).
4. **Primary parameter source:** Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH (2004) PNAS 101:7287–7292.

---

## Documented Corner Cases and Failure Modes

### From Zuker & Stiegler (1981) and NNDB Turner 2004

1. **Minimum hairpin loop = 3 nt:** a hairpin must enclose ≥3 unpaired nucleotides; pairs (i,j) with `j − i − 1 < 3` cannot close a hairpin (nearest-neighbor rules prohibit shorter loops). A sequence shorter than `minLoopSize + 2` cannot form any pair.
2. **No structure ⇒ ΔG = 0:** an unfoldable sequence (e.g. a homopolymer with no complementary bases) has MFE = 0 — the empty (open-chain) structure is always available and has ΔG° = 0 (the optimum is never positive).
3. **Empty / null input:** `null` or empty sequence has MFE = 0 (no structure).
4. **Intramolecular ⇒ no helix-initiation constant:** for a unimolecular fold the intermolecular helix-initiation term is **not** added (NNDB hairpin-example-1 note).

### From Ward et al. (2017)

5. **Multiloop model choice is correctness-affecting:** the affine (linear) multiloop model `a + b·unpaired + c·branches` gives O(n³); a logarithmic model changes both the optimum and the complexity. Seqeron uses the linear model with the sourced Turner 2004 DP values a = 9.3, b = 0 (per unpaired), c = −0.9 (per branch, closing pair included) — ViennaRNA 2.7.2 `default.c`; b = 0 is the published value, not a simplification.

---

## Test Datasets

### Dataset: NNDB Hairpin Example 1 — full MFE of a single hairpin

**Source:** NNDB Turner 2004 hairpin-example-1.html (Mathews et al. 2004); Zuker–Stiegler decomposition.

Constructed sequence whose unique optimal fold is exactly NNDB Example 1's stem-loop:
5'-`CACA` (stem) + `AAAAAA` (6-nt loop, first=last=A) + `UGUG` (stem)-3' = **`CACAAAAAAAUGUG`** (length 14).
Pairs C-G, A-U, C-G, A-U (outer→inner), closing pair A-U, loop A…A.

| Parameter | Value |
|-----------|-------|
| Sequence | `CACAAAAAAAUGUG` |
| Structure (dot-bracket) | `((((......))))` |
| Stacks (CG/AU, AU/CG, CG/AU) | −2.11, −2.24, −2.11 |
| AU end penalty | +0.45 |
| Terminal mismatch (AU·AA) | −0.8 |
| Hairpin initiation(6) | +5.4 |
| **MFE ΔG°37** | **−1.41 kcal/mol** (NNDB rounds to −1.4) |

### Dataset: NNDB Hairpin Example 2 — full MFE (GG first mismatch)

**Source:** NNDB Turner 2004 hairpin-example-2.html.

Constructed sequence: 5'-`CACA` + `GAAAG` (5-nt loop, first=last=G) + `UGUG`-3' = **`CACAGAAAGUGUG`** (length 13).

| Parameter | Value |
|-----------|-------|
| Sequence | `CACAGAAAGUGUG` |
| Stacks | −2.11, −2.24, −2.11 |
| AU end penalty | +0.45 |
| Terminal mismatch (AU·GG) + GG bonus | −0.8 + −0.8 |
| Hairpin initiation(5) | +5.7 |
| **MFE ΔG°37** | **−1.91 kcal/mol** (NNDB rounds to −1.9) |

### Dataset: Unfoldable / boundary sequences

**Source:** Zuker & Stiegler (1981) decomposition; nearest-neighbor minimum-loop rule.

| Input | Expected MFE | Why |
|-------|--------------|-----|
| `""` / `null` | 0 | no sequence |
| `AAAAAAAA` (homopolymer, no pairs) | 0 | no base pair possible; open chain ΔG = 0 |
| `GC` / any length < 5 (`minLoopSize+2`) | 0 | too short to enclose a 3-nt hairpin |

### Dataset: Algorithmic invariants (property-based)

**Source:** Zuker & Stiegler (1981) — the MFE optimum sums non-positive loop contributions over the optimal structure.

| Invariant | Statement |
|-----------|-----------|
| INV-01 | MFE ≤ 0 for every input (the empty structure with ΔG = 0 is always available; the optimum cannot be positive). |
| INV-02 | Adding nucleotides cannot raise the optimum: `MFE(s)` ≤ `MFE(prefix of s)` — extending the sequence only adds folding options (monotonic non-increase under suffix extension). |
| INV-03 | The MFE is deterministic (repeated evaluation gives the identical value and structure). The internal classic baseline uses a different per-pair model and is NOT expected to match (corrected 2026-11-25; the earlier "two engines agree" claim was false). |

---

## Assumptions

<!-- Behaviour not directly confirmed by an authoritative worked example. -->

1. ~~ASSUMPTION: multiloop per-unpaired cost c = 0, a = 9.25, helix term −0.63~~ — **withdrawn (2026-11-25)**: this described the pre-RNA-STRUCT-001 model. The DP now uses the sourced linear multiloop model **a = 9.3, b = 0 per unpaired nucleotide, c = −0.9 per branch** (closing pair included), ViennaRNA 2.7.2 `params/default.c` `ML_closing37 = 930`, `ML_BASE37 = 0`, `ML_intern37 = −90` (Turner 2004 DP set, Mathews et al. 2004; Ward et al. 2017 for the model family; also `rna_turner2004.par` `# ML_params` 0 0 930 3000 −90 −220 and `RNA.param(RNA.md())` MLclosing/MLbase/MLintern = 930/0/−90). The NNDB efn2 values 9.25 / 0.91·asymmetry / −0.63 + strain are not decomposable by the O(n³) DP and are used only by `CalculateMultibranchLoopEnergy(IReadOnlyList<int>, …)` (RNA-ENERGY-001 F15). Multiloop folds are cross-checked against `RNA.fold` (below).
2. **ASSUMPTION: rounding** — NNDB tabulates final ΔG°37 to one decimal; the implementation rounds to two decimals (`Math.Round(.., 2)`). Tests assert the exact two-decimal arithmetic sum of the cited per-term parameters with `.Within(1e-9)`, and additionally that `Math.Round(mfe, 1)` equals the one-decimal NNDB total, so the choice changes no source-defined parameter.

---

## Review 2026-11-25 (B12, Stage A + Stage B vs ViennaRNA 2.7.2)

**Oracle:** Python `RNA` 2.7.2, `fold_compound(seq, md).mfe()` / `eval_structure` with `md = RNA.md()` (`dangles = 2`, `noLP = 0`,
`min_loop_size` as tested); parameter files (a) default, (b) NNDB-precision (Evidence RNA-STRUCT-001 §12.2), loaded with
`RNA.params_load` and passed through an explicit `md` (gotcha in RNA-HAIRPIN-001 Evidence); sanity check GGGAAAACCC `(((....)))` = −2.02 (b).
For each sequence: (i) ours = RNA.fold MFE, (ii) identical structure, (iii) scalar = traceback energy = `CalculateFreeEnergy(our structure)`,
(iv) model match `RNA.eval_structure(our structure)` + NNDB-only terms (special GU closure, all-C, GGUC 3-stack, n = 1 bulge
degeneracy/special C, hairpin > 30 rounding) = ours, (v) optimality ours ≤ `CalculateFreeEnergy(RNA.fold structure)`.

| Set (random ACGU, `random.Random(seed)`) | file | (i) exact | (ii) same structure | (iii) | (iv) | (v) |
|---|---|---|---|---|---|---|
| 300 × 20–200 nt, no repeated neighbours, seed 20261126 | b | 300 | 299 (1 tie) | 300 | 300 | 300 |
| 200 × 20–200, no rep., seed 1 | b | 200 | 200 | 200 | 200 | 200 |
| 300 × 20–200, unrestricted, seed 20261128 | b | 202 (rest: NNDB-only terms) | 239 | 300 | 300 | 300 |
| 150, minLoopSize 4, seed 5 | b | 87 (NNDB-only terms) | 106 | 150 | 150 | 150 |
| 150, minLoopSize 6, no rep., seed 6 | b | 150 | 150 | 150 | 150 | 150 |
| 150, minLoopSize 0, no rep., seed 7 | b | 150 | 149 (tie) | 150 | 150 | 150 |
| 150 over ACGUN, seed 8 | b | 138 | 141 | 150 | 150 | 150 |
| 20 × 400–800 nt, no rep., seed 9 | b | 20 | 20 | 20 | 20 | 20 |
| 300 × 20–200, no rep., seed 20261127 | a | 14 (max/mean \|Δ\| 0.84/0.146: parameter rounding) | 253 | 300 | — | 300 |
| 300 × 20–200, unrestricted, seed 20261128 | a | 9 (max/mean 2.59/0.365) | 213 | 300 | — | 300 |

Contracts (all = `RNA.fold`, file b): empty → 0 (ViennaRNA rejects ""); lower case / T / N / IUPAC / `-` / `.` / X: same MFE and structure
(e.g. GGGNNNNCCC −1.42, GGG-AAA-CCC −1.32, GGGXAAAACCC −1.92); `md.min_loop_size` 0–2 behave as 3 in ViennaRNA too (hairpins < 3 are INF),
4/5/6/8 forbid shorter hairpins exactly as ours (GGGGAAACCCC: ml4 −4.22, ml5 −2.32); MAXLOOP 30: a 30-nt bulge forms (−12.07), 31-nt does not
(−8.46); lonely pairs allowed (61/300 MFE structures contain one, 299/300 identical to RNA.fold, `noLP = 0`); G·U helix ends via terminal
GU penalty + d2 mismatch (ACGUCGAUGCAGAUAU −0.96); ties: deterministic traceback, valid structure, eval = MFE (GUGCGCGAUAGCUGAGAUACUAU −1.52).

**Performance** (Release, this container, random ACGU): n = 200 / 500 / 1000 / 2000 → 0.33 / 0.82 / 3.6 / 17.7 s (Debug: 0.6 / 3.4 / 15.3 / 71 s);
ViennaRNA `RNA.fold` 0.03 / 0.26 / 1.07 / 4.6 s. Growth ×4.4–4.7 per doubling = the O(n²·MAXLOOP²) interior term with the O(n³) multiloop
split still sub-dominant; not pathological (≈ 3× RNAfold). Memory: two n×n double matrices (n = 1000: 16 MB).

**PredictStructure (Stage B finding).** It is not the MFE: it assembles `FindStemLoops` candidates (perfect stems ≥ minStemLength, loops
minLoopSize..maxLoopSize). Before: greedy by isolated stem-loop energy, accepting positive-energy stem-loops, and `MinimumFreeEnergy` = sum of
isolated energies (no exterior d2 terms). 304 random 20–200-nt sequences: `MinimumFreeEnergy` = eval of its own structure 1/304; eval > 0 in
51; ACAAAUACACAGGUUAAU `...(((......)))...` reported +2.91 (RNA.eval_structure +1.91; RNA.fold: open chain 0); GGGAAACCCAGGGAAACCC reported
−2.24 vs eval/RNA.fold −4.14. Fix: Zuker–Stiegler exterior-loop recursion W(j) over the candidates with in-context energies (under d2 the exterior
term depends on the sequence only, so disjoint stem-loops add) — exact optimum over every candidate combination (brute force on 60 sequences,
test R10); energy = `CalculateFreeEnergy` of the returned structure (304/304), never > 0, ≥ MFE. A heuristic stem-loop assembly is kept as a
separate contract (stem/loop bounds, `StemLoops` output) and documented as not the MFE fold; the MFE structure was not exposed via MCP — added
`rna_mfe_structure` (→ `PredictStructureMfe`).

---

## Recommendations for Test Coverage

1. **MUST Test:** `CalculateMinimumFreeEnergy("CACAAAAAAAUGUG")` = −1.41 (NNDB Example 1 full stem-loop) — Evidence: hairpin-example-1.html.
2. **MUST Test:** `CalculateMinimumFreeEnergy("CACAGAAAGUGUG")` = −1.91 (NNDB Example 2, GG first mismatch) — Evidence: hairpin-example-2.html.
3. **MUST Test:** unfoldable homopolymer `AAAAAAAA` → 0; empty/null → 0 — Evidence: Zuker–Stiegler decomposition (no pair possible ⇒ open chain ΔG = 0).
4. **MUST Test:** sequence shorter than `minLoopSize + 2` (e.g. `GCGC`) → 0 — Evidence: minimum-hairpin-loop rule (≥3 nt).
5. **MUST Test (property/invariant INV-01):** MFE ≤ 0 for arbitrary deterministic sequences — Evidence: Zuker–Stiegler optimum is over a set that includes the 0-energy open chain.
6. **MUST Test (property/invariant INV-02):** monotonic non-increase under suffix extension — Evidence: extending sequence only adds folding options.
7. **MUST Test:** `PredictStructure("CACAAAAAAAUGUG")` yields dot-bracket `((((......))))` with 4 base pairs C-G/A-U/C-G/A-U — Evidence: hairpin-example-1.html structure.
8. **SHOULD Test:** `PredictStructure("")` → empty result (no pairs, empty dot-bracket); homopolymer → all dots, 0 pairs — Rationale: documented empty/no-structure modes.
9. **COULD Test (INV-03 / performance baseline):** optimized DP equals classic baseline across lengths and record timing — Rationale: O(n³) algorithm requires a performance baseline (covered by `RnaSecondaryStructure_MFE_Benchmark`).

---

## References

1. Zuker M, Stiegler P (1981). Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information. *Nucleic Acids Research* 9(1):133–148. https://doi.org/10.1093/nar/9.1.133 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC326673/)
2. Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH (2004). Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure. *PNAS* 101:7287–7292. https://doi.org/10.1073/pnas.0401799101
3. Lorenz R, Bernhart SH, Höner zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology* 6:26. https://doi.org/10.1186/1748-7188-6-26 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3319429/)
4. Ward M, Datta A, Wise M, Mathews DH (2017). Advanced multi-loop algorithms for RNA secondary structure prediction reveal that the simplest model is best. *Nucleic Acids Research* 45(14):8541–8552. https://doi.org/10.1093/nar/gkx512 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC5737859/)
5. NNDB Turner 2004 Hairpin Loop Example 1. https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-1.html (accessed 2026-06-14 via Wayback snapshot 20240709061712)
6. NNDB Turner 2004 Hairpin Loop Example 2. https://rna.urmc.rochester.edu/NNDB/turner04/hairpin-example-2.html (accessed 2026-06-14 via Wayback snapshot 20240709061712)
7. NNDB Turner 2004 Watson-Crick Helix Example 2. https://rna.urmc.rochester.edu/NNDB/turner04/wc-nsc-example.html (accessed 2026-06-14 via Wayback snapshot 20240709061712)

---

## Change History

- **2026-06-14**: Initial documentation (RNA-MFE-001).
- **2026-11-25**: B12 review — multiloop assumption withdrawn (sourced 9.3/0/−0.9), INV-03 corrected, hairpin initiation 10–30 from loop.txt, ViennaRNA contract/edge-case cross-check, performance, PredictStructure exterior-DP fix, MCP `rna_mfe_structure`.
