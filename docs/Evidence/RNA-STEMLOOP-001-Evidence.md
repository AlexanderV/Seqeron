# RNA-STEMLOOP-001 Evidence Document

**Test Unit:** RNA-STEMLOOP-001  
**Title:** Stem-Loop Detection  
**Area:** RnaStructure  
**Created:** 2026-02-05  
**Status:** Complete

---

## 1. Overview

This document provides evidence for testing the stem-loop (hairpin) detection algorithms in RNA sequences. Stem-loops are fundamental secondary structure motifs that form when a single-stranded RNA folds back on itself via Watson-Crick and wobble base pairing.

---

## 2. Authoritative Sources

### 2.1 Primary References

| Source | Type | Key Information | URL/Citation |
|--------|------|-----------------|--------------|
| Wikipedia - Stem-loop | Encyclopedia | Formation, stability, structural contexts | https://en.wikipedia.org/wiki/Stem-loop |
| Wikipedia - Tetraloop | Encyclopedia | GNRA, UNCG, CUUG special loops | https://en.wikipedia.org/wiki/Tetraloop |
| Wikipedia - Pseudoknot | Encyclopedia | Crossing base pairs, detection difficulty | https://en.wikipedia.org/wiki/Pseudoknot |
| Woese et al. (1990) | Peer-reviewed | GNRA/UNCG tetraloops make up 70% in 16S rRNA | PNAS 87(21):8467-8471 |
| Heus & Pardi (1991) | Peer-reviewed | GNRA tetraloop stability mechanisms | Science 253(5016):191-194 |
| Antao et al. (1991) | Peer-reviewed | UUCG is most stable tetraloop | Nucleic Acids Res 19(21):5901-5905 |
| Svoboda & Cara (2006) | Review | Hairpin RNA biology and importance | Cell Mol Life Sci 63(7):901-908 |
| Rivas & Eddy (1999) | Peer-reviewed | Pseudoknot detection algorithm | J Mol Biol 285(5):2053-2068 |

### 2.2 Key Evidence from Sources

**From Wikipedia - Stem-loop:**
> "Optimal loop length tends to be about 4-8 bases long; loops that are fewer than three bases long are sterically impossible and thus do not form"

**From Wikipedia - Tetraloop:**
> "UUCG and GNRA tetraloops make up 70% of all tetraloops in 16S-rRNA"
> "The UUCG tetraloop is the most stable tetraloop"

**From Wikipedia - Pseudoknot:**
> "A pseudoknot occurs when base pairs cross each other, i.e., for pairs (i,j) and (k,l), we have i < k < j < l"

---

## 3. Algorithm Specification

### 3.1 Stem-Loop (Hairpin) Definition

A stem-loop consists of:
1. **Stem:** Double-stranded helical region formed by antiparallel Watson-Crick/wobble base pairing
2. **Loop:** Single-stranded region at the end of the stem (minimum 3 nucleotides due to steric constraints)

**Structure:**
```
      Loop (≥3 nt)
    5'-G-A-A-A-3'
       |     |
  5'-G-C     G-C-3'
     | |     | |
     G-C     G-C
     | |     | |
     C-G     C-G
     Stem    Stem
```

### 3.2 Detection Algorithm

**Input:** RNA sequence, minStemLength, minLoopSize, maxLoopSize, allowWobble

**Approach (one candidate per hairpin loop):**
1. Scan sequence for potential loop positions and sizes (minLoopSize..maxLoopSize, closed by a canonical pair)
2. Perfect stem (default): extend the helix outward while the next pair can form (WC or wobble)
3. Opt-in (`maxInternalLoopSize = M > 0`): the stem is the MFE unbranched chain (stacks, bulges, internal loops with n1 + n2 ≤ M) — Zuker–Stiegler V(i,j) without multiloops (§11)
4. Return stem-loops with ≥ minStemLength pairs; energy = isolated-motif evaluation (§11)

**Complexity:** perfect stems O(n · L · s) (s = motif span); unbranched mode O(n² · M²) time, O(n²) memory

### 3.3 Special Cases

**Tetraloops (4-nucleotide loops):**
- GNRA (N=any, R=purine): GAAA, GCAA, GGAA, GUAA
- UNCG: UACG, UCCG, UGCG, UUCG
- CUUG: CUUG, CCUG

**Energies (corrected 2026-11, B12):** the former "~3.0 kcal/mol stability bonus for GNRA/UNCG" had no source. Turner 2004
(NNDB special hairpin loops) gives *measured total loop free energies* for 18 tetraloops with their closing pair (e.g. CUUCGG +3.7,
CCUCGG +2.5), 2 triloops and 4 hexaloops; GNRA loops have no special entry and get only the generic GA first-mismatch bonus (−0.9).

### 3.4 Pseudoknot Detection

Pseudoknots occur when base pairs "cross" each other:
- For pairs (i,j) and (k,l): pseudoknot if i < k < j < l

**Note:** Standard secondary structure prediction cannot detect pseudoknots efficiently (NP-complete for general case).

---

## 4. Constraints and Invariants

### 4.1 Biological Constraints

| Constraint | Value | Source |
|------------|-------|--------|
| Minimum loop size | 3 nucleotides | Wikipedia (steric impossibility for <3) |
| Optimal loop size | 4-8 nucleotides | Wikipedia |
| Minimum stem length | 2-3 bp typically | Implementation choice |
| Valid base pairs | A-U, U-A, G-C, C-G, G-U, U-G | IUPAC, Wikipedia |

### 4.2 Algorithm Invariants

1. Within one stem-loop each base pairs at most once and pairs are properly nested (candidates of the enumeration may overlap each other — the method returns all of them)
2. Each base can participate in at most one base pair of a given candidate
3. Loop must be contiguous
4. Stem regions must be antiparallel (5'→3' pairs with 3'→5')

---

## 5. Test Methodology

### 5.1 Methods Under Test

| Method | Type | Priority |
|--------|------|----------|
| `FindStemLoops(sequence, minStem, minLoop, maxLoop, allowWobble)` | Canonical | Must |
| `FindHairpins(sequence, params)` | Variant | Should |
| `FindPseudoknots(sequence)` | Structural | Should |

### 5.2 Test Categories

**Must Tests:**
1. Basic hairpin detection (GGGAAAACCC)
2. Minimum loop size constraint (3 nt)
3. Wobble pair inclusion/exclusion
4. Empty/invalid input handling
5. Case insensitivity
6. Tetraloop detection
7. Dot-bracket notation generation

**Should Tests:**
1. Multiple stem-loops in sequence
2. Pseudoknot detection from base pairs
3. Energy calculation for stems
4. tRNA-like complex structures

---

## 6. Reference Test Data

### 6.1 Simple Hairpin

**Sequence:** `GGGAAAACCC`
**Expected:**
- Stem: 3 bp (GGG:CCC)
- Loop: 4 nt (AAAA)
- Structure: `(((....)))` (3 opening brackets, 4 dots for 4nt loop, 3 closing brackets)

**Source:** Basic Watson-Crick complementarity

### 6.2 GNRA Tetraloop

**Sequence:** `GGGCGAAAGCCC`
**Expected:**
- Stem: 4 bp (GGG-C paired with G-CCC; closing pair C-G)
- Loop: GAAA (4 nt) — canonical GNRA tetraloop (G=G, N=A, R=A, A=A)
- GA first mismatch bonus (−0.9 kcal/mol) applies via standard Turner 2004 model
- No NNDB special loop entry exists for GNRA with C-G closing; stability
  comes from GA mismatch bonus and structural tertiary interactions

**Source:** Wikipedia Tetraloop ("GNRA = G + any N + purine R + A"),
Heus & Pardi (1991), NNDB Turner 2004 hairpin-special-parameters.html
(confirms all special tetraloops are UNCG-family with C-G closing)

### 6.3 No Structure Possible

**Sequence:** `AAAAAAAAAAAAAAA`
**Expected:** No stem-loops (no complementary bases)

### 6.4 Minimum Loop Size

**Sequence with <3 loop:** Cannot form (steric constraint)
**Sequence:** `GCAUC` → Too short

### 6.5 Pseudoknot Example

**Base pairs:** (0,6) and (3,9)
- Positions: 0 < 3 < 6 < 9 → Crossing detected
- Result: Pseudoknot identified

**Source:** Wikipedia Pseudoknot

---

## 7. Edge Cases

| Case | Input | Expected | Source |
|------|-------|----------|--------|
| Empty string | `""` | Empty result | Implementation |
| Null | `null` | Empty result or exception | Implementation |
| Too short | `"GC"` | No stem-loops | Constraint |
| No complement | `"AAAA"` | Empty | Biology |
| Lowercase | `"gggaaaaccc"` | Same as uppercase | Convention |
| Minimum stem | 2 bp | Depends on minStemLength param | Param |
| Maximum loop | >30 nt | Higher energy penalty | Turner 2004 |
| All-C loop | `"CCCC"` | Energy penalty | Turner 2004 |
| Wobble only | Stem with only G-U | Valid if allowWobble=true | Biology |

---

## 8. Known Limitations

1. **No pseudoknot prediction from sequence** in this method (see `PredictStructurePseudoknot`, a separate unit)
2. **Energy model:** identical to ViennaRNA 2.7.2 `eval_structure` of the isolated motif plus sourced NNDB-only terms (§11); hairpin loops > 30 nt differ by ≤ 0.01 kcal/mol (ViennaRNA truncates the Jacobson–Stockmayer term to integer dcal/mol)
3. **Bulges / internal loops:** opt-in (`maxInternalLoopSize`), MFE unbranched stem; multiloops are out of scope (not a single stem-loop)
4. **Isolated-motif energy:** no exterior dangles (the motif has no neighbours)

---

## 9. Implementation Notes

The implementation uses a scanning approach:
1. Iterate over potential loop start positions
2. For each loop size in [minLoopSize, maxLoopSize]
3. Extend stem by checking base pair validity
4. Collect results meeting minStemLength requirement

**Key design decisions:**
- Allow wobble pairs by default (biologically valid)
- Default minimum loop size = 3 (steric constraint)
- Default minimum stem length = 3 (stability consideration)

---

## 10. Sign-off Checklist

- [x] Primary sources identified and cited
- [x] Algorithm specification documented
- [x] Test data with expected values
- [x] Edge cases enumerated
- [x] Invariants defined
- [x] Limitations acknowledged

---

## 11. Validation 2026-11 (B12 review, RNA-STEMLOOP-001)

### 11.1 Sources opened
- NNDB Turner 2004 hairpin loops (rna.urmc.rochester.edu/NNDB/turner04/hairpin.html): ΔG°hairpin(>3 nt) = initiation(n) + terminal
  mismatch + UU/GA first mismatch (not AG) + GG first mismatch + special GU closure (GU closing pair, not UG, preceded by two Gs) + all-C
  penalty; 3-nt loops: initiation(3) + C3 penalty only; "hairpin loops with fewer than 3 nucleotides are prohibited"; worked examples §5.3
  (−1.4, −1.9, −4.1 CcgagG, +2.1 all-C, −4.8 special GU) with ΔG = ΔG(helix incl. per-AU/GU-end 0.45) + ΔG(hairpin loop); special
  hairpin table (turner_2004_special_hairpin, hexaloop_dg.txt).
- ViennaRNA 2.7.2 sdist: `eval/eval_hairpin.c` (`vrna_E_hairpin` returns `Tetraloop_E` / `Hexaloop_E` / `Triloop_E` directly; size > 30:
  `hairpin[30] + (int)(lxc·ln(n/30))`), `params/default.c` (`mismatchH37`, `mismatchExt37`, `mismatchM37`, `dangle5/3_37`: every row/column 0 = 'N'
  equals the **maximum over A/C/G/U** — verified programmatically, 54/54 + 54/54 + 54/54 + 12/12 entries), `sequences/alphabet.c` (T ≡ U, other
  symbols → 0 = N).
- Zuker M, Stiegler P (1981) Nucleic Acids Res 9:133 — V(i,j) = min{hairpin, min interior/bulge/stack + V(p,q), multiloop}; the unbranched
  stem is V without the multiloop branch.

### 11.2 Oracle and method
ViennaRNA 2.7.2 Python `RNA`, `md.dangles = 2`, NNDB-precision parameter file (RNA-STRUCT-001 Evidence §12.2, file b). Every candidate of
`FindStemLoops` was re-evaluated as `RNA.fold_compound(seq[Start..End]).eval_structure(DotBracketNotation)` and NNDB-only terms (special GU
closure on ≥ 4-nt loops, all-C loops, GGUC/CUGG 3-stack) were added in Python. For the unbranched mode an independent Python DP computes, per
hairpin (p,q), min over unbranched chains of 0.45·[AU/GU outer pair] + Σ `fc.eval_int_loop` + `fc.eval_hp_loop` on sequences without identical
neighbours (no NNDB-only term can occur), and checks both the optimum and `eval_structure` of the reported structure. Scripts: scratchpad
`sl/cmp2.py`, `sl/cls.py`, `sl/ub.py`, `sl/fixed.py` (not in the repo).

### 11.3 Findings
- **SL-F1 A-U-closed special hexaloops double-penalised in `TotalFreeEnergy`.** `TotalFreeEnergy` = `CalculateStemEnergy` (NNDB helix, inner
  AU end +0.45) + measured hexaloop total (which already contains the closure; RNA-STRUCT-001 F5). `GGGACAGUGUUCCC` `((((......))))` −6.62 →
  **−7.07**; `CGACAGUACUCG` `(((......)))` −1.46 → **−1.91** (ViennaRNA eval −7.07 / −1.91 file b; −7.2 / −2.0 default file). Fix:
  `TotalFreeEnergy` = shared `EvaluateStructureEnergy` of the isolated motif (no second energy implementation).
- **SL-F2 Unknown base (N) as hairpin first mismatch scored 0.** `GGGAAANCCC` `(((....)))` −0.92 → **−2.02**; `GGGNAAANCCC` `(((.....)))`
  −0.82 → **−1.32**; `GGGCNNNNGCCC` `((((....))))` −4.34 → **−5.14**; `CalculateHairpinLoopEnergy("AAAN",'G','C')` 4.50 (= `eval_hp_loop`
  450). Fix: ViennaRNA N convention (max over A/C/G/U of mismatch + bonuses) in `CalculateHairpinLoopEnergy` (shared → also the MFE engines).
- **SL-F3 `CalculateHairpinLoopEnergy` ignored DNA T.** `("TTCG",'C','G')` 5.6 → **3.7** (special CUUCGG; `eval_hp_loop("CUUCGG")` 370).
- **SL-F4 `CalculateStemEnergy` (public, MCP `stem_energy`) read case-sensitive U-only keys and trusted `BasePair.Type`.** lower-case g·c,g·c
  0 → **−3.26**; A·T,A·T 0 → **−0.03** (AA/UU −0.93 + 2 × 0.45); a G·U pair labelled WatsonCrick lost its +0.45. Fix: bases normalised
  (upper case, T→U), AU/GU class from the bases (`Turner2004Parameters.IsTerminalPenaltyPair`).
- **SL-F5 Perfect-stem-only (simplification).** Stems with bulges/internal loops are sourced (Zuker–Stiegler V without multiloops) and
  implementable → added as opt-in `maxInternalLoopSize` (default 0 keeps the previous enumeration bit for bit). Example
  `GGGACAGGAAACUGCCC`: perfect `(((....)))` −1.09; `maxInternalLoopSize = 1` → `(((.(((....))))))` −7.23 = `RNA.fold` MFE (file b; −7.30 default).
- Doc claims without source removed/corrected: "~3.0 kcal/mol GNRA/UNCG bonus" (§3.3), "stem-loops must not overlap" (§4.2), "simplified energy
  model" (§8); TestSpec methods `FindHairpins` / `FindPseudoknots(sequence)` do not exist.

### 11.4 Numbers
| Check | Candidates | Mismatches before | after |
|---|---|---|---|
| perfect stems, random ACGU 15–60 nt, minStem 2, loops 3–10 (seeds 1, 7, 8, 9, 10) vs eval + NNDB-only | 9,147 + 31,829 + 31,864 + 62,901 + 61,521 | A-U-closed hexaloop candidates (SL-F1) | **0** |
| perfect stems, ACGUN, loops 3–12 (seed 12) | 15,745 | N-mismatch cases | **0** |
| perfect stems, `acgtNRY`, no wobble (seed 13) | 1,868 | — | **0** |
| perfect stems, minStem 1, loops 3–45, 15–90 nt (seed 11) | 144,437 | — | 9,615 differ, all loops > 30 nt, max \|Δ\| 0.010 (ViennaRNA integer truncation of lxc·ln(n/30); NNDB formula kept) |
| unbranched, M = 4 / 8 / 30 (200 / 200 / 100 seqs, noadj) vs Python DP over ViennaRNA loop functions | 13,836 / 13,587 / 7,732 hairpins | — | **0** (optimum and eval of reported structure) |
| `GCAUGCUACGUAGCAUCGAUCAGCUAGCUAGC`, M = 8 | 75 | — | 0 (5 rows locked in tests) |

