---
type: source
title: "Evidence: RNA-STEMLOOP-001 (stem-loop / hairpin enumeration)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-STEMLOOP-001-Evidence.md
sources:
  - docs/Evidence/RNA-STEMLOOP-001-Evidence.md
source_commit: 05292f4bc746f5b7f5f6637a2953864d096833cc
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-STEMLOOP-001

The validation-evidence artifact for test unit **RNA-STEMLOOP-001** — **Stem-Loop (Hairpin)
Detection** (area `RnaStructure`). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the synthesizing concept is
[[rna-stem-loop-enumeration]]. [[test-unit-registry]] tracks the unit.

This unit is a **sequence-scanning enumerator** of stem-loops — it *finds* every hairpin a sequence
can form under size constraints. That makes it distinct from its RNA-family neighbours: it is **not**
the hairpin *energy* calculator [[rna-hairpin-001-evidence|RNA-HAIRPIN-001]] (which scores one given
hairpin), **not** the single-optimal-fold [[rna-minimum-free-energy-folding|MFE folder]] (RNA-MFE-001,
DP over structure space), and **not** the miRNA-specific [[pre-mirna-hairpin-detection|pre-miRNA
precursor detector]]. It builds on the shared [[rna-base-pairing]] {A-U, G-C} + G-U wobble primitive
and emits [[rna-dot-bracket-notation|dot-bracket]] structures.

## What this file records

- **Methods under test:** `FindStemLoops(sequence, minStem, minLoop, maxLoop, allowWobble)` (canonical,
  Must), `FindHairpins(sequence, params)` (variant, Should), `FindPseudoknots(sequence)` (structural,
  Should).
- **Definition.** A stem-loop = a **stem** (antiparallel Watson-Crick/wobble double-stranded helix) +
  a single-stranded **loop** at its apex. Loop must be **≥ 3 nt** (loops < 3 are sterically
  impossible), optimal **4-8 nt** (Wikipedia). Valid pairs: A-U, U-A, G-C, C-G, G-U, U-G.
- **Enumeration algorithm.** Scan candidate loop positions; for each loop size in
  `[minLoopSize, maxLoopSize]`, extend the stem outward on both sides while base pairs remain valid
  (WC or, if `allowWobble`, wobble); collect stems meeting `minStemLength`. **Complexity O(n² · L)**
  (n = sequence length, L = max loop size).
- **Tetraloops (special 4-nt loops).** `GNRA` (N=any, R=purine — GAAA/GCAA/GGAA/GUAA), `UNCG`
  (UACG/UCCG/UGCG/UUCG), `CUUG`. **~70 %** of 16S-rRNA tetraloops are UUCG/GNRA (Woese 1990);
  **UUCG is the most stable** tetraloop (Antao 1991); a **~3.0 kcal/mol** stability bonus applies to
  GNRA/UNCG. GA first-mismatch bonus (−0.9 kcal/mol) via the Turner 2004 model.
- **Pseudoknot detection.** `FindPseudoknots` flags crossing base pairs — pairs (i,j),(k,l) cross when
  `i < k < j < l` — the same crossing primitive as [[rna-pseudoknot-detection]]. Prediction of
  pseudoknots from *sequence* is NP-complete for the general case and out of scope here.

### Authoritative sources

Wikipedia (Stem-loop, Tetraloop, Pseudoknot); Woese et al. 1990 (PNAS 87:8467); Heus & Pardi 1991
(Science 253:191); Antao et al. 1991 (NAR 19:5901, UUCG most stable); Svoboda & Cara 2006 (CMLS
63:901); Rivas & Eddy 1999 (JMB 285:2053, pseudoknot algorithm). Tetraloop energetics reconciled
against NNDB Turner 2004 hairpin-special-parameters.

## Datasets / oracles

| Sequence | Expected |
|----------|----------|
| `GGGAAAACCC` | 3-bp stem (GGG:CCC) + 4-nt loop (AAAA) → `(((....)))` |
| `GGGCGAAAGCCC` | 4-bp stem + GAAA GNRA tetraloop; GA first-mismatch bonus (−0.9); no NNDB special-loop entry for GNRA with C-G closing (stability from GA mismatch + tertiary interactions) |
| `AAAAAAAAAAAAAAA` | no stem-loops (no complementary bases) |
| `GCAUC` / `GC` | too short — cannot form |
| base pairs (0,6)+(3,9) | `0<3<6<9` → crossing → pseudoknot |

## Corner cases

- **Empty `""` / null** → empty result (or exception for null) — no throw contract for empty.
- **No complement (`AAAA`)** → empty. **Too short (`GC`)** → none.
- **Lowercase (`gggaaaaccc`)** → same as uppercase (case-insensitive).
- **Wobble-only stem (all G-U)** → valid iff `allowWobble = true`.
- **Loops > 30 nt / all-C loops** → higher Turner 2004 energy penalty.

## Invariants

1. Stem-loops must **not overlap** (simple structure). 2. Each base participates in **≤ 1** base pair.
3. Loop is **contiguous**. 4. Stem regions are **antiparallel** (5'→3' pairs with 3'→5').

## Known limitations (from the artifact)

1. **No pseudoknot prediction from sequence** — only detection from known base pairs.
2. **Simplified energy model** — less accurate than ViennaRNA.
3. **No internal loops / bulges** — hairpin loops only.
4. **Single structure** — does not enumerate suboptimal structures.

**Defaults:** wobble allowed, `minLoopSize = 3` (steric), `minStemLength = 3` (stability). No source
contradictions — Wikipedia, Woese 1990, Antao 1991, Heus & Pardi 1991, and NNDB Turner 2004 are
mutually consistent.
</content>
</invoke>
