---
type: concept
title: "Pre-miRNA precursor hairpin detection"
tags: [mirna, algorithm]
mcp_tools:
  - find_pre_mirna_hairpins
sources:
  - docs/Evidence/MIRNA-PRECURSOR-001-Evidence.md
  - docs/algorithms/MiRNA/Pre_miRNA_Detection.md
  - docs/Evidence/RNA-STEMLOOP-001-Evidence.md
source_commit: 05292f4bc746f5b7f5f6637a2953864d096833cc
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: mirna-precursor-001-evidence
      evidence: "Test Unit ID: MIRNA-PRECURSOR-001 ... Algorithm: Pre-miRNA Hairpin Detection ... Algorithm Group: MiRNA"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:rna-base-pairing
      source: mirna-precursor-001-evidence
      evidence: "uninterrupted pairing between mirrored positions from the two ends inward under the pairing set {A-U, U-A, G-C, C-G, G-U, U-G} ... allowing Watson-Crick and G:U wobble pairs"
      confidence: high
      status: current
---

# Pre-miRNA precursor hairpin detection

Finding **stem-loop precursor hairpins** in RNA — the ~60-120 nt pre-miRNA intermediates from which
Dicer excises a ~22 nt mature miRNA:miRNA* duplex (Bartel 2004; Ambros 2003). This is the **second
ingested unit of the MiRNA family** (test unit **MIRNA-PRECURSOR-001**, `MiRnaAnalyzer`); the record
is [[mirna-precursor-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern. It **builds on the shared RNA
base-pairing primitive** [[rna-base-pairing]] — a candidate stem is scored as consecutive
{A-U, G-C} + **G-U wobble** pairs between mirrored positions. This unit is the precursor-hairpin
sibling of the miRNA-target duplex; the **general** RNA secondary-structure stem-loop enumerator is
its neighbour [[rna-stem-loop-enumeration]] (RNA-STEMLOOP-001), not this page.

The unit has one **default heuristic** path plus three **opt-in** production paths — the default is
never changed by the opt-ins.

## 1. Default heuristic — `FindPreMiRnaHairpins`

Slides candidate windows (default length **55-120 nt**) across the T→U-normalised sequence and, for
each, counts **uninterrupted complementary base pairs from both ends inward**, breaking at the first
non-pairing mirrored position. A window `S = 5'stem + loop + 3'stem` is accepted only when:

- **stem length ≥ 18 bp** (Krol 2004; the effective stem is capped at `n/2 − 5`),
- **loop length = n − 2·stem, with 3 ≤ loop ≤ 25 nt** (Bartel 2004),
- the window is within the caller's length bounds.

It then extracts the **mature** strand from the 5' arm and the **star (miRNA\*)** strand from the
mirrored 3' arm (both length `min(matureLength, stemLength)`), builds a **dot-bracket** structure
(`(`×stem, `.`×loop, `)`×stem — always balanced, `|Structure| = |Sequence|`; the notation is parsed
and validated by [[rna-dot-bracket-notation]]), and computes a
**Turner 2004** nearest-neighbor hairpin `FreeEnergy` ([[rna-free-energy-turner-model|stacking +
loop-initiation + terminal-mismatch + a fixed 0.45 kcal/mol AU/GU-end penalty]]). Coordinates are
0-based inclusive.

**Documented limitation (accepted, not a bug):** the consecutive-pairing model is stricter than real
pre-miRNA structure — natural precursors carry internal mismatches, asymmetric bulges, and G:U
wobbles that offset the pairing frame. Real miRBase precursors are therefore **rejected** by the
default: **hsa-mir-21** yields only 16 consecutive end-pairs (< 18) and **hsa-let-7a-1** only 5
(tests M18/M19). Full folding (Zuker/Nussinov) is needed to catch them — supplied by the opt-in MFE
path below. Arm switching / opposite-arm mature products are also not modeled (ASM-03).

## 2. Opt-in MFE-structure assessment — `AssessHairpinByMfe` / `FindPreMiRnaHairpinsByMfe`

Replaces the consecutive-pairing scan with the **real minimum-free-energy structure** from the
validated Zuker–Stiegler folder (`RnaSecondaryStructure.CalculateMfeStructure`, RNA-STRUCT-001,
Turner 2004 NN model — reused, not re-derived). Hairpin features are read from the **actual MFE
dot-bracket**, so internal bulges/loops are tolerated. Acceptance:

- a **single dominant hairpin** — exactly one apical loop, nested stem, **no multibranch** (no `(`
  after a `)`) — the fold-back single-arm duplex of the annotation criteria (Ambros 2003; Meyers 2008);
- **stem base pairs ≥ 16** (Ambros 2003 uniform-system ≥16 complementary bases);
- **terminal loop ∈ [3, 25] nt** (Bartel 2004);
- **MFEI ≥ 0.85** (Zhang 2006), where `AMFE = 100·|ΔG°|/n` and `MFEI = AMFE / (G+C)%`.

`ΔG°` is `CalculateMinimumFreeEnergy` verbatim (negative); MFEI uses `|ΔG°|` so the published
"MFEI > 0.85" cutoff applies directly. Because genuine pre-miRNAs fold far more stably than random
sequence (Bonnet 2004), this path **detects the natural precursors the heuristic rejects**.

**Oracles (folded this session, exact expected values):**

| Candidate | Len | GC% | ΔG° | Stem bp | Loop | MFEI | Verdict |
|-----------|-----|-----|-----|---------|------|------|---------|
| `ValidHairpin57` | 57 | 43.86 | −48.48 | 27 | 3 | 1.9392 | ACCEPT |
| `hsa-mir-21` (MI0000077) | 72 | 48.61 | −35.13 | 32 | 3 | 1.0037 | **ACCEPT** (heuristic rejects) |
| `hsa-let-7a-1` (MI0000060) | 80 | 42.50 | −34.31 | 32 | 4 | 1.0091 | **ACCEPT** (heuristic rejects) |
| `5S-rRNA-like` (multibranch) | 120 | 64.17 | −47.04 | — | — | — | **REJECT** (not single hairpin) |
| `NoComplementarity` | 70 | 50.00 | 0.00 | 0 | — | 0 | REJECT |

The `5S-rRNA-like` case proves acceptance rests on **structure**, not merely a strong ΔG°: its
ΔG° (−47.04) beats `ValidHairpin57`'s yet it is rejected as multibranch.

## 3. Opt-in Drosha/Dicer cleavage-site prediction — `PredictDroshaDicerCleavage`

Predicts excision coordinates from the **published measuring ("ruler") rules only** — no trained
model:

- **Drosha:** cuts ~**11 bp** (≈ one helical turn) from the basal ssRNA–dsRNA junction (Han 2006);
  `DroshaCut5' = basalJunction + 11` = the 5' end of the 5p mature.
- **Dicer:** cuts ~**22 nt** from the Drosha-generated 5' end (5'-counting rule, Park 2011), fixing
  mature length at 22 nt.
- **2-nt 3' overhang:** each RNase III cut (Drosha, Dicer) leaves a 2-nt 3' overhang (Lee 2003).
- **CNNC motif:** a C-N-N-C 16-18 nt 3' of the Drosha cut (Auyeung 2013) sets an **optional**
  `HasCnncMotif` confidence flag — reported, not required.

Cross-checked against miRBase **hsa-miR-21-5p**: feeding a pri-miRNA whose 11-bp lower stem places
the +11 Drosha cut at the annotated 5p start reproduces `UAGCUUAUCAGACUGAUGUUGA` (22 nt) exactly.

## 4. Opt-in trained classifier — `ClassifyPreMiRna`

A **logistic-regression** natural-vs-background classifier over the published feature family
`[FreeEnergy, AMFE, MFEI, GC, PairedFraction]` (microPred/Xue 2005 lineage), standardised by
train-set mean/std. **Positives:** 13 public-domain miRBase precursors; **negatives:**
`DinucleotideShuffle` (Altschul–Erickson 1985 Eulerian-walk, dinucleotide-preserving) of the
positives — the Bonnet 2004 shuffled-background convention (4 shuffles/positive → 52 negatives).
Fixed-seed 70/30 split, batch gradient ascent (η=0.1, 20k epochs, L2 λ=1e-3); **held-out accuracy =
AUC = 1.0**. Bundled coefficients ship in `MiRnaAnalyzer.cs`. **No GPL miRDeep2 code or weights** are
used (its method was consulted only).

## Scope and limitations

A [[scientific-rigor|research-grade]] hairpin-detection reference. **Not implemented:** the
**read-stacking** (small-RNA-seq pileup) signal of miRDeep2 — it needs the caller's sequencing reads
and cannot be derived from sequence/structure alone (use miRDeep2 with your own reads). Pseudoknotted
precursors are out of scope (the RNA-STRUCT-001 energy-model floor). No source contradictions —
Bartel 2004/2009, Ambros 2003, Krol 2004, Bonnet 2004, Zhang 2006, Han 2006, Park 2011, and Auyeung
2013 are mutually consistent; the two recorded items are accepted **assumptions/deviations** (5'-arm
mature extraction; uninterrupted-stem strictness of the default path), both mitigated by the opt-in
MFE fold.
