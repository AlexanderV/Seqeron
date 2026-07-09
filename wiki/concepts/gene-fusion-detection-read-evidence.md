---
type: concept
title: "Gene-fusion detection from split + spanning reads (with in-frame codon-phase check)"
tags: [oncology, structural-variant, algorithm]
sources:
  - docs/Evidence/ONCO-FUSION-001-Evidence.md
source_commit: ea13dcc183c950560fe910068244e507f45a455f
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-fusion-001-evidence
      evidence: "Test Unit ID: ONCO-FUSION-001, Algorithm: Fusion Gene Detection (candidate fusion calling from breakpoint-supporting reads)"
      confidence: high
      status: current
---

# Gene-fusion detection from split + spanning reads

The **fourteenth ingested Oncology unit** (ONCO-FUSION-001) and the wiki's **first
gene-fusion / read-evidence structural-rearrangement method**. It decides whether a candidate
gene fusion (e.g. **BCR–ABL1**, **EML4–ALK**, **TMPRSS2–ERG**) is **detected** from its
supporting-read counts, and whether the fused product is **in-frame**. Validated under test
unit **ONCO-FUSION-001** ([[onco-fusion-001-evidence]]); [[test-unit-registry]] tracks the
unit and [[algorithm-validation-evidence]] describes the artifact pattern. The decision rule is
**exact** with respect to the cited STAR-Fusion / Arriba thresholds.

The evidence model is the standard **STAR-Fusion** (Haas 2017/2019) and **Arriba** (Uhrig 2021)
split-read + discordant-pair + minimum-support paradigm — see [[onco-fusion-001-evidence]] for
the source-by-source trace.

## Two evidence classes and total support

A fusion is supported by two classes of read evidence bracketing the breakpoint:

- **Junction / split reads** — a single read whose two segments align noncontiguously across the
  breakpoint. Arriba splits these into **`split_reads1` / `split_reads2`** by **anchor** = the
  gene the read's **longer** aligned segment maps to.
- **Discordant / spanning mates** (a.k.a. bridge reads) — paired-end reads whose mates align in a
  nonlinear way, bracketing the breakpoint. Arriba's **`discordant_mates`**.

```
total support = split_reads1 + split_reads2 + discordant_mates
junction reads = split_reads1 + split_reads2
```

## Detection rule (STAR-Fusion defaults)

A candidate is **DETECTED** iff it passes both the support thresholds and the distinct-gene
invariant:

- **With junction reads present** (`junction ≥ MIN_JUNCTION_READS = 1`): require
  **`total ≥ MIN_SUM_FRAGS = 2`**. So junction = 1, total = 1 → **REJECTED** (fails min-sum);
  junction = 1, total = 2 → **DETECTED**.
- **With zero junction reads**: the min-sum rule does not apply; require
  **`discordant ≥ MIN_SPANNING_FRAGS_ONLY = 5`**. So 4 discordant with no split reads →
  **REJECTED**; 5 → **DETECTED**.
- **Distinct-gene invariant**: `gene5p ≠ gene3p` — a "fusion" of a gene with itself is not a gene
  fusion → **REJECTED**.

Candidates are reported **ordered by descending total support** (most-supported fusions first),
reflecting STAR-Fusion's scoring by "the abundance of fusion-supporting reads". (`MIN_FFPM = 0.1`,
the abundance-normalized fragments-per-million filter, is recorded but not part of the count-based
rule. The thresholds are the CLI-exposed defaults and can be overridden.)

## In-frame check (exon-phase / modulo-3)

Separately from detection, the fused transcript is **in-frame** iff the coding bases the 5'
partner contributes before the breakpoint keep the 3' partner's codons in phase:

```
in-frame  ⇔  (5' coding bases before breakpoint − 3' partner coding-start phase) mod 3 == 0
```

An in-frame fusion reads through into a protein made of parts of both partners; **out-of-frame**
shifts the downstream reading frame and is unlikely to yield a functional protein. This uses
codon **phase** only — it does **not** scan the spliced transcript for a premature stop codon
(Arriba's `stop-codon` reading_frame value), which needs full transcript reconstruction — that is
the separate protein-consequence unit
[[fusion-breakpoint-frame-and-protein-prediction]] (ONCO-FUSION-003), which elaborates this binary
in-frame check into a four-state `BreakpointFrameStatus` plus `PredictFusionProtein`. Naming a detected fusion — the canonical HGNC `5′::3′` designation and
directional known-fusion matching — is the separate annotation unit
[[gene-fusion-nomenclature-known-fusion-lookup]] (ONCO-FUSION-002).

## Invariants and edge cases

- **Read-through false positives.** Adjacent same-strand neighbouring genes produce common
  read-through chimeras; guarded against by the support thresholds and the distinct-gene rule.
- **Stop codon before junction.** Even when numerically in-frame, an upstream stop means the 3'
  partner is not translated — out of scope here (transcript reconstruction), handled by
  [[fusion-breakpoint-frame-and-protein-prediction]] (ONCO-FUSION-003).
- **Input validation.** Null input → `ArgumentNullException`; negative counts → `ArgumentException`.

Worked oracles (from [[onco-fusion-001-evidence]], `split1,split2,discordant`): **EML4-ALK**
(3,2,4) → DETECTED; **TMPRSS2-ERG** (1,0,1) → DETECTED; **CD74-ROS1** (0,0,5) → DETECTED;
**NCOA4-RET** (0,0,4) → REJECTED (span < 5); **KIF5B-RET** (1,0,0) → REJECTED (sum < 2);
**ALK-ALK** (5,5,5) → REJECTED (same gene). Frame: 300/phase-0 → in-frame, 301/phase-0 →
out-of-frame, 301/phase-1 → in-frame.

## Scope and limitations

A [[research-grade-limitations|research-grade]] method, **not for clinical use**. It operates on
**already-grouped candidate breakpoints** with per-class supporting-read counts (the Arriba output
schema), **not** raw BAM records — extracting chimeric reads from alignments is a separate
`FindChimericReads` step out of this unit's canonical-threshold scope. It is the **read-evidence
fusion / rearrangement** Oncology unit, orthogonal to the copy-number / clonal-structure ONCO units
(e.g. [[copy-number-alteration-classification]], [[cancer-cell-fraction-clonal-clustering]]) and to
the clinical-interpretation ONCO units ([[clinical-actionability-oncokb-levels]],
[[cancer-variant-tier-classification-amp-asco-cap]]) that would consume a called fusion. Distinct
from the signed-permutation, gene-order [[genome-rearrangement-breakpoint-distance]] (which measures
rearrangement *distance* between whole gene orders, not read-evidence fusion calling). No source
contradictions.
