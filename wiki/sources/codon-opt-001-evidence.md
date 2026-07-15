---
type: source
title: "Evidence: CODON-OPT-001 (Codon optimization — OptimizeSequence)"
tags: [validation, annotation]
doc_path: docs/Evidence/CODON-OPT-001-Evidence.md
sources:
  - docs/Evidence/CODON-OPT-001-Evidence.md
source_commit: 8d1b85e321fa52d6dea20205e2d0a4d2f28d1dbc
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: CODON-OPT-001

The validation-evidence artifact for test unit **CODON-OPT-001** (Sequence Optimization —
`OptimizeSequence`). One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm itself is summarized in
[[codon-optimization]], the rewriting operation of the codon-usage family whose measures are
[[relative-synonymous-codon-usage]] / [[codon-adaptation-index]] / [[effective-number-of-codons]].
See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** — Wikipedia "Codon usage bias" (optimization adjusts codons to host tRNA
  abundances; strategies = local mRNA folding, codon-pair bias, codon ramp, codon harmonization; rare
  codons deplete ribosomes; translation rate affects cotranslational folding); Wikipedia "Codon
  Adaptation Index" for the `w_i = f_i/max(f_j)` + `CAI = (∏ w_i)^{1/L}` scoring; **Sharp & Li (1987)**
  (PMID 3547335) for the CAI definition; **Plotkin & Kudla (2011)**, *Nat. Rev. Genet.* 12(1):32–42
  as the strategy review; **Mignon et al. (2018)**, *FEBS Lett.* 592(9):1554–1564 for the
  HarmonizeExpression / codon-harmonization strategy; Kazusa Codon Usage Database as the
  organism-table source.
- **Five strategies** — MaximizeCAI (most-frequent codon per aa), BalancedOptimization (CAI vs
  40–60% GC), HarmonizeExpression (match host distribution — Mignon 2018), AvoidRareCodons (replace
  only sub-threshold-frequency codons), MinimizeSecondary (avoid mRNA secondary structure).
- **Invariants** — protein preservation across all strategies (synonymous substitution only);
  CAI ∈ (0,1]; Met/AUG and Trp/UGG unchanged (unique codons); stop codons preserved.
- **Datasets** — organism preferred-codon tables with Kazusa species IDs and CDS counts: E. coli K12
  (316407, W3110, 4332 CDS; Leu CUG 0.50, Arg CGC/CGU 0.40/0.38 with AGA/AGG rare, Pro CCG 0.53);
  S. cerevisiae (4932, 14411 CDS; Leu UUA/UUG, Arg AGA 0.48, Pro CCA 0.42); H. sapiens (9606, 93487
  CDS; Leu CUG 0.40, Arg AGA/AGG 0.21, Pro CCC 0.32).
- **Edge cases** — empty → empty result; incomplete codons trimmed to `length % 3 == 0`; DNA (T)
  converted to RNA (U); lowercase handled case-insensitively; stop codons preserved not optimized.
- **Known failure modes** — invalid/unknown codon → `X` or error; non-RNA characters handled
  gracefully; a mid-sequence stop codon may terminate the protein prematurely.

## Implementation notes (from the artifact)

1. Implementation uses **RNA notation** (U not T), auto-converts T→U, and trims to complete codons.
2. **BalancedOptimization** targets 40–60% GC and **rebuilds the `Changes` list** after GC balancing
   so all modifications are reflected.
3. **Zero-frequency codons clamped to `1e-6`** in the CAI calculation (per Sharp & Li prescription) —
   the same guard documented on [[codon-adaptation-index]].
4. **MinimizeSecondary** delegates to BalancedOptimization for codon selection; a dedicated
   `ReduceSecondaryStructure` method handles secondary-structure reduction separately.

**Contradictions:** none — Wikipedia's strategy catalogue, Sharp & Li 1987 (CAI), Plotkin & Kudla
2011 (review) and Mignon 2018 (harmonization) agree; the five strategies each cite a distinct source
point. Several behaviours (strategies, invariants, edge cases) are recorded as "from theory" /
"from implementation" rather than a single external worked oracle — the correctness anchor is the
protein-preservation invariant plus the CAI formula, both source-backed.
