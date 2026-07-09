---
type: source
title: "Evidence: ONCO-PHYLO-001 (tumor phylogeny reconstruction — clonal tree from CCF clusters)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-PHYLO-001-Evidence.md
sources:
  - docs/Evidence/ONCO-PHYLO-001-Evidence.md
source_commit: ea992b89032ebee5bf103593a140cf59a8d032d8
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-PHYLO-001

The validation-evidence artifact for test unit **ONCO-PHYLO-001** — **Tumor Phylogeny Reconstruction:
a clonal-evolution tree from CCF clusters governed by the sum rule + lineage-precedence rule**. The
**twenty-sixth ingested unit of the Oncology family** and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is
synthesized in [[tumor-phylogeny-clonal-tree-reconstruction]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (mutually consistent — LICHeE VAF form and PICTograph CCF form state the same two
  constraints, no contradictions):**
  - **Popic et al. (2015) LICHeE** *Genome Biology* 16:91 (DOI 10.1186/s13059-015-0647-8; arXiv
    preprint 1412.8574), rank 1 — the **perfect-phylogeny ordering constraints** (verbatim): (1) a
    mutation in a sample set cannot be a successor of one present in a smaller subset; (2) a mutation
    cannot have a VAF higher than its predecessor (except via CNVs); (3) the summed VAFs of disjoint
    subclones cannot exceed a common predecessor's VAF. The **ancestor ≥ descendant edge rule** (Eq. 2):
    edge `(u→v)` added only if for all samples i `u.VAFᵢ ≥ v.VAFᵢ − ϵ_uv` and `u.VAFᵢ = 0 ⇒ v.VAFᵢ = 0`.
    The **sum rule** (Eq. 5): a valid lineage tree is a spanning tree where `∀u ∀i :
    Σ_{v:(u→v)} v.VAFᵢ ≤ u.VAFᵢ + ϵ` — children's VAF centroids may not exceed the parent's.
    **Inequality not equality** because not all true branches need be observed.
  - **Zheng et al. (2022) PICTograph** *Bioinformatics* 38(15):3677 (DOI 10.1093/bioinformatics/btac367,
    PMC9344857), rank 1 — the **Sum Condition** in CCF form (verbatim): "the CCF of an ancestral clone
    must be greater than or equal to the sum of CCFs of its descendants" (trees exceeding by > ε₂ = 0.2
    excluded); **Lineage Precedence** (verbatim): "the CCF of any mutation cannot exceed the CCF of its
    ancestor" (descendant present in a subset of the ancestor's samples, CCF at most ε₁ = 0.1 greater).
    The Sum Condition **generalizes the pigeonhole principle** per node.

- **Documented corner cases / failure modes:** **insufficient parent budget** — three sibling groups
  whose CCFs sum > parent CCF cannot all descend from the parent; conflicting branches must be
  removed/re-placed; **under-determined private mutations** — single-sample groups are
  under-constrained (multiple valid ancestors), needing a deterministic tie-break; **noise margin** —
  equality relaxed to `≤ u + ϵ` to tolerate CCF measurement noise; **presence-pattern constraint** — a
  descendant cluster is present only in (a subset of) the ancestor's samples.

- **Datasets (deterministic, hand-derived from Eq. 2 + Eq. 5):**
  - **Linear chain (single sample):** Normal 1.0 / A 1.0 / B 0.6 / C 0.3 → edges Normal→A→B→C, Trunk
    {A}, Branches {B, C}.
  - **Branching (two samples):** Normal (1,1) / A (1,1) / B (0.6, 0.0) / C (0.0, 0.7) → Normal→A, A→B,
    A→C; sum rule under A per sample s1 = 0.6 ≤ 1.0, s2 = 0.7 ≤ 1.0; B, C mutually non-ancestral
    (each absent in the other's sample).
  - **Sum-rule forces a chain:** Normal 1.0 / A 1.0 / B 0.6 / C 0.6 → B and C cannot both be children
    of A (0.6 + 0.6 = 1.2 > 1.0); C nests under B → chain Normal→A→B→C.

- **Coverage recommendations (10 items):** MUST — linear chain from descending single-sample CCFs;
  branching tree from two private single-sample clusters; sum rule rejects two equal-CCF siblings and
  forces a chain; trunk = CCF-≈1 root-path clusters / branch = the rest; **invariant** ancestor CCF ≥
  descendant CCF on every edge (Eq. 2); **invariant** per-node sum rule holds (Eq. 5). SHOULD — empty
  input → root-only tree; single cluster → child of root, is the trunk; null / CCF ∉ [0,1] / NaN /
  inconsistent sample count → exceptions. COULD — tolerance ε > 0 admits a near-violating edge that
  ε = 0 rejects.

## Deviations and assumptions

- **ASSUMPTION — deterministic tie-break for under-constrained placement.** Popic et al. note
  private / under-constrained clusters admit multiple valid ancestors. This implementation attaches
  each cluster to its **deepest valid ancestor** (smallest-total-CCF candidate parent whose per-sample
  sum-rule budget still admits the child), ties broken by ascending cluster id — the most-recent common
  ancestor consistent with all cited constraints. Changes only **which single valid tree is returned**,
  not the valid-tree set.
- **ASSUMPTION — noise margin ε = 0.** The sources relax the inequalities by a configurable ε (LICHeE
  ϵ; PICTograph ε₁ = 0.1, ε₂ = 0.2). Because this unit consumes already-clustered CCF point estimates
  (the noise model lives in ONCO-CCF-001), the default uses ε = 0 (strict inequalities), exposed as an
  optional tolerance. Setting ε > 0 only widens admissibility — never flips a strictly-satisfied
  relationship.

No source contradictions — LICHeE (VAF form) and PICTograph (CCF form) state the identical two
constraints, and both assumptions are source-consistent rather than source-contradicting.
