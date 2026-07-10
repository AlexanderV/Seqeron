---
type: source
title: "Evidence: PHYLO-BOOT-001 (Phylogenetic Bootstrap Analysis — Felsenstein's bootstrap proportions)"
tags: [validation, phylogenetics]
doc_path: docs/Evidence/PHYLO-BOOT-001-Evidence.md
sources:
  - docs/Evidence/PHYLO-BOOT-001-Evidence.md
source_commit: fd10c2dd29b1edd788a25d3aed9b67cdab80685d
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PHYLO-BOOT-001

The validation-evidence artifact for test unit **PHYLO-BOOT-001** — **Phylogenetic Bootstrap Analysis
(Felsenstein's bootstrap proportions, FBP)**: resample the **alignment columns (sites) with
replacement** into same-length pseudo-alignments, rebuild a tree from each, and score every internal
clade of the reference tree by the **fraction of replicate trees containing that clade** (its bootstrap
support). This is the **first phylogenetics-family (`PHYLO-*`) Evidence file** and one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the procedure,
clade-matching rule, worked oracles, and the two units/representation assumptions are summarized in the
dedicated anchor concept [[phylogenetic-bootstrap-support]]. See [[test-unit-registry]] for how units
are tracked.

## What this file records

- **Online sources:**
  - **Felsenstein, J. (1985)** *"Confidence Limits on Phylogenies: An Approach Using the Bootstrap"*
    (Evolution 39(4):783–791, authority 1 — origin of the method; abstract fetched via the OSTI
    bibliographic record after the Wiley DOI returned HTTP 403) — bootstrap "resampl[es] points from
    one's own data, with replacement, to create bootstrap samples of the same size"; "keep all of the
    original species while sampling characters with replacement"; *P* for a group = fraction of
    bootstrap samples in which that group appears; "if a group shows up 95% of the time or more, the
    evidence for it is taken to be statistically significant" (descriptive threshold, not a
    computation parameter).
  - **Lemoine et al. (2018)** *"Renewing Felsenstein's Phylogenetic Bootstrap in the Era of Big Data"*
    (Nature 556:452–456, PMC6030568, authority 1) — "resample, with replacement, the sites of the
    alignment to obtain pseudo-alignments of the same length"; support of every branch = "proportion
    of pseudo-trees containing that branch" (binary per-replicate); "any branch of an X-tree defines a
    bipartition of X"; the values are "Felsenstein's bootstrap proportions (FBPs)", practitioners use
    thresholds such as 70%.
  - **Biopython `Bio.Phylo.Consensus`** (`bootstrap`, `bootstrap_trees`, `get_support`; master branch,
    authority 3) — column resampling `for j in range(length): col = random.randint(0, length - 1)`
    (resampled columns = original length); support counting `c.confidence = (t + 1) * 100.0 / size`;
    clades compared **by terminal (leaf) name set**, not branch length or internal labels; only
    non-terminal clades (`find_clades(terminal=False)`) are scored.
- **Corner cases / failure modes:** reference tree provides the entity set (clades only in replicates
  are not reported); binary per-replicate scoring (exact match or nothing, no partial credit);
  same-length resampling required; clade equality by leaf-name set (order/branch-length irrelevant);
  terminal single-leaf clades excluded.
- **Datasets (documented oracles):**
  - *Two-group deterministic alignment* — taxa A,B,C,D with `A=B=AAAAAAAAAA`, `C=D=GGGGGGGGGG`
    (length 10), UPGMA + JukesCantor, 100 replicates, seed 42. Every column is all-`A` or all-`G`, so
    any resample leaves distances unchanged (`d(A,B)=d(C,D)=0`, cross-pairs JC-saturated `+∞`) → every
    replicate reproduces the same topology → `support({A,B}) = support({C,D}) = 100/100 = 1.0` for all
    seeds.
  - *All-identical alignment* — taxa A,B,C each `ACGTACGT` → identical zero-distance matrix every
    replicate → every reported clade has support 1.0.
- **Coverage recommendations:** MUST-test two-group → support 1.0 for {A,B}/{C,D}; all support values
  ∈ [0,1]; determinism (same seed → identical output); returned keys = exactly the non-trivial clades
  of the reference tree; support = `k/replicates` (denominator = replicate count). SHOULD-test
  all-identical → all clades 1.0; input validation (null / <2 sequences / `replicates < 1` throw).
  COULD-test resampled alignment length equals the original length.

## Deviations and assumptions

No deviations from the literature; **two documented assumptions**, both source-consistent. (1)
**Rooted-clade scoring rather than unrooted bipartitions** — Felsenstein/Lemoine describe bipartitions
of unrooted trees while this unit scores rooted clades (subtree leaf-sets) of the UPGMA/NJ tree,
matching Biopython's `get_support` (compares clades by terminal set); for rooted ultrametric UPGMA
trees this is the conventional, consistent representation. (2) **Support reported as a proportion in
[0,1] rather than a percentage** — the literature/Biopython express *P* as a percent (×100), this unit
returns raw `count/replicates ∈ [0,1]` (multiply by 100 to recover the published percentage); a
units/labeling choice that changes neither which clades are reported nor their ranking. No source
contradictions.
