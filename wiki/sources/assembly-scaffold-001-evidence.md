---
type: source
title: "Evidence: ASSEMBLY-SCAFFOLD-001 (Scaffolding)"
tags: [validation, assembly]
doc_path: docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md
sources:
  - docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md
source_commit: 3fda6bc37e2510c7c559dcc3b40efb62ecf722bd
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ASSEMBLY-SCAFFOLD-001

The validation-evidence artifact for test unit **ASSEMBLY-SCAFFOLD-001** — scaffolding: joining
ordered contigs into a scaffold in which adjacent contigs are separated by runs of the character
`N` whose length equals the estimated inter-contig distance. One instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the construction rule,
the non-positive-gap → 100-N default, oracles and the scoping assumption are summarized in
[[scaffolding]], the anchor for the assembly SCAFFOLD family. See [[test-unit-registry]] for how
units are tracked.

## What this file records

- **Online sources** (all accessed 2026-06-13 via WebFetch):
  - **Jackman et al., ABySS 2.0 — "Scaffolding a genome sequence assembly using ABySS"** (rank 1;
    *Genome Research* 2017, 27:768–777) — the canonical scaffold-construction rule ("sequences of
    the vertices in a path are concatenated, interspersed with gaps represented by a run of the
    character N, whose length corresponds to the estimate of the distance"), gap length = distance
    estimate, ML distance estimator is upstream (this unit consumes a supplied estimate), and
    **negative estimate = overlap → contigs merged if the overlap is found**.
  - **NCBI AGP Specification v2.1** (rank 2; official NCBI/INSDC file format) — gap component type
    `N` (specified size) vs `U` (unknown size); "Gap lengths must be positive. Negative gaps and
    gap lines with zero length are not valid."; **negative/unknown gaps → `U` with size 100**, the
    GenBank/EMBL/DDBJ standard unknown-gap length.
  - **Sahlin et al. (2012) — "Improved gap size estimation for scaffolding algorithms"** (rank 1;
    *Bioinformatics* 28(17):2215–2222) — gap = the unknown distance `d` between two contigs; the
    **negative-gap case frequently occurs** because a de Bruijn assembler splits contigs leaving a
    one-k-mer overlap, confirming non-positive gaps are a real expected input class.
  - **Pop, Kosack & Salzberg (2004) — "Hierarchical Scaffolding With Bambus"** (via Wikipedia
    primary citation, rank 4 → primary *Genome Research* 14(1):149–159) — scaffolding links a
    non-contiguous series of sequences separated by gaps of known length; Bambus greedy scaffolder
    "joins together contigs with the most links first" (here links supplied pre-ordered).
- **Datasets (published oracles):**
  - contigs `["ACGT","TTGG","CCAA"]`, links `[(0,1,3),(1,2,2)]`, default gap `N` →
    `ACGTNNNTTGGNNCCAA` (length 17), **1** scaffold — ABySS construction rule.
  - contigs `["AAAA","TTTT"]`, links `[(0,1,-5)]` → `AAAA`+(`N`×100)+`TTTT` (length 108) — the
    AGP unknown-size default for a non-positive estimate.
- **Corner cases / failure modes** — zero-length gap is not a valid AGP gap line (→ unknown-gap
  default 100, never "0 N"); a negative gap is not a valid positive run (→ 100 N unless overlap is
  resolved); a negative distance estimate indicates the contigs should overlap (ABySS merges them
  if the overlap is found, otherwise no positive `N`-run represents it).
- **Recommended coverage** — MUST: positive gaps concatenate to one scaffold with exact `N`-runs
  (`ACGTNNNTTGGNNCCAA`); a gap of size `g` emits exactly `g` fill chars and scaffold length =
  Σ|contig| + Σgap; a non-positive (zero/negative) estimate emits exactly 100 fill chars; a custom
  gap character is used verbatim instead of `N`; each contig appears in at most one scaffold (a link
  to an already-placed contig is skipped); an unlinked contig becomes its own single-contig
  scaffold. SHOULD: out-of-range / self link indices are ignored; null `contigs` / null `links` throw
  `ArgumentNullException` (mirrors `MergeContigs`). COULD: empty contig list → empty result.

## Assumptions (from the artifact)

One assumption record, a **scoping** decision rather than an invented value: **the unresolved-overlap
placeholder uses the AGP unknown-gap length (100)**. ABySS merges contigs when a negative estimate's
overlap is actually *found*; this unit performs no overlap resolution, so a non-positive estimate is
emitted as a gap of unknown size. The chosen length (100 `N`) is the GenBank/EMBL/DDBJ standard for
unknown-size gaps per NCBI AGP v2.1 — so the **constant is source-backed**; only the decision to fall
back to it (rather than resolve the overlap) is the assumption. No numeric value is invented.

No contradictions among the sources — ABySS 2.0, the NCBI AGP spec, Sahlin et al. and Bambus give the
same "ordered contigs + sized `N`-gaps" scaffold model. The AGP 100-N unknown-size default and the
ABySS negative-gap = overlap rule are complementary, covering the same non-positive-gap case from the
file-format side and the assembler side respectively.
