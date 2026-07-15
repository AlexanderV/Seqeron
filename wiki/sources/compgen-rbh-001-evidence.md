---
type: source
title: "Evidence: COMPGEN-RBH-001 (Reciprocal Best Hits ortholog identification)"
tags: [validation, comparative-genomics]
doc_path: docs/Evidence/COMPGEN-RBH-001-Evidence.md
sources:
  - docs/Evidence/COMPGEN-RBH-001-Evidence.md
source_commit: fbf303b3f1b14f641169d98548fbaf753faed752
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: COMPGEN-RBH-001

The validation-evidence artifact for test unit **COMPGEN-RBH-001** — **ortholog identification by
Reciprocal Best Hits (RBH / bidirectional best hits)**. This is a **Comparative-genomics** family
Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm, its rules, invariants,
worked oracles, and corner cases are summarized in the dedicated concept
[[ortholog-detection-reciprocal-best-hits]], which was **deliberately scoped as the shared RBH
anchor this unit reuses** rather than re-deriving. COMPGEN-RBH-001 is the **RBH-only** slice: it
covers the between-genome ortholog rule but not the within-genome in-paralog rule that the broader
sibling artifact [[compgen-ortho-001-evidence|COMPGEN-ORTHO-001]] adds. RBH is also the shared-gene
sub-unit behind the [[genome-comparison-core-dispensable]] pipeline (whose conserved/core set *is*
these RBH pairs). Its sibling COMPGEN units are [[average-nucleotide-identity]],
[[synteny-and-rearrangement-detection]], [[conserved-gene-clusters-common-intervals]], and
[[dot-plot-word-match]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (both authority rank 1):**
  - **Moreno-Hagelsieb & Latimer 2008** (*Bioinformatics* 24:319–324, DOI 10.1093/bioinformatics/btm585)
    — the **operational RBH definition**, quoted verbatim from the fetched OUP article: two genes in
    two genomes are orthologs iff their products *"find each other as the best hit in the opposite
    genome."* Best hit = BLASTP hits *"sorted … from highest to lowest bit-score, then, if the
    bit-scores were identical, from smallest to highest E-values"* ⇒ **maximum-similarity candidate
    with a deterministic tie-break**. Gates: **≥ 50 % coverage** of *"any of the protein sequences in
    the alignments"* and a **maximum E-value of 1×10⁻⁶** (a significance gate → minimum-similarity
    gate in alignment-free space).
  - **Tatusov, Koonin & Lipman 1997** (*Science* 278:631–637, DOI 10.1126/science.278.5338.631;
    method taken from the NCBI Handbook COG page NBK21090 after the *Science* PDF returned HTTP 404)
    — COGs built from **genome-specific best hits (BeTs)**: proteins *"more similar to each other than
    they are to any other proteins from the same genomes"* form an orthologous set; COGs are
    *"triangles of mutually consistent … BeTs"*, each COG *"a group of three or more proteins …
    inferred to be orthologs."* The pairwise special case of a mutually consistent BeT is exactly the
    reciprocal best hit.

- **Contradiction resolved in-file (data-quality note, not a source conflict):** a search-engine
  summary claiming a **60 %** coverage threshold was **contradicted by the fetched article body
  (50 %) and rejected**. The wiki records the source-backed value of ≥ 50 %.

- **Documented corner cases / failure modes:** tie in best hit → deterministic tie-break (smaller
  E-value) so the best hit is unique and matching is order-independent; coverage filter rejects a
  high-scoring local match over a short shared region (< 50 % coverage); **non-reciprocity**
  (A's best hit is B but B's best hit is C ≠ A) is *not* an ortholog pair — the guarded defect
  class; a gene with no above-threshold hit yields no pair; an empty genome yields no between-genome
  pair.

- **Datasets (documented oracles):**
  - *Reciprocity* — G1 {a1=`ACGTACGTACGTAC`, a2=`TTTTGGGGCCCCAA`}, G2 {b1, b2} identical by row
    → RBH {a1↔b1, a2↔b2} (identical sequences ⇒ maximal similarity each way; identity = coverage =
    1.0).
  - *Non-reciprocity* — a1 identical to b1; b2 = a1's superstring (`ACGTACGTACGTACGTACGT`, shares all
    of a1's 5-mers). a1→b1 (Jaccard 1.0 > b2), b1→a1, b2→a1 but a1↛b2 → RBH {a1↔b1} only,
    **RBH count = 1** (b2 excluded).
  - *Coverage / minimum-identity gate* — a pair below `minIdentity` or `minCoverage` yields no RBH
    pair even when it is the mutual top candidate; `minIdentity = 1.0` rejects every non-identical
    pair.

- **Test-coverage recommendations:** MUST — mutual best hits returned as RBH pairs; one-directional
  best hit NOT returned (count = 1 for the superstring set); returned pair carries the actual hit
  identity/coverage (1.0 for identical, not hardcoded); sub-threshold pairs rejected by the
  min-identity/coverage gate; empty genome → no pairs; null → `ArgumentNullException`. SHOULD — RBH
  is a partial matching (each gene at most once); empty-sequence genes skipped. COULD — determinism
  across input orderings.

## Deviations and assumptions

**One ASSUMPTION**, source-backed, not a correctness gap:

- **Alignment-free similarity as the best-hit ranking metric.** Moreno-Hagelsieb ranks best hits by
  BLAST bit-score, but the `ComparativeGenomics` class does not reference the
  `Seqeron.Genomics.Alignment` project, so candidates are ranked by **5-mer Jaccard similarity** — a
  monotone alignment-free measure (cf. Mash, Ondov et al. 2016). This affects *which* pair wins only
  among near-identical ties; the correctness-critical parts — the **RBH reciprocity rule**, the
  deterministic tie-break, the coverage gate (mapped to "shared k-mers ≥ 50 % of the smaller k-mer
  set"), and the minimum-similarity gate — are source-backed. Identical sequences score 1.0, so the
  metric is order-preserving for the datasets above.

No contradictions among the two sources — Tatusov (symmetrical/mutually-consistent best hits behind
COGs) and Moreno-Hagelsieb (operational RBH definition + coverage/E-value thresholds) are mutually
consistent, the latter being the pairwise operationalization of the former. The alignment-free metric
is the sole documented substitution; the reciprocity, coverage, and threshold semantics are unchanged.
