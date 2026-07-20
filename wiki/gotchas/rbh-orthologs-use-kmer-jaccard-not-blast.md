---
type: gotcha
title: "find_orthologs ranks best hits by 5-mer Jaccard, not BLAST bit-score — RBH results won't match a BLAST pipeline"
tags: [comparative-genomics, gotcha]
mcp_tools:
  - find_orthologs
  - find_reciprocal_best_hits
sources:
  - docs/algorithms/Comparative_Genomics/Reciprocal_Best_Hits.md
  - docs/algorithms/Comparative_Genomics/Ortholog_Identification.md
source_commit: 87f349f06987e26172ea87e365e9f25f38330997
created: 2026-07-20
updated: 2026-07-20
---

# find_orthologs ranks best hits by 5-mer Jaccard, not BLAST bit-score — RBH results won't match a BLAST pipeline

**The trap.** Canonical reciprocal-best-hits (Moreno-Hagelsieb 2008) ranks best hits by **BLAST
bit-score**. Seqeron's `ComparativeGenomics` has no alignment back end, so it ranks candidates by
**5-mer Jaccard similarity** — a monotone **alignment-free** metric (Mash-style). The coverage gate
is likewise expressed as "shared k-mers ≥ 50 % of the smaller k-mer set", not alignment coverage.

**Why it bites.** The **RBH reciprocity rule, coverage gate, and minimum-similarity threshold are
source-backed**, but the *ranking metric* is substituted. Among **near-identical** candidates (close
paralogs, recent duplicates) the pair that "wins" as the best hit can differ from a BLAST-RBH
pipeline, so your ortholog set won't reproduce a BLAST/DIAMOND result exactly. Sequence-identity
edge cases where k-mer similarity and alignment score diverge are where this shows.

**What to rely on instead.** Use `find_orthologs` / `find_reciprocal_best_hits` as a fast
alignment-free ortholog screen; for a BLAST-grade set, confirm ambiguous best-hit pairs with an
alignment-based reciprocal search. Identical sequences still score 1.0, so unambiguous orthologs are
unaffected. Full model: [[ortholog-detection-reciprocal-best-hits]]; the k-mer metric is
[[kmer-jaccard-similarity]]. Pipeline view: [[comparative-genomics-pipeline-silent-traps]].
