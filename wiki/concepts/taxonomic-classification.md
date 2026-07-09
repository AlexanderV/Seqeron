---
type: concept
title: "Taxonomic classification (Kraken k-mer / LCA / RTL read assignment)"
tags: [metagenomics, algorithm]
sources:
  - docs/Evidence/META-CLASS-001-Evidence.md
  - docs/algorithms/Metagenomics/Taxonomic_Classification.md
source_commit: a4a2f861b85acaf767a0e8992a75bb41cbc27227
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: meta-class-001-evidence
      evidence: "Test Unit ID: META-CLASS-001, Area: Metagenomics, Methods ClassifyReads/BuildKmerDatabase/TaxonomyTree.Lca"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:metagenomic-binning
      source: meta-class-001-evidence
      evidence: "Both assign community sequences to taxa: classification is per-read reference-database (Kraken k-mer/LCA) assignment; binning is unsupervised per-contig grouping into MAGs. Wikipedia Metagenomics frames both under sequence-to-organism assignment."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:alpha-diversity
      source: meta-class-001-evidence
      evidence: "Per-read classification produces the taxon→abundance profile that the within-sample diversity indices (alpha) and between-sample dissimilarities (beta) summarize."
      confidence: medium
      status: current
---

# Taxonomic classification (Kraken k-mer / LCA / RTL read assignment)

**Taxonomic classification** assigns a taxonomic label to each query read/contig by matching its
constituent **k-mers** against a reference database — the standard fast, alignment-free metagenomic read
assignment strategy. Seqeron implements the **faithful Kraken algorithm** (Wood & Salzberg 2014): a
taxonomy tree with a lowest-common-ancestor (LCA) operation, a canonical-k-mer → taxon database in which
a k-mer shared by several taxa is stored as their **LCA**, and a per-read **root-to-leaf (RTL)
maximum-weight-path** classifier. This is the **fourth ingested unit of the Metagenomics family**; its
siblings are [[metagenomic-binning]] (unsupervised per-contig grouping into MAGs) and the diversity
pair [[alpha-diversity]] / [[beta-diversity]] (which summarize the abundance profile classification
produces), and [[functional-prediction]] (which transfers *functional* annotation to genes and tests
pathway over-representation — the "what can they do" complement to this "who is there" assignment).
Validated under test unit **META-CLASS-001**; the record is [[meta-class-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the artifact
pattern.

This concept is deliberately scoped to the **per-read classification / LCA assignment**. Abundance
*profiling* (community composition / relative-abundance estimation) is the separate
[[taxonomic-profile]] unit (META-PROF-001), which aggregates *this* unit's per-read output
(`GenerateTaxonomicProfile` takes an `IEnumerable<TaxonomicClassification>`) into normalized
per-taxon relative abundances at four ranks.

## Database build (`BuildKmerDatabase`)

Each **canonical** k-mer stores the **LCA of all reference taxa whose genomes contain it** — Kraken's
"a k-mer and the LCA of all organisms whose genomes contain that k-mer." As references are processed, a
k-mer already owned by a taxon folds in the next owner via LCA: "if a k-mer … has had its LCA value
previously set, then the LCA of the stored value and the current sequence's taxon is calculated."

```
canonical(kmer) = min(kmer, reverse_complement(kmer))     (lexicographic; DnaSequence.GetReverseComplementString)
```

Canonical form is used **both** at build and query time. Default **k = 31** (Kraken standard).
Taxonomy strings are semicolon/pipe-delimited hierarchies
(`Bacteria|Proteobacteria|Gamma|Escherichia|coli`) parsed into rank levels.

## Per-read classification (`ClassifyReads`)

1. Extract all k-mers from the read; **skip ambiguous** (non-ACGT) k-mers — they are **not counted in
   Q** (`TotalKmers`).
2. Canonicalize each remaining k-mer and query the database.
3. Tally per-taxon hit counts. **No hits ⇒ the read is unclassified** (reported as the taxonomy root,
   `TaxonomyTree.RootId`, `C = 0`).
4. The hit taxa **plus their ancestors** form the **classification tree**; each node is weighted by its
   k-mer hit count.
5. Score every **root-to-leaf (RTL) path** as the sum of node weights along it. The **maximum-scoring
   path** is the classification path and its **leaf** is the assigned taxon.
6. **Tie-break:** if several paths share the maximum score, assign the **LCA of their leaves**
   (`TaxonomyTree.Lca`).
7. **Confidence = C / Q** (Kraken 2), where `C` = k-mers mapped to a taxon in the **clade rooted at the
   assigned label** and `Q` = non-ambiguous k-mers queried.

`TaxonomyTree.Lca` is the shared primitive: it supplies both the DB-build fold and the RTL tie-break.
**LCA correctness:** siblings → parent; ancestor/descendant → the ancestor; self → self; disjoint
lineages → root.

## Invariants and edge cases

- **Output count:** `|output| = |input reads|` (one classification per read).
- **Confidence range:** `0 ≤ Confidence = C/Q ≤ 1`.
- **Q (TotalKmers):** count of non-ambiguous k-mers (`= len − k + 1` for an all-ACGT read).
- **C (MatchedKmers) ≤ Q:** matched clade k-mers never exceed queried k-mers.
- **Canonical uniqueness:** exactly one canonical form per k-mer.
- **Unclassified path:** empty read, read shorter than `k` (no k-mers extractable), or zero database
  hits → **Unclassified** (root, `C = 0`). Empty read list → empty output; empty database → every read
  Unclassified.
- **Robustness:** mixed-case input is upper-cased internally; all-N / low-complexity runs are handled
  gracefully (production use should DUST-filter to avoid false-positive hits, Morgulis 2006).

Worked references (from [[meta-class-001-evidence]]): a k-mer shared by several taxa in the DB is stored
as their LCA; multiple taxon matches on a read → the leaf of the max-scoring root-to-leaf path, with
equal-score paths collapsing to the LCA of their leaves.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the Kraken exact-k-mer LCA/RTL
classifier. The Evidence file lists **no open questions** and **no literature deviations** — the
behaviour is fully pinned by the Kraken references (Wood & Salzberg 2014, Kraken 1/2 manuals). Scope is
**per-read classification only**: it does not perform abundance profiling (a separate unit), does not
ship a prebuilt reference database (the caller supplies labeled references to `BuildKmerDatabase`), and
its speed/sensitivity trade-off is the intrinsic exact-k-mer one (no gapped/minimizer/spaced-seed
extension). No source contradictions.
