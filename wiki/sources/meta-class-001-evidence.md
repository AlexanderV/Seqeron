---
type: source
title: "Evidence: META-CLASS-001 (taxonomic classification — Kraken k-mer / LCA / RTL)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-CLASS-001-Evidence.md
sources:
  - docs/Evidence/META-CLASS-001-Evidence.md
source_commit: a4a2f861b85acaf767a0e8992a75bb41cbc27227
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-CLASS-001

The validation-evidence artifact for test unit **META-CLASS-001** — **k-mer based taxonomic
classification**, the per-read assignment of taxonomic labels computed by
`MetagenomicsAnalyzer.ClassifyReads` over a `BuildKmerDatabase` index and a `TaxonomyTree` (+ `Lca`).
Fourth ingested unit of the Metagenomics family, alongside the diversity siblings
[[meta-alpha-001-evidence|META-ALPHA-001]] / [[meta-beta-001-evidence|META-BETA-001]] and the binning
sibling [[meta-bin-001-evidence|META-BIN-001]]. One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The method is synthesized in its own
concept, [[taxonomic-classification]]; [[test-unit-registry]] tracks the unit. See
`docs/Evidence/META-CLASS-001-Evidence.md`.

## What this file records

The unit implements the **faithful Kraken algorithm** (Wood & Salzberg 2014) — a lowest-common-ancestor
(LCA) database + a root-to-leaf (RTL) maximum-weight-path classifier. The Evidence file notes that a
**C1 enhancement** replaced an earlier flat best-hit classifier (no LCA); that obsolete wording is
superseded.

- **Online sources (mutually consistent):**
  - **Wood & Salzberg (2014), Genome Biology 15:R46 (DOI 10.1186/gb-2014-15-3-r46)** — the primary
    reference describing exact-k-mer metagenomic classification.
  - **Kraken 1 manual (CCB, JHU)** — canonical implementation: default **k = 31**, canonical k-mers in
    both DB build and query, ambiguous (non-ACGT) k-mers not queried, confidence **C/Q**.
  - **Kraken 2 manual (GitHub wiki)** — confirms the score fraction **C/Q**: C = k-mers mapped to LCA
    values in the clade rooted at the assigned label, Q = k-mers in the sequence lacking an ambiguous
    nucleotide (i.e. actually queried).
  - **Wikipedia — Metagenomics** — defines composition- vs similarity-based binning and taxonomic
    profiling context.

## Algorithm (from the Evidence file)

**Database build (`BuildKmerDatabase`).** Each canonical k-mer stores the **LCA of all reference taxa
that contain it**: "if a k-mer … has had its LCA value previously set, then the LCA of the stored value
and the current sequence's taxon is calculated." Canonical form `= min(kmer, revcomp(kmer))`.

**Per-read classification (`ClassifyReads`).**
1. Extract all k-mers; skip those with ambiguous (non-ACGT) nucleotides (**not counted in Q**).
2. Canonicalize each k-mer and query the database.
3. Tally per-taxon hit counts; **no hits ⇒ unclassified (taxonomy root)**.
4. Hit taxa + their ancestors form the **classification tree**, each node weighted by its k-mer count.
5. Score every **root-to-leaf (RTL)** path as the sum of node weights; the **max-scoring path's leaf**
   is the assigned taxon.
6. **Tie-break:** equal-max paths → assign the **LCA of their leaves** (`TaxonomyTree.Lca`).
7. **Confidence = C/Q** (Kraken 2).

Taxonomy strings are semicolon/pipe-delimited hierarchies
(`Bacteria|Proteobacteria|Gamma|Escherichia|coli`) parsed into rank levels; sequences are upper-cased
internally.

## Source-verified invariants and oracles

**Invariants (from the Evidence file):**
- **Output count:** `|output| = |input reads|`.
- **Confidence range:** `0 ≤ C/Q ≤ 1`; `Confidence = C/Q`.
- **Q (TotalKmers):** count of non-ambiguous k-mers (`= len − k + 1` for all-ACGT reads).
- **C (MatchedKmers) ≤ Q:** C = k-mers in the clade rooted at the assigned label.
- **Unclassified:** no k-mer hits ⇒ `TaxonId = TaxonomyTree.RootId`, `C = 0`.
- **Canonical uniqueness:** exactly one canonical form per k-mer.
- **DB-build LCA:** a canonical k-mer owned by several taxa is stored as their LCA.
- **RTL assignment:** assigned taxon = leaf of the max-scoring root-to-leaf path; ties → LCA of tied
  leaves.
- **LCA correctness:** siblings → parent; ancestor/descendant → ancestor; self → self; disjoint → root.

**Documented test cases / edge cases:** empty sequence → Unclassified; read shorter than k → no k-mers
→ Unclassified; no matching k-mers → Unclassified; all-N / low-complexity → handled gracefully (DUST
filtering recommended in production, Morgulis 2006); empty read list → empty output; empty database →
all reads Unclassified; mixed-case input upper-cased.

**Test datasets:** Kraken accuracy set (three FASTA files, 10,000 simulated 100 bp reads each, 2.1% SNP
/ 1.1% indel). Testing methodology from the literature = precision/sensitivity at genus level over
simulated metagenomes of known composition.

## Deviations

The Evidence file lists **no open questions** and no literature deviations — behaviour is well-defined
by the Kraken references. The only historical note is the superseded pre-C1 flat best-hit wording (now
obsolete; the current unit is the faithful LCA/RTL Kraken algorithm).
