---
type: concept
title: "Antibiotic-resistance gene detection (ResFinder-style identity/coverage dual-threshold best-match screen)"
tags: [metagenomics, algorithm]
sources:
  - docs/Evidence/META-RESIST-001-Evidence.md
  - docs/algorithms/Metagenomics/Antibiotic_Resistance_Detection.md
source_commit: c81ef58a4d2d8fd4a9ceb1e322d2c6a1ee237cfc
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: meta-resist-001-evidence
      evidence: "Test Unit ID META-RESIST-001, Area Metagenomics, Method MetagenomicsAnalyzer.FindAntibioticResistanceGenes"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:functional-prediction
      source: meta-resist-001-evidence
      evidence: "Both are BLAST-style homology screens against a caller-supplied reference DB with a best-hit rule; AMR detection scores a nucleotide identity/coverage sliding match, functional-prediction a BLOSUM62 protein bit-score/E-value — shared homology-search machinery, distinct scoring layer (Zankari 2012 best-matching gene; Heng Li 2018 identity)"
      confidence: high
      status: current
---

# Antibiotic-resistance gene detection (ResFinder-style)

**Antibiotic-resistance (AMR) gene detection** answers "does this community carry a known resistance
gene?": it screens assembled contigs against a caller-supplied reference database of acquired resistance
genes (CARD / ResFinder-style) and reports, per contig, the single **best-matching** reference gene whose
**BLAST percent identity** and **coverage of the reference gene** both clear user-selectable thresholds —
the ResFinder acquired-gene detection method (`MetagenomicsAnalyzer.FindAntibioticResistanceGenes`).

This is the **seventh ingested unit of the Metagenomics family**. Where [[taxonomic-classification]] and
[[metagenomic-binning]] ask **who is there**, [[functional-prediction]] asks **what functions they
encode**, and the [[alpha-diversity]] / [[beta-diversity]] pair asks **how diverse** the community is, AMR
detection asks the sharp clinical-epidemiology question **which specific known resistance genes are
present**. Validated under test unit **META-RESIST-001**; the record is [[meta-resist-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the artifact
pattern.

The numerical core (identity, coverage, the dual-threshold reporting rule, best-match selection) is
**exact** with respect to the cited ResFinder / BLAST definitions; only the *alignment* step is simplified
(gapless ungapped sliding match, not gapped BLAST).

## Core model

For a contig *C* (length *n*) and a reference gene *R* (length *m*):

- **Percent identity** (BLAST definition, Heng Li 2018): identical positions ÷ alignment columns. For a
  gapless alignment there are no gap columns, so the
  denominator is the aligned **window** length *w*: `identity = matches / w`. This is the same
  denominator-includes-all-columns rule captured in [[alignment-statistics]], specialized to the gapless
  case.
- **Coverage** — measured against the **reference** gene, not the contig: `coverage = w / m`, the fraction
  of the reference gene length spanned by the alignment. Because coverage is referenced to *m*, a gene
  truncated at a contig edge or split across two contigs is still scored and detectable (the reason the
  floor is < 1).
- **Reporting rule (INV-02):** report *R* for *C* only if `identity ≥ idThreshold` **AND**
  `coverage ≥ covThreshold`.
- **Best-matching gene (INV-03):** when several references pass, only the **best match per contig** is
  returned — max identity, ties broken by greater coverage. Mirrors ResFinder's "best-matching gene"
  single-output convention and CARD RGI's best-hit ranking.

### The `BestUngappedMatch` scan

For each (contig, reference) pair the reference is slid across the contig at **every offset from
`-(m-1)` to `n-1`** — overhanging both ends so contig-edge / truncated genes are scored against the full
reference length. At each offset the identical positions over the overlapping window are counted, and the
offset with the **most matches** is kept; ties are broken toward the **shorter** (higher-identity) window,
so the chosen alignment is never padded with flanking mismatches — padding would dilute identity and could
spuriously fail the identity threshold even when a perfect high-scoring segment exists (mirroring BLAST
reporting the best-scoring HSP). Cost is `O(n·m)` per pair, `O(c·d·n·m)` over *c* contigs × *d* references.

## Thresholds

- **Defaults:** `identityThreshold = 0.90`, `coverageThreshold = 0.60`
  (`DefaultResistanceIdentityThreshold` / `DefaultResistanceCoverageThreshold`), both named constants and
  user-selectable.
- **Provenance / operating points:** the ResFinder GitHub service ships `-t 0.80` / `-l 0.60`; Zankari
  (2012) *selected* **98% identity** for their study because lower thresholds "give too much noise (e.g.
  fragments of genes)", with a **2/5 (later 60%) coverage** floor so edge/split genes are not missed
  (Sci Rep 2023, JAC 2016). The shipped 0.90 default is a mid-point between the 0.80 service default and
  the 0.98 study threshold — see the threshold-provenance note in [[meta-resist-001-evidence]].

## Relationship to functional prediction — shared machinery, different scoring layer

AMR detection and [[functional-prediction]] are **both** BLAST-style homology screens: each takes a
caller-supplied reference/signature database, finds where a database entry occurs in the query, and applies
a **best-hit rule**. But the scoring layers differ:

| | AMR detection (this unit) | [[functional-prediction]] `PredictFunctions` |
|---|---|---|
| Query / reference | **nucleotide** contig vs resistance-gene DB | **protein** gene vs signature DB |
| Match model | mismatch-tolerant gapless sliding window | exact-substring (`string.Contains`) occurrence |
| Score | percent **identity** + **coverage** of reference | BLOSUM62 self-score → **bit score** + **E-value** |
| Report gate | dual threshold (`id ≥ t` AND `cov ≥ t`) | lowest E-value (highest bit score) |
| Best hit | max identity, tie → max coverage | min E-value |

So AMR is the **similarity/coverage** flavour of the homology screen, functional prediction the
**significance-statistic** flavour. Both are distinct from the post-alignment percent-identity/similarity
metrics of [[alignment-statistics]], though AMR's identity denominator is the same rule.

## Invariants and edge cases

- **INV-01:** `0 ≤ identity ≤ 1` and `0 ≤ coverage ≤ 1` (matches ≤ window ≤ m).
- **INV-02:** a hit is reported only if identity ≥ idThreshold AND coverage ≥ covThreshold.
- **INV-03:** at most one hit per contig (best match: max identity, tie → max coverage).
- **INV-04:** exact full-length match ⇒ identity = 1.0, coverage = 1.0 (CARD "Perfect" match).
- **INV-05:** default thresholds = 0.90 identity, 0.60 coverage.
- **Robustness:** empty contig sequences and empty-sequence reference genes are skipped; a contig where no
  reference passes yields no hit; comparison is **case-sensitive**, nucleotide-only, no T↔U conversion
  (caller normalizes). Null `contigs`/`referenceGenes` → `ArgumentNullException`; a threshold outside
  `[0,1]` → `ArgumentOutOfRangeException`.

Worked oracles (from [[meta-resist-001-evidence]]): `CGTACGT` (m=7) in `AAACGTACGT` → identity 7/7 = 1.0,
coverage 7/7 = 1.0; `CGTTCGT` vs `CGTACGT` → identity 6/7 ≈ 0.857142857, coverage 1.0; contig ending
`CGTA` vs `CGTACGT` → identity 4/4 = 1.0, coverage 4/7 ≈ 0.571428571.

## Scope and limitations

A [[research-grade-limitations|research-grade]] (*Simplified*) detector: the identity / coverage / dual
threshold / best-match numerics are exact and source-backed, but the alignment is **gapless ungapped**
sliding rather than full gapped BLAST — the single accepted deviation (**ASM-01**), which under-scores
genes whose true alignment needs indels while scoring substitution divergence and contig-edge truncation
exactly. No SNP / protein-homology models (CARD's homolog / variant models are out of scope), no bundled
curated gene catalogue (the caller supplies a ResFinder/CARD-derived reference set — curated tables are not
fabricated), nucleotide-only case-sensitive comparison. A legacy `FindResistanceGenes` motif-containment
stub is retained separately for an existing MCP tool and is **not** part of this unit. No source
contradictions in the algorithm itself; a threshold-default provenance note (evidence records 0.80/0.98,
implementation ships 0.90) is flagged on [[meta-resist-001-evidence]].
