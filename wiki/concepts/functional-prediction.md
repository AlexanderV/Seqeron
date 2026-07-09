---
type: concept
title: "Functional prediction (homology-based annotation transfer + hypergeometric pathway ORA)"
tags: [metagenomics, algorithm]
sources:
  - docs/Evidence/META-FUNC-001-Evidence.md
  - docs/algorithms/Metagenomics/Functional_Prediction.md
source_commit: ab6548a3b59cfe5105eb5aa7cb252e63ffed5de7
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: meta-func-001-evidence
      evidence: "Test Unit ID: META-FUNC-001, Area: Metagenomics, Methods MetagenomicsAnalyzer.PredictFunctions + FindPathwayEnrichment"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:pathway-enrichment-ora
      source: meta-func-001-evidence
      evidence: "Component B FindPathwayEnrichment / HypergeometricUpperTail is the shared ORA machinery owned by pathway-enrichment-ora (META-PATHWAY-001, its dedicated unit)"
      confidence: high
      status: current
---

# Functional prediction (homology-based annotation transfer + pathway ORA)

**Functional prediction** answers "what can this community *do*?": it assigns biological function to a
metagenome's predicted genes/proteins and then asks which metabolic **pathways** are statistically
over-represented among them — the PICRUSt / KO-pathway style of turning sequence into functional
capability. Seqeron implements two specification-driven, deterministic pieces:

1. **`PredictFunctions`** — **homology-based annotation transfer**: match a query protein against a
   signature database and transfer the annotation (Function, Pathway, KO) of the **best hit**, reporting
   significance as a **BLAST bit score** and **E-value** from Karlin-Altschul statistics.
2. **`FindPathwayEnrichment`** — **pathway over-representation analysis (ORA)**: score each pathway with
   the **hypergeometric** right-tail test.

This is the **fifth ingested unit of the Metagenomics family**; its siblings are the diversity pair
[[alpha-diversity]] / [[beta-diversity]], the per-read reference assignment [[taxonomic-classification]],
and the genome-reconstruction unit [[metagenomic-binning]]. Where classification and binning ask **who is
there**, functional prediction asks **what functions they encode**. Validated under test unit
**META-FUNC-001**; the record is [[meta-func-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

The numerical core (bit score, E-value, hypergeometric p-value) is **exact** with respect to the cited
formulas; only the *matching* step is simplified (exact-signature occurrence, not gapped alignment).

## A — Homology-based annotation transfer (`PredictFunctions`)

For an ungapped alignment of raw score `S`, the expected number of high-scoring segment pairs (HSPs)
with score ≥ S is the **Karlin-Altschul E-value**, and the normalized **bit score** `S'` and its
E-value form follow (Altschul et al., NCBI BLAST statistics tutorial):

```
E  = K · m · n · e^(−λ·S)          expected #HSPs ≥ S; m, n = sequence lengths
S' = (λ·S − ln K) / ln 2           normalized bit score (standard units)
E  = m · n · 2^(−S')               bit-score → E-value (K, λ drop out)
```

- **Ungapped BLOSUM62 parameters** (NCBI `blast_stat.c` `blosum62_values`, gap-open = gap-extend =
  INT2_MAX row): **`λ = 0.3176`, `K = 0.134`** (H = 0.4012). Using these ungapped constants is exactly
  why the matching step is restricted to ungapped exact matches (ASM-01) — no gap model is needed.
- **Raw score** of a matched segment = the sum of the **BLOSUM62 diagonal** (self-match) scores over its
  residues (W 11, C 9, A 4, …), upper-casing residues and treating unknown residues as score 0.
- **Best-hit rule (INV-06):** among all database signatures occurring exactly in a protein, transfer the
  annotation of the one with the **lowest E-value** (equivalently the **highest bit score**). E-value
  ranks significance.

**Matching step (simplified, ASM-01):** a database signature is a hit iff it occurs as an **exact,
ordinal `string.Contains` substring** of the query protein. The repository suffix tree was *evaluated and
not used* — each gene is scanned once against the database (a short search per gene × signature), so the
operative cost is the per-residue BLOSUM62 scoring, not occurrence enumeration. Consequence: only verbatim
sub-sequence hits are found; **divergent homologs a full Smith-Waterman/BLAST search would detect are
missed**. `CogCategory` is a keyword-inferred label (a convenience, not a scoring input; does not affect
the bit score, E-value, or which hit is selected).

Note this BLAST **significance** model (Karlin-Altschul E-value / bit score of a local match) is a
different layer from the post-alignment percent-identity/similarity metrics captured in
[[alignment-statistics]].

## B — Pathway over-representation (`FindPathwayEnrichment`)

Given a **background universe** of `N` genes containing `M` members of a pathway, and a **query** of `n`
genes of which `x` fall in the pathway, the overlap `X` is hypergeometric (sampling `n` without
replacement from `N`) and the right-tail over-representation p-value is (PNNL ORA §8.2):

```
P(X ≥ x) = 1 − Σ_{i=0}^{x−1}  C(M,i) · C(N−M, n−i) / C(N, n)
```

- The tail is summed in **log-space** via a Lanczos **log-Gamma** to avoid factorial overflow (validated
  to `N = 8000`); the final p is clamped to `[0,1]`; pathways are returned **sorted ascending by p-value**.
- **Background / query handling:** an explicit `backgroundGenes` sets `N`; when null/empty, `N` defaults
  to the **union of all pathway members**. Pathway members are **intersected with the background** before
  counting (so `M` and `x` are measured against the actual sampled universe), and the **query is always
  unioned into the background**.
- **Not implemented:** multiple-testing correction (BH/Bonferroni) — apply your own over the returned
  p-values.

> **This ORA half has its own dedicated unit, now ingested.** `FindPathwayEnrichment` /
> `HypergeometricUpperTail` is shared with test unit **META-PATHWAY-001**
> (`docs/algorithms/Metagenomics/Pathway_Enrichment_ORA.md`), which frames ORA standalone and **owns** the
> statistic on [[pathway-enrichment-ora]] — see that page for the M↔n symmetry invariant, the exact
> rational oracles (1/252, 5/6, 251/252), and the GSEA/FDR scope notes. META-FUNC-001 validates only
> `Functional_Prediction.md`; the pathway ORA details above are the shared machinery synthesized in full
> on [[pathway-enrichment-ora]].

## Invariants and edge cases

- **INV-01:** bit score `S'` is strictly increasing in raw score `S` (slope `λ/ln2 > 0`).
- **INV-02:** `K·m·n·e^(−λS) = m·n·2^(−S')` (the two E-value forms are algebraically identical).
- **INV-03:** E-value is strictly decreasing in `S` (monotonicity: larger score ⇒ larger bit ⇒ smaller E).
- **INV-04:** hypergeometric p-value ∈ `[0,1]`.
- **INV-05 / INV-06:** ORA `p = 1` when `x = 0`, `M = 0`, or `n = 0` (empty sum / degenerate margins);
  best-hit = lowest-E-value annotation transferred.
- **Robustness:** empty/whitespace protein sequences and empty signatures are skipped; a gene matching no
  signature yields no annotation; empty/null function database → no annotations; duplicate query genes are
  counted once. Null `proteins`/`functionDatabase`/`queryGenes`/`pathwayDatabase` → `ArgumentNullException`.

Worked oracles (from [[meta-func-001-evidence]]): signature `"WWW"` → `S = 33`, bit `S' =
18.0202932787533`, `E = 3.3852730346546 × 10⁻⁵` (m = n = 3, both forms agree); ORA `N=8000, M=400,
n=100, x=20` → `P(X ≥ 20) = 7.88 × 10⁻⁸`.

## Scope and limitations

A [[research-grade-limitations|research-grade]] (*Simplified*) functional-prediction unit: the bit-score /
E-value / hypergeometric numerics are exact and source-backed, but annotation transfer uses **exact
substring** signature matching rather than gapped local alignment — no divergent-homolog recovery, no
gapped `λ/K`, no composition-adjusted statistics, no HMM/Pfam domain modeling, and no FDR control across
pathway p-values. For production-scale annotation and multiple-testing control, external tools (NCBI
BLAST+, DIAMOND, clusterProfiler) apply. The bit-score model assumes the ungapped BLOSUM62 scoring system;
another matrix would need its own `λ, K`. The single accepted deviation (ASM-01) affects which hits are
found, not the significance formulas. No source contradictions.
</content>
