---
type: concept
title: "Ortholog & in-paralog detection (Reciprocal Best Hits)"
tags: [comparative-genomics, algorithm]
mcp_tools:
  - find_orthologs
  - find_reciprocal_best_hits
sources:
  - docs/Evidence/COMPGEN-ORTHO-001-Evidence.md
  - docs/Evidence/COMPGEN-RBH-001-Evidence.md
  - docs/Validation/reports/COMPGEN-ORTHO-001.md
  - docs/Validation/reports/COMPGEN-RBH-001.md
  - docs/algorithms/Comparative_Genomics/Ortholog_Identification.md
  - docs/algorithms/Comparative_Genomics/Reciprocal_Best_Hits.md
source_commit: 00c5ea423943b75acdee59bd9ce88cc801bb2a37
created: 2026-07-09
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: compgen-ortho-001-evidence
      evidence: "Test Unit ID: COMPGEN-ORTHO-001 ... Algorithm: Ortholog identification by Reciprocal Best Hits (RBH); paralog (in-paralog) identification by within-genome best hits"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:genome-comparison-core-dispensable
      source: compgen-ortho-001-evidence
      evidence: "CompareGenomes' conserved/core set is exactly the RBH ortholog pairs of this unit; COMPGEN-COMPARE-001 composes RBH ortholog detection as its shared-gene sub-unit"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:average-nucleotide-identity
      source: compgen-ortho-001-evidence
      evidence: "Sibling Comparative-genomics unit quantifying two-genome relatedness — ANI at nucleotide identity, RBH at the gene-ortholog level; both are alignment-free-scored in Seqeron"
      confidence: medium
      status: current
---

# Ortholog & in-paralog detection (Reciprocal Best Hits)

The **Comparative-genomics** family's homology-classification unit: it decides which genes of two
genomes are **orthologs** (diverged by *speciation*) versus **paralogs** (diverged by *gene
duplication*), using the operational **Reciprocal Best Hits (RBH)** criterion for orthologs and a
**within-genome best-hit** rule for recent (in-)paralogs. It is the shared sub-unit behind the
[[genome-comparison-core-dispensable|genome-comparison pipeline]]'s core/dispensable partition and
is the shared anchor the **COMPGEN-RBH-001** unit reuses rather than re-deriving. Its
Comparative-genomics siblings are [[average-nucleotide-identity]] (nucleotide-level
relatedness), [[synteny-and-rearrangement-detection]] (gene-order conservation),
[[conserved-gene-clusters-common-intervals]] (order-free gene sets), and [[dot-plot-word-match]]
(visual comparison). Validated under **two test units**: the RBH-plus-in-paralog
**COMPGEN-ORTHO-001** (record [[compgen-ortho-001-evidence]]) and the RBH-only
**COMPGEN-RBH-001** (record [[compgen-rbh-001-evidence]]) — the latter is the between-genome
ortholog slice without the within-genome in-paralog rule. [[test-unit-registry]] tracks the units
and [[algorithm-validation-evidence]] describes the artifact pattern. The independent two-stage
re-validation of COMPGEN-ORTHO-001 ([[compgen-ortho-001-report]]) graded it **Stage A/B
PASS-WITH-NOTES · End state ✅ CLEAN** — no code defect (the historical non-reciprocal-`FindOrthologs`
bug is already fixed), the notes being the honestly-documented alignment-free simplifications below;
three test-quality improvements were applied in-session (including correcting a comment that wrongly
called TtBlock/GcBlock "Jaccard 0.0" when they share 8/12 5-mers = 0.5). The RBH-only unit's own
re-validation ([[compgen-rbh-001-report]]) graded it **Stage A PASS-WITH-NOTES · Stage B PASS · End
state ✅ CLEAN** — again no code defect; the only in-session work was closing two test-coverage gaps on
the `FindReciprocalBestHits` primitive (a missing `minCoverage`-gate test M7 with exact 6/16 = 0.375,
6/11 = 0.5455, and a missing `< k` short-sequence-edge test S4).

## Orthology vs paralogy (Fitch 1970)

Walter Fitch (1970) drew the distinction still used today:

- **Orthologs** — homologs whose split reflects a **speciation** event, so the gene history
  mirrors the species history (his example: hemoglobin in man and mouse). Orthologs sit in
  *different* genomes.
- **Paralogs** — homologs whose split reflects a **gene duplication**, both copies descending
  side-by-side within a lineage (his example: α- and β-hemoglobin). Paralogs sit *within* a
  genome/lineage.
- **In-paralogs vs out-paralogs** (Remm, Storm & Sonnhammer 2001) — in-paralogs duplicated
  *after* the speciation of interest (they are the recent, within-species best matches);
  out-paralogs predate it and are excluded.

## The RBH ortholog rule

The operational definition (Moreno-Hagelsieb & Latimer 2008; the symmetrical-best-hit rule of
Tatusov, Koonin & Lipman 1997 behind COGs): **two genes in two different genomes are orthologs iff
each is the other's best hit** in the opposite genome.

- **Best hit** = the candidate with the **maximum similarity score**, ties broken
  deterministically (Moreno-Hagelsieb ranks by descending BLAST bit-score, then ascending E-value;
  the implementation applies an analogous deterministic tie-break so the best hit is unique).
- **Reciprocity is required** — if A's best hit is B but B's best hit is C ≠ A, then A–B is **not**
  an ortholog pair. A one-directional best hit is not an ortholog. This is the defect class the
  unit guards against.
- **Gates** — Moreno-Hagelsieb requires **≥ 50 % coverage** of one of the aligned sequences (a
  short high-scoring local match over a shared domain is rejected) and a maximum **E-value 1×10⁻⁶**
  significance gate. In the alignment-free implementation these map to a shared-k-mer coverage gate
  and a minimum-similarity gate.

## The in-paralog rule (within-genome best hits)

Recent paralogs are found by the same best-hit machinery turned **inward**: a within-genome pair of
mutual best hits (Remm et al. 2001: within-species hits that are reciprocally better than the
between-species hit to the seed ortholog). The within-genome best-hit rule alone identifies the
closest within-genome relatives (recent paralogs) and is the documented operational proxy;
out-paralogs (pre-speciation duplicates) are not in-paralogs.

## Outputs and invariants

- Ortholog output is a **partial matching** between the two genomes — each gene is paired at most
  once per direction (RBH yields a matching, not a many-to-many map); every pair is symmetric.
- A gene with **no above-threshold hit** in the other genome yields no ortholog pair.
- **Determinism** — the same input yields an identical pair set across runs (order-independent).
- Genes with **empty/absent sequence** are skipped (similarity undefined without sequence);
  null inputs throw `ArgumentNullException` (repository contract).

## Documented oracles

- **Reciprocity** — G1 {a1, a2}, G2 {b1, b2} with a1=b1 and a2=b2 (identical sequences) →
  orthologs {a1↔b1, a2↔b2} (mutual best hits, maximal similarity each way).
- **Non-reciprocity** — a1 (G1) identical to b1 (G2); b2 (G2) is a1's superstring sharing all its
  5-mers. a1's best hit is b1 (Jaccard 1.0 > b2's), b1's best hit is a1, but b2's best hit is a1
  (not reciprocated) → orthologs {a1↔b1} only, **RBH count = 1** (b2 excluded).
- **In-paralog** — G1 {p1, p2 (duplicate of p1), q1 (unrelated)} → within-genome paralogs {p1↔p2}
  (mutual best hits, similarity 1.0); q1 has no qualifying within-genome partner.
- **Empty / single-gene** — an empty genome yields no orthologs (no between-genome pair possible);
  a single-gene genome yields no paralogs (no within-genome pair possible).

## Assumption (source-backed, not a correctness gap)

**Alignment-free similarity replaces the BLAST bit-score ranking.** Moreno-Hagelsieb (2008) ranks
best hits by BLAST bit-score, but the `ComparativeGenomics` class does not reference the alignment
project, so candidates are ranked by **5-mer Jaccard similarity** — a monotone alignment-free
similarity (cf. Mash, Ondov et al. 2016). This affects *which* pair wins only among near-identical
sequences; the correctness-critical parts — the **RBH reciprocity rule**, the coverage gate
(mapped to "shared k-mers ≥ 50 % of the smaller k-mer set"), and the minimum-similarity threshold —
are source-backed. The metric is order-preserving on the test datasets (identical sequences score
1.0). **Deviations: none** beyond this documented metric substitution.

## Reference tools

Definitions trace to **Fitch 1970** (orthology/paralogy, *Systematic Zoology* 19:99–113),
**Tatusov, Koonin & Lipman 1997** (COG symmetrical best hits, *Science* 278:631–637),
**Moreno-Hagelsieb & Latimer 2008** (operational RBH definition + coverage/E-value thresholds,
*Bioinformatics* 24:319–324), **Remm, Storm & Sonnhammer 2001** (InParanoid RBH seed + in-paralog
rule, *J. Mol. Biol.* 314:1041–1052) with OrthoMCL corroboration (Li et al. 2003), and **Ondov et
al. 2016** (Mash, the alignment-free k-mer-similarity basis). No source contradictions — the four
orthology sources are mutually consistent, each governing a distinct part of the rule.

**Sharp edge:** [[rbh-orthologs-use-kmer-jaccard-not-blast]] — best hits ranked by **5-mer Jaccard**, not BLAST bit-score — RBH won't match a BLAST pipeline.
