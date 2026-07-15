---
type: source
title: "Evidence: COMPGEN-ORTHO-001 (Ortholog detection by Reciprocal Best Hits + in-paralog within-genome best hits)"
tags: [validation, comparative-genomics]
doc_path: docs/Evidence/COMPGEN-ORTHO-001-Evidence.md
sources:
  - docs/Evidence/COMPGEN-ORTHO-001-Evidence.md
source_commit: 29ce869c75abab6c57f5482a2e2ec51750607d3c
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: COMPGEN-ORTHO-001

The validation-evidence artifact for test unit **COMPGEN-ORTHO-001** — **ortholog identification by
Reciprocal Best Hits (RBH)** and **in-paralog identification by within-genome best hits**. This is a
**Comparative-genomics** family Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm, its rules, invariants,
worked oracles, and corner cases are summarized in the dedicated concept
[[ortholog-detection-reciprocal-best-hits]]. This unit is the shared RBH/ortholog sub-unit behind
the [[genome-comparison-core-dispensable]] pipeline (whose conserved/core set *is* these RBH pairs)
and the anchor the future **COMPGEN-RBH-001** unit reuses. Its sibling COMPGEN units are
[[average-nucleotide-identity]], [[synteny-and-rearrangement-detection]],
[[conserved-gene-clusters-common-intervals]], and [[dot-plot-word-match]]. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (all authority rank 1 unless noted):**
  - **Fitch 1970** (*Systematic Zoology* 19:99–113, quoted verbatim via the Koonin 2011 PMC3178060
    obituary) — the founding **orthology/paralogy** definitions: orthology = homology from
    *speciation* (gene history mirrors species history); paralogy = homology from *gene
    duplication* (copies descend side-by-side within a lineage). Consequence: orthologs split
    between genomes, paralogs split within a genome.
  - **Tatusov, Koonin & Lipman 1997** (*Science* 278:631–637) — COGs built from **reciprocal /
    symmetrical best BLAST hits**; each COG is orthologous proteins *or orthologous sets of
    paralogs* across ≥ 3 lineages (recent paralogs co-clustered with the ortholog). Full text was
    403/404-blocked; the method statement is from the search summary + a scirp citation page, DOI
    10.1126/science.278.5338.631 confirmed.
  - **Moreno-Hagelsieb & Latimer 2008** (*Bioinformatics* 24:319–324) — the **operational RBH
    definition**: two genes in two genomes are orthologs iff their products *"find each other as
    the best hit in the opposite genome."* Best hit = descending bit-score, ties by ascending
    E-value ⇒ **maximum-similarity candidate, deterministic tie-break**. Gates: **≥ 50 % coverage**
    of one aligned sequence, **max E-value 1×10⁻⁶** (a significance gate → minimum-similarity gate
    in alignment-free space).
  - **Remm, Storm & Sonnhammer 2001** (*J. Mol. Biol.* 314:1041–1052, InParanoid; Frontiers review
    PMC5674930 corroboration) — an ortholog cluster is **seeded by a between-genome RBH**;
    **in-paralogs** are within-species hits reciprocally better than between-species hits (recent,
    post-speciation duplicates), out-paralogs excluded.

- **Documented corner cases / failure modes:** tie in best hit → deterministic tie-break so the
  best hit is unique; coverage filter rejects a high-scoring match over a short shared region;
  **non-reciprocity** (A→B but B→C≠A) is *not* an ortholog pair (the guarded defect class); a gene
  with no above-threshold hit yields no pair; empty genome → no between-genome pair, single-gene
  genome → no within-genome paralog; out-paralogs (pre-speciation within-genome pairs) are not
  in-paralogs.

- **Datasets (documented oracles):**
  - *Reciprocity* — G1 {a1=`ACGTACGTACGTAC`, a2=`TTTTGGGGCCCCAA`}, G2 {b1, b2} identical by column
    → orthologs {a1↔b1, a2↔b2} (mutual best hits).
  - *Non-reciprocity* — a1 identical to b1; b2 = a1's superstring (`ACGTACGTACGTACGTACGT`, shares
    all a1's 5-mers). a1→b1 (Jaccard 1.0 > b2), b1→a1, b2→a1 but a1↛b2 → orthologs {a1↔b1} only,
    **RBH count = 1**.
  - *In-paralog* — G1 {p1=`GGGGCCCCAAAATT`, p2 = duplicate of p1, q1=`ACGTACGTACGTAC` unrelated}
    → within-genome paralogs {p1↔p2} (similarity 1.0); q1 has no partner.

- **Test-coverage recommendations:** MUST — RBH mutual best hits returned; one-directional best hit
  NOT returned; within-genome mutual best hits returned as paralogs (unrelated gene excluded);
  coverage/min-identity threshold rejects sub-threshold pairs; empty→no orthologs / single-gene→no
  paralogs; null → `ArgumentNullException`. SHOULD — pairs symmetric / no gene paired twice one
  direction; empty-sequence genes skipped. COULD — determinism across runs.

## Deviations and assumptions

**One ASSUMPTION**, source-backed, not a correctness gap:

- **Alignment-free similarity as the best-hit ranking metric.** Moreno-Hagelsieb ranks by BLAST
  bit-score, but `ComparativeGenomics` does not reference the `Seqeron.Genomics.Alignment` project,
  so candidates are ranked by **5-mer Jaccard similarity** — a monotone alignment-free measure (cf.
  Mash, Ondov et al. 2016). This affects *which* pair wins ties only among near-identical
  sequences; the correctness-critical parts — the **RBH reciprocity rule**, the coverage gate
  (mapped to "shared k-mers ≥ 50 % of the smaller k-mer set"), and the minimum-similarity gate —
  are source-backed, and the metric is order-preserving for the datasets above (identical sequences
  score 1.0).

No contradictions among sources — Fitch (orthology/paralogy definitions), Tatusov (symmetrical best
hits), Moreno-Hagelsieb (operational RBH + thresholds), and Remm et al. (RBH seed + in-paralog rule)
are mutually consistent, each governing a distinct part of the rule. The alignment-free metric is
the sole documented substitution; the paper reciprocity/coverage/threshold semantics are unchanged.
