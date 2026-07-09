---
type: source
title: "Evidence: EPIGEN-METHYL-001 (methylation context classification CpG/CHG/CHH + profile)"
tags: [validation, epigenetics]
doc_path: docs/Evidence/EPIGEN-METHYL-001-Evidence.md
sources:
  - docs/Evidence/EPIGEN-METHYL-001-Evidence.md
source_commit: f59358fda234a26aa74779bcfa3844f8fbdae7f8
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: EPIGEN-METHYL-001

The validation-evidence artifact for test unit **EPIGEN-METHYL-001** — **methylation site detection**,
**sequence-context classification** into **CpG / CHG / CHH**, and the per-context **methylation
profile**. This is the **sixth (final) ingested unit of the Epigenetics family** and one instance of
the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct
context-classification method is synthesized in its own concept,
[[methylation-context-classification]]; the shared Schultz weighted profile is covered by
[[bisulfite-methylation-calling]]. [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (all mutually consistent, no contradictions):**
  - **Cornish-Bowden A (1985)** IUPAC-IUB nucleotide nomenclature, *Nucleic Acids Research*
    13:3021–3030 (rank 2, via the Los Alamos HIV DB table) — the ambiguity code **H = A or C or T**
    ("not G"), the exact predicate that separates CHG/CHH from CpG.
  - **Krueger F & Andrews SR (2011)** Bismark, *Bioinformatics* 27(11):1571–1572 (rank 3, reference
    implementation) — methylation calls "discriminate between cytosines in CpG, CHG and CHH context",
    "H can be either A, T or C"; **CpG and CHG are symmetric, CHH asymmetric**; strand-specific output
    for hemi-/CHH methylation.
  - **Lister R et al. (2009)** Human DNA methylomes at base resolution, *Nature* 462:315–322 (rank 1,
    primary) — defines non-CG contexts **mCHG and mCHH (H = A, C or T)**; **prevalence realism**: H1 ES
    cells ~25% of methylation is non-CG, IMR90 somatic cells **99.98% CG-context** (non-CG is
    stem-cell/plant, essentially absent in differentiated cells).
  - **Schultz MD, Schmitz RJ & Ecker JR (2012)** 'Leveling' the playing field, *Trends in Genetics*
    28(12):583–585 (rank 1, defines the metric) — **weighted methylation level = Σ(methylated reads) /
    Σ(total reads)** per sequence context, computed separately for CG/CHG/CHH (reproduced in
    PMC6686912, rank 3). Per-cytosine level = methylated / total reads ∈ [0, 1].

- **Documented corner cases / failure modes:** H excludes G — `CGG` middle is G → **CpG not CHG**;
  CHH is **asymmetric** (no symmetric complement, unlike CpG/CHG); **incomplete trailing context** —
  a terminal C with a truncated 3' window cannot be CHG-vs-CHH classified, but a terminal `CG` is still
  an unambiguous CpG (only two bases needed); **non-ACGT / ambiguous** base in the H or third position
  → context undetermined, cytosine not classified.

- **Datasets (documented oracles):**
  - **Context classification** — `CGACAGCAA`: C@0 `CG…`→**CpG**, C@3 `CAG`→**CHG** (H=A then G),
    C@6 `CAA`→**CHH** (H=A, H=A); 3 classifiable C sites, 0-based positions.
  - **Weighted methylation level** — Site A (meth 8 / total 10 → 0.8) + Site B (meth 2 / total 10 →
    0.2): weighted CpG = (8+2)/(10+10) = **0.5**; here equal to the unweighted mean of fractions
    (0.8+0.2)/2 = 0.5 (they differ only under unequal coverage).
  - **Biological prevalence** — IMR90 99.98% CG-context; H1 ES ~25% non-CG (context-realism check).

## Deviations and assumptions

- **ASSUMPTION 1 — sequence-only level defaults to 0.** `FindMethylationSites(sequence)` has no
  bisulfite read data, so each site's `MethylationLevel` and `Coverage` are 0 until measured — a
  representational default (the site is *potentially* methylatable), **not** a claim that the level is
  zero. Real levels come from `CalculateMethylationFromBisulfite`
  ([[bisulfite-methylation-calling]], out of this unit's scope) or caller-supplied sites.
- **ASSUMPTION 2 — 0.5 cutoff for the `MethylatedCpGSites` count.** Schultz (2012) recommends a
  binomial test rather than a fixed fraction; 0.5 is a descriptive convenience for the count field and
  does **not** affect the continuous methylation-level outputs, which follow Schultz (2012) exactly.
- **Deviation (from the algorithm doc, not the Evidence file):** the profile method is named
  `GenerateMethylationProfile` (Registry name `CalculateMethylationProfile`) — naming only, same
  behaviour, kept for API stability.

No source contradictions — IUPAC (Cornish-Bowden 1985), Bismark (Krueger & Andrews 2011), Lister
(2009), and Schultz (2012) are mutually consistent.
</content>
</invoke>
