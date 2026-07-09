---
type: source
title: "Evidence: ONCO-ARTIFACT-001 (sequencing artifact detection — OxoG / FFPE deamination + strand bias)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-ARTIFACT-001-Evidence.md
sources:
  - docs/Evidence/ONCO-ARTIFACT-001-Evidence.md
source_commit: d4ef2c36c5c292c694f25a2fba12074d63939467
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ONCO-ARTIFACT-001

The validation-evidence artifact for test unit **ONCO-ARTIFACT-001** — **Sequencing Artifact Detection**
(`FilterArtifacts`): OxoG / FFPE deamination substitution classification plus strand-orientation bias.
This is the **third ingested unit of the Oncology family** and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct rule-based
artifact-classification method is synthesized in its own concept,
[[sequencing-artifact-detection]]; the two clinical-significance siblings are
[[clinical-actionability-oncokb-levels]] and [[cancer-variant-tier-classification-amp-asco-cap]];
[[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (four, each covering a disjoint signal; mutually consistent):**
  - **Chen L. et al. (2017)** "DNA damage is a pervasive cause of sequencing errors" — *Science*
    355(6326):752–756 (rank 1, primary; abstract + Nature Methods write-up + PMC7350422 comment).
    **OxoG signature** = **G>T on read 1 / C>A on read 2** (a single G:C>T:A class with read-orientation
    imbalance). **GIV (Global Imbalance Value)** = per-substitution-type ratio of R1 to R2 variant counts,
    `GIV_G_T = count(G>T in R1) / count(G>T in R2)`; **GIV = 1 = undamaged, GIV > 1.5 = damaged**
    (Nature Methods summary), GIV = 2 = twice the non-8-oxoG error rate, GIV > 5 = significant damage.
  - **Ettwiller `Damage-estimator`** (rank 3, the paper authors' reference implementation) — damage
    estimation is based on the systematic R1-vs-R2 mutation-rate difference; `estimate_damage.pl` outputs
    a GIV score per substitution type; standard acoustic shearing gives **GIV_G_T ≈ 2** (confirms the
    neutral baseline is 1).
  - **Do & Dobrovic (2015)** "Deamination Effects in FFPE Tissue Samples" — *Clin Chem* / ScienceDirect
    (rank 1, review; cross-checked with the Oxford NAR review). **FFPE signature** = **C>T / G>A**
    (collectively C:G>T:A), caused by cytosine → uracil deamination (uracil pairs with A). Explicitly
    disjoint from oxidation: FFPE = C>T/G>A (deamination) vs C>A/G>T (oxidation).
  - **GATK FisherStrand / StrandBiasTest (Broad)** (rank 3, canonical variant-calling toolkit, source
    read). **2×2 table** `[ref_fwd, ref_rev, alt_fwd, alt_rev]` (StrandBiasTest ARRAY_DIM=2,
    ARRAY_SIZE=4); **FS = -10·log10(two-sided Fisher exact p)**, floor `MIN_PVALUE = 1e-320`. No bias
    ⇒ p ≈ 1 ⇒ FS ≈ 0; higher FS ⇒ stronger strand bias ⇒ more likely artifact.

- **Documented corner cases / failure modes:**
  - **FisherStrand:** no strand bias (even distribution) → p ≈ 1, FS ≈ 0; perfect segregation → minimal
    p, large FS; empty / single-orientation (zero-margin) table → p = 1 (no bias evidence).
  - **GIV:** zero R2 G>T denominator handled (both zero → no imbalance; only R2 zero → maximal
    imbalance); balanced R1 = R2 → GIV = 1 (undamaged).
  - **Substitution class:** a class that is neither C>T/G>A nor G>T/C>A (e.g. A>G) is neither a
    deamination nor an oxidation artifact.

- **Datasets (documented oracles):**
  - **GIV:** R1 G>T = 200, R2 G>T = 100 → GIV = 2.0 (damaged, > 1.5); R1 = R2 = 100 → GIV = 1.0
    (undamaged).
  - **FisherStrand:** balanced `[10,10,10,10]` → p = 1.0 → FS = 0.000; segregated `[20,0,0,20]` →
    exact hypergeometric two-sided p → FS large (> 0).
  - **Substitution class:** C>T, G>A → FFPE deamination; G>T, C>A → OxoG oxidation; A>G → neither.

## Deviations and assumptions

- **ASSUMPTION 1 — no BAM parser.** The checklist signature `FilterArtifacts(variants, bamFile)` implies a
  BAM reader, but the repository has none; the read-orientation / per-strand alt/ref evidence a BAM would
  supply is passed **directly on the variant observation record** instead. This is an API-shape decision,
  not correctness-affecting — the substitution-class, GIV-ratio, and Fisher-strand rules are unchanged.
- **ASSUMPTION 2 — GIV neutral / decision thresholds.** GIV = 1 (undamaged) and GIV > 1.5 (damaged) are
  taken verbatim from the Nature Methods summary of Chen 2017; the underlying R1 G>T / R2 G>T ratio is
  from the paper and Damage-estimator. The 1.5 cutoff is a documented operational threshold, not invented.
- **Coverage recommendations:** MUST-test each substitution-class mapping (C>T/G>A → FFPE, G>T/C>A → OxoG,
  A>G → not-an-artifact); GIV = R1/R2 (200/100 → 2.0 damaged, R1 = R2 → 1.0 undamaged);
  FS = -10·log10(two-sided Fisher p) (balanced → 0.0, segregated → > 0); `FilterArtifacts` removes flagged
  artifacts, keeps real variants (result ⊆ input). SHOULD-test GIV zero-R2 handling (no division error),
  null/empty → throw / empty. COULD-test strand-bias FS monotonicity (more segregation ⇒ higher FS).

No source contradictions — Chen 2017 / Damage-estimator (oxidation), Do & Dobrovic 2015 (FFPE deamination),
and GATK FisherStrand (strand bias) each describe a disjoint artifact signal and are mutually consistent.
</content>
