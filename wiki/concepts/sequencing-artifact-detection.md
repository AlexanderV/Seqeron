---
type: concept
title: "Sequencing artifact detection (OxoG / FFPE deamination + strand bias)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-ARTIFACT-001-Evidence.md
source_commit: d4ef2c36c5c292c694f25a2fba12074d63939467
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-artifact-001-evidence
      evidence: "Test Unit ID: ONCO-ARTIFACT-001 ... Algorithm: Sequencing Artifact Detection (OxoG / FFPE deamination substitution classification + strand-orientation bias)"
      confidence: high
      status: current
---

# Sequencing artifact detection (OxoG / FFPE deamination + strand bias)

Distinguishing **technical sequencing artifacts** from true somatic variants (`FilterArtifacts`) — the
**third ingested unit of the Oncology family** and the upstream quality-control step that precedes
clinical interpretation. Unlike the two clinical-significance siblings — the therapeutic ranking
[[clinical-actionability-oncokb-levels]] and the four-tier [[cancer-variant-tier-classification-amp-asco-cap]]
— this unit does not judge biological significance; it removes **false-positive calls caused by DNA
damage and mapping bias** before any variant reaches those classifiers. `FilterArtifacts` removes flagged
artifacts and keeps real variants, so the result is a **subset of the input**. Validated under test unit
**ONCO-ARTIFACT-001**; the record is [[onco-artifact-001-evidence]], [[test-unit-registry]] tracks the
unit, and [[algorithm-validation-evidence]] describes the artifact pattern. Its **biological-origin
counterpart** is [[clonal-hematopoiesis-cfdna-filtering]] — same pre-interpretation QC stage, but
removing clonal-hematopoiesis (CHIP) blood-clone false positives rather than technical DNA-damage
artifacts.

## Three independent artifact signals

The filter composes three rule layers, each sourced from a canonical reference.

### 1. Substitution-class classification

Two damage chemistries produce **disjoint** substitution classes; a third catch-all covers everything
else. The class alone (independent of counts) identifies the *candidate* artifact type:

| ref→alt | Class | Chemistry |
|---------|-------|-----------|
| **G>T** / **C>A** | **OxoG oxidation** | 8-oxo-guanine (8-oxo-dG) mispairs, G:C→T:A (Chen 2017) |
| **C>T** / **G>A** | **FFPE deamination** | cytosine → uracil → pairs with A, C:G→T:A (Do & Dobrovic 2015) |
| e.g. **A>G** | **neither** (not an artifact class) | not a known damage signature |

The two classes are disjoint by substitution type: **C>A/G>T is oxidation**, **C>T/G>A is deamination**.
G>T on read 1 appears as its reverse complement C>A on read 2 (same G:C>T:A event, read-orientation
imbalance) — the read-orientation asymmetry is what the GIV score quantifies.

### 2. GIV (Global Imbalance Value) — read-pair orientation imbalance

**Oxidative (OxoG) damage** during library prep (e.g. acoustic shearing) creates a **G>T excess on read 1**
and the reverse-complement **C>A excess on read 2**. The **GIV score** (Chen et al. 2017; `Damage-estimator`,
Ettwiller lab) is computed **per substitution type** as the ratio of R1 to R2 variant counts:

```
GIV_G_T = (count of G>T variants in R1) / (count of G>T variants in R2)
```

- **GIV = 1** ⇒ no DNA damage (balanced R1/R2) — **undamaged**.
- **GIV > 1.5** ⇒ **damaged DNA** (Nature Methods summary of Chen 2017, verbatim operational threshold).
- Standard acoustic-shearing protocol gives GIV_G_T ≈ **2** (Damage-estimator README); GIV_G_T = 2 means
  the G>T error rate in the 8-oxoG mode is twice the non-8-oxoG rate; GIV > 5 = significant damage.

The ratio is neutral at 1 (not 0). A separate GIV is defined for each of the 12 possible mutation types.

### 3. FisherStrand (FS) — strand-orientation bias

**GATK FisherStrand** (Broad Institute) flags variants whose alt reads segregate onto one strand — a
mapping/damage artifact signature. It builds a **2×2 contingency table** and reports the **Phred-scaled
p-value of a two-sided Fisher's exact test**:

```
table = [ ref_fwd, ref_rev, alt_fwd, alt_rev ]   (StrandBiasTest cell order)
FS = -10 · log10( max(p, MIN_PVALUE) )            MIN_PVALUE = 1e-320
```

- **Null hypothesis:** ref/alt are evenly distributed across strands (no bias). No bias ⇒ p ≈ 1 ⇒ **FS ≈ 0**.
- **Higher FS ⇒ stronger strand bias ⇒ more likely artifact.** Perfect segregation (all alt on one strand)
  ⇒ minimal p ⇒ large FS. The `MIN_PVALUE = 1e-320` floor caps FS to avoid infinity.
- A table with a **zero margin** (e.g. no reverse reads at all) yields p = 1 (no evidence of bias) under
  the two-sided test.

## Worked oracles

- **GIV:** R1 G>T = 200, R2 G>T = 100 → GIV = **2.0** (damaged, > 1.5); balanced R1 = R2 = 100 → GIV = **1.0**
  (undamaged).
- **FisherStrand:** balanced table `[10,10,10,10]` → two-sided p = 1.0 → FS = **0.000**; segregated
  `[20,0,0,20]` → exact hypergeometric two-sided p → FS **large (> 0)**.
- **Substitution class:** C>T, G>A → FFPE deamination; G>T, C>A → OxoG oxidation; A>G → **neither**.

## Invariants and edge cases

- **INV:** `FilterArtifacts` output ⊆ input (removes flagged artifacts, keeps real variants).
- **GIV zero denominator:** 0 G>T in R2 must be handled — both zero ⇒ no imbalance evidence; only R2 zero
  ⇒ maximal imbalance (no division error).
- **FS empty / single-orientation table:** a zero-margin table ⇒ p = 1 ⇒ FS ≈ 0 (no bias evidence).
- **Non-artifact substitution:** a class that is neither C>T/G>A nor G>T/C>A (e.g. A>G) is not flagged by
  the substitution-class rule.
- **Strand-bias monotonicity:** more segregation ⇒ higher FS (invariant check).
- Null / empty inputs → throw / return empty (API contract per sibling methods).

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the **rule-based artifact
classification** only. **No BAM parser:** the checklist signature is `FilterArtifacts(variants, bamFile)`,
but the repository has no BAM reader; the read-orientation and per-strand alt/ref evidence a BAM would
supply is passed **directly on the variant observation record** instead of parsed from a file. This is an
API-shape decision — the classification rules (substitution class, GIV ratio, Fisher strand p) are
unchanged. The GIV thresholds (neutral 1, damaged > 1.5) are documented operational cutoffs taken verbatim
from the Nature Methods summary of Chen 2017, not invented. **Not for clinical or diagnostic use.** No
source contradictions — Chen 2017 / Damage-estimator (oxidation), Do & Dobrovic 2015 (FFPE deamination),
and GATK FisherStrand (strand bias) each cover a disjoint signal and are mutually consistent.
</content>
</invoke>
