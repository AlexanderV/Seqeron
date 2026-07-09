---
type: source
title: "Evidence: ONCO-SIG-004 (Mutational process classification)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-SIG-004-Evidence.md
sources:
  - docs/Evidence/ONCO-SIG-004-Evidence.md
source_commit: 3c5a975fe365264937f55c3a66b25eeff9d0bb0f
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-SIG-004

The validation-evidence artifact for test unit **ONCO-SIG-004** — **mutational process classification**, the
aetiology-annotation layer that turns fitted SBS signature exposures into a named set of active mutational
processes. It is the **thirty-second ingested unit of the Oncology family** and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in
its own concept, [[mutational-process-classification]]; [[test-unit-registry]] tracks the unit.

This unit answers "**which biological processes are active?**" — not a new decomposition but a
normalize → cutoff → map → aggregate pipeline over the per-signature exposures produced by the ONCO-SIG-002 NNLS
refit ([[mutational-signature-fitting-and-extraction]]).

## What this file records

Three sources, one hand-derived dataset:

- **deconstructSigs — Rosenthal et al. 2016** (*Genome Biology* 17:31; rank 1 paper / rank 3 reference source).
  - **Normalized relative contributions:** "the weights W are normalized between 0 and 1" — each signature's
    reported activity is its fraction of the total reconstructed burden.
  - **6% presence cutoff:** "any signature with Wᵢ < 6% is excluded"; the reference `whichSignatures.R` declares
    the default `signature.cutoff = 0.06` and applies `weights[weights < signature.cutoff] <- 0` (both returned
    **verbatim** from the raw source). The comparison is **strict `<`**, so exactly 0.06 is retained.
  - **Calibration:** the 6% threshold "only resulted in 38 instances where a signature was incorrectly excluded
    for a false negative rate of 1.4%" across 2,646 simulated signatures — the source-recommended,
    empirically-calibrated absence threshold.
- **COSMIC Mutational Signatures — SBS** (Wellcome Sanger; rank 5, underpinned by Alexandrov 2020). Proposed
  aetiology strings quoted per signature, giving the label→process map used: **SBS1** "Spontaneous deamination
  of 5-methylcytosine (clock-like)" and **SBS5** "Unknown (clock-like)" → **Aging**; **SBS2/SBS13** "Activity of
  APOBEC family of cytidine deaminases" → **APOBEC**; **SBS4** "Tobacco smoking" → **Tobacco**; **SBS7a–d**
  "Ultraviolet light exposure" → **UV**; **SBS6/SBS15/SBS26** "Defective DNA mismatch repair" and **SBS20**
  "Concurrent POLD1 mutations and defective DNA mismatch repair" → **MMR deficiency**.
- **Alexandrov et al. 2020** (*Nature* 578(7793):94–101; rank 1) — the primary aetiology reference: 4,645 WGS +
  19,184 exomes, **81 SBS/DBS/indel signatures**, the published basis behind the COSMIC aetiology assignments.

## Corner cases / failure modes

- **Sub-cutoff signatures dropped:** any signature below 6% is forced to zero, so surviving contributions can
  sum to **< 1** (missing mass → "unknown"). A process is declared active only from surviving contributions.
- **Multiple active processes:** a tumour can show several processes simultaneously (APOBEC SBS2+SBS13 plus
  Aging SBS1); classification reports the **full active set**, not a single label.
- **Unknown/unmapped aetiology:** a signature label absent from the process map contributes to no named process;
  SBS5's aetiology is literally "Unknown (clock-like)" yet is still grouped as clock-like Aging by the map.
- **All-zero exposures (Σ = 0):** no normalization possible → no active processes, no dominant process.

## Dataset (deterministic worked oracle)

Hand-derived from the 6% cutoff rule + COSMIC label→process map (reference **profiles** are caller-supplied, not
fabricated; only exposures and labels are used). Raw exposures `SBS2=50, SBS13=30, SBS1=15, SBS4=5`, Σ = 100:

| Signature | Normalized | ≥ 0.06? | Process |
|-----------|------------|---------|---------|
| SBS2  | 0.50 | yes | APOBEC |
| SBS13 | 0.30 | yes | APOBEC |
| SBS1  | 0.15 | yes | Aging |
| SBS4  | 0.05 | **no** | Tobacco (excluded) |

Per-process: APOBEC = 0.80, Aging = 0.15, Tobacco = 0. Active set = {APOBEC, Aging}; **dominant = APOBEC (0.80)**.

## Deviations and assumptions

- **ASSUMPTION — per-process aggregation by summation.** When several signatures map to one process (SBS2+SBS13
  → APOBEC; SBS6/15/20/26 → MMR), the process contribution is the **sum** of member surviving contributions.
  COSMIC defines per-signature aetiologies but no aggregation rule; summation is the natural reading of additive
  relative contributions (deconstructSigs weights are additive fractions of one reconstruction). Affects only the
  per-process total, not the per-signature cutoff decision.
- **ASSUMPTION — cutoff applied per-signature, then grouped.** Following deconstructSigs exactly, the 6% cutoff
  is applied to each individual signature's normalized contribution; surviving contributions are aggregated by
  process afterwards (deconstructSigs operates per signature; processes are a downstream grouping).

No source contradictions — deconstructSigs (contribution normalization + cutoff), COSMIC (per-signature
aetiologies), and Alexandrov 2020 (the catalogue) are mutually consistent; the two assumptions fill gaps COSMIC
leaves open rather than departing from any source.
</content>
