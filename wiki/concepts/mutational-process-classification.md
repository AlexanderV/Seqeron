---
type: concept
title: "Mutational process classification (SBS exposure → active aetiology processes)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-SIG-004-Evidence.md
  - docs/algorithms/Oncology/Mutational_Process_Classification.md
source_commit: 7783b8d65bc2a9ec09d6958e86147e92a7e07908
created: 2026-07-10
updated: 2026-07-14
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-sig-004-evidence
      evidence: "Test Unit ID: ONCO-SIG-004 ... Algorithm: Mutational Process Classification (SBS exposure → active mutational processes)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:mutational-signature-fitting-and-extraction
      source: onco-sig-004-evidence
      evidence: "only the per-signature exposures (activities) and their COSMIC labels are used ... normalized relative contribution = exposureᵢ / Σ exposure — the exposures are the output of the ONCO-SIG-002 NNLS refit that this unit classifies."
      confidence: high
      status: current
---

# Mutational process classification

The Oncology family's **aetiology-annotation** unit (**ONCO-SIG-004**): the interpretation step that turns the
per-signature **exposures** produced by the fit ([[mutational-signature-fitting-and-extraction]], ONCO-SIG-002)
into a human-readable set of **active mutational processes** — APOBEC, Aging (clock-like), Tobacco smoking, UV,
MMR deficiency, etc. It adds **no new decomposition**: it normalizes exposures to relative contributions,
applies a presence cutoff, maps each surviving COSMIC signature label to its **proposed aetiology**, and
aggregates per process. The literature-traced record is [[onco-sig-004-evidence]]; [[test-unit-registry]] tracks
the unit and [[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## The pipeline: normalize → cutoff → map → aggregate

Four steps, each pinned to a source:

1. **Normalize to relative contributions.** Each signature's activity is its fraction of the total
   reconstructed burden: `Wᵢ = exposureᵢ / Σⱼ exposureⱼ` (deconstructSigs "weights W are normalized between 0
   and 1"). This is the same proportion form `x / Σx` reported by ONCO-SIG-002.
2. **Presence cutoff = 6%.** Any signature with `Wᵢ < 0.06` is set to zero (declared **absent**); `Wᵢ ≥ 0.06`
   is retained (declared **present**). Verbatim from the deconstructSigs reference implementation:
   `signature.cutoff = 0.06` default, applied as `weights[weights < signature.cutoff] <- 0`. The comparison is
   **strict `<`**, so a signature at **exactly 0.06 is retained**. The 6% threshold is empirically calibrated —
   a **1.4% false-negative rate** (38 wrong exclusions across 2,646 simulated signatures).
3. **Map label → process.** Each surviving COSMIC SBS label is mapped to its **proposed aetiology** (the
   caller-supplied label→process map). COSMIC assignments used: **SBS1/SBS5 → Aging (clock-like)**,
   **SBS2/SBS13 → APOBEC**, **SBS4 → Tobacco smoking**, **SBS7a–d → UV**, **SBS6/SBS15/SBS20/SBS26 → MMR
   deficiency** (Alexandrov 2020 / COSMIC SBS aetiology strings).
4. **Aggregate per process.** A process's contribution is the **sum** of its surviving member signatures'
   normalized contributions. The **active-process set** is every process with a non-zero total; the **dominant
   process** is the `argmax` per-process contribution.

## Worked oracle (hand-derived)

Raw exposures `SBS2=50, SBS13=30, SBS1=15, SBS4=5` (Σ = 100):

| Signature | Normalized `Wᵢ` | `≥ 0.06`? | Process |
|-----------|-----------------|-----------|---------|
| SBS2  | 0.50 | yes | APOBEC |
| SBS13 | 0.30 | yes | APOBEC |
| SBS1  | 0.15 | yes | Aging |
| SBS4  | 0.05 | **no** (0.05 < 0.06) | Tobacco (excluded) |

Per-process totals: **APOBEC = 0.50 + 0.30 = 0.80**, **Aging = 0.15**, **Tobacco = 0** (SBS4 below cutoff).
Active processes = `{APOBEC, Aging}`; **dominant = APOBEC (0.80)**.

## Corner cases and invariants

- **Sub-cutoff mass is dropped.** Surviving contributions can sum to **< 1** — the missing mass is attributed to
  "unknown". A process is active only from contributions that survive the 6% cutoff.
- **Multiple simultaneous processes.** A tumour can show several active processes at once (APOBEC + Aging in the
  oracle); the classifier reports the **full active set**, not a single label.
- **Unmapped / unknown-aetiology labels.** A signature absent from the label→process map (e.g. a COSMIC
  "Unknown"-aetiology signature) contributes to **no** named process. Note SBS5's COSMIC aetiology is literally
  "Unknown (clock-like signature)" yet is still grouped as a clock-like **Aging** process by the supplied map.
- **All-zero exposures (Σ = 0).** No normalization is possible → no active processes, no dominant process.
- **Custom cutoff.** The 0.06 default is a parameter (deconstructSigs exposes `signature.cutoff`); overriding it
  changes which signatures survive.

Two source-aligned **assumptions** carry the ambiguities COSMIC leaves open: (1) **per-process aggregation by
summation** — COSMIC defines per-signature aetiologies but no aggregation rule; summation is the natural reading
of additive relative contributions; (2) **the cutoff is applied per-signature, then grouped** — following
deconstructSigs exactly (it operates per signature; processes are a downstream grouping), not to the per-process
total.

## Implementation contract (ONCO-SIG-004 API)

The algorithm spec pins the concrete C# surface in
`OncologyAnalyzer.cs`. Two entry points:

- `OncologyAnalyzer.ClassifyMutationalProcess(exposures, contributionCutoff = 0.06)` — takes
  `IReadOnlyList<(string Signature, double Exposure)>`, returns `ActiveProcesses`
  (`IReadOnlyList<ProcessActivity>`, each a process + its summed contribution ∈ [0,1], ordered by
  **descending contribution then process enum**) and `DominantProcess` (a `MutationalProcess` enum,
  `Unknown` when no process is active).
- `OncologyAnalyzer.GetMutationalProcess(signatureLabel)` — the O(1) dictionary lookup that maps a COSMIC
  SBS label to its `MutationalProcess`. **Label matching is case-insensitive**; labels outside the map resolve
  to `Unknown` and contribute to no named process.

**Validation contract:** null `exposures` or a null label → `ArgumentNullException`; a negative or `NaN`
exposure → `ArgumentException`; `contributionCutoff` that is `NaN` or outside `[0, 1)` →
`ArgumentOutOfRangeException`. An empty list or a zero-total list yields an empty active set and `Unknown`
dominant. **Complexity** is `O(k log k)` time / `O(k)` space for `k` signatures — the `log k` is only the final
ordering of the ≤ 5 processes, so it is effectively `O(k)`. This is not a search/matching unit, so the
repository suffix tree does not apply. Confidence-based presence (a bootstrap-CI lower bound above zero rather
than a point cutoff) is **not implemented** here — pair this cutoff rule with
[[signature-exposure-bootstrap-confidence-intervals]] (ONCO-SIG-003) for interval estimates.

## Relation to the oncology family

This is the **interpretation layer directly above** [[mutational-signature-fitting-and-extraction]]
(ONCO-SIG-002): it `depends_on` that unit's exposure vector as its input and annotates it with aetiology, three
steps above the [[sbs96-mutational-signature-catalog]] (ONCO-SIG-001) 96-channel spectrum. Where ONCO-SIG-002
answers "**how much of each signature?**" and ONCO-SIG-003
([[signature-exposure-bootstrap-confidence-intervals]]) answers "**how confident is that fit?**", this unit
answers "**which biological processes are active?**" — the label a clinician reads (an APOBEC-driven vs
tobacco-driven vs MMR-deficient tumour). The named-process output is a somatic-mutation-process biomarker
orthogonal to the copy-number-scar [[homologous-recombination-deficiency-score]] and the mismatch-repair
[[microsatellite-instability-detection]] (note MMR-deficiency signatures SBS6/15/20/26 here are a *signature-
based* echo of that separate MSI biomarker), and it feeds the clinical-interpretation units
([[cancer-variant-tier-classification-amp-asco-cap]], [[clinical-actionability-oncokb-levels]]) with a
process label rather than an anonymous exposure vector.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for mutational-process classification. The models are
deconstructSigs (Rosenthal 2016 — normalized relative contributions + the 6% presence cutoff, verbatim
`signature.cutoff = 0.06` / `weights[weights < signature.cutoff] <- 0`, 1.4% false-negative calibration), COSMIC
SBS proposed-aetiology strings, and Alexandrov 2020 (*Nature* 578:94, the 81-signature catalogue behind the
COSMIC aetiologies). The **per-process summation** and **per-signature cutoff-then-group** are the two source-
aligned assumptions; the caller supplies the reference signature profiles and the label→process map (neither is
fabricated here). **Not for clinical or diagnostic use.** No source contradictions.
</content>
</invoke>
