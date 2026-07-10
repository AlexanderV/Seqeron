---
type: concept
title: "CRISPR guide RNA design and on-target efficacy scoring"
tags: [primer, algorithm, validation]
sources:
  - docs/Validation/reports/CRISPR-GUIDE-001.md
  - docs/algorithms/MolTools/Guide_RNA_Design.md
source_commit: 8ab783ae77cc9e5a6a05c3daefb18a3a4ad4a52d
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: crispr-guide-001-report
      evidence: "Test Unit ID CRISPR-GUIDE-001 (MolTools), doc docs/algorithms/MolTools/Guide_RNA_Design.md; the two-stage report validates CrisprDesigner.CalculateOnTargetDoench2014 / CalculateOnTargetRuleSet2."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:primer-dimer-thermodynamics-tm
      source: crispr-guide-001-report
      evidence: "Sibling MolTools reagent-design unit in the CrisprDesigner/MolTools family, alongside the primer and probe design units."
      confidence: high
      status: current
---

# CRISPR guide RNA design and on-target efficacy scoring

The **CRISPR guide RNA (gRNA/sgRNA) design** surface (test unit **CRISPR-GUIDE-001**, MolTools) — a
distinct reagent-design algorithm in the `CrisprDesigner` family, sibling to the primer units
[[primer-dimer-thermodynamics-tm]] / [[primer3-weighted-penalty-objective]] and the probe units
[[taqman-probe-design-rules]] / [[probe-offtarget-specificity-scan]]. It has **two clearly separated
layers**: (1) a **composition-based heuristic** that extracts candidate guides at PAM sites and ranks
them, and (2) **two experimentally-calibrated on-target efficacy models** (Doench 2014 Rule Set 1 and
Doench 2016 Rule Set 2 / Azimuth) that predict editing activity from a fixed 30-nt context. The
independent two-stage validation verdict is [[crispr-guide-001-report]]; [[test-unit-registry]] tracks
the unit and [[algorithm-validation-evidence]] describes the artifact pattern. Implementation:
`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs`.

## Layer 1 — heuristic guide extraction and ranking

`DesignGuideRnas(sequence, regionStart, regionEnd, systemType, parameters)` resolves the CRISPR system
into its PAM pattern / guide length / orientation, finds PAM sites in the region, extracts each
PAM-adjacent candidate spacer, and scores it with `EvaluateGuideRna`. The score is **heuristic, not
learned** — a base of `100` with deductions:

| Rule | Deduction |
|------|-----------|
| GC below `MinGcContent` (default 40) | `(MinGcContent − actual) × 2` |
| GC above `MaxGcContent` (default 70) | `(actual − MaxGcContent) × 2` |
| Contains `TTTT` (poly-T) | `20` |
| Self-complementarity > 0.3 | `selfComplementarity × 30` |
| Seed-region GC outside 30–80 % | `5` |
| Common restriction site present | `5` |

Score is clamped `≥ 0`; the **seed region** is the 10 nt at the PAM-proximal end (last 10 for
PAM-after-target systems, first 10 for PAM-before). `DesignGuideRnas` yields only candidates with
`Score ≥ MinScore` (default 50); `FullGuideRna` appends a fixed scaffold to the spacer. This layer
**ranks/filters** candidates; it does **not** reproduce a learned activity model — for that, call the
Doench methods below (the algorithm doc explicitly directs experimentally-calibrated efficacy to them).

**Documented simplifications (by design, not defects):** the boolean parameters `AvoidPolyT` and
`CheckSelfComplementarity` are stored but **not honored** by the scoring path (the penalties always
apply); the seed is fixed at 10 nt (upper bound of the cited 8–10 range); the internal helper only
remaps `SpCas9`/`SaCas9`/`Cas12a` metadata, so other named systems (`SpCas9NAG`, `As/LbCas12a`, `CasX`)
reuse `SpCas9` metadata in the returned candidate even though extraction still follows the selected PAM
geometry; `Position` is copied from `pamSite.TargetStart` and not remapped to a forward-source
coordinate for reverse-strand designs.

## Layer 2 — on-target efficacy: Doench 2014 Rule Set 1

`CalculateOnTargetDoench2014(string context30Mer)` is a **logistic-regression linear model** over a
fixed **30-nt context** = `[4 nt 5′] + [20 nt protospacer] + [3 nt PAM] + [3 nt 3′]`. The model is
grounded byte-for-byte against the **CRISPOR reference `doenchScore.py`** (Haeussler 2016), which
implements Doench et al. 2014 (*Nat Biotechnol* 32:1262, PMID 25184501):

- `score = intercept` (**0.59763615**);
- **GC term** over `seq[4..24)`: `score += abs(10 − gcCount) · gcWeight`, where `gcWeight = gcLow
  (−0.2026259)` if `gcCount ≤ 10` else `gcHigh (−0.1665878)` — boundary is `≤ 10 → gcLow`;
- **feature loop:** for each `(position, subsequence, weight)` in a **70-entry** coefficient table, add
  `weight` when `seq[pos:pos+len] == subsequence` (`CompareOrdinal`);
- **output:** `1/(1+e^−score)` ∈ (0,1), the repo scales ×100.

The 70-entry table matches the reference exactly, including its **intentional quirks** — `(24,'AG'/'CG'/
'TG')` reuse the `(24,'A'/'C'/'T')` single-base weights, and `(26,'GT') = (27,'T') = 0.11787758`.

**Input contract (stricter than the guard-free reference):** wrong length / null / empty / non-ACGT
throw; lowercase is upper-cased first (identical value); an **NGG PAM guard** (context offsets 25–26 ==
`GG`) enforces SpCas9 specificity. The guard is an input check, **not** a scoring-model change.

**Worked oracles** (reference formula reproduced from scratch, ≤ 2e-8 agreement — deltas are float-print
precision in the reference's quoted literals):

- `TATAGCTGCGATCTGAGGTAGGGAGGGACC` → 0.7130893 (ref 0.713089368)
- `TCCGCACCTGTCACGGTCGGGGCTTGGCGC` → 0.0189838 (ref 0.0189838464)
- `AAAAAAAAAAAAAAAAAAAAAAAAAGGAAA` ×100 → 4.4338168085 (test oracle M-003)

## Layer 2 — on-target efficacy: Doench 2016 Rule Set 2 / Azimuth

`CalculateOnTargetRuleSet2(...)` is a trained scikit-learn **gradient-boosted regression tree (GBRT)**,
**not a coefficient table** — it cannot be reproduced from published numbers, only from the shipped
Microsoft Research **Azimuth** (BSD-3-Clause) pickles. The repository reproduces it **sklearn-free**:
`AzimuthRuleSet2.cs` (model reader + featurizer + tree traversal), embedded pickles
`Resources/azimuth_rs2_{nopos,full}.bin`, extractor `scripts/azimuth/extract_azimuth_model.py`. Test
oracles are **externally derived** — the CSVs `scripts/azimuth/oracle/{nopos,full}_oracle.csv` (947
rows) carry both the verified reference `ref_score` **and** the upstream Azimuth `upstream` prediction
with an `agrees` flag (documented ~38 % upstream-fixture drift), so tests assert against the reference,
never against the C# output. The `nopos`/`full` wrappers delegate to one engine.

## Validation status

Independently re-validated 2026-06-24: **Stage A PASS · Stage B PASS · State CLEAN**, 54/54 CRISPR tests
green, **no code or test change**. This session re-downloaded the CRISPOR reference and re-confirmed the
intercept, gcLow/gcHigh, the 4+20+3+3 layout, the GC-term boundary, and the full 70-entry table as
byte-identical, reproducing three worked-example scores to ≤ 2e-8; Rule Set 2 provenance and oracle
externality re-confirmed. Full detail and the test-quality (green-washing) audit are in
[[crispr-guide-001-report]]. Research-grade, not for clinical use.

## Assumptions and contract

- **Two distinct scoring surfaces.** The Layer-1 `EvaluateGuideRna` heuristic (composition penalties) is
  **not** an activity predictor; learned on-target efficacy is the separate Doench Rule Set 1 / Rule Set
  2 path over the 30-nt context. Callers must not read the heuristic score as an efficacy estimate.
- **Rule Set 1 requires a valid 30-nt SpCas9 context** (4+20+3+3, NGG PAM at offsets 25–26); malformed
  input throws rather than silently scoring.
- **Rule Set 2 is only reproducible from the shipped pickles**, not from published constants; its tests
  are anchored to externally-derived oracles.

No source contradictions: the Doench 2014 model matches the CRISPOR reference exactly, and the Azimuth
provenance is consistent with the prior multi-session report.
