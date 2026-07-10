---
type: concept
title: "CRISPR guide RNA design, on-target efficacy, and off-target scoring"
tags: [primer, algorithm, validation]
sources:
  - docs/Validation/reports/CRISPR-GUIDE-001.md
  - docs/Validation/reports/CRISPR-OFF-001.md
  - docs/algorithms/MolTools/Guide_RNA_Design.md
source_commit: c763b50cd147a56f659298a8626d2ede865297f1
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
      object: concept:test-unit-registry
      source: crispr-off-001-report
      evidence: "Test Unit ID CRISPR-OFF-001 (MolTools); the two-stage report validates the off-target scoring surface CrisprDesigner.CalculateMitHitScore / CalculateMitSpecificityScore (MIT/Hsu 2013) and CalculateCfdScore (Doench 2016)."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:primer-dimer-thermodynamics-tm
      source: crispr-guide-001-report
      evidence: "Sibling MolTools reagent-design unit in the CrisprDesigner/MolTools family, alongside the primer and probe design units."
      confidence: high
      status: current
---

# CRISPR guide RNA design, on-target efficacy, and off-target scoring

The **CRISPR guide RNA (gRNA/sgRNA) design** surface — a distinct reagent-design algorithm in the
`CrisprDesigner` family, sibling to the primer units [[primer-dimer-thermodynamics-tm]] /
[[primer3-weighted-penalty-objective]] and the probe units [[taqman-probe-design-rules]] /
[[probe-offtarget-specificity-scan]]. It has **three clearly separated layers**: (1) a
**composition-based heuristic** that extracts candidate guides at PAM sites and ranks them
(test unit **CRISPR-GUIDE-001**); (2) **two experimentally-calibrated on-target efficacy models**
(Doench 2014 Rule Set 1 and Doench 2016 Rule Set 2 / Azimuth) that predict editing activity from a
fixed 30-nt context (CRISPR-GUIDE-001); and (3) **off-target risk scoring** — scanning a genome for
near-matches to a guide and scoring each hit's cutting risk, via the **MIT/Hsu-Zhang 2013**
specificity model and the **CFD (Doench 2016)** score (test unit **CRISPR-OFF-001**). The two
independent two-stage validation verdicts are [[crispr-guide-001-report]] (on-target) and
[[crispr-off-001-report]] (off-target); [[test-unit-registry]] tracks both units and
[[algorithm-validation-evidence]] describes the artifact pattern. Implementation:
`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs`.

On-target and off-target scoring are the **two complementary halves of guide selection**: on-target
predicts *how well* a guide cuts its intended site, off-target predicts *how much* it also cuts elsewhere.
Both are learned/experimentally-calibrated models grounded byte-for-byte against the CRISPOR reference
implementation; both are research-grade, not for clinical use.

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

## Layer 3 — off-target risk scoring (CRISPR-OFF-001)

Off-target scoring answers the complementary question: given a designed guide, **scan a genome for
near-matches and score each hit's cutting risk**. Two independent models live on `CrisprDesigner`, both
grounded against the CRISPOR reference and validated in [[crispr-off-001-report]] (Stage A/B PASS, CLEAN;
71/71 off-target+CFD tests). The pre-existing honest-heuristic `FindOffTargets` /
`CalculateSpecificityScore` / `CalculateOffTargetScore` are unchanged; the scored models below were added
on top of them.

**MIT / Hsu-Zhang 2013 specificity** — `CalculateMitHitScore(guide20, offTarget20)` returns a single-hit
score in [0,100]; `CalculateMitSpecificityScore(...)` aggregates per-hit scores into a guide specificity
(with a genome-scanning overload that reuses `FindOffTargets`). Realised verbatim from CRISPOR
`calcHitScore` / `calcMitGuideScore` over a **20-element mismatch-penalty vector `W`** (byte-identical to
the reference `hitScoreM`):

- `score1 = Π over mismatched positions i of (1 − W[i])`;
- `score2 = 1` if `nmm < 2` else `1/(((19 − meanInterMismatchDist)/19)·4 + 1)` (`maxDist = 19`);
- `score3 = 1` if `nmm == 0` else `1/nmm²`;
- `hitScore = score1·score2·score3·100`; aggregate `= 100/(100 + Σ hitScores)·100`.

**Orientation** (pinned from the source): `W` index 0 = PAM-distal 5' (low/zero weight), index 19 =
PAM-proximal seed (max `W[13]=0.851`) — mismatches in the seed are least tolerated. Worked oracles:
perfect→100, mm@0 (W=0)→100, mm@5→60.5, mm@13→14.9, mm@19→41.7, two-mm{5,15}→0.8987, aggregate
{60.5}→62.30529595. A committed orientation-guard test fails if `W` is reversed.

**CFD (Cutting Frequency Determination, Doench 2016)** — `CalculateCfdScore(sgRna20, offTarget20,
offTargetPam)` returns a score in [0,1]. `score = 1`; guide/off `T→U`; per position `i` (0..19, 5'→3'),
matched → ×1, else × `mm_scores["r{guide[i]}:d{complement(off[i])},{i+1}"]`; finally × `pam_scores[last-2
PAM-nt]`. The **240-entry mismatch matrix + 16-entry PAM table** (Doench et al. 2016, *Nat Biotechnol*
34:184) are cross-checked **240/240 + 16/16 identical across two independent sources** (CRISPOR + iGWOS)
and reproduced at full `double` precision. PAM table: `GG=1.0`, `AG=0.259259`, `CG=0.107143`,
`GA=0.069444`, `TG=0.038961`, `GC=0.022222`, `GT=0.016129`, all others `0.0`. Same orientation convention
(index 0 = PAM-distal); key `rX` = guide base T→U, `dY` = complement of the off-target base. Oracles:
perfect+GG→1.0, published iGWOS doctests 0.4635989007074176 / 0.5140384614450001, single mm rC:dT,20→0.5.
This closed finding **C7** (CFD was its last off-target residual; Rule Set 2 / Azimuth cleared separately).

## Validation status

Independently re-validated 2026-06-24 (on-target, CRISPR-GUIDE-001): **Stage A PASS · Stage B PASS · State
CLEAN**, 54/54 CRISPR tests green, **no code or test change**. Off-target (CRISPR-OFF-001) likewise
**Stage A/B PASS · CLEAN**, 71/71 off-target+CFD tests, no production-code change (one in-tree MIT
orientation-guard test verified source-correct and committed); full suite 6812 passed / 0 failed. This session re-downloaded the CRISPOR reference and re-confirmed the
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
- **Off-target scoring is a separate surface from on-target efficacy.** MIT/Hsu (`CalculateMitHitScore`,
  [0,100]) and CFD (`CalculateCfdScore`, [0,1]) score a guide against a **specific 20-nt off-target
  protospacer** (CFD also takes the off-target PAM); both require 20-nt ACGT inputs and throw otherwise.
  Their **orientation is index 0 = PAM-distal, seed = PAM-proximal** — reversing it is the classic bug
  and is caught by dedicated guard tests. A HIGH MIT/CFD hit score means MORE off-target cutting risk;
  the aggregate MIT specificity inverts this (high = specific guide).

No source contradictions: the Doench 2014 model, the MIT/Hsu `W` vector, and the CFD matrices all match
the CRISPOR reference exactly (CFD additionally 240/240 + 16/16 identical across CRISPOR and iGWOS), and
the Azimuth provenance is consistent with the prior multi-session report.
