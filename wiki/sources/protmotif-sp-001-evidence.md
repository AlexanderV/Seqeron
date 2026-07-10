---
type: source
title: "Evidence: PROTMOTIF-SP-001 (Signal-peptide cleavage-site prediction — von Heijne 1986 weight matrix / EMBOSS sigcleave)"
tags: [validation, protein, motif]
doc_path: docs/Evidence/PROTMOTIF-SP-001-Evidence.md
sources:
  - docs/Evidence/PROTMOTIF-SP-001-Evidence.md
source_commit: 1513221ea012107594feec09dd5b4e850e3d5f37
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: supersedes
      object: source:protmotif-domain-001-evidence
      source: protmotif-sp-001-evidence
      evidence: "Change History / spec §7: the ProteinMotifFinder.PredictSignalPeptide record was redesigned to the honest von Heijne (1986) weight-matrix model (CleavagePosition/Score/SignalSequence/WindowSequence/IsLikelySignalPeptide), REPLACING the prior fabricated n/h/c tripartite fields recorded in PROTMOTIF-DOMAIN-001. Scoped to the signal-peptide method; DOMAIN-001 still covers domain prediction."
      confidence: high
      status: current
---

# Evidence: PROTMOTIF-SP-001

The validation-evidence artifact for test unit **PROTMOTIF-SP-001** — **signal-peptide
cleavage-site prediction** on `ProteinMotifFinder.PredictSignalPeptide` (same method as
[[protmotif-domain-001-evidence|PROTMOTIF-DOMAIN-001]], now on a **redesigned** algorithm). One
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the model, matrix, contract, invariants and worked oracle are synthesized in
[[protein-domain-and-signal-peptide-prediction]]. See [[test-unit-registry]] for how units are tracked.

## Supersession — the important finding

This unit documents the **von Heijne (1986) weight-matrix** method (EMBOSS `sigcleave`), which
**replaces the earlier tripartite n/h/c + −1,−3 model** that PROTMOTIF-DOMAIN-001 recorded for the
same `PredictSignalPeptide` method. The prior tripartite scoring (constants 0.95/0.825, fields
`NRegion`/`Probability`, score range [0,1]) was found to be **fabricated** and was removed; the
record now returns `CleavagePosition`, `Score`, `SignalSequence`, `WindowSequence`,
`IsLikelySignalPeptide`. The current code (`ProteinMotifFinder.PredictSignalPeptide`) implements
the weight-matrix method described here. The DOMAIN-001 evidence's signal-peptide description is
therefore **historical** for this method.

## What this file records

- **Method (rank 1/3):** cleavage-site prediction by the **von Heijne (1986) log-odds weight
  matrix**, as implemented in **EMBOSS 6.6.0 `sigcleave`** (doc + `sigcleave.c` source + `Esig.euk`/
  `Esig.pro` data files). "The EMBOSS implementation of the weight matrix method (von Heijne 1986)."
- **Scoring:** at each candidate site the score is `Σ ln(count/expect)` over positions **−13..+2**
  (15 matrix columns; the `Expect` column is not scored). Natural log. Zero counts are replaced by
  `1.0e-10` at the conserved columns **−3 and −1** (a strong penalty) and by `1.0` elsewhere before
  the log (`sigcleave_readSig`).
- **Selection:** the single prediction is the **argmax** of the weight over all positions
  (`maxweight`/`maxsite`); cleavage is between `−1` and `+1`; `+1` is the first residue of the mature
  protein → `CleavagePosition` (1-based mature start).
- **Threshold:** `-minweight` default **3.5** → `IsLikelySignalPeptide ⇔ Score ≥ 3.5`; at 3.5 the
  method identifies ~95% of signal peptides / rejects ~95% of non-signal peptides; cleavage site
  correct in only **75–80%** of cases (heuristic, not exact).
- **Matrices:** eukaryotic 20×16 counts from **161** aligned sequences (default); prokaryotic from
  **36** (via `prokaryote:true`), both from von Heijne (1986), copied verbatim into
  `EukaryoticCounts`/`EukaryoticExpect` / `ProkaryoticCounts`/`ProkaryoticExpect`; each position column
  sums to the sample size (EMBOSS sanity check).
- **Worked oracle (rank 1 reproduction):** **ACH2_DROME (UniProt P17644**, 576 aa) → **maximum score
  13.739**, mature-protein start residue **42**, window `LLVLLLLCETVQA` (positions −13..−1, residues
  29–41). Independently re-derived in Python (natural-log transform + pseudocount rule) → 13.739 exactly.
- **Corner cases:** a best site is always reported for in-window sequences (`-minweight` only flags
  likelihood, no intrinsic cutoff); sequences shorter than one full **15-residue window** return
  `null` (ASSUMPTION-1); non-standard residues (X/B/Z/*) contribute 0 and do not throw; case-insensitive
  (input upper-cased).

## Deviations and assumptions

One assumption: **minimum input length = one full 15-residue window** (shorter inputs → `null`); EMBOSS
scores any length by skipping off-window columns, but a truncated-window score is not meaningful and no
in-scope signal peptide is < 15 aa. Two documented pitfalls the tests guard: using log base 10 instead
of natural log changes every score; forgetting the −3/−1 zero-count `1.0e-10` penalty inflates scores at
sites missing conserved residues. No source contradictions within this unit.

## Recommended coverage

MUST: ACH2_DROME → `CleavagePosition == 42`, `Score == 13.739` (±1e-3); `Score` = argmax of the
Σ ln(count/expect) log-odds over −13..+2 (hand-computed window); a known runner-up scores strictly lower;
`IsLikelySignalPeptide ⇔ Score ≥ 3.5`; null/empty/short (< 15 aa) → `null`. SHOULD: case-insensitivity;
prokaryotic matrix (`prokaryote:true`) scores against `Esig.pro` and may differ. COULD: non-standard
residue `X` contributes 0 without throwing.
