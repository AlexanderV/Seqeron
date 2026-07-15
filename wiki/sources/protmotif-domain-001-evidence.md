---
type: source
title: "Evidence: PROTMOTIF-DOMAIN-001 (Protein domain prediction — PROSITE pattern + Plan7 profile-HMM — & signal-peptide prediction)"
tags: [validation, protein, motif]
doc_path: docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md
sources:
  - docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md
source_commit: d7b350cc67d3684eb247a261ab387137e91e7fba
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROTMOTIF-DOMAIN-001

The validation-evidence artifact for test unit **PROTMOTIF-DOMAIN-001** — **domain prediction &
signal-peptide prediction** on `ProteinMotifFinder`: `FindDomains` (deterministic PROSITE-pattern
domain scan), the opt-in `FindDomainsByHmm` / `Plan7ProfileHmm` profile-HMM engine, and
`PredictSignalPeptide` (von Heijne tripartite + −1,−3 rule). One instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the model, patterns,
engine, contract, invariants and worked oracles are synthesized in
[[protein-domain-and-signal-peptide-prediction]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Domain PROSITE patterns (rank 2, official spec):** PS00028 C2H2 zinc finger
  `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` (PF00096); PS00017 P-loop/Walker A `[AG]-x(4)-G-K-[ST]`
  (PF00069); PS00678 WD40 `[LIVMSTAC]-[LIVMFYWSTAGC]-[LIMSTAG]-[LIVMSTAGC]-x(2)-[DN]-x-{P}-[LIVMWSTAC]-{DP}-[LIVMFSTAG]-W-[DEN]-[LIVMFSTAGCN]`
  (14 elements / 15 residues, PF00400). ScanProsite → regex rules (`x`→`.`, `{X}`→`[^X]`, `x(n)`→`{n}`).
  Real WD40 positive: GBB1_HUMAN (UniProt P62873) matched at 0-based 69/156/284.
- **Profile-only families (honest residual):** SH3 (PS50002) and PDZ (PS50106) are PROSITE **PROFILEs**,
  not patterns — no deterministic regex exists, so they are NOT found by `FindDomains` (unsourced ad-hoc
  regexes removed). Covered instead by the Plan7 engine.
- **Plan7 profile-HMM engine (opt-in, 5 addenda 2026-06-25):** reproduces HMMER3 over 3 bundled **CC0
  Pfam** HMMs — PF00018.35 SH3 (LENG 48), PF00595.30 PDZ (LENG 81), PF00400.39 WD40 (LENG 39). Verified
  against **pyhmmer 0.12.1** ground truth and a from-scratch Python re-derivation:
  - Viterbi/Forward log-odds (Durbin §5.4; insert log-odds ≈ 0) — exact to 1e-9 on a hand-built HMM
    (0.5187937934… nats).
  - hmmsearch-parity local-multihit Forward `pre_score` scored vs Swiss-Prot `bg->f` (not `COMPO`):
    SH3 68.7097 / PDZ 84.8629 / WD40 213.4120 bits, ~1e-5-bit parity.
  - null2 biased-composition correction (`omega=1/256`); E/P-value layer (Gumbel MSV/Viterbi,
    exponential Forward, λ≈log2, `E=P·Z`, from `STATS LOCAL`).
  - Multi-domain envelope decomposition (`p7_domaindef` rt1=0.25/rt2=0.10/rt3=0.20) → GBB1/PF00400 = **7**
    envelopes (exact coords); stochastic-traceback clustering (Easel LCG seed 42) for closely-overlapping
    tandem domains (`FindDomains(seq, clusterOverlapping:true)`).
- **Signal-peptide model (von Heijne, rank 1):** N-terminal 16–30 aa, tripartite n (1–5, K/R+) / h (7–15
  hydrophobic α-helix) / c (3–7 polar) regions; **−1,−3 rule** small neutral **{A,G,S}**. Scoring
  `(nScore+2·hScore+cScore)/4`; constraints `nScore>0` & `hScore≥0.5` (no arbitrary threshold);
  `Probability=Score`.
- **Datasets / oracles:** C2H2 `AAAACXXCXXXLXXXXXXXXHXXXHAAA`→4..24; P-loop `AAAAGXXXXGKSAAAA`→4..11;
  signal peptide `MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF`→detected, cleavage ≈25.
- **Corner cases:** empty/null → no domains; tandem repeats (multiple zinc fingers / 7 WD40 blades) all
  detected; sequences < ~15 aa or charged/cytoplasmic → no signal peptide; case-insensitive.

## Deviations and assumptions

Six items, all documented as **resolved design decisions** (not open assumptions): signal-peptide 1:2:1
region weights (h double, von Heijne 1985); evidence-based detection constraints replacing the old 0.4
threshold; `Probability=Score` (dropped the ×1.2 scaling); strict {A,G,S} −1,−3 set (dropped T/C);
PROSITE-derived regex patterns as a deliberate scope boundary; `FindDomains` naming (pattern search, not
probabilistic prediction). **Honest residuals:** SH3/PDZ have no deterministic PROSITE pattern; the Plan7
engine bundles only 3 CC0 profiles and does not reimplement hmmsearch MSV/bias prefilters or exact-RNG
trace-ensemble bit parity.

## Recommended coverage

MUST: FindDomains detects C2H2 (PS00028) and P-loop (PS00017) at the expected coords with correct
name/accession/start/end; empty/null → empty; PredictSignalPeptide detects the tripartite structure,
enforces the −1,−3 rule, returns the three regions, and returns null for charged/too-short sequences;
case-insensitive. SHOULD: FindDomains detects WD40 (PF00400) and (via the HMM engine) SH3/PDZ;
no-match → empty. COULD: multiple tandem domain instances. No source contradictions.
