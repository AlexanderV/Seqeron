# Validation Report: PROTMOTIF-HMM-001 — Plan7 Profile-HMM Domain Search

- **Validated:** 2026-06-25   **Area:** ProteinMotif
- **Canonical method(s):** `Plan7ProfileHmm.{Viterbi,Forward,LocalForward}Score`, `HmmSearchBitScore`,
  `LocalForwardBitScore`, `Null2BiasBits`, `FindDomains`, the Gumbel/exponential survival + P/E-value
  helpers; `ProteinMotifFinder.{FindDomainsByHmm, ScoreDomainHmm, FindDomainHitsByHmm,
  ScoreDomainHmmEValue, FindDomainEnvelopes}`.
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Stage A — Description

### Sources opened (this session)
- **Durbin, Eddy, Krogh & Mitchison (1998), *Biological Sequence Analysis* §5.4** — profile-HMM Plan7
  architecture (M/I/D per node, mute B/E), the Viterbi and Forward log-odds recurrences. The code's
  glocal recurrence (`RunGlocal`) is the §5.4 recurrence with emissions scored as log-odds against the
  null `q` (COMPO); Forward = same with max → log-sum-exp. Confirmed.
- **Eddy (2011), *PLoS Comput Biol* 7:e1002195 "Accelerated Profile HMM Searches"** — the HMMER3 local
  N/B/M/I/D/E/J/C special-state architecture, occupancy-weighted local entry `t(B→Mk)=occ[k]/Z`, local
  exit `esc=0`, length model `pmove=(2+nj)/(L+2+nj)`, multihit E→{J,C} split `−ln2`, and the
  biased-composition (null2) correction. The code's `RunLocalForward`/`RunForwardBackward`/`LocalEntryLn`
  realise exactly this. Confirmed.
- **HMMER User's Guide v3.4 (Eddy, Aug 2023), HMM file format** — parameters stored as negative
  natural-log probabilities, `*` = zero; alphabet `ACDEFGHIKLMNPQRSTVWY` (K=20); COMPO line;
  `STATS LOCAL {MSV|VITERBI|FORWARD} <loc> <slope>` with "all three or none", λ>0, scores in bits.
  The parser (`Parse`) matches verbatim. Confirmed.
- **EddyRivasLab/hmmer source semantics** — `p7_bg.c` (p1=L/(L+1) null1; omega=1/256), `generic_null2.c`
  `p7_GNull2_ByExpectation`, `p7_domaindef.c` (rt1=0.25, rt2=0.10, rt3=0.20; nsamples=200; min_overlap=0.8;
  max_diagdiff=4; min_posterior=0.25; min_endpointp=0.02; pipeline seed 42), `p7_spensemble.c` clustering,
  `esl_random.c` LCG (`esl_mix3`, x←69069x+1), `hmmer.c` `p7_AminoFrequencies` (Swiss-Prot background).
  Each constant in the code is annotated with its source and matches.
- **Eddy (2008) *PLoS Comput Biol* 4:e1000069** — Viterbi/MSV bit scores Gumbel-distributed; Forward
  high-scoring tail exponential; λ≈log2. Survival functions (`esl_gumbel_surv`, `esl_exp_surv`) ported
  verbatim, including the `eslSMALLX1=5e-9` tail branch. Confirmed.
- **Pfam PF00018 (SH3_1), PF00595 (PDZ), PF00400 (WD40)** — CC0; bundled HMMER3/f files retrieved
  2026-06-25 from EMBL-EBI InterPro; headers (NAME/ACC/LENG/GA/STATS) parsed correctly.

### Formula / convention checks
- Glocal Viterbi/Forward log-odds, emissions against COMPO null — Durbin §5.4. ✓
- Local-multihit Forward (hmmsearch parity) scores match emissions against HMMER's standard Swiss-Prot
  background `bg->f` (not COMPO), per `p7_bg.c`; the code uses a separate `EmissionLogOddsHmmer`. ✓
- null1 `nullsc = L·ln(p1)+ln(1−p1)`, p1=L/(L+1). ✓  null2 seqbias `logsumexp(0, ln(omega)+Σ ln null2)`. ✓
- Per-domain score `(envsc + (n−Ld)·ln(n/(n+3)) − (nullsc+dombias))/ln2`; i-Evalue `Z·exp(−λ(S−τ))`
  clamped below τ — `p7_pipeline.c`. ✓
- Coordinates: envelopes 1-based inclusive (HMMER env from/to); `FindDomainsByHmm` reports 0-based
  whole-sequence span (documented; the glocal score is a whole-sequence quantity). ✓
- Edge semantics: empty sequence → no domains / −∞ bit score / 0 bias; null → throw; unknown accession
  → throw; out-of-alphabet residue → background odds 1 (HMMER degenerate-residue convention). All defined.

### Independent cross-check (numbers, this session)
**Hand-derivation (Python, full double precision):**
- H1 glocal Viterbi of "AC" through B→M1(A)→M2(C)→E = `ln(0.7/0.6)+ln(0.9)+ln(0.8/0.4)+ln(0.8)`
  = **0.5187937934151676 nats** — matches the test constant exactly.
- H18 local-toy Forward of "A" (1-node, B→M1=1) = **1.272400756045032 nats**, pre-bits
  **3.835686260769536** — matches exactly.
- Viterbi Gumbel P@40 (μ=−8.2932,λ=0.71923) = **8.227179545686635e-16**; Forward exp P@40
  (τ=−4.5735) = **1.1943390031599535e-14**; MSV Gumbel P@40 (μ=−8.1284) = **9.262484858441732e-16**
  — all match.

**pyhmmer 0.12.1 ground-truth hmmsearch (Z=1, domZ=1, seed=42)** vs C# engine:

| Profile / target | hmmsearch (pyhmmer) | C# engine | Δ |
|---|---|---|---|
| SH3/PF00018 vs SRC core | pre 68.709740, dom env 3-50 sc 68.540695, iE 1.452861e-23, bias 0.025574 | matches within 1e-4 bit | ✓ |
| PDZ/PF00595 vs PSD-95 PDZ1 | pre 84.862930, dom env 5-89 sc 84.650108 | matches within 1e-4 bit | ✓ |
| WD40/PF00400 vs GBB1 (L=340) | pre 213.411926, **7** envelopes 45-83/87-125/133-170/174-212/216-254/259-298/303-340, scores 31.139/19.004/25.054/35.552/40.454/23.443/27.824, iE 1.21e-11…1.36e-10 | 7 envelopes, coords exact, scores ≤1e-2 bit, iE ≤5% | ✓ |
| Overlapping tandem SH3 (rt3 ensemble, seed 42) trim 4/12/16 | trim12: env 1-37 sc 48.0472 + 37-84 sc 66.6785; trim4: 1-46/45-92; trim16: 1-33/33-80 | coords exact, scores ≤0.1 bit | ✓ |
| Well-separated tandem SH3 | env 1-48 + 49-96, sc 66.1928 each | exact | ✓ |

Every numeric constant hard-coded in the test fixture was reproduced from pyhmmer/hand-derivation in
this session — they are genuine ground truth, not code echoes.

### Stage A findings
None. The description (XML docs, TestSpec contract) is biologically and mathematically correct and
matches the cited primary sources exactly. The two documented residuals (MSV-stage prefilter; the
glocal `FindDomainsByHmm` reporting a whole-sequence span) are honestly stated, not defects.

## Stage B — Implementation

### Code path reviewed
- `Plan7ProfileHmm.cs` — `Parse` (header/COMPO/node parsing, neg-log→ln), `RunGlocal` (Viterbi/Forward),
  `RunLocalForward`/`RunForwardBackward` (local N..C Forward+Backward), `LocalEntryLn` (occupancy),
  `ComputeNull2`/`Null2ByExpectation`, `FindDomains` (region id rt1/rt2 + rt3 ensemble), `RescoreEnvelope`,
  `RegionTraceEnsemble`/`StochasticTrace`/`ClusterEnsemble`/`SingleLinkage`/`EslLcgRandom`, the
  Gumbel/exp survival + P/E-value helpers.
- `ProteinMotifFinder.cs` lines 1414-1628 — `FindDomainsByHmm`, `ScoreDomainHmm`, `FindDomainHitsByHmm`,
  `ScoreDomainHmmEValue`, `FindDomainEnvelopes` (×2 overloads), bundled-profile registry.

### Formula realised correctly?
Yes. The DP recurrences, length model, local entry/exit, null1/null2, domain decomposition, and the
LCG/ensemble clustering each match the cited HMMER/Easel source. Verified numerically by the
pyhmmer parity table above (the fixture's hard-asserted exact values all reproduce).

### Variant / delegate consistency
- `ScoreDomainHmm`/`FindDomainsByHmm`/`FindDomainHitsByHmm` all route through `ViterbiScore`/`ScoreBits`. ✓
- `HmmSearchBitScore = LocalForwardBitScore − Null2BiasBits` (pipeline identity) — test H18 pins it. ✓
- `ForwardPValue` delegates to `ExponentialSurvival`; `Viterbi/MsvPValue` to `GumbelSurvival`;
  `*EValue = EValue(P,Z) = P·Z`. ✓
- The two `FindDomainEnvelopes` overloads (all-families / by-accession) agree. ✓

### Test quality audit
The fixture (`ProteinMotifFinder_FindDomainsByHmm_Tests`) asserts exact pyhmmer/hand-derived values
with source-justified tolerances (1e-9 hand DP; 1e-4 bit float-vs-double parity; ≤1e-2 bit + 5% iE for
the multi-domain decomposition). No green-washing: tolerances are physically motivated (HMMER single
precision), no skips, the opt-out test (H20b) and cross-family test (H7) guard against degenerate
"always-N" implementations.

**Gap found & fixed (Stage-B test defect):** the original fixture left three public surfaces / Stage-A
edge cases uncovered. Added region **H21** (10 tests):
- `MsvPValue` — was never invoked (only the parsed `MsvMu` constant was checked). Added hand-derived
  Gumbel P@40 = 9.262484858441732e-16, a distinct-from-Viterbi guard, and the uncalibrated-throw guard.
- `ExponentialSurvival` (static) — never called directly. Added the exact hand-derived value and the
  below-τ clamp.
- `EValue` (static) positive path — only the negative-Z throw was tested. Added `EValue(0.5,10)=5.0`.
- All-X / out-of-alphabet residues (Stage-A edge case, `EmissionLogOdds` `residueIndex<0` branch):
  added Viterbi (finite, ≤0, deterministic), `ScoreDomainHmm`/`FindDomainsByHmm` (below threshold,
  no domain), and `HmmSearchBitScore` (finite, ≈0.48 bits, below threshold — hmmsearch itself reports
  no hit for all-X). *Note:* an early draft asserted the per-seq bit score ≤0; pyhmmer/measurement showed
  the local-multihit Forward retains a small positive structural mass (~0.48 bits), so the assertion was
  corrected to the genuine sourced property (near-zero, below the detection threshold) — not a fudge.

### Stage B findings
No code defect. One test-coverage defect (above) fully fixed in this session.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.**
- Full unfiltered `dotnet test Seqeron.sln -c Debug`: **Failed: 0** (Genomics 18737 passed; all other
  projects green). 0 warnings / 0 errors on the changed test project.
- No open defects. Documented (non-defect) residuals: HMMER's MSV/bias prefilter stages are not run
  (the engine computes the full Forward directly, which is a superset, not a divergence), and
  `FindDomainsByHmm` reports a whole-sequence span (use `FindDomainEnvelopes` for HMMER env from/to).
