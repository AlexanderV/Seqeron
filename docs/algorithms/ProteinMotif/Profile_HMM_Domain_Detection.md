# Profile-HMM Domain Detection (Plan7)

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinMotif |
| Test Unit ID | PROTMOTIF-DOMAIN-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-25 |

## 1. Overview

A faithful Plan7 profile hidden Markov model (the HMMER3 architecture) scorer that detects
protein domains lacking any deterministic PROSITE pattern — SH3 (Pfam PF00018), PDZ (PF00595)
and WD40 (PF00400) — by scoring a query against bundled, CC0-licensed Pfam profile HMMs. It is
an **opt-in** path: the existing exact PROSITE-pattern `FindDomains` and its defaults are
unchanged. Scoring is the log-odds Viterbi (optimal-path) and Forward (full-likelihood)
recurrences in log space against a null (background) model; it is exact for the modelled
glocal path. An additional **opt-in** layer reproduces HMMER's `hmmsearch` bit-score pipeline —
the local-multihit Forward score and the null2 biased-composition correction — verified against
pyhmmer/HMMER3 ground truth to single-precision rounding (see §5.3, §6.2).
An **opt-in** P-value / E-value layer reads the profile's `STATS LOCAL` calibration lines and
applies the HMMER significance statistics — a Gumbel survival function for MSV/Viterbi bit scores
and an exponential-tail survival function for Forward bit scores, with `E = P·Z` — exactly as the
Easel/HMMER survival functions compute them (see §2.3).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Many protein domains are defined by a curated multiple alignment summarised as a **profile HMM**,
not a short consensus pattern [1][4]. A profile HMM has, per alignment column (node), a Match (M)
state with position-specific residue-emission probabilities, an Insert (I) state for residues
between consensus columns, and a Delete (D) state for skipped columns; a mute Begin (B) and End
(E) flank the chain [2][3]. Pfam builds such HMMs with HMMER `hmmbuild` [4].

### 2.2 Core Model

**HMMER3/f file storage** [2, "HMMER profile HMM files"]: all probability parameters are stored
as **negative natural-log probabilities** to five decimals (e.g. probability 0.25 → −ln 0.25 =
1.38629); a zero probability is stored as `*`. For `ALPH amino` the alphabet size K = 20 in order
`ACDEFGHIKLMNPQRSTVWY`. The `COMPO` line holds the model's mean match-state composition, used as
the background residue composition of the null model. Each node has three lines: match emissions
(prefixed by node number), insert emissions, and 7 state transitions in the order
`Mk→Mk+1, Mk→Ik, Mk→Dk+1, Ik→Mk+1, Ik→Ik, Dk→Mk+1, Dk→Dk+1` [2, "State transition line"]. The
first two body lines after the optional `COMPO` are the BEGIN node: insert-0 emissions, then the
transitions `B→M1, B→I0, B→D1, I0→M1, I0→I0, 0.0, *` [2, "main model section"].

**Bit score** [2, "Null model"]: a HMMER bit score is the log of the ratio of the sequence's
probability under the profile (homology hypothesis) over its probability under the null model
(non-homology). The null model is a one-state i.i.d. background; emission probabilities are turned
into odds ratios against it. Match log-odds for residue `x` at state `k` is `log(p_k(x) / f_x)`
where `f_x` is the background frequency [3, "Scoring system for Profile HMM"].

**Viterbi log-odds recurrence** (Durbin et al. 1998, §5.4; reproduced verbatim in Stanford CS273
lecture 7 [3]), with emissions as log-odds against background `q`:

```
V^M_j(i) = log(e_Mj(x_i)/q_xi) + max{ V^M_{j-1}(i-1) + log a_{M(j-1)M(j)},
                                       V^I_{j-1}(i-1) + log a_{I(j-1)M(j)},
                                       V^D_{j-1}(i-1) + log a_{D(j-1)M(j)} }
V^I_j(i) = log(e_Ij(x_i)/q_xi) + max{ V^M_j(i-1) + log a_{M(j)I(j)},
                                       V^I_j(i-1) + log a_{I(j)I(j)},
                                       V^D_j(i-1) + log a_{D(j)I(j)} }
V^D_j(i) =                       max{ V^M_{j-1}(i)  + log a_{M(j-1)D(j)},
                                       V^I_{j-1}(i)  + log a_{I(j-1)D(j)},
                                       V^D_{j-1}(i)  + log a_{D(j-1)D(j)} }
```

The Forward recurrence is identical with `max` replaced by log-sum-exp [3, Durbin §3.6/§5.4].

### 2.3 E-value / P-value statistics

**STATS calibration lines** [2, HMM file format]: `STATS <s1> <s2> <f1> <f2>` — "`<f1>` and `<f2>`
are two real-valued parameters controlling location and slope of each distribution, respectively;
µ and λ for Gumbel distributions for MSV and Viterbi scores, and τ and λ for exponential tails for
Forward scores. λ values must be positive. All three lines or none of them must be present; when
all three are present, the model is considered to be calibrated for E-value statistics." All
parameters are in **bits** (HMMER reports scores in bits; λ ≈ log 2) [2][5].

**Distributions** [5, Eddy 2008]: Viterbi (and MSV) bit scores are **Gumbel** (Type-I extreme
value) distributed; the high-scoring tail of **Forward** scores is **exponentially** distributed,
both with parametric λ = log 2. The survival functions (P-value = P(S ≥ x)) are taken verbatim from
the Easel library used by HMMER:

```
Gumbel       P(S ≥ x) = 1 − exp(−exp(−λ(x − μ)))          [esl_gumbel_surv; with the |ey|<5e-9
                                                            tail branch returning exp(−λ(x−μ))]
Exponential  P(S ≥ x) = exp(−λ(x − τ))   (= 1 for x < τ)  [esl_exp_surv]
```

**E-value** [2, tutorial]: `E = P · Z`, where Z is the number of target sequences searched — a hit
"would be expected to happen Z times as often" (worked example: `1.2e-16 × 539165 = 6.47e-11`).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-HMM-01 | Forward score ≥ Viterbi score (same units) | Forward sums (log-sum-exp) over all paths including the single optimal one [3]. |
| INV-HMM-02 | A true-domain sequence scores strictly above a non-domain sequence of similar length | Log-odds is positive when the profile explains the sequence better than the background [2]. |
| INV-HMM-03 | Scoring is deterministic | The DP has no randomness; identical inputs give identical scores. |
| INV-HMM-04 | A `*` (zero-probability) parameter maps to −∞ in log space and propagates (the path is forbidden) | −ln 0 = +∞ stored; negated to −∞ log-prob [2]. |
| INV-HMM-05 | E-value is monotone decreasing in bit score and exactly linear in Z (`E = P·Z`) | Gumbel/exponential survival is strictly decreasing in x; the E-value multiplies P by Z [2][5]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `proteinSequence` | `string` | required | Amino-acid sequence | Case-insensitive; non-amino characters treated as background-odds (neutral) |
| `minBitScore` | `double` | 10.0 | Minimum Viterbi bit score to report a domain | bits |
| `accession` | `string` | required (`ScoreDomainHmm`) | Bundled Pfam accession | `PF00018`, `PF00595`, or `PF00400` |
| `databaseSize` (Z) | `double` | 1.0 | Number of target sequences searched (E-value scaling) | ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `ProteinDomain.Score` | `double` | Viterbi log-odds score in **bits** (nats / ln 2) |
| `ProteinDomain.Start/End` | `int` | 0-based, whole-sequence span (glocal score is a whole-sequence quantity) |
| `ProteinDomainHit.EValue` | `double` | Viterbi E-value (`E = P·Z`) from the profile's Gumbel `STATS` |
| `ScoreDomainHmm` return | `double` | Viterbi bit score for the one profile |
| `ScoreDomainHmmEValue` return | `(double BitScore, double EValue)` | bit score + Viterbi E-value |
| `Plan7ProfileHmm.Statistics` | `ScoreStatistics?` | parsed `STATS LOCAL` µ/λ/τ (null if uncalibrated) |

### 3.3 Preconditions and Validation

Null/empty sequence → empty result (`FindDomainsByHmm`) or `ArgumentNullException`
(`ScoreDomainHmm`/`Plan7ProfileHmm.Parse`). Unknown accession → `ArgumentException`. A profile
that is not HMMER3 or not `ALPH amino`, or is truncated, → `FormatException`. Residues are matched
case-insensitively; characters outside the 20 canonical amino acids contribute neutral (zero)
log-odds.

## 4. Algorithm

### 4.1 High-Level Steps

1. Parse the HMMER3/f profile: header (NAME/ACC/LENG/ALPH/GA), `COMPO` background, BEGIN node,
   then per-node match/insert emissions and 7 transitions. Convert stored −ln probabilities to ln.
2. For the query, run the glocal Viterbi (or Forward) DP over nodes 1..M and residues 1..n with the
   recurrence in §2.2, using match/insert log-odds against the `COMPO` background.
3. The score at E (from M_M / D_M) is the natural-log log-odds; divide by ln 2 for bits.
4. `FindDomainsByHmm` scores all three bundled profiles and reports those at/above `minBitScore`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Viterbi / Forward | O(n·M) | O(M) | n = sequence length, M = profile length; two-row DP |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [Plan7ProfileHmm.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/Plan7ProfileHmm.cs),
[ProteinMotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)

- `Plan7ProfileHmm.Parse(...)`: parses a HMMER3/f ASCII profile.
- `Plan7ProfileHmm.ViterbiScore(seq)` / `ForwardScore(seq)`: glocal log-odds score in nats (default mode).
- `Plan7ProfileHmm.LocalForwardScore(seq)`: **opt-in** HMMER local-multihit Forward score (nats).
- `Plan7ProfileHmm.LocalForwardBitScore(seq)`: local-multihit Forward bit score = HMMER `pre_score`.
- `Plan7ProfileHmm.Null2BiasBits(seq)`: the null2 biased-composition correction in bits (HMMER `bias`).
- `Plan7ProfileHmm.HmmSearchBitScore(seq)`: null2-corrected per-sequence bit score = `pre_score − bias`.
- `Plan7ProfileHmm.ViterbiPValue/MsvPValue/ForwardPValue(bits)`, `ViterbiEValue/ForwardEValue(bits, Z)`,
  and static `GumbelSurvival/ExponentialSurvival/EValue`: STATS-based significance (opt-in).
- `ProteinMotifFinder.FindDomainsByHmm(seq, minBitScore)`: detect SH3/PDZ/WD40 (opt-in, bit-score only).
- `ProteinMotifFinder.FindDomainHitsByHmm(seq, Z, minBitScore)`: detect with bit score **and** E-value.
- `ProteinMotifFinder.ScoreDomainHmm(seq, accession)`: Viterbi bit score for one bundled profile.
- `ProteinMotifFinder.ScoreDomainHmmEValue(seq, accession, Z)`: bit score + Viterbi E-value.
- `Plan7ProfileHmm.FindDomains(seq[, clusterOverlapping])`: **opt-in** HMMER `p7_domaindef`
  multi-domain envelope decomposition — returns one `DomainEnvelope` (env coords + null2-corrected
  bit score + i-Evalue) per domain. `clusterOverlapping` (default `true`) resolves a region flagged
  multi-domain by the `rt3` test via stochastic-traceback clustering (`region_trace_ensemble`); set
  `false` to keep the prior behaviour (such a region emitted as a single envelope).
- `ProteinMotifFinder.FindDomainEnvelopes(seq[, minBitScore])` / `FindDomainEnvelopes(seq, accession)`:
  decompose a protein into per-domain envelope hits against the bundled / one named profile.

Bundled CC0 profiles: [Resources/](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/Resources/) (PF00018/PF00595/PF00400).

### 5.2 Current Behavior

Profiles are loaded lazily once and cached. Insert emissions are scored with their stored
log-odds (HMMER sets these ≈ background, so the contribution is ≈ 0). Forward uses log-sum-exp for
numerical stability. The suffix tree was **not** used: this is a scoring-based DP alignment, not an
exact-substring search, so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- HMMER3/f ASCII parsing: −ln-probability storage, `*` → −∞, K=20 `ACDEFGHIKLMNPQRSTVWY`, COMPO
  background, BEGIN node, 7-transition node layout [2].
- Plan7 glocal Viterbi and Forward log-odds recurrences with match/insert/delete states [3].
- Bit score = log-odds (nats) / ln 2 against the null background [2][3].
- `STATS LOCAL MSV/VITERBI/FORWARD` parsing into µ/λ/τ; Gumbel survival for MSV/Viterbi and
  exponential-tail survival for Forward, verbatim from Easel's `esl_gumbel_surv` / `esl_exp_surv`
  (including the `|ey|<5e-9` Gumbel tail branch), and `E = P·Z` [2][5][6][7].
- **Opt-in HMMER `hmmsearch`-parity layer** (defaults unchanged): the **local-multihit** Forward
  score (occupancy-weighted local entry `t(B→M_k)=occ[k]/Σ occ[i]·(M−i+1)`, local exit `esc=0` for
  all k, multihit E→{C,J} split `−ln 2`, length config `pmove=(2+nj)/(L+2+nj)`), scored against
  HMMER's standard amino background (`p7_AminoFrequencies`, not COMPO), and the **null2
  biased-composition correction** by posterior expectation (Forward/Backward/decoding +
  `p7_GNull2_ByExpectation`), giving `score = (fwd − (nullsc + seqbias))/ln 2` with
  `seqbias = logsumexp(0, ln(1/256) + Σ ln null2[x_i])` [8][9][10][11][12][13][14]. **Verified to
  hmmsearch (pyhmmer 0.12.1):** local-multihit `pre_score` reproduces hmmsearch to ~1e-5 bits for
  PF00018/PF00595/PF00400 (68.709740 / 84.862930 / 213.411926); the null2 `bias` reproduces
  hmmsearch's reported single-domain bias to 3e-5 bits (SH3 envelope 0.025574).
- **Opt-in HMMER multi-domain envelope decomposition** (`p7_domaindef`; defaults unchanged):
  multihit Forward+Backward + posterior decoding to `mocc`/`btot`/`etot` (`p7_GDomainDecoding`:
  `btot[i]=btot[i-1]+P(B@i-1)`, `etot[i]=etot[i-1]+P(E@i)`, `mocc[i]=1−(N/J/C residue posteriors)`);
  region identification by the `rt1=0.25` trigger / `rt2=0.10` flank bound; the `is_multidomain_region`
  `rt3=0.20` test; per-envelope rescore in **unihit** mode at the full length n with per-domain bit
  score `(envsc + (n−Ld)·ln(n/(n+3)) − (nullsc + dombias))/ln 2`,
  `dombias = logsumexp(0, ln(1/256) + Σ ln null2[x])`, and i-Evalue `= Z·exp(−λ(score−τ))`
  [31][32][13][22]. **Verified to hmmsearch (pyhmmer 0.12.1):** the GBB1_HUMAN 7-blade WD40
  β-propeller (PF00400, L=340) decomposes into the **same 7 envelopes** (env 45-83, 87-125, 133-170,
  174-212, 216-254, 259-298, 303-340) with per-domain scores matching to ≈1e-3 bits and i-Evalues to
  ≥3 sig figs; a single SH3 domain gives one envelope (3-50, 68.54 bits) — HMMER uses float32, this
  engine float64.
- **Opt-in stochastic-traceback clustering of closely-overlapping multi-domain regions**
  (`region_trace_ensemble` → `p7_spensemble_Cluster`; default `clusterOverlapping=true`): for a region
  the `rt3` test flags multi-domain, sample `nsamples=200` stochastic tracebacks of the region's
  multihit Forward matrix (`p7_GStochasticTrace`: backward walk normalising predecessor log-scores with
  `esl_vec_FLogNorm` and sampling via `esl_rnd_FChoose` over a verbatim port of Easel's fixed-seed LCG
  RNG, `--seed 42`), index each into `B..E` segment pairs (`p7_trace_Index`), single-linkage-cluster the
  ensemble (`link_spsamples`: ≥0.8 overlap of the smaller segment in seq+hmm, start/end within 4
  diagonals), keep clusters with posterior ≥0.25, take consensus endpoints (widest with frequency ≥
  `ceil(ninc·0.02)`), drop dominated clusters (≥0.8 mutual overlap), and rescore each consensus envelope
  by the same null2-corrected per-domain path [31][34][35][36][37][38][39]. **Verified to hmmsearch
  (pyhmmer 0.12.1):** a closely-overlapping tandem-SH3 construct (two SH3 cores with the first truncated)
  splits into the **same overlapping envelopes** — trim=4 → 1-46 / 45-92; trim=12 → 1-37 / 37-84;
  trim=16 → 1-33 / 33-80 — envelope coordinates **exact**, per-domain scores within ~0.06 bits; the
  ensemble is deterministic (fixed-seed LCG, reproducible across runs).

**Intentionally simplified:**

- **Default** alignment mode remains **glocal full-profile** (`ViterbiScore`/`ForwardScore`,
  B→1..M→E spanning the whole query); the local-multihit parity path is opt-in via the
  `LocalForward*`/`HmmSearchBitScore` methods. The glocal `FindDomainsByHmm`/`ScoreDomainHmm` path is
  unchanged: it reports whole-sequence spans, not sub-envelopes.

**Not implemented:**

- **Exact RNG-bit parity** of the stochastic-traceback ensemble: the clustering reproduces the
  domain **count** and envelope **coordinates** of `hmmsearch` exactly with a verbatim port of Easel's
  fixed-seed LCG, but bit-for-bit identity of every sampled trace would additionally require HMMER's
  single-precision (`float`) Forward matrix and `esl_vec_FLogNorm`/`esl_rnd_FChoose` float arithmetic
  (this engine computes the Forward matrix in `double`). The consensus endpoints (a tally over 200
  samples) are robust to that difference, so coordinates match; per-sample trace-by-trace bit parity
  is the residual.
- **MSV / bias prefilters** are not reimplemented: they only gate which sequences reach the Forward
  stage; they do not change a reported hit's bit score, so they are not needed for score parity.
- The full Pfam library beyond the three bundled (caller-supplied `.hmm`) profiles is out of scope.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Default mode is glocal; local-multihit is opt-in | Design | `FindDomainsByHmm` spans whole query; `LocalForward*` give hmmsearch-mode scores | accepted | §5.3 |
| 2 | Multi-domain envelope decomposition is opt-in (`FindDomains`); stochastic clustering of overlapping-domain regions is implemented (default `clusterOverlapping=true`) but not RNG-bit-exact | Design | Well-separated AND closely-overlapping domains decompose with hmmsearch envelope-coordinate parity; only per-sample trace-by-trace bit parity (float32 RNG/DP) is unreached | accepted | §5.3 "Implemented"/"Not implemented"; §6.2 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null / empty sequence | empty result / `ArgumentNullException` | contract §3.3 |
| Non-domain sequence | bit score well below threshold; not reported | log-odds negative [2] |
| Unknown accession | `ArgumentException` | only 3 profiles bundled |
| `*` parameter on the optimal path | path forbidden (−∞) | −ln 0 [2] |
| Lower-case / mixed-case input | matched case-insensitively | implementation contract |

### 6.2 Limitations

Only three Pfam domains are bundled (SH3, PDZ, WD40); the full Pfam library is not embedded (it is
caller-supplied via the `.hmm` parser). The HMMER `hmmsearch` **bit-score pipeline is now reproduced
(opt-in) and verified against pyhmmer 0.12.1**: the local-multihit Forward `pre_score` matches
hmmsearch to ~1e-5 bits, and the null2 biased-composition `bias` matches hmmsearch's reported
single-domain value to 3e-5 bits. HMMER's automatic **multi-domain envelope decomposition**
(`p7_domaindef`) is now reproduced (opt-in `FindDomains`/`FindDomainEnvelopes`) and **verified against
pyhmmer 0.12.1**: the GBB1 7-blade WD40 β-propeller decomposes into the same 7 envelopes with
per-domain scores to ≈1e-3 bits and i-Evalues to ≥3 sig figs; a single SH3 gives one matching
envelope. The **stochastic-traceback clustering** of *closely-overlapping* domains
(`region_trace_ensemble` → `p7_spensemble_Cluster`, the path a region the `rt3` test flags as
multi-domain takes) is now also reproduced (opt-in, default on) and **verified against pyhmmer
0.12.1**: a closely-overlapping tandem-SH3 construct splits into the same overlapping envelopes
(coordinates exact, scores within ~0.06 bits, deterministic). The only remaining residual is exact
RNG-bit parity of the per-sample trace ensemble (HMMER samples in float32 with a fixed-seed LCG; this
engine ports the LCG verbatim but computes the Forward matrix in float64 — the consensus endpoints are
robust to this, so envelope coordinates match). The MSV/bias prefilters are not reimplemented (they gate which sequences reach
Forward; they do not change a reported hit's score). The default `FindDomainsByHmm`/glocal path
reports whole-sequence spans, not per-domain envelopes.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
foreach (var d in ProteinMotifFinder.FindDomainsByHmm(proteinSequence))
    Console.WriteLine($"{d.Name} {d.Accession}: {d.Score:F1} bits");

double sh3Bits = ProteinMotifFinder.ScoreDomainHmm(sh3Sequence, "PF00018");
var (bits, eValue) = ProteinMotifFinder.ScoreDomainHmmEValue(sh3Sequence, "PF00018", databaseSize: 1000);
```

**Numerical walk-through (DP arithmetic pin):** a tiny 2-match-state HMM over {A,C} with
background `qA=0.6, qC=0.4`, match emissions `M1: A=0.7,C=0.3` / `M2: A=0.2,C=0.8`, and
transitions `B→M1=0.9, M1→M2=0.8`, scoring "AC" by the path B→M1(A)→M2(C)→E:
`ln(0.7/0.6) + ln 0.9 + ln(0.8/0.4) + ln 0.8 = 0.5187937934151676` nats. The engine reproduces
this exactly (tested to 1e-9).

**E-value pin (STATS statistics):** with the bundled `PF00018_SH3_1.hmm` Viterbi STATS
`µ=−8.2932, λ=0.71923` and Forward STATS `τ=−4.5735, λ=0.71923`, at a bit score `S = 40`:
Gumbel `P = 1 − exp(−exp(−λ(S−µ))) = 8.227179545686635e-16` (Easel tail branch) → `E(Z=1000) =
8.227179545686635e-13`; exponential `P = exp(−λ(S−τ)) = 1.1943390031599535e-14` → `E(Z=1000) =
1.1943390031599535e-11`. The engine reproduces both to 1e-9 relative.

**hmmsearch local-multihit + null2 parity pin (pyhmmer 0.12.1 ground truth):** for SRC_HUMAN SH3 vs
`PF00018`, `LocalForwardBitScore = 68.7097` bits (hmmsearch `pre_score` 68.709740); the null2 bias
over the domain envelope (positions 3–50) `= 0.02554` bits (hmmsearch `bias` 0.025574); PDZ vs
`PF00595` pre = 84.8629 (hmmsearch 84.862930); WD40 vs `PF00400` pre = 213.4120 (hmmsearch
213.411926). All within single-precision rounding. A 1-node hand HMM emitting A (B→M1=1) gives
`LocalForwardScore("A") = 1.272400756045032` nats exactly.

**Multi-domain decomposition pin (pyhmmer 0.12.1 ground truth):** the GBB1_HUMAN 7-blade WD40
β-propeller (P62873, L=340) vs `PF00400` decomposes into **7** envelopes — `env_from..env_to`
(1-based) `45-83 / 87-125 / 133-170 / 174-212 / 216-254 / 259-298 / 303-340` with per-domain bit
scores `31.139467 / 19.004278 / 25.053679 / 35.552242 / 40.454269 / 23.443121 / 27.824228` and
i-Evalues `1.21e-11 / 8.41e-08 / 1.02e-09 / 4.85e-13 / 1.36e-14 / 3.31e-09 / 1.36e-10`.
`Plan7ProfileHmm.FindDomains` reproduces the envelope bounds **exactly** and the scores/i-Evalues to
single precision (HMMER float32 vs this engine's float64). A single SH3 (P12931, L=55) vs `PF00018`
gives **1** envelope `3-50`, score `68.540695`, i-Evalue `1.45e-23`.

```csharp
foreach (var d in ProteinMotifFinder.FindDomainEnvelopes(gbb1Sequence))
    Console.WriteLine($"{d.Name} env {d.EnvelopeStart}-{d.EnvelopeEnd}: {d.BitScore:F1} bits, E={d.IndependentEValue:E2}");
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProteinMotifFinder_FindDomainsByHmm_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs) — covers `INV-HMM-01`..`INV-HMM-04`, the H18 local-mode + null2 hmmsearch-parity cases, the H19 multi-domain decomposition (GBB1 7-domain / SH3 1-domain pyhmmer parity), and the H20 stochastic-traceback clustering of closely-overlapping domains (tandem-SH3 pyhmmer parity)
- Evidence: [PROTMOTIF-DOMAIN-001-Evidence.md](../../../docs/Evidence/PROTMOTIF-DOMAIN-001-Evidence.md)
- Related algorithms: [Domain_Prediction](./Domain_Prediction.md) (exact PROSITE patterns)

## 8. References

1. Mistry J, et al. 2021. Pfam: The protein families database in 2021. *Nucleic Acids Res* 49:D412–D419. https://doi.org/10.1093/nar/gkaa913
2. Eddy SR & the HMMER team. 2023. HMMER User's Guide, v3.4. http://hmmer.org (PDF: http://eddylab.org/software/hmmer/Userguide.pdf) — "HMMER profile HMM files"; "The HMMER profile/sequence comparison pipeline".
3. Durbin R, Eddy SR, Krogh A, Mitchison G. 1998. *Biological Sequence Analysis*, Ch. 5.4 (profile-HMM Viterbi/Forward recurrences); reproduced in Stanford CS273 Lecture 7. https://web.stanford.edu/class/cs273/scribing/scribe7.pdf
4. Eddy SR. 2011. Accelerated Profile HMM Searches. *PLoS Comput Biol* 7:e1002195. https://doi.org/10.1371/journal.pcbi.1002195
5. Eddy SR. 2008. A Probabilistic Model of Local Sequence Alignment That Simplifies Statistical Significance Estimation. *PLoS Comput Biol* 4:e1000069. https://doi.org/10.1371/journal.pcbi.1000069 (PMC2396288) — Viterbi/MSV Gumbel (λ = log 2), Forward exponential tail.
6. Eddy SR & contributors. Easel `esl_gumbel.c` — `esl_gumbel_surv(x, µ, λ) = 1 − exp(−exp(−λ(x−µ)))` (survivor function). https://github.com/EddyRivasLab/easel/blob/master/esl_gumbel.c
7. Eddy SR & contributors. Easel `esl_exponential.c` — `esl_exp_surv(x, µ, λ) = exp(−λ(x−µ))`, `=1` for `x<µ`. https://github.com/EddyRivasLab/easel/blob/master/esl_exponential.c
8. HMMER `modelconfig.c` — `p7_ProfileConfig` local entry/exit (`occ[k]/Z`, `esc=0`), `p7_ReconfigLength` (`pmove=(2+nj)/(L+2+nj)`). https://github.com/EddyRivasLab/hmmer/blob/master/src/modelconfig.c
9. HMMER `p7_hmm.c` — `p7_hmm_CalculateOccupancy`. https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_hmm.c
10. HMMER `generic_fwdback.c` — `p7_GForward` / `p7_GBackward` (local `esc=0`). https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_fwdback.c
11. HMMER `generic_decoding.c` — `p7_GDecoding` (posterior probabilities). https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_decoding.c
12. HMMER `generic_null2.c` — `p7_GNull2_ByExpectation` (posterior-expectation null2). https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_null2.c
13. HMMER `p7_domaindef.c` / `p7_pipeline.c` — per-domain `domcorrection`, `seqbias = logsum(0, ln ω + Σn2sc)`, `pre_score`, `bias`. https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_pipeline.c
14. HMMER `p7_bg.c` (`omega=1/256`, `p7_bg_NullOne`, `p1=L/(L+1)`) and `hmmer.c` (`p7_AminoFrequencies`); ground-truth scores via pyhmmer 0.12.1. https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_bg.c
22. HMMER `modelconfig.c` — `p7_ReconfigUnihit` / `p7_ReconfigMultihit` / `p7_ReconfigLength` (envelope rescore length config). https://github.com/EddyRivasLab/hmmer/blob/master/src/modelconfig.c
31. HMMER `p7_domaindef.c` — `p7_domaindef_ByPosteriorHeuristics` (region identification `rt1=0.25`/`rt2=0.10`, `is_multidomain_region` `rt3=0.20`, `region_trace_ensemble`, `rescore_isolated_domain`). https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_domaindef.c
32. HMMER `generic_decoding.c` — `p7_GDomainDecoding` (`btot`/`etot`/`mocc`). https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_decoding.c
33. pyhmmer 0.12.1 — HMMER3 Cython binding; multi-domain + overlapping-domain `hmmsearch` ground truth. https://pyhmmer.readthedocs.io/
34. HMMER `generic_stotrace.c` — `p7_GStochasticTrace` (stochastic traceback of the Forward matrix). https://github.com/EddyRivasLab/hmmer/blob/master/src/generic_stotrace.c
35. HMMER `p7_spensemble.c` — `p7_spensemble_Add` / `p7_spensemble_Cluster` / `link_spsamples` (consensus envelope clustering). https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_spensemble.c
36. HMMER `p7_trace.c` — `p7_trace_Index` (split a trace into `B..E` domains). https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_trace.c
37. HMMER `p7_pipeline.c` — pipeline RNG (`esl_randomness_CreateFast(42)`, `do_reseeding`). https://github.com/EddyRivasLab/hmmer/blob/master/src/p7_pipeline.c
38. Easel `esl_cluster.c` — `esl_cluster_SingleLinkage` (generalized single-linkage clustering). https://github.com/EddyRivasLab/easel/blob/master/esl_cluster.c
39. Easel `esl_random.c` (LCG/knuth, `esl_rnd_FChoose`) + `easel.c` (`esl_mix3`). https://github.com/EddyRivasLab/easel/blob/master/esl_random.c
