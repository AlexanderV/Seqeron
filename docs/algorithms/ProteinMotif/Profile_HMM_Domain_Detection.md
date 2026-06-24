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
glocal path but does not reproduce HMMER's full `hmmsearch` bit-score pipeline (see §6.2).
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
- `Plan7ProfileHmm.ViterbiScore(seq)` / `ForwardScore(seq)`: glocal log-odds score in nats.
- `Plan7ProfileHmm.ViterbiPValue/MsvPValue/ForwardPValue(bits)`, `ViterbiEValue/ForwardEValue(bits, Z)`,
  and static `GumbelSurvival/ExponentialSurvival/EValue`: STATS-based significance (opt-in).
- `ProteinMotifFinder.FindDomainsByHmm(seq, minBitScore)`: detect SH3/PDZ/WD40 (opt-in, bit-score only).
- `ProteinMotifFinder.FindDomainHitsByHmm(seq, Z, minBitScore)`: detect with bit score **and** E-value.
- `ProteinMotifFinder.ScoreDomainHmm(seq, accession)`: Viterbi bit score for one bundled profile.
- `ProteinMotifFinder.ScoreDomainHmmEValue(seq, accession, Z)`: bit score + Viterbi E-value.

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

**Intentionally simplified:**

- Alignment mode: **glocal full-profile** path (B→1..M→E spanning the whole query); **consequence:**
  the score is a whole-sequence quantity, not HMMER's local multihit (N/C/J flanks, per-domain
  envelopes). A domain embedded in a long protein is still detected by its strong log-odds, but
  Start/End are reported as the whole sequence rather than a sub-envelope.

**Not implemented:**

- Exact `hmmsearch`-reported E-value **pipeline parity**: the P-value/E-value formulas are exact
  given a bit score, but HMMER applies them to its *local-multihit* sequence bit score after the
  MSV/bias prefilters and the **null2 biased-composition correction**, which this glocal scorer does
  not compute. **Consequence:** for the same query the absolute bit score (and therefore the
  absolute E-value) differs from `hmmsearch`; ranking and detection are faithful. **Users should
  rely on:** HMMER `hmmsearch` for exact reported-E-value parity. The full Pfam library beyond the
  three bundled (caller-supplied `.hmm`) profiles is likewise out of scope.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Glocal vs local-multihit alignment | Deviation | Start/End span whole query; no sub-envelopes | accepted | §5.3 "Intentionally simplified" |
| 2 | Absolute bit score ≠ `hmmsearch` (no null2 / MSV-bias pipeline) | Deviation | Ranking correct; reported E-value not pipeline-calibrated | accepted | honest residual; §6.2 |

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
caller-supplied via the `.hmm` parser). The Gumbel/exponential P-value and `E = P·Z` formulas are
implemented exactly (and read the profile's own `STATS LOCAL` calibration). What is **not**
reproduced is exact `hmmsearch`-reported E-value *pipeline* parity: HMMER applies those formulas to
its local-multihit sequence bit score after the MSV/bias prefilters and the null2 biased-composition
correction, none of which this glocal scorer computes — so the absolute bit score (and hence the
absolute reported E-value) differs from `hmmsearch` even though ranking and detection are faithful.
Glocal scoring reports whole-sequence spans, not per-domain envelopes.

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

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ProteinMotifFinder_FindDomainsByHmm_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs) — covers `INV-HMM-01`..`INV-HMM-04`
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
