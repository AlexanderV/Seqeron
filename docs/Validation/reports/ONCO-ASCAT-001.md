# Validation Report: ONCO-ASCAT-001 — Allele-specific copy number (ASCAT fit + ASPCF segmentation + sub-clonal CN)

- **Validated:** 2026-06-24   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.SegmentAlleleSpecific`, `FitPurityPloidy`,
  `DeriveMultiplicity`, `SegmentAlleleSpecificAspcf`, `FitSubclonalCopyNumber`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

This is a NEW unit created during the limitations campaign (commits `940995dd` greedy
segmentation + ASCAT fit + multiplicity; `31cfe817` ASPCF + sub-clonal). It had never been
validated under the protocol. Both commits add NEW methods; no pre-existing oncology algorithm
was modified.

---

## Stage A — Description

### Sources opened (live, this session)

| Source | What it confirmed |
|--------|-------------------|
| `ascat.runAscat.R` (VanLoo-lab/ascat, master, raw fetch) | nA/nB equations, GoF distance `d`, `TheoretMaxdist`, `goodnessOfFit`, round/clamp |
| Nilsen et al. 2012, PMC3582591 (full text) | PCF criterion `L(S\|y,γ)=ΣΣ(y−ȳ)²+γ\|S\|`, DP recurrence `e_k=min_j(d_jk+e_{j−1}+γ)`, joint `L(y1)+L(y2)`, default γ=40 |
| DeCiFering, PMC8542635 (full text) | CCF closed form `c=vF/(ρM)`, `F=ρN_tot+2(1−ρ)` |
| McGranahan 2016 (Science) + PICTograph/DeCiFer cross-refs | n_mut = VAF·F/ρ; M = round(n_mut) at CCF=1, clamp to [1, major CN] |

### Formula check (verbatim vs source)

**ASCAT nA/nB** — `AscatRawCopyNumbers` (line 7382) vs `ascat.runAscat.R`:
```
R:    nA = (rho-1 - (b-1)*2^(r/gamma) * ((1-rho)*2+rho*psi))/rho
      nB = (rho-1 +  b   *2^(r/gamma) * ((1-rho)*2+rho*psi))/rho
code: scaledTotal = 2^(r/gamma)*((1-rho)*2+rho*psi)
      nA = (rho-1 - (b-1)*scaledTotal)/rho ;  nB = (rho-1 + b*scaledTotal)/rho
```
Byte-for-byte identical. ✅

**GoF distance + percentage** (lines 7472, 7510) vs R:
```
R:    d = Σ |nMinor − pmax(round(nMinor),0)|² · length · (b==0.5 ? 0.05 : 1)
      TheoretMaxdist = Σ 0.25 · length · (b==0.5 ? 0.05 : 1)
      goodnessOfFit = (1 − m/TheoretMaxdist)·100
code: minorDistance += minorDev² · weights[i];  weights = length · (balanced? 0.05 : 1)
      theoreticalMaxDistance += 0.25 · weights[i]
      GoF = (1 − bestMinorDistance/theoreticalMaxDistance)·100
```
Matches. ✅ (γ=1 for sequencing per ASCAT README — `AscatSequencingGamma = 1`.)

**PCF/ASPCF** — `SegmentChromosomeAspcf` (line 7679): `e[0]=0`,
`e[k]=min_j (e[j−1] + SSE_logR(j..k) + SSE_BAF(j..k) + γ)`, with `SSE = Σx² − (Σx)²/m`. Exactly
the Nilsen 2012 DP with the joint two-track cost (γ charged once per segment, separate per-track
means via `BuildSegmentSummary`). DP runs per chromosome, so no segment crosses a contig. ✅

**Multiplicity** — `DeriveMultiplicity` (line 7574): `n_mut = VAF·(ρ·N_T+2(1−ρ))/ρ`, round
away-from-zero, `Clamp(·,1,majorCN)`. = McGranahan n_mut at CCF=1 / DeCiFer F·v/(ρM) with M=n_mut. ✅

**Sub-clonal** — `FitSubclonalCopyNumber` (line 7825): integer-within-tolerance ⇒ clonal single
state f=1; else `n_obs = f·ceil + (1−f)·floor` with a single shared least-squares fraction over
both alleles, two ceil/floor pairings tried, smaller-residual kept; primary = larger-fraction
state (Battenberg frac1≥frac2, frac1+frac2=1). Matches the Battenberg two-population algebra. ✅

### Edge-case semantics

- Balanced BAF=0.5 segments down-weighted ×0.05 in GoF (sourced, line 7438). ✅
- 2n/4n degeneracy: tie broken by lower ploidy (Van Loo parsimony) + a both-allele selection
  distance (line 7476/7487). Sourced convention. ✅
- Multiplicity clamp to [1, majorCN]: an observed variant sits on ≥1 and ≤major copies. ✅
- γ→large ⇒ 1 segment; small γ ⇒ each level; chromosome boundary never crossed. ✅
- Integer n_obs ⇒ single clonal state. ✅

### Independent cross-checks (hand-computed, then matched to code/tests)

1. **ASCAT recovery, planted ρ=0.80 ψ=2.2, segment (nA=2,nB=0):**
   denom=0.8·2+0.4=2.0, D=0.8·2.2+0.4=2.16, r=log2(2.0/2.16)=−0.1110, b=0.2/2.0=0.1.
   scaledTotal=(2.0/2.16)·2.16=2.0 ⇒ nA=(−0.2−(−0.9)·2.0)/0.8=1.6/0.8=**2.0**,
   nB=(−0.2+0.1·2.0)/0.8=**0.0**. Exact (2,0). ✅
2. **Multiplicity (M6):** VAF=4/7, ρ=0.8, N_T=3 ⇒ F=2.8, n_mut=(4/7)·2.8/0.8=**2.0** ⇒ M=2. ✅
3. **Sub-clonal (M-SUB-1):** n_obs major=1.4, minor=0.6. Anti-monotone pairing
   hi=(2,0), lo=(1,1): da=1, db=−1 ⇒ f=(0.4+0.4)/2=**0.4**, residual 0 (beats co-monotone f=0.5,
   residual 0.02). ⇒ primary (1,1) at 0.6, secondary (2,0) at 0.4. Exactly the planted mixture. ✅

### Findings / divergences (Stage A)

None material. One simplification: ASCAT R selects `nMinor` **globally** (whichever of nA/nB has
the smaller total sum across segments) so the GoF uses one consistent allele as "minor"; the code
uses the **per-segment** `min(nA,nB)`. Because the GoF term is the squared distance of the *minor*
allele to its nearest integer and the forward model keeps nB ≤ nA (BAF mirrored to ≤0.5), the two
conventions coincide on every planted dataset and the squared-distance objective is symmetric to a
major/minor swap. This is benign (no defect). The documented residual (multi-sample asmultipcf,
3+-population mixtures, WGD refit search) is by-design and stays in scope.

---

## Stage B — Implementation

### Code path reviewed
`OncologyAnalyzer.cs`: `SegmentAlleleSpecific` 7282, `BuildSegmentSummary` 7354,
`AscatRawCopyNumbers` 7382, `FitPurityPloidy`+`ValidateGrid` 7411/7518, `DeriveMultiplicity` 7574,
`SegmentAlleleSpecificAspcf`+`SegmentChromosomeAspcf`+`AspcfSegmentSse` 7636/7679/7750,
`FitSubclonalCopyNumber`+`SolveSharedFraction` 7825/7918. Constants: `NormalDiploidCopyNumber=2`,
`AscatSequencingGamma=1`, `AscatWorstCaseIntegerDistance=0.25`, `AscatBalancedSegmentWeight=0.05`,
`BalancedBaf=0.5`, `AspcfDefaultPenalty=40`, `SubclonalIntegerTolerance=0.05`.

### Formula realised correctly?
Yes — every formula above maps line-for-line to the source. SSE uses the numerically stable prefix-
sum form `Σx²−(Σx)²/m` with a round-off guard; copy numbers round away-from-zero and clamp at 0.

### Cross-verification vs code (tests run)
`OncologyAnalyzer_AscatDerivation_Tests` (23 tests): **23 passed, 0 failed** (8 ms).
Covers M1–M12, S1–S3, C1, M-ASPCF-1..3, S-ASPCF-1..2, C-ASPCF-1, M-SUB-1..2, C-SUB-1. The planted-
truth values (ρ=0.80, ψ=2.2 and 3.0; M=1 and 2; CCF=1.0; ASPCF breakpoint at index 10; sub-clonal
f0=0.4) all match my hand computations above.

### Variant/delegate consistency
`EstimateCcf` delegate (ONCO-CCF-001) drives the end-to-end M9 path: fit → DeriveMultiplicity →
EstimateCcf ⇒ CCF=1.0 on the planted clonal mutation. `SegmentAlleleSpecific` (greedy) and
`SegmentAlleleSpecificAspcf` (DP global optimum) agree on clean data; M-ASPCF-2 verifies the DP
cost ≤ greedy cost on a noisy track (global-optimality invariant INV-7).

### Numerical robustness
Grid loops use `+1e-12` epsilon guards on the upper bounds; SSE guarded against tiny negative
round-off; div-by-zero avoided (denom checks in `SolveSharedFraction`, `m≤1 ⇒ SSE 0`). Validation
guards (`ValidateGrid`, arg checks) cover NaN, non-positive, and out-of-range inputs.

### Test quality audit
Tests assert exact sourced values (integer CN, GoF≈100, M, CCF, breakpoint index, mirrored BAF
1.0, f≈0.4) — not "no throw" tautologies. Inputs are synthesised by the exact algebraic inverse of
the cited ASCAT equations, independent of the implementation. Deterministic (fixed jitter array for
the noisy ASPCF case). Edge/contract cases (M10–M12, C-ASPCF-1, C-SUB-1) covered.

### Findings / defects (Stage B)
None.

---

## Verdict & follow-ups

- **Stage A: PASS** — all formulas (ASCAT nA/nB + GoF, Nilsen PCF/ASPCF DP, McGranahan/DeCiFer
  multiplicity & CCF, Battenberg two-state mixture) verified verbatim against live primary sources
  and reference code; three independent hand cross-checks reproduce the planted values.
- **Stage B: PASS** — code realises each formula line-for-line; 23/23 unit tests green and matching
  hand computations.
- **End state: CLEAN.** The only gaps are the by-design documented residuals (multi-sample
  asmultipcf, 3+-population per-segment mixtures, WGD refit search, odd-ploidy TAI branch); the
  per-segment vs global `nMinor` selection is a benign, GoF-invariant simplification. No code change
  required.

CodeChanged: no. FullSuite: N/A (no code modified; the unit's 23 tests run green).
