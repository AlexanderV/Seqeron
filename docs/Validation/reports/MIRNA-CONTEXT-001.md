# Validation Report: MIRNA-CONTEXT-001 — TargetScan context++ Scoring

- **Validated:** 2026-06-25   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus` (+ SA wiring via the McCaskill partition function; PCT via Friedman 2009 Bls + sigmoid)
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Authoritative sources opened THIS session

| # | Source | Retrieval |
|---|--------|-----------|
| 1 | Agarwal V, Bell GW, Nam JW, Bartel DP (2015) eLife 4:e05005, "Predicting effective microRNA target sites in mammalian mRNAs" (the context++ model: one MLR per site type; features site-type intercept, Local_AU, Min_dist, 3P_score, SA, sRNA1/sRNA8, PCT …). doi:10.7554/eLife.05005 | peer-reviewed model |
| 2 | `Agarwal_2015_parameters.txt` (TargetScan distribution) — the fitted coefficients + min/max per site type. `https://raw.githubusercontent.com/nsoranzo/targetscan/main/Agarwal_2015_parameters.txt` | curl, verbatim |
| 3 | `targetscan_70_context_scores.pl` — reference implementation: `getAgarwalContribution`, `getLocalAU_contribution`, `get_sRNA1_8_contributions`, `getSite8_contribution`, `getSA_contribution`, `get3primePairingContribution`, `getMinDist_weighted_contribution`, `get_len3UTR_weighted_contribution`, `getOffset6mer_weighted_contribution`/`getOffset6merSites`, `runRNAplfold_all_UTRs`. `https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_context_scores.pl` | curl, verbatim |
| 4 | `targetscan_70_BL_PCT.pl` — `calculatePCTthisBL` (Bls→PCT logistic). Same mirror. | curl, verbatim |
| 5 | Friedman RC, Farh KK, Burge CB, Bartel DP (2009) Genome Res 19:92 — branch-length score (Bls) definition. | model |

## Stage A — Description

**The context++ model.** Agarwal 2015 fits a separate multiple-linear-regression model for each
of the four canonical seed-match site types (8mer / 7mer-m8 / 7mer-A1 / 6mer). The score is the
SUM of an intercept plus per-feature `coeff × value` terms. `targetscan_70_context_scores.pl`
`getAgarwalContribution` confirms the exact arithmetic: features in
`{TA_3UTR, SPS, Local_AU, 3P_score, SA, Len_ORF, Len_3UTR, Min_dist, PCT}` are min-max scaled
`scaled = (raw − min)/(max − min)` then multiplied by the per-site-type coefficient; the indicator
features (`sRNA1A/C/G`, `sRNA8A/C/G`, `Site8A/C/G`) and the count features (`Off6m`, `ORF8m`) are
used **raw** (not scaled). Scaled values are **not** clamped to [0,1].

**Feature definitions confirmed against the perl:**
- *Intercept* — raw per-site-type term.
- *Local_AU* — `getLocalAU_contribution`: position-weighted A/U fraction over the 30 nt up- and
  downstream of the site; upstream weight `1/(i+1)` for 8mer(3)/7mer-m8(2), else `1/(i+2)`;
  downstream weight `1/(i+2)` for 8mer(3)/7mer-A1(1), else `1/(i+1)`; then min-max scaled.
- *sRNA1/sRNA8* — `get_sRNA1_8_contributions`: only scored when miRNA nt1 (resp. nt8) ≠ U; the
  indicator for the actual A/C/G is the coefficient (others 0).
- *Site8* — `getSite8_contribution`: only for 7mer-A1 and 6mer; identity of the target base
  opposite miRNA position 8 (`substr(subseqForAlignment,14,1)`; 0 when U).
- *SA* — `getSA_contribution`: `log10` of the unpaired probability of the 14-nt RNAplfold window
  (`-u` column 14) of the row 7 nt 3′ of the seed-match start → the 14-nt region centred on the
  match to miRNA nt 7–8 (eLife Fig 4A), then min-max scaled. Missing/zero plfold ⇒ contribution 0.
- *3P_score* — `get3primePairingContribution`: best single-gap-offset 3′-supplementary pairing run
  score (1.0 per pair in offset-adjusted positions 4–7, else 0.5; only runs ≥2; minus
  `max(0,(offset−2)/2)`), min-max scaled (min 1, max 3.5).
- *Min_dist* — `log10` of distance to the nearest 3′UTR end, min-max scaled.
- *Len_3UTR* — `log10` of the 3′UTR length, min-max scaled.
- *Off6m* — count of `first-6-of-revcomp(seed nt2-8)` occurrences in the 3′UTR, used raw.
- *SPS / TA_3UTR / Len_ORF / ORF8m* — require a seed-stability table / transcriptome / the ORF,
  i.e. data the library cannot derive from the miRNA+3′UTR alone. Caller-supplied (documented
  boundary), each computed FAITHFULLY (verbatim coeff × scaling) when supplied.
- *PCT* — Friedman 2009 branch-length score (Bls) of the minimal subtree connecting the conserved
  species → PCT via `calculatePCTthisBL` `pct = b0 + b1/(1 + e^((0−b2)·BL + b3))` (truncate at 0),
  then min-max scaled by the Agarwal PCT coefficient. The per-family b0..b3 are in TargetScan's
  compiled data tables (not published as numbers in Friedman 2009) ⇒ caller-supplied.

**Coefficient byte-check (every row of `Agarwal_2015_parameters.txt` vs the C# constants):** all
22 rows match exactly — Intercept (−0.589 / −0.224 / −0.195 / −0.079), 3P_score, Len_3UTR, Len_ORF,
Local_AU, Min_dist, Off6m, ORF8m, PCT (−0.103/−0.048/−0.048/0.005, min 0, max 0.816/0.364/0.449/0.193),
SA (−0.115/−0.134/−0.077/−0.028, mins −4.356/−5.218/−4.23/−5.082, maxs −0.661/−0.725/−0.588/−0.666),
Site8A/C/G (7mer-A1 & 6mer only), SPS, sRNA1A/C/G, sRNA8A/C/G, TA_3UTR. Column order 8mer / 7mer-m8 /
7mer-A1 / 6mer in both. No coefficient transcription error.

**Independent cross-check (hand computation, recorded numbers):**
- 8mer all-G flanks: Local_AU = −0.254×((0−0.308)/(0.814−0.308)) = **0.154608695652174**;
  Min_dist (min 5) = 0.118×((log10 5 −1.415)/(3.113−1.415)) = **−0.049759446106213065**;
  Len_3UTR (len 18) = 0.310×((log10 18 −2.392)/(3.637−2.392)) = **−0.2830405810586145**;
  3P raw 0 = −0.040×((0−1)/2.5) = **+0.016**; Off6m 1× = **−0.020**; sRNA8G(8mer) = **+0.015**;
  intercept −0.589 ⇒ partial = **−0.7561913315126536**.
- 6mer (miR-21 fixture) downstream-only A/U: fraction = **0.357142857142857**, Local_AU = **−0.005211141060198**.
- Supplied inputs (8mer): SPS(−8) = **0.11716577540106952**, TA(3.5) = **0.11424734042553189**,
  Len_ORF(1000) = **0.04503626943005184**, ORF8m(2) raw = **−0.236**.
- PCT: Bls{A,B}=3.0 → PCT = 1/(1+e⁻³) = **0.952574126822433**; 8mer contribution
  −0.103×(PCT/0.816) = **−0.120239136106263**. 7mer-m8 Bls{A,C}=7.0 → PCT = **0.999088948805599**,
  contribution −0.048×(PCT/0.364) = **−0.131747993249090**.

Every hand-computed number reproduced the corresponding test expectation EXACTLY. The min-max
scaling and "scale-these-only" feature set in `getAgarwalContribution` match the C#
`ScaledContribution`. The sigmoid in `calculatePCTthisBL` matches `PctFromBranchLength`.

Stage A **PASS** — feature definitions, coefficients, scaling, and site-type handling are all
confirmed against the eLife model + the parameter file + the perl reference.

## Stage B — Implementation

Code path: `src/.../MiRnaAnalyzer.cs` lines 486–1527 (`ScoreTargetSiteContextPlusPlus` and the
per-feature helpers).

- **Coefficients** (lines 619–709) match `Agarwal_2015_parameters.txt` byte-exact (verified above).
- **`ScaledContribution`** (1049) = `coeff×((raw−min)/(max−min))` (or `coeff×raw` when min==max),
  unclamped — matches the perl `getAgarwalContribution` scaled branch; the non-scaled features
  (sRNA1/8, Site8, Off6m, ORF8m) pass the value raw, matching the perl `else` branch.
- **Local_AU** (872): upstream/downstream walk, A/U test, the `1/(i+1)` vs `1/(i+2)` weight split
  per site type, division by `maxRaw`, then min-max scale — faithful port; 0/1-based coordinate
  mapping verified (utrUpStart = utrStart−31 ⇒ idx siteStart−1−i; utrDownStart = utrEnd ⇒ idx siteEnd+1+i).
- **sRNA1/sRNA8** (929/947): U-gate then A/C/G indicator → coefficient — matches the perl `ne "U"` gate.
- **Site8** (968): only 7mer-A1/6mer; reads the base opposite miRNA pos 8 (`mrna[Start−1]`); U-gate.
- **SA** (1000): window-fit check; **plfold = `RnaSecondaryStructure.CalculateRegionUnpairedProbability`**,
  which is `Z_open/Z` from the Turner-2004 **McCaskill partition function** (`FillPartitionDp` with the
  window's bases forbidden from pairing) — **not** an MFE approximation (verified by reading
  RnaSecondaryStructure.cs §2645–2719). `log10`, then min-max scale. Window mapping (utrStart+7 row,
  L=14) matches `getSA_contribution`; out-of-fit ⇒ 0 + reported omitted (perl: missing plfold ⇒ 0).
- **3P_score** (1063): full port of `extractSubseqForAlignment`/`modifySubseqForAlignment`/
  `get3primePairingContribution` (reverse, base-code product Watson-Crick test, run scoring,
  offset penalty), min-max scaled. Two raw scores (4.5, 6) and the 8mer raw=6 fixture are
  documented as cross-checked against the perl run; the contribution arithmetic verified.
- **Min_dist / Len_3UTR / Off6m** (1277/1304/1324): single-isoform (AIR=1) ports; log10/raw + scale.
  Off6m pattern = first 6 nt of `reverse_complement(seed nt2-8)` (T→U) counted over the UTR — matches
  `getOffset6merSites`.
- **PCT** (1435/1506/1517): `ComputeBranchLengthScore` (Friedman 2009 connecting-subtree branch-length
  sum), `PctFromBranchLength` (verbatim sigmoid, truncate at 0), `PctContribution` (min-max scale).
- **Partial sum** (818) sums exactly the 15 realised contributions; `BuildOmittedFeatures` reports
  the honest residual (SA when window doesn't fit; PCT/SPS/TA/Len_ORF/ORF8m when not supplied).
- **Contract**: non-seed site types throw `ArgumentException` (770).

**Test-quality audit** (`MiRnaAnalyzer_TargetPrediction_Tests.cs` CTX-001…CTX-PCT-006, ~15 cases):
expected values are hand-derived from the parameter file + perl logic (per-test derivation comments),
not code echoes — I independently reproduced CTX-001/004/011 and all PCT cases and they match to the
last digit. Coverage: all four site types; intercept; Local_AU (zero-fraction and mixed-flank);
sRNA1 non-U branch; sRNA8 U-gate; Site8 (7mer-A1/6mer only); 3P_score (perl-verified raws);
Min_dist; Len_3UTR; Off6m (raw count, 1× and 2×); caller-supplied SPS/TA/Len_ORF/ORF8m and their
removal from the residual; SA computed from McCaskill (independent recompute, asserts ≠0 and not in
residual, and presence in the partial sum — so it would fail a wrong coefficient/window or an MFE
implementation); PCT (Bls on a worked tree, sigmoid, negative-truncation, per-site-type params,
removal from residual); invalid-site-type throw; honest-residual reporting. No green-washing
(no bare "no-throw"/tautology assertions on the numeric paths).

**Cross-verification:** see Stage-A numbers; all recomputed values equal the test expectations and
the parameter file.

**Findings / defects:** none.

## Verdict & follow-ups

- **Stage A: ✅ PASS · Stage B: ✅ PASS · State: ✅ CLEAN.**
- The score is, by design, a PARTIAL context++ (the locally-computable subset + faithfully-computed
  caller-supplied features). Full-transcript inputs (TA / SPS / Len_ORF / ORF8m) and conservation
  (PCT) are **caller-supplied** because they require a seed-stability table / transcriptome / the
  ORF / a multi-species alignment+tree that the library cannot derive from a single miRNA+3′UTR.
  This is an explicit, documented boundary (each computed verbatim when supplied; reported in
  `OmittedFeatures` otherwise) — **NOT** a limitation. The computable-subset coefficients + feature
  math verify byte-exact against Agarwal 2015 / the perl reference.
- SA derives from the validated McCaskill partition function (RNA-ACCESS-001), confirmed by source
  read — not an MFE proxy.
- Full unfiltered `dotnet test Seqeron.sln -c Debug` green: Seqeron.Genomics.Tests 18762 passed /
  0 failed (all other projects pass; the empty Annotation.Tests project reports "no test available").
- No defect logged.
