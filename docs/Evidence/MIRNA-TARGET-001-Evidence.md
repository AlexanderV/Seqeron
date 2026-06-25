# MIRNA-TARGET-001 Evidence: Target Site Prediction

## Sources

| Source | Key Finding | URL / Reference |
|--------|-------------|-----------------|
| Bartel (2009) | Defines seed site hierarchy: 8mer > 7mer-m8 > 7mer-A1 > 6mer; seed positions 2–8; antiparallel binding | PMID: 19167326, doi:10.1016/j.cell.2009.01.002 |
| Lewis et al. (2005) | "Conserved seed pairing, often flanked by adenosines" — adenosine at position 1 opposite correlates with efficacy | PMID: 15652477, doi:10.1016/j.cell.2004.12.035 |
| Grimson et al. (2007) | Context scores: site-type efficacy weights (8mer=0.31, 7mer-m8=0.161, 7mer-A1=0.099); AU-rich context favors targeting | PMID: 17612493, doi:10.1016/j.molcel.2007.06.017 |
| Agarwal et al. (2015) | TargetScan 7: context++ model; 6mer has minimal but detectable efficacy | PMID: 26267216, doi:10.7554/eLife.05005 |
| TargetScan 8.0 FAQ | Centered sites removed (not reliably functional); 3 canonical + offset 6mer | https://www.targetscan.org/vert_80/docs/help.html |
| Wikipedia: MicroRNA | Seed region = positions 2–7/8; target recognition in 3' UTR; antiparallel Watson-Crick+ G:U wobble | https://en.wikipedia.org/wiki/MicroRNA |
| miRBase | Authoritative miRNA sequences for test data | https://mirbase.org/ |
| Friedman et al. (2009) | Probability of conserved targeting (PCT); most mammalian mRNAs are miRNA targets | PMID: 18955434 |

## Seed Site Type Definitions (Bartel 2009, TargetScan)

Target sites are on the mRNA and pair **antiparallel** to the miRNA. The seed region is positions 2–8 (1-indexed from 5' end of miRNA). The reverse complement of the seed is sought in the mRNA (5'→3').

### Site Types (decreasing efficacy)

| Type | Definition | miRNA Positions | Additional Requirement |
|------|-----------|-----------------|----------------------|
| **8mer** | Match to miRNA positions 2-8 + A opposite position 1 | 2–8 + pos 1 | 'A' on mRNA at position opposite miRNA pos 1 |
| **7mer-m8** | Match to miRNA positions 2-8 | 2–8 | None |
| **7mer-A1** | Match to miRNA positions 2-7 + A opposite position 1 | 2–7 + pos 1 | 'A' on mRNA at position opposite miRNA pos 1 |
| **6mer** | Match to miRNA positions 2-7 | 2–7 | None |
| **Offset 6mer** | Match to miRNA positions 3-8 | 3–8 | None (marginal efficacy) |

### Orientation Details
- miRNA binds 3'→5' to mRNA 5'→3' (antiparallel)
- The "reverse complement" of the miRNA seed is what appears on the mRNA target
- Position 1 of miRNA is at the 3' end of the target site on mRNA

### Scoring Hierarchy (Grimson 2007)
From experimental data in mammalian cells:
- 8mer most effective (largest median fold repression)
- 7mer-m8 intermediate-high
- 7mer-A1 intermediate
- 6mer weak but detectable
- Offset 6mer marginal

## Worked Example: hsa-let-7a-5p

**miRNA**: `UGAGGUAGUAGGUUGUAUAGUU`
**Seed (pos 2-8)**: `GAGGUAG`
**Reverse complement of seed**: `CUACCUC`

For an 8mer site on mRNA: `...CUACCUCA...` (the 'A' at the 3' end corresponds to position 1 of miRNA)
For a 7mer-m8 site: `...CUACCUC...` (match to pos 2-8 only)
For a 7mer-A1 site: `...UACCUCA...` (match to pos 2-7 + 'A' at pos 1 position)
For a 6mer site: `...UACCUC...` (match to pos 2-7 only)

## Alignment Conventions
- miRNA aligns from 5'→3' (left to right)
- Target aligns from 3'→5' (left to right, i.e., reverse of mRNA direction)
- Watson-Crick pairs: A:U, G:C → '|'
- Wobble pairs: G:U → ':'
- Mismatches: ' '

## Edge Cases & Boundary Conditions

1. **Empty/null inputs**: FindTargetSites returns empty collection
2. **mRNA shorter than seed**: No possible matches, return empty
3. **Multiple sites in same mRNA**: All should be found independently
4. **Overlapping sites**: Same seed RC at positions that overlap — each should be reported
5. **DNA input (T instead of U)**: Implementation converts T→U before matching
6. **Case insensitivity**: Normalized to uppercase
7. **minScore threshold**: Sites below threshold are filtered out
8. **Score range**: All scores should be in [0.0, 1.0]
9. **8mer vs 7mer-m8 at same position**: 8mer should be reported (higher rank)

## Known Numerical Properties
- Score for 8mer ≥ Score for 7mer-m8 ≥ Score for 7mer-A1 ≥ Score for 6mer ≥ Score for offset 6mer
- FreeEnergy for well-paired duplex should be negative
- SeedMatchLength: 8 for 8mer, 7 for 7mer types, 6 for 6mer types

---

## TargetScan context++ Scoring (Agarwal et al. 2015) — opt-in

### Online Sources (retrieved verbatim this session)

| Source | How retrieved | Key fact extracted |
|--------|---------------|--------------------|
| `Agarwal_2015_parameters.txt` (TargetScan distribution) | `curl https://raw.githubusercontent.com/nsoranzo/targetscan/main/Agarwal_2015_parameters.txt` (2026-06-24) | Plaintext fitted coefficients (coeff/min/max) for 21 feature rows + Intercept, one column **per site type** (column order: 8mer / 7mer-m8 / 7mer-A1 / 6mer). |
| `targetscan_70_context_scores.pl` (TargetScan distribution) | `curl https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_context_scores.pl` (2026-06-24) | `getAgarwalContribution` (l.1678), `getLocalAU_contribution` (l.1081), `get_sRNA1_8_contributions`, `getSite8_contribution`, `get3primePairingContribution` (l.1234) + `extractSubseqForAlignment` (l.936) + `modifySubseqForAlignment` (l.1026), `getMinDist_weighted_contribution` (l.1182), `get_len3UTR_weighted_contribution` (l.1912), `getOffset6merSites` (l.2185) + `getOffset6mer_weighted_contribution` (l.2046), CS sum (l.404–411). Defines feature computation + scaling. 3P_score raw values reproduced by running these subs on test fixtures. |
| Agarwal V, Bell GW, Nam JW, Bartel DP (2015) eLife 4:e05005 | WebFetch https://pmc.ncbi.nlm.nih.gov/articles/PMC4532895/ (2026-06-24) | "Using the 14 robustly selected features, we trained multiple linear regression models on all of the data. The resulting models, **one for each of the four site types**, were collectively called the context++ model." Continuous features scaled to comparable ranges; nucleotide-identity features categorical. doi:10.7554/eLife.05005 |
| hsa-miR-122-5p | WebFetch https://mirbase.org/mature/MIMAT0000421 (2026-06-24) | Mature sequence `UGGAGUGUGACAAUGGUGUUUG` (cross-check of nt1/nt8 handling). |

### context++ model (Agarwal 2015)

The context++ score (CS) of a site is the **sum of per-feature contributions** of a multiple-linear-regression model fit **separately per site type** (8mer, 7mer-m8, 7mer-A1, 6mer). Each feature contribution = `coeff × scaledScore` (`getAgarwalContribution`):

- **Continuous features** (`TA_3UTR, SPS, Local_AU, 3P_score, SA, Len_ORF, Len_3UTR, Min_dist, PCT`) are **min-max scaled**: `scaled = (raw − min) / (max − min)`, using the per-site-type min/max from the parameter file.
- **Binary indicator features** (`sRNA1A/C/G, sRNA8A/C/G, Site8A/C/G, Off6m, ORF8m`) are **used raw** (the indicator is 1 for the matching nucleotide, else 0).
- **Intercept** (the site-type term) is the raw coefficient, added once.

### Fitted coefficients (verbatim, parameter-file column order 8mer / 7mer-m8 / 7mer-A1 / 6mer)

| Feature | 8mer | 7mer-m8 | 7mer-A1 | 6mer | min (per type) | max (per type) |
|---------|------|---------|---------|------|----------------|----------------|
| Intercept | -0.589 | -0.224 | -0.195 | -0.079 | — | — |
| Local_AU | -0.254 | -0.177 | -0.075 | -0.040 | 0.308 / 0.277 / 0.342 / 0.295 | 0.814 / 0.782 / 0.801 / 0.772 |
| sRNA1A | -0.018 | 0.010 | -0.025 | -0.002 | 0 | 1 |
| sRNA1C | -0.021 | 0.014 | -0.021 | 0.004 | 0 | 1 |
| sRNA1G | 0.060 | 0.062 | 0.030 | 0.018 | 0 | 1 |
| sRNA8A | 0.022 | 0.004 | -0.049 | -0.015 | 0 | 1 |
| sRNA8C | 0.012 | -0.031 | 0.033 | 0.016 | 0 | 1 |
| sRNA8G | 0.015 | -0.008 | -0.017 | 0.006 | 0 | 1 |
| Site8A | (n/a) | (n/a) | 0.000 | -0.002 | 0 | 1 |
| Site8C | (n/a) | (n/a) | 0.036 | 0.015 | 0 | 1 |
| Site8G | (n/a) | (n/a) | 0.015 | 0.012 | 0 | 1 |
| 3P_score | -0.040 | -0.055 | -0.060 | -0.024 | 1 / 1 / 1 / 1 | 3.5 / 3.5 / 3.5 / 3.5 |
| Min_dist | 0.118 | 0.056 | 0.045 | 0.036 | 1.415 / 1.491 / 1.431 / 1.477 | 3.113 / 3.096 / 3.117 / 3.106 |
| Len_3UTR | 0.310 | 0.154 | 0.129 | 0.045 | 2.392 / 2.409 / 2.413 / 2.405 | 3.637 / 3.615 / 3.630 / 3.620 |
| Off6m | -0.020 | -0.011 | -0.020 | -0.010 | 0 (used raw) | — |
| SPS | 0.210 | 0.135 | 0.095 | 0.035 | -11.13 / -11.13 / -8.41 / -8.57 | -5.52 / -5.49 / -3.33 / -3.33 |
| TA_3UTR | 0.222 | 0.139 | 0.117 | 0.058 | 3.113 / 3.067 / 3.145 / 3.113 | 3.865 / 3.887 / 3.887 / 3.887 |
| Len_ORF | 0.205 | 0.100 | 0.063 | 0.029 | 2.788 / 2.773 / 2.773 / 2.775 | 3.753 / 3.729 / 3.730 / 3.731 |
| ORF8m | -0.118 | -0.044 | -0.058 | -0.060 | 0 (used raw) | — |
| SA | -0.115 | -0.134 | -0.077 | -0.028 | -4.356 / -5.218 / -4.230 / -5.082 | -0.661 / -0.725 / -0.588 / -0.666 |
| PCT | -0.103 | -0.048 | -0.048 | 0.005 | 0 | 0.816 / 0.364 / 0.449 / 0.193 |

Site8 features are computed by the perl **only for 7mer-A1 (siteType 1) and 6mer (siteType 4)**; the parameter file leaves the 8mer/7mer-m8 cells blank. (Verbatim from `Agarwal_2015_parameters.txt`, column order 8mer / 7mer-m8 / 7mer-A1 / 6mer.)

### 3P_score — 3' supplementary pairing (verbatim from `get3primePairingContribution`)

The UTR subsequence and reversed mature miRNA are built by `extractSubseqForAlignment` (UTR window from `utrStart−16` to `utrEnd(+1)`, padded to `$DESIRED_UTR_ALIGNMENT_LENGTH = 23` with N's) and `modifySubseqForAlignment`. `get3primePairingContribution` then reverses both, drops the seed-paired prefix using per-site-type `seedinfo{utrstart}` / `mirnastart` / `overhang`, and for every single-gap offset (in both the UTR and the miRNA strand) sums, over each contiguous run of ≥2 base pairs, `+1.0` when the offset-adjusted position `(i+offset−overhang)` (top) or `(i−overhang)` (bottom) is in `4..7`, else `+0.5`; base pairs are A:U/U:A (`code product = 2`) or G:C/C:G (`= 12`). The per-offset run score has `max(0,(offset−2)/2)` subtracted; the raw 3P score is the maximum over all offsets/orientations. Then `getAgarwalContribution(type,"3P_score",raw)` min-max scales with min=1, max=3.5.

### Min_dist (verbatim from `getMinDist_weighted_contribution`)

Single-isoform: `distTo5primeEndOfUTR = siteStart − 1`, `distTo3primeEndOfUTR = AIR_end − siteEnd` (AIR_end = UTR length), `distToNearestEndOfUTR = min(...)`; transformed `log10(dist)` (0 when dist=0), then min-max scaled by the per-site-type Min_dist coeff.

### Len_3UTR (verbatim from `get_len3UTR_weighted_contribution`)

Single-isoform: `utrLength = AIR_end` (UTR length); `log10(utrLength)`, then min-max scaled by the per-site-type Len_3UTR coeff.

### Off6m (verbatim from `getOffset6merSites` + `getOffset6mer_weighted_contribution`)

`siteToFind = substr(reverse_complement(seedRegion), 0, 6)` where `seedRegion = miRNA nt 2–8`; Off6m = number of (case-insensitive) occurrences of that 6mer in the UTR (single-isoform). Used **raw** (Off6m is not in the min-max regex): contribution = `coeff × count`.

### Caller-supplied (data-blocked) features — SPS / TA_3UTR / Len_ORF / ORF8m

These cannot be derived from a single 3'UTR sequence: **SPS** is a lookup in the Garcia et al. (2011) `TA_SPS_by_seed_region.txt` table (4096 seed regions; not a formula), **TA_3UTR** is a transcriptome-wide site count, **Len_ORF** needs the transcript ORF length, **ORF8m** needs a precomputed ORF 8mer count. The implementation accepts each as an OPTIONAL caller-supplied value and applies the exact `getAgarwalContribution` math (SPS/TA/Len_ORF min-max scaled — Len_ORF after `log10` — and ORF8m used raw); when not supplied they stay residual.

### SA — structural accessibility (now computed via the Turner-2004 McCaskill partition function)

**Definition (verbatim, retrieved this session).** `getSA_contribution($transcriptID,$speciesID,$utrStart,$siteType)` in `targetscan_70_context_scores.pl`:

```perl
if ($siteType == 1 || $siteType == 5) { $utrStart--; }
$plfold = (split /\t/, `grep -A7 -P '^$utrStart\t' $RNAplfold_outfile | tail -1 | cut -f 2-`)[13];
...
if (isNonzeroNumber($plfold)) { $log10_plfold = log($plfold) / log(10); } else { $log10_plfold = 0; }
$SA_contribution = getAgarwalContribution($siteType, "SA", $log10_plfold);
```

- Site-type integers (verbatim): `7mer-1a → 1`, `7mer-m8 → 2`, `8mer-1a → 3`, `6mer → 4`. So `$utrStart--` applies to **7mer-A1 only** (the `|| == 5` branch is dead legacy code; type 5 does not exist).
- `$utrStart` = the **1-based 5' start of the seed match** in the UTR. The `_lunp` row read is **7 nt downstream** (`grep -A7 … tail -1`), and `cut -f 2-` drops the position column so **index 13 = the 14th unpaired-probability column = L = 14**: the probability that the **14-nt stretch ending at that row is unpaired** (RNAplfold man page; the `_lunp` row `i`, column `L`, is "the probability that `[i−L+1..i]` is unpaired").
- The `_lunp` file is generated by `runRNAplfold_all_UTRs`: **`RNAplfold -L 40 -W 80 -u 20`** (window `W=80`, max base-pair span `L=40`, unpaired lengths up to 20).
- Then `log10(plfold)` (0 if not a nonzero number), passed to `getAgarwalContribution(type,"SA",log10_plfold)` → min-max scaled `(log10_plfold − min)/(max − min)` × `coeff` (SA row above).
- **Agarwal et al. (2015) eLife Fig 4A (verbatim):** SA is "scored as proportional to the **log10 value of the unpaired probability for a 14-nt region centered on the match to miRNA nucleotides 7 and 8**." This matches `row = utrStart+7`, `L = 14`.

**McCaskill partition-function sources (retrieved this session).**

| Source | How retrieved | Key fact |
|--------|---------------|----------|
| McCaskill JS (1990) Biopolymers 29:1105–1119 | WebSearch (PMID 1695107; ViennaRNA refman cites it) | Partition function `Z = Σ_P exp(−E(P)/RT)` over nested structures; base-pair probability `P^bp_{i,j} = Σ_{P∋(i,j)} exp(−E(P)/RT) / Z`; O(N³)/O(N²) inside/outside DP. |
| Lorenz et al. (2011) ViennaRNA Package 2.0, Algorithms Mol Biol 6:26 | WebSearch + ViennaRNA refman `pf_fold.html` (WebFetch 2026-06-25) | Boltzmann ensemble `Z = Σ_{s∈Ω} e^{−βE(s)}`, `p(s) = e^{−βE(s)}/Z`, β = 1/(kT), k ≈ 1.987e−3 kcal/(mol·K); per-base unpaired probability `p_unpaired(i) = 1 − Σ_j p_ij`; partition-function methods consider the entire equilibrium ensemble (ensemble FE `−RT ln Z` ≤ MFE). |
| RNAplfold man page (`tbi.univie.ac.at/RNA/RNAplfold.1.html`) | WebFetch (2026-06-25) | `-u` "mean probability that regions of length 1 to a given length are unpaired … saved in a `_lunp` file"; row `i`, column `L` = P(`[i−L+1..i]` unpaired). `-W` window-averaging size, `-L` max base-pair span. |

**Implementation.** The library now computes the 14-nt window unpaired probability EXACTLY from the **Turner-2004 McCaskill partition function** (`RnaSecondaryStructure.CalculateRegionUnpairedProbability`): `P(window unpaired) = Z_open / Z`, where `Z_open` is the partition function over structures in which every base of the window is unpaired (computed by forbidding any pair incident to the window — exactly Σ over window-unpaired structures). The fold is performed over a local context of up to `W = 80` nt centred on the window (mirroring RNAplfold's local-folding intent; a base can pair only within ±`L = 40` nt). The result is `log10`-transformed and min-max scaled by the SA row coefficients exactly as `getSA_contribution`/`getAgarwalContribution` do. SA is reported in `OmittedFeatures` only when the 14-nt window cannot be placed within the 3'UTR (matching the perl's "missing `_lunp` row → plfold = 0 → return 0"). **PCT** (probability of conserved targeting) remains an honest residual (needs a multiple-species alignment / branch length).

**Analytic verification of the partition function (the solid pin).** For `GAAAC` (the only admissible pair is G(0)–C(4), closing a 3-nt hairpin "AAA"), the Turner hairpin energy is `CalculateHairpinLoopEnergy("AAA",'G','C',false) = 5.4` kcal/mol (G–C, no terminal-AU penalty). With `RT = 1.987·310.15/1000 = 0.61626805` kcal/mol and `w = exp(−5.4/RT)`:
- `Z = 1 + w = 1.0001565052764922` (engine: exact match).
- `P(0,4) = w/Z = 0.00015648078642340854` (engine: exact match).
- `p_unpaired(0) = p_unpaired(4) = 1 − P(0,4) = 0.9998435192135765`; `p_unpaired(2) = 1` (interior loop base, never pairs).
- Ensemble FE `= −RT·ln Z = −9.64416549414892e−05` ≤ MFE = 0.
- `CalculateRegionUnpairedProbability("GAAAC", windowEnd=4, length=5) = 1/Z = 0.9998435192135765` (only the open chain leaves the whole window unpaired, `Z_open = 1`).

**Worked SA contribution (hand-derived, used as a test expectation).** For an 8mer let-7a site at `Start=20` in the 48-nt UTR `GGGGCCCCGGGGCCCCGGGG`+`CUACCUCA`+`GGGGCCCCGGGGCCCCGGGG` (8mer ⇒ no `utrStart--`), `windowEnd0 = Start+7 = 27`, window `[14..27]`. With the SA coeff/min/max for 8mer (`−0.115 / −4.356 / −0.661`): `SaContribution = −0.115·((log10 pu) − (−4.356))/((−0.661) − (−4.356))`, which the engine reproduces from `pu = CalculateRegionUnpairedProbability(localContext, …, 14)` to `1e−12` (`SaContribution == coeff·(log10 pu − min)/(max − min)`).

### Local_AU computation (verbatim from `getLocalAU_contribution`)

- `utrUp` = up to 30 nt immediately 5' of the site; `utrDown` = up to 30 nt immediately 3' of the site.
- Upstream position weight: `1/(i+1)` for 8mer & 7mer-m8, else `1/(i+2)` (i=0 at the base adjacent to the site).
- Downstream position weight: `1/(i+2)` for 8mer & 7mer-A1, else `1/(i+1)`.
- For each position, A or U contributes its weight to the score; every position contributes its weight to the denominator. `fraction = score/denominator`; then min-max scaled × Local_AU coeff.
- (The perl rounds intermediate contributions to 3 decimals for output; the implementation keeps full precision — a more-accurate, documented choice.)

### sRNA1 / sRNA8 / Site8 (verbatim)

- `sRNA1_nt = miRNA[1]`, `sRNA8_nt = miRNA[8]` (1-based). sRNA1 contributions are scored **only when nt1 ≠ U**; sRNA8 only when nt8 ≠ U.
- Site8 contributions are scored **only when the target site-position-8 base ≠ U**, and only for 7mer-A1 / 6mer; the relevant base is the nucleotide opposite miRNA position 8 (in this repo's geometry: the base immediately 5' of the 6mer core, `mRNA[site.Start − 1]`).

### Worked examples (hand-derived; used as exact test expectations)

| Case | miRNA / layout | Realised contributions | Partial CS |
|------|----------------|------------------------|-----------|
| 8mer | let-7a (nt1=U,nt8=G); `GGGGG+CUACCUCA+GGGGG`; all-G flanks ⇒ fraction 0 | Int=-0.589; Local_AU=-0.254×((0-0.308)/(0.814-0.308))=+0.154608695652174; sRNA8G=+0.015 | **-0.419391304347826** |
| 7mer-m8 | let-7a; `GGGGG+CUACCUC+GGGG`; fraction 0 | Int=-0.224; Local_AU=-0.177×((0-0.277)/(0.782-0.277))=+0.097087128712871; sRNA8G=-0.008 | **-0.134912871287129** |
| 7mer-A1 | let-7a; `GGGGG+UACCUC+A+GGGG`; Site8 base 'G'; fraction 0 | Int=-0.195; Local_AU=-0.075×((0-0.342)/(0.801-0.342))=+0.055882352941176; sRNA8G=-0.017; Site8G=+0.015 | **-0.141117647058824** |
| 6mer | miR-21 (nt1=U,nt8=U); `GGGGC+UAAGCU+UGAGG`; Site8 base 'C'; fraction 0.357142857142857 | Int=-0.079; Local_AU=-0.040×((0.357142857142857-0.295)/(0.772-0.295))=-0.005211141060198; Site8C=+0.015 | **-0.069211141060198** |

## PCT — probability of conserved targeting (Friedman et al. 2009)

### Online Sources (retrieved verbatim this session)

- **Friedman RC, Farh KKH, Burge CB, Bartel DP. 2009. "Most mammalian mRNAs are conserved targets of microRNAs." Genome Research 19(1):92-105. doi:10.1101/gr.082701.108.** Retrieved via `WebFetch` of `https://pmc.ncbi.nlm.nih.gov/articles/PMC2612969/`. Methods (verbatim): *"The conservation of a given sequence (e.g., an 8mer miR-1 site in a particular 3′UTR) was then assessed by summing the total branch length in the phylogenetic tree connecting the subset of species having the sequence perfectly aligned…"* and *"the P_CT was defined as E[(S − B)/S] where B, the background estimate, is a constant, and S is a random variable."*
- **TargetScan PCT documentation.** Retrieved via `WebFetch` of `https://www.targetscan.org/docs/pct.html`: PCT *"is approximately equal to (S/B − 1)/(S/B) (or near zero, for sites with S/B < 1)"*; aggregate PCT `= 1 − Π(1 − PCT_site_i)`.
- **`targetscan_70_BL_PCT.pl`** (`calculatePCTthisBL`). Retrieved via `WebFetch` of `https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_BL_PCT.pl`. Verbatim: `my $pct = $b0 + ( $b1 / (1 + $eConstant ** ( (0 - $b2) * $BL + $b3)));` with `$eConstant = 2.71828182845904523536;` and `if ($pct < 0) { $pct = "0.0"; }`. Branch length ≤ 0 ⇒ PCT 0; the four coefficients are retrieved **per miRNA family and site type** from the compiled `8mer_PCT_parameters.txt` / `7mer_m8_PCT_parameters.txt` / `7mer_1a_PCT_parameters.txt` (column 1 = miRNA family, columns 2–5 = `b0,b1,b2,b3`).
- **`targetscan_70_context_scores.pl`** (`getPCT_contribution` → `getAgarwalContribution`). Retrieved this session: PCT matches the min-max-scaling regex branch `/^TA_3UTR$|^SPS$|^Local_AU$|^3P_score$|^SA$|^Len_ORF$|^Len_3UTR$|^Min_dist$|^PCT$/`, so its context++ contribution is `coeff × (pct − min)/(max − min)` exactly like the other scaled features.

### Bls definition used in the implementation

**Bls = total branch length of the minimal subtree of the phylogenetic tree connecting the species in which the seed match is perfectly aligned** (Friedman 2009). An edge belongs to that connecting subtree iff at least one conserved species lies on each side of the edge. A single conserved species (or none) ⇒ Bls = 0 (consistent with the perl's `branchLength == 0 ⇒ PCT 0`).

### PCT formula (parameterised; parameters caller-supplied)

`PCT(Bls) = b0 + b1 / (1 + e^(−b2·Bls + b3))`, truncated at 0. **Licence/sourcing note:** Friedman 2009 publishes the Bls definition and `PCT = E[(S−B)/S]` but does **not** publish the fitted `b0..b3` as numbers; those live in TargetScan's compiled, citation-required `*_PCT_parameters.txt` tables (per miRNA family). Per the stop-rule, the `b0..b3` are therefore **caller-supplied** (`PctSigmoidParameters`) — not bundled, not invented. Only the published **equation** and the Agarwal **PCT coefficient** (a permissive scientific-fact constant, like the other context++ coefficients) are bundled.

### Worked examples (hand-derived; used as exact test expectations)

Worked tree (Newick): `((A:1.0,B:2.0):0.5,(C:1.5,D:3.0):4.0);`

| Conserved species | Bls (hand-derived) | Reasoning |
|-------------------|--------------------|-----------|
| {A, B} | **3.0** | A(1.0)+B(2.0); the (A,B)→root edge (0.5) is not counted (no conserved species outside the {A,B} subtree) |
| {A, C} | **7.0** | A(1.0)+(A,B)→root(0.5)+(C,D)→root(4.0)+C(1.5) |
| {A, B, C, D} | **12.0** | all four leaf edges (1+2+1.5+3) + both internal edges (0.5+4.0) |
| {A} | **0.0** | a single conserved species — no connecting subtree |

PCT with the simple sigmoid `b0=0, b1=1, b2=1, b3=0` (so PCT = 1/(1+e^(−Bls))):
- `PCT(3.0) = 1/(1+e^-3) = 0.952574126822433`.
- `PCT(7.0) = 1/(1+e^-7) = 0.999088948805599`.
- Negative-truncation case `b0=-0.5, b1=0.3, b2=1, b3=5`, `Bls=0`: raw `= -0.5 + 0.3/(1+e^5) ≈ -0.498 < 0 ⇒ 0`.

context++ PCT contribution (verbatim PCT row of `Agarwal_2015_parameters.txt`, min 0):
- 8mer (coeff -0.103, max 0.816), PCT(3.0): `-0.103 × (0.952574126822433/0.816) = -0.120239136106263`.
- 7mer-m8 (coeff -0.048, max 0.364), PCT(7.0): `-0.048 × (0.999088948805599/0.364) = -0.131747993249090`.
