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
| `targetscan_70_context_scores.pl` (TargetScan distribution) | `curl https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_context_scores.pl` (2026-06-24) | `getAgarwalContribution`, `getLocalAU_contribution`, `get_sRNA1_8_contributions`, `getSite8_contribution`, `readAgarwalParameters`, CS sum (lines ~404–411). Defines feature computation + scaling. |
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

Site8 features are computed by the perl **only for 7mer-A1 (siteType 1) and 6mer (siteType 4)**; the parameter file leaves the 8mer/7mer-m8 cells blank.

### Locally-computable subset realised (and the honest residual)

Only features depending solely on the **miRNA** and the **local 3'UTR** are realised: **Intercept, Local_AU, sRNA1{A,C,G}, sRNA8{A,C,G}, Site8{A,C,G}**. The remaining full-transcript features are **NOT computed** and are reported as the residual: `3P_score` (3' supplementary pairing), `SPS`, `TA_3UTR`, `Min_dist`, `SA` (structural accessibility), `PCT`, `Len_3UTR`, `Len_ORF`, `ORF8m`, `Off6m`. Because most omitted coefficients are negative, the partial CS is an **upper bound** on the full context++ score.

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
