# Target Site Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | MiRNA |
| Test Unit ID | MIRNA-TARGET-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-24 |

## 1. Overview

Target-site prediction identifies mRNA segments that can pair to the canonical miRNA seed and classifies them into the standard site hierarchy used in animal miRNA biology [1][2][3][5]. The repository implements exact seed-based scanning in `MiRnaAnalyzer.FindTargetSites`, then builds a full miRNA-target duplex alignment and a heuristic score for each retained hit. Canonical 8mer, 7mer-m8, 7mer-A1, 6mer, and offset-6mer sites are supported, and higher-priority seed classes suppress overlapping offset-6mer calls at the same positions. The implementation is deterministic and useful for exact seed-complement discovery, but it does not model conservation, site accessibility, or noncanonical site discovery in the main finder [1][3][4].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Animal miRNAs typically recognize targets through antiparallel pairing between the miRNA seed and sites in the mRNA, especially in 3' UTRs [1][2]. Canonical site classes are ordered by efficacy, with 8mer sites generally strongest, followed by 7mer-m8, 7mer-A1, and 6mer classes [1][2][3]. Offset 6mer sites, based on miRNA positions 3-8 rather than 2-7, are weaker but still used in some computational frameworks [4][5].

### 2.2 Core Model

Let the canonical seed of a miRNA be $m_2\ldots m_8$. The finder first constructs the reverse complement of that 7-mer seed:

$$
\operatorname{seedRC} = \operatorname{revcomp}(m_2\ldots m_8)
$$

and then derives two exact seed-matching patterns:

$$
\mathrm{6mer\ core} = \operatorname{seedRC}[2\ldots 7], \qquad
\mathrm{offset\ 6mer} = \operatorname{seedRC}[1\ldots 6]
$$

using biological positions in the underlying miRNA [1][2]. Site classes follow the standard seed taxonomy [1][2][3]:

| Site Type | Pairing Rule |
|-----------|--------------|
| 8mer | miRNA positions 2-8 match, plus an `A` opposite miRNA position 1 |
| 7mer-m8 | miRNA positions 2-8 match |
| 7mer-A1 | miRNA positions 2-7 match, plus an `A` opposite miRNA position 1 |
| 6mer | miRNA positions 2-7 match |
| Offset 6mer | miRNA positions 3-8 match |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported site has `0 <= Score <= 1`. | `CalculateTargetScore` clamps the final score into `[0, 1]`. |
| INV-02 | `SeedMatchLength` is `8`, `7`, or `6` according to the emitted site class. | `FindTargetSites` assigns the seed-match length when classifying each site. |
| INV-03 | `Start` and `End` are zero-based inclusive target coordinates. | `CreateTargetSite` sets `End = Start + length - 1`. |
| INV-04 | Pass-1 canonical seed sites take priority over overlapping offset-6mer sites. | Covered positions from pass 1 are tracked and used to suppress pass-2 offset-6mer calls. |
| INV-05 | Duplex alignment is antiparallel. | `AlignMiRnaToTarget` reads the target string in reverse when comparing it to the miRNA. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `mRnaSequence` | `string` | required | RNA or DNA sequence to scan for seed-complement sites. | `null` or empty input yields no results; `T` is normalized to `U`. |
| `miRna` | `MiRna` | required | miRNA record containing the normalized sequence and seed. | Empty stored sequence yields no results. |
| `minScore` | `double` | `0.5` | Minimum heuristic score required for a site to be emitted. | Not explicitly validated or clamped before filtering. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | Zero-based inclusive start coordinate of the site on the target sequence. |
| `End` | `int` | Zero-based inclusive end coordinate of the canonical site window. |
| `TargetSequence` | `string` | Target subsequence extended up to the miRNA length or the remaining mRNA tail. |
| `MiRnaName` | `string` | Name copied from the supplied `MiRna` record. |
| `Type` | `TargetSiteType` | Site class: `Seed8mer`, `Seed7merM8`, `Seed7merA1`, `Seed6mer`, or `Offset6mer` in current emitted results. |
| `SeedMatchLength` | `int` | Nominal number of matched seed positions used for the site class. |
| `Score` | `double` | Heuristic site score in `[0, 1]`. |
| `FreeEnergy` | `double` | Heuristic duplex energy derived from alignment features. |
| `Alignment` | `string` | Alignment string using `|` for Watson-Crick matches, `:` for G:U wobble pairs, and space for mismatches. |

### 3.3 Preconditions and Validation

`FindTargetSites` is non-throwing for empty sequence inputs: if `mRnaSequence` is `null` or empty, or if `miRna.Sequence` is empty, the method yields no results. Both mRNA and miRNA sequences are normalized to uppercase RNA by replacing `T` with `U`. If the seed reverse complement is shorter than 7 nt, the method yields no sites. Coordinates are zero-based and inclusive. `minScore` is compared exactly as supplied, so values above `1.0` suppress all results and values below `0.0` admit every scored candidate.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the mRNA input to uppercase RNA and read the normalized miRNA sequence from the supplied `MiRna` record.
2. Compute the reverse complement of the miRNA seed, then derive the canonical 6mer core and the offset-6mer pattern.
3. Pass 1: scan the mRNA for exact 6mer-core matches and upgrade each hit to 8mer, 7mer-m8, 7mer-A1, or 6mer depending on the upstream position-8 base and downstream `A1` base.
4. For each pass-1 hit, build a full duplex alignment, compute its heuristic score, and retain the site if `Score >= minScore`, marking its covered coordinates.
5. Pass 2: scan for offset-6mer matches, skipping sites that overlap already retained pass-1 sites or that would become a full seed match at the same location.
6. Emit retained `TargetSite` records in scan order.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Current site classification is entirely rule-based:

- If the 6mer core matches and both the upstream position-8 base and downstream `A1` base are present, emit `Seed8mer`.
- If only the upstream position-8 base is present, emit `Seed7merM8`.
- If only the downstream `A1` base is present, emit `Seed7merA1`.
- Otherwise emit `Seed6mer`.
- Offset-6mer sites are discovered only in pass 2 and only when they do not overlap coordinates already covered by higher-priority pass-1 sites.

Alignment allows Watson-Crick pairs (`A-U`, `U-A`, `G-C`, `C-G`) and G:U wobble pairs (`G-U`, `U-G`).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Seed-based site scan | `O(n + hL)` | `O(k + L)` | `n` is mRNA length, `h` is the number of candidate sites that require full-duplex scoring, `L` is miRNA length, and `k` is the number of covered target coordinates retained from pass 1. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MiRnaAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs)

- `MiRnaAnalyzer.FindTargetSites(string, MiRna, double)`: Canonical seed-based target finder.
- `MiRnaAnalyzer.GetReverseComplement(string)`: Builds the RNA reverse complement used in seed matching.
- `MiRnaAnalyzer.AlignMiRnaToTarget(string, string)`: Generates the antiparallel duplex summary.
- `MiRnaAnalyzer.AnalyzeTargetContext(string, int, int, int)`: Separate context annotation helper, not part of site discovery.
- `MiRnaAnalyzer.CalculateSiteAccessibility(string, int, int)`: Separate accessibility heuristic, not part of site discovery.
- `MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(string, MiRna, TargetSite, ContextPlusPlusInputs)`: Opt-in TargetScan context++ scorer (Agarwal et al. 2015 [4]); realises the miRNA- and 3'UTR-computable feature subset (site-type intercept, local-AU, sRNA position-1/8 identity, target site-position-8 identity, **3' supplementary pairing (3P_score), minimum distance to the 3'UTR end (Min_dist), 3'UTR length (Len_3UTR), offset-6mer count (Off6m)**) and, when the caller supplies them via the optional `ContextPlusPlusInputs`, the data-blocked features **SPS, TA_3UTR, Len_ORF, ORF8m**. Returns a `ContextPlusPlusScore` breakdown with the still-residual features listed. Does not replace the default `Score`.

### 5.2 Current Behavior

The implementation performs target discovery in two passes, not through a single generic seed-matching helper. Pass 1 discovers canonical 6mer-core-derived sites and upgrades them to `Seed8mer`, `Seed7merM8`, `Seed7merA1`, or `Seed6mer`; pass 2 discovers `Offset6mer` sites only when they do not overlap already retained higher-priority sites. `CreateTargetSite` extends `TargetSequence` to the shorter of the miRNA length or the remaining mRNA tail, so the stored target segment can be longer than the canonical seed window. Duplex energy is heuristic rather than thermodynamic: stacked Watson-Crick matches contribute `-2.0`, isolated Watson-Crick matches `-1.0`, G:U wobble pairs `-0.5`, and mismatches `+0.5`, but the current loop bounds omit the final alignment position from that running total. Base site scores are currently `1.0` (8mer), `0.52` (7mer-m8), `0.32` (7mer-A1), `0.15` (6mer), `0.10` (offset 6mer), with a `+0.05` bonus when the duplex has more than 10 matches and a `-0.01` penalty per mismatch before clamping to `[0, 1]` [3][4]. Although the `TargetSiteType` enum also contains `Supplementary` and `Centered`, `FindTargetSites` does not emit those classes.

Separately, `ScoreTargetSiteContextPlusPlus` provides an **opt-in** TargetScan context++ score (Agarwal et al. 2015 [4]) for a discovered seed-match site, leaving the default `Score` unchanged. context++ is a multiple-linear-regression model with a distinct coefficient set **per site type**, where the score is the sum of per-feature contributions (`contribution = coeff × scaledScore`). Per `getAgarwalContribution`, the continuous features `TA_3UTR, SPS, Local_AU, 3P_score, SA, Len_ORF, Len_3UTR, Min_dist, PCT` are min-max scaled `(raw − min)/(max − min)` (NOT clamped to `[0,1]`), while the nucleotide-identity indicators and the `Off6m`/`ORF8m` counts are used raw [7]. This implementation realises every feature derivable from the miRNA and the supplied 3'UTR:

- **Intercept**, **Local_AU** (position-weighted A/U fraction of the 30 nt up/downstream, scaled), the **sRNA1{A,C,G}** / **sRNA8{A,C,G}** miRNA position-1/8 identity indicators (skipped when that nucleotide is U), and the **Site8{A,C,G}** target site-position-8 indicators (defined only for 7mer-A1 and 6mer).
- **3P_score** — 3' supplementary pairing, a faithful port of `get3primePairingContribution` (fed by `extractSubseqForAlignment` + `modifySubseqForAlignment`): the maximum, over single-gap offsets in both strands, of summed run-of-≥2 base pairs (1.0 in the offset-adjusted window 4..7, else 0.5) minus `max(0,(offset−2)/2)`; then min-max scaled (min=1, max=3.5). Raw scores were cross-checked against the reference perl.
- **Min_dist** — `log10(min(distTo5′, distTo3′))` of the site to the nearest 3'UTR end (single-isoform: the UTR length), min-max scaled (`getMinDist_weighted_contribution`).
- **Len_3UTR** — `log10(3'UTR length)`, min-max scaled (`get_len3UTR_weighted_contribution`).
- **Off6m** — count of the offset-6mer pattern (first 6 nt of `reverse_complement(seed 2-8)`) in the 3'UTR, used raw (`getOffset6merSites` + `getOffset6mer_weighted_contribution`).

The features the library cannot derive from sequence alone — **SPS** (Garcia et al. 2011 seed-region table value), **TA_3UTR** (transcriptome-wide site abundance), **Len_ORF** and **ORF8m** (require the transcript ORF) — are computed faithfully **only when the caller supplies them** via `ContextPlusPlusInputs`; otherwise they contribute 0 and are listed in `OmittedFeatures`. **SA** (RNAplfold partition-function unpaired probability of a 14-nt window) and **PCT** (probability of conserved targeting) are always honest residuals (a partition-function accessibility / multi-species alignment is out of scope). All coefficients are the verbatim values from `Agarwal_2015_parameters.txt` [6] and every feature/scaling formula is ported from `targetscan_70_context_scores.pl` [7]. The returned `ContextScorePartial` is therefore a partial context++ score, not the published headline CS.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Canonical 8mer, 7mer-m8, 7mer-A1, 6mer, and offset-6mer seed classes [1][2][3][5].
- Antiparallel seed reverse-complement matching between miRNA and target [1][2].
- Priority ordering that prefers stronger canonical seed classes over weaker offset-6mer interpretations at overlapping positions [1][3].
- (Opt-in `ScoreTargetSiteContextPlusPlus`) The context++ features computable from the miRNA and the supplied 3'UTR — site-type Intercept, Local_AU, sRNA position-1/8 and target site-position-8 nucleotide-identity indicators, **3' supplementary pairing (3P_score), minimum distance to the 3'UTR end (Min_dist), 3'UTR length (Len_3UTR), and offset-6mer count (Off6m)** — plus the data-blocked features **SPS, TA_3UTR, Len_ORF, ORF8m** when supplied by the caller; all with the **verbatim fitted coefficients** from `Agarwal_2015_parameters.txt` [6] and the per-site-type min-max scaling / DP / indicator logic ported exactly from `targetscan_70_context_scores.pl` [7] (3P_score raw values cross-checked against the reference perl).

**Intentionally simplified:**

- The default `Score` is a heuristic normalization of site class, mismatch burden, and extra pairing rather than a fitted repression model; **consequence:** `Score` ranks hits within this implementation but is not a calibrated repression probability. (The opt-in `ScoreTargetSiteContextPlusPlus` provides the source-fitted context++ alternative.)
- (Opt-in context++) The score still omits the genuinely data-/method-blocked features: **SA** (RNAplfold partition-function unpaired probability — the library has only MFE folding, so SA is left as an honest residual rather than approximated by MFE), **PCT** (requires a multi-species alignment / branch length), and **SPS / TA_3UTR / Len_ORF / ORF8m** unless the caller supplies them via `ContextPlusPlusInputs`; **consequence:** `ContextScorePartial` is a partial context++ score, reported alongside the explicit `OmittedFeatures` list of whatever remains residual for that call.
- Duplex free energy is derived from simple match, wobble, and mismatch rules over an approximate alignment scan; **consequence:** the reported `FreeEnergy` is useful comparatively inside the repository but is not a nearest-neighbor folding energy and should not be interpreted as a full per-position thermodynamic score.
- Context and accessibility are exposed only as separate helper methods; **consequence:** the main finder does not integrate AU context, conservation, or accessibility into site discovery or ranking.

**Not implemented:**

- Discovery of centered sites, supplementary-only sites, and bulged seed sites; **users should rely on:** `AnalyzeTargetContext` and `CalculateSiteAccessibility` only for separate annotations, or use an external miRNA target-prediction tool for broader site classes.
- The full context++ headline score: the still-residual features (`SA`, `PCT`, and — unless supplied — `SPS`, `TA_3UTR`, `Len_ORF`, `ORF8m`) require RNA partition-function structure, multi-species conservation, the Garcia et al. 2011 SPS table, a transcriptome, and the transcript ORF — none derivable from a single 3'UTR sequence; **users should rely on:** the official TargetScan distribution for the full CS, supply the data-blocked features via `ContextPlusPlusInputs`, or treat `ContextScorePartial` as a partial score.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `TargetSiteType.Centered` and `TargetSiteType.Supplementary` exist in the enum but are not emitted by `FindTargetSites`. | Deviation | The public type system advertises more site classes than the current finder can actually discover. | accepted | Discovery remains limited to canonical seed-driven classes and offset 6mers. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or `null` mRNA | Returns no sites. | The finder exits early for empty input. |
| Empty `miRna.Sequence` | Returns no sites. | The finder exits before seed matching if the miRNA sequence is empty. |
| Seed shorter than 7 nt | Returns no sites. | Canonical seed reverse-complement matching requires the 7-mer seed. |
| DNA target sequence | `T` is treated as `U`. | mRNA input is normalized with `Replace('T', 'U')`. |
| `minScore > 1.0` | Returns no sites. | All heuristic scores are clamped to `[0, 1]` before comparison. |
| `minScore < 0.0` | All candidate sites that survive class-specific filtering are emitted. | Scores are compared directly against the caller-supplied threshold. |

### 6.2 Limitations

The finder reports only exact canonical seed-complement sites and offset-6mer sites. It does not model bulges, gaps, centered sites, supplementary-only pairing, evolutionary conservation, transcript isoforms, or experimentally fitted repression models. Its heuristic `FreeEnergy` and `Score` are implementation-level ranking aids rather than biophysical or transcriptome-scale confidence measures.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
using Seqeron.Genomics.Annotation;

var let7a = MiRnaAnalyzer.CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU");
string mrna = "GGGGGCUACCUCAGGGGG";

var sites = MiRnaAnalyzer.FindTargetSites(mrna, let7a, minScore: 0.1).ToList();

// sites[0].Type == TargetSiteType.Seed8mer
// sites[0].SeedMatchLength == 8

// Opt-in TargetScan context++ score (Agarwal et al. 2015) for the discovered site:
var ctx = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna, let7a, sites[0]);
// ctx.ContextScorePartial    -> partial CS (Intercept + Local_AU + sRNA1/8 + Site8
//                               + 3P_score + Min_dist + Len_3UTR + Off6m)
// ctx.OmittedFeatures        -> still-residual features (SA, PCT, + SPS/TA/Len_ORF/ORF8m if not supplied)

// Supply the data-blocked features to fold them in (else they stay residual):
var withData = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(
    mrna, let7a, sites[0],
    new MiRnaAnalyzer.ContextPlusPlusInputs(Sps: -8.0, Ta: 3.5, OrfLength: 1000, Orf8mCount: 2));
```

**Numerical walk-through (8mer, let-7a, all-G flanks):** with `mRNA = GGGGG + CUACCUCA + GGGGG` (length 18), the 8mer occupies `Start=5,End=12`; both 30 nt flanks are all-G so the local-AU fraction is 0. The realised contributions are: `Intercept(8mer)=-0.589`; `Local_AU = -0.254 × ((0 − 0.308)/(0.814 − 0.308)) = +0.154608695652174`; `sRNA8G(8mer)=+0.015` (nt1=U ⇒ sRNA1 unscored; Site8 undefined for 8mer); `3P_score`: raw=0 ⇒ `-0.040 × ((0 − 1)/(3.5 − 1)) = +0.016`; `Min_dist`: nearest end = min(5, 18−13) = 5 ⇒ `0.118 × ((log10 5 − 1.415)/(3.113 − 1.415)) = -0.049759446106213065`; `Len_3UTR`: `0.310 × ((log10 18 − 2.392)/(3.637 − 2.392)) = -0.2830405810586145`; `Off6m`: pattern `CUACCU` occurs once ⇒ `-0.020 × 1 = -0.020`. Sum = **-0.7561913315126536**. All numbers are the verbatim Agarwal_2015_parameters.txt coefficients.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MiRnaAnalyzer_TargetPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_TargetPrediction_Tests.cs)
- Test spec: [MIRNA-TARGET-001.md](../../../tests/TestSpecs/MIRNA-TARGET-001.md)
- Evidence: [MIRNA-TARGET-001-Evidence.md](../../../docs/Evidence/MIRNA-TARGET-001-Evidence.md)
- Related property tests: [MiRnaProperties.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Properties/MiRnaProperties.cs)
- Related snapshots: [MiRnaSnapshotTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Snapshots/MiRnaSnapshotTests.cs)

## 8. References

1. Bartel DP. 2009. MicroRNAs: target recognition and regulatory functions. Cell. 136(2):215-233.
2. Lewis BP, Burge CB, Bartel DP. 2005. Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. Cell. 120(1):15-20.
3. Grimson A, Farh KKH, Johnston WK, Garrett-Engele P, Lim LP, Bartel DP. 2007. MicroRNA targeting specificity in mammals: determinants beyond seed pairing. Molecular Cell. 27(1):91-105.
4. Agarwal V, Bell GW, Nam JW, Bartel DP. 2015. Predicting effective microRNA target sites in mammalian mRNAs. eLife. 4:e05005. doi:10.7554/eLife.05005. https://pmc.ncbi.nlm.nih.gov/articles/PMC4532895/
5. TargetScan Human 8.0. 2021. Canonical site classes and ranking for animal miRNA target prediction. https://www.targetscan.org/
6. TargetScan distribution. `Agarwal_2015_parameters.txt` (context++ fitted coefficients per site type). https://raw.githubusercontent.com/nsoranzo/targetscan/main/Agarwal_2015_parameters.txt
7. TargetScan distribution. `targetscan_70_context_scores.pl` (context++ feature computation and scaling: `getAgarwalContribution`, `getLocalAU_contribution`, `get_sRNA1_8_contributions`, `getSite8_contribution`). https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_context_scores.pl
