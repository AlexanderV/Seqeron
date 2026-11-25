# Target Site Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | MiRNA |
| Test Unit ID | MIRNA-TARGET-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-09-28 |

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
| `minScore` | `double` | `0.0` | Minimum predicted repression (`Score`) for a site to be emitted; the default returns every site. | Not validated; negative values admit every site. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | Zero-based inclusive start coordinate of the site on the target sequence. |
| `End` | `int` | Zero-based inclusive end coordinate of the canonical site window. |
| `TargetSequence` | `string` | mRNA segment antiparallel to the full miRNA: `mrna[max(0, p1 − L + 1) .. p1]`, where `p1` is the base opposite miRNA nt 1 (just 3' of the 6mer core; 3' of the offset-6mer match + 1). The miRNA 3' end pairs **upstream** of the seed match [1][7]. Truncated at the mRNA ends. |
| `MiRnaName` | `string` | Name copied from the supplied `MiRna` record. |
| `Type` | `TargetSiteType` | Site class: `Seed8mer`, `Seed7merM8`, `Seed7merA1`, `Seed6mer`, or `Offset6mer` in current emitted results. |
| `SeedMatchLength` | `int` | Nominal number of matched seed positions used for the site class. |
| `Score` | `double` | Predicted repression `max(0, −ContextScore)` in log2 units (TargetScan 5.0 / Grimson 2007 context score; 0 for 6mer, offset 6mer and sites starting within the first 15 nt). The signed score is `TargetSite.ContextScore` (`double?`, `null` when undefined). |
| `FreeEnergy` | `double` | Turner 2004 nearest-neighbour ΔG°37 of the ungapped duplex (`AlignMiRnaToTarget`, MIRNA-PAIR-001: initiation, stacks, internal loops, AU/GU ends, dangling ends / terminal mismatches). |
| `Alignment` | `string` | Per-miRNA-position string (index 0 = nt 1) using `|` for Watson-Crick matches, `:` for G:U wobble pairs, and space for mismatches / no partner. Seed positions of a canonical site are always `|`. |

### 3.3 Preconditions and Validation

`FindTargetSites` is non-throwing for empty sequence inputs: if `mRnaSequence` is `null` or empty, or if `miRna.Sequence` is empty, the method yields no results. Both mRNA and miRNA sequences are normalized to uppercase RNA by replacing `T` with `U`. If the seed reverse complement is shorter than 7 nt, the method yields no sites. Coordinates are zero-based and inclusive. `minScore` is compared exactly as supplied, so values above `1.0` suppress all results and values below `0.0` admit every scored candidate.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the mRNA input to uppercase RNA and read the normalized miRNA sequence from the supplied `MiRna` record.
2. Compute the reverse complement of the miRNA seed, then derive the canonical 6mer core and the offset-6mer pattern.
3. Pass 1: scan the mRNA for exact 6mer-core matches and upgrade each hit to 8mer, 7mer-m8, 7mer-A1, or 6mer depending on the upstream position-8 base and downstream `A1` base.
4. For each pass-1 hit, build a full duplex alignment (`TargetSequence` / `Alignment` / `FreeEnergy`), compute the Grimson (2007) / TargetScan 5.0 context score (`ScoreTargetSiteGrimson`: site type + 3' pairing + local AU + position), and retain the site if `Score >= minScore`. Every pass-1 hit marks its coordinates as covered whether or not it is retained, so the offset-6mer suppression does not depend on `minScore`.
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
- `MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(string, MiRna, TargetSite, ContextPlusPlusInputs)`: Opt-in TargetScan context++ scorer (Agarwal et al. 2015 [4]); realises the miRNA- and 3'UTR-computable feature subset (site-type intercept, local-AU, sRNA position-1/8 identity, target site-position-8 identity, **3' supplementary pairing (3P_score), minimum distance to the 3'UTR end (Min_dist), 3'UTR length (Len_3UTR), offset-6mer count (Off6m)**) **SPS and TA_3UTR from the bundled TargetScan 7.0 `TA_SPS_by_seed_region.txt`** (caller values override), **Len_ORF and ORF8m from `ContextPlusPlusInputs.OrfSequence`** (as `targetscan_count_8mers.pl`; or caller numbers), PCT from a conservation input, and SA from a supplied RNAplfold value or the local estimate. Returns a `ContextPlusPlusScore` breakdown with the still-residual features listed. Does not replace the default `Score`.
- `MiRnaAnalyzer.TryGetTargetScanTaSps(string seedRegion, TargetSiteType, out double ta, out double sps)` / `TargetScanTaSpsSeedRegionCount` (16 384): the bundled TargetScan 7.0 TA/SPS table, keyed by the 7-nt seed region (miRNA nt 2–8).
- `MiRnaAnalyzer.FindTargetScanSiteGroups(string seedOrMature, IReadOnlyDictionary<string,string> alignedUtrs, IReadOnlyCollection<string>? speciesWithMiRna)` (B13 F35): port of `targetscan_70.pl` — seed-match sites (8mer / 7mer-m8 / 7mer-A1 / 6mer, gap-tolerant) in every species of a multi-species 3'UTR alignment and their cross-species groups (alignment + ungapped UTR coordinates, group type, per-site-type species lists, ORF-overlap flag) [12]. `AssignTargetScanUtrBin` (`targetscan_70_BL_bins.pl`), `GetTargetScanBinTree(Newick)` / `GetTargetScanGenericTree` / `TargetScanUtrBinThresholds` (bundled TargetScan 7.0 trees, single-line Newick) and the end-to-end `ScoreTargetScanConservation` (alignment → groups → bin → Bls → PCT → Conserved, `targetscan_70_BL_PCT.pl`) (B13 F36). Implementation: `MiRnaAnalyzer.TargetScanConservation.cs`, trees in the generated `MiRnaAnalyzer.TargetScanTrees.cs`.
- `MiRnaAnalyzer.CountTargetScanOrf8merSites(string orfSequence, string seedRegion)`: ORF 8mer count as `targetscan_count_8mers.pl` (gaps removed, T→U, case-insensitive, non-overlapping `revcomp(seed)+A`).
- `MiRnaAnalyzer.ComputeTa3Utr(MiRna, IEnumerable<string>)` / `CountSeedSites3Utr(MiRna, IEnumerable<string>)`: Opt-in computation of the context++ feature **`TA_3UTR`** (transcriptome-wide target-site abundance) from a supplied 3'UTR set, per Garcia et al. (2011) [10] / Agarwal et al. (2015) [4]: `TA = log10(N)`, with `N` = total non-overlapping **8mer + 7mer-m8 + 7mer-A1** sites of the seed across the 3'UTRs (`CountSeedSites3Utr` exposes the raw `N`). The result feeds `ScoreTargetSiteContextPlusPlus` via `ContextPlusPlusInputs.Ta`. Converts TA from an opaque caller-supplied number into a value computed from a 3'UTR set per the published definition.

### 5.2 Current Behavior

The implementation performs target discovery in two passes, not through a single generic seed-matching helper. Pass 1 discovers canonical 6mer-core-derived sites and upgrades them to `Seed8mer`, `Seed7merM8`, `Seed7merA1`, or `Seed6mer`; pass 2 discovers `Offset6mer` sites only when they do not overlap already retained higher-priority sites. `CreateTargetSite` aligns the full miRNA in the register fixed by the seed match: miRNA nt *k* faces mRNA index `p1 − (k − 1)`, where `p1` is the base opposite nt 1, so the duplex window `TargetSequence = mrna[max(0, p1 − L + 1) .. p1]` lies **upstream** (5') of and including the seed match — the miRNA 3' end pairs 5' of the site (Bartel 2009 Fig. 1 [1]; TargetScan `extractSubseqForAlignment` takes `utrStart − 16 .. utrEnd` [7]). If `p1` is past the mRNA 3' end (a 6mer / offset 6mer flush with the end), nt 1 is reported unpaired. (Fixed 2026-09-28: the window was previously taken *downstream* from `Start`, mis-registering the seed — e.g. the let-7a 8mer worked example aligned nt 1 to a flanking G and reported the seed unpaired, making the mismatch penalty, `Alignment` and `FreeEnergy` meaningless.) Duplex `FreeEnergy` sums Turner-2004 stacking energies over consecutive paired positions (no loop terms). Base site scores are currently `1.0` (8mer), `0.52` (7mer-m8), `0.32` (7mer-A1), `0.15` (6mer), `0.10` (offset 6mer), with a `+0.05` bonus when the duplex has more than 10 matches and a `-0.01` penalty per mismatch before clamping to `[0, 1]` [3][4]. Although the `TargetSiteType` enum also contains `Supplementary` and `Centered`, `FindTargetSites` does not emit those classes.

Separately, `ScoreTargetSiteContextPlusPlus` provides an **opt-in** TargetScan context++ score (Agarwal et al. 2015 [4]) for a discovered seed-match site, leaving the default `Score` unchanged. context++ is a multiple-linear-regression model with a distinct coefficient set **per site type**, where the score is the sum of per-feature contributions (`contribution = coeff × scaledScore`). Per `getAgarwalContribution`, the continuous features `TA_3UTR, SPS, Local_AU, 3P_score, SA, Len_ORF, Len_3UTR, Min_dist, PCT` are min-max scaled `(raw − min)/(max − min)` (NOT clamped to `[0,1]`), while the nucleotide-identity indicators and the `Off6m`/`ORF8m` counts are used raw [7]. This implementation realises every feature derivable from the miRNA and the supplied 3'UTR:

- **Intercept**, **Local_AU** (position-weighted A/U fraction of the 30 nt up/downstream, rounded to 3 decimals as the perl does, then scaled), the **sRNA1{A,C,G}** / **sRNA8{A,C,G}** miRNA position-1/8 identity indicators (skipped when that nucleotide is U), and the **Site8{A,C,G}** target site-position-8 indicators (defined only for 7mer-A1 and 6mer).
- **3P_score** — 3' supplementary pairing, a faithful port of `get3primePairingContribution` (fed by `extractSubseqForAlignment` + `modifySubseqForAlignment`): the maximum, over single-gap offsets in both strands, of summed run-of-≥2 base pairs (1.0 in the offset-adjusted window 4..7, else 0.5) minus `max(0,(offset−2)/2)`; then min-max scaled (min=1, max=3.5). Raw scores were cross-checked against the reference perl.
- **Min_dist** — `log10(min(distTo5′, distTo3′))` of the site to the nearest 3'UTR end (single-isoform: the UTR length), min-max scaled (`getMinDist_weighted_contribution`).
- **Len_3UTR** — `log10(3'UTR length)`, min-max scaled (`get_len3UTR_weighted_contribution`).
- **Off6m** — count of the offset-6mer pattern (first 6 nt of `reverse_complement(seed 2-8)`) in the 3'UTR, non-overlapping like the perl's `/g` regex, used raw (`getOffset6merSites` + `getOffset6mer_weighted_contribution`). The public perl matches against `$utr_seq{$transcriptID}` (a hash reference, not the sequence) and stores the match list in scalar context, so it reports Off6m = 1 for every site; the library counts the occurrences, as Agarwal 2015 defines the feature.
- **Reported score** — `ContextScore` follows `readTargets`: each contribution rounded to 3 decimals (`getAgarwalContribution`), the sum rounded to 3 decimals and capped at `$maxContextScore` (8mer −0.03, 7mer-m8 −0.02, 7mer-A1 −0.01, 6mer 0); `null` for a site whose 1-based start is < 15 (`too_close`). `ContextScorePartial` stays the unrounded, uncapped sum.
- **Reference run (2026-11-25)** — `targetscan_70.pl` + `targetscan_70_context_scores.pl` run with perl 5.38 on a 181-nt UTR (let-7a; ORF 1500 nt, one ORF 8mer, TA/SPS from `TA_SPS_by_seed_region.txt`, RNAplfold values from ViennaRNA 2.8 `pfl_fold_up(seq, 20, 80, 40)`): with the RNAplfold value supplied, every printed contribution and the context++ score agree (8mer −0.742, 6mer −0.118, 7mer-A1 −0.257; test `MiRnaAnalyzer_ContextPlusPlusTargetScanParity_Tests`). Two perl defects are not reproduced: Off6m (above), and Site8 for 6mer sites — the 6mer alignment subsequence is 22 nt, padded with two leading N's, so `substr($subseq, 14, 1)` reads the base two positions 5' of the one opposite miRNA nt 8.

The **`TA_3UTR`** (transcriptome-wide target-site abundance) feature is now computable directly from a 3'UTR set via `ComputeTa3Utr`: per Garcia et al. (2011) [10] (Online Methods: *"the number of non-overlapping 3′UTR 8mer, 7mer-m8, and 7mer-A1 sites in the reference mRNAs"*) and the TargetScan precompute (`TA_SPS_by_seed_region.txt`, which stores `log10(count)` and feeds it as-is to `getAgarwalContribution` [7]), `TA = log10(N)` where `N` is the total non-overlapping 8mer + 7mer-m8 + 7mer-A1 site count of the seed across the supplied 3'UTRs. The computed value can be passed through `ContextPlusPlusInputs.Ta`. **Default TA_3UTR and SPS (B13 F33):** TargetScan itself does not count TA per run — `targetscan_70_context_scores.pl` reads TA and SPS from `TA_SPS_by_seed_region.txt` ("from Garcia et al., Supp Data 5"; one row for each of the 16 384 seed regions; column 2 = SPS of 8mer/7mer-m8, column 3 = SPS of 7mer-1a/6mer, column 4 = TA) keyed by `substr($matureMiRNA, 1, 7)`. The file (identical in the official `TargetScan7_context_scores.zip` and the nsoranzo mirror, SHA-256 d39a045c…8620) is bundled verbatim in `MiRnaAnalyzer.TargetScanTaSps.cs` (generated, one string per column, parsed once; header: source, rows, SHA-256, citations) and is used whenever the caller supplies no `Sps` / `Ta`; a caller value overrides that feature only. A seed region that is not 7 nt over ACGU has no row (TargetScan: contribution 0) and is listed as omitted. On the kit's own output (`targetscan_70_context_scores_output.txt`, 675 rows) every TA and SPS contribution is reproduced from the bundled table. **Len_ORF and ORF8m (B13 F34)** are derived from `ContextPlusPlusInputs.OrfSequence` exactly as the kit's `targetscan_count_8mers.pl` writes the two files the context script reads: gap characters `-` and `.` removed (README: "The ORF sequences can be aligned (with gaps) or not, but the alignment is ignored; gaps are removed"), length = remaining characters (Len_ORF = scaled log10, 0 for length 0); ORF8m = number of case-insensitive, non-overlapping (perl `/g`) occurrences of `reverse_complement(seed region) + "A"` on the given strand (T read as U), used raw. Caller `OrfLength` / `Orf8mCount` override. On the kit's ORF sample the lengths and 8mer counts equal `ORF_Sequences_sample.lengths.txt` / `ORF_8mer_counts_sample.txt` for every family, and the Len_ORF / ORF8m columns of the kit output are reproduced; on the worked example, only the ORF sequence + RNAplfold value (no TA/SPS/ORF numbers) reproduces every printed column of the perl run. Without an ORF these two features contribute 0 and are listed in `OmittedFeatures`. **PCT** (probability of conserved targeting) is computed when the caller supplies a `Conservation` input (a phylogenetic tree + the species in which the site is conserved; the per-site-type sigmoid parameters default to the bundled TargetScan 7.0 tables, a caller row overrides): the library derives the **Friedman et al. (2009) branch-length score (Bls)** — the total branch length of the minimal subtree connecting the conserved species [8] — then maps it to a PCT value via the published logistic `PCT(Bls) = b0 + b1/(1 + e^(−b2·Bls + b3))` (`targetscan_70_BL_PCT.pl`, `calculatePCTthisBL` [9]), and applies the bundled Agarwal PCT coefficient with the verbatim min-max scaling. TargetScan's reporting rules are reproduced: Bls is `sprintf("%.4f")` (`getBranchLength`) and the sigmoid is evaluated on that value; PCT is `sprintf("%.4f")` then truncated at 0; `getPCT` gives PCT 0 for Bls = 0 and for 6mer sites, and no PCT ("NA", contribution 0) for a negative Bls. Species that are not leaves of the tree are ignored. `IsConservedByBranchLength` gives the script's "Conserved" call (Bls ≥ 1.8 / 2.8 / 3.6 for 8mer / 7mer-m8 / 7mer-1a on the TargetScanHuman 7 trees; mouse thresholds 0.6 / 1.8 / 2.5; 6mer never). The per-seed-family `b0..b3` of TargetScan's BL/PCT kit (`targetscan_70_BL_PCT.zip`, `PCT_parameters/{8mer,7mer_m8,7mer_1a}_PCT_parameters.txt`, 111 seed families × 3 site types; TargetScan FAQ: free to use for research with a TargetScan citation — Friedman et al. 2009 [8] for BL/PCT) are **bundled verbatim** in `MiRnaAnalyzer.TargetScanPctParameters.cs` (generated from the kit files; file names, row counts and SHA-256 of the sources in its header). `TryGetTargetScanPctParameters(seedRegion, siteType, out p)` returns the row keyed as `targetscan_70_BL_PCT.pl` keys it — the 7-nt seed region (miRNA nt 2–8, `s/T/U/gi`; lower case is also accepted), 8mer / 7mer-m8 / 7mer-A1 only — and `IsTargetScanPctSeedFamily` tells whether a seed has PCT at all. When `PctConservation.SigmoidParameters` is `null`, context++ uses the bundled row for the miRNA's seed; a seed outside the 111 families gets PCT "NA" (`PctIsNa`, contribution 0, as `getPCT` + `readTargets`; not an omitted feature — it is TargetScan's own full score). **PCT evidence from an alignment (B13 F35/F36):** the species-with-site list and the bin tree no longer have to come from the caller. `FindTargetScanSiteGroups` ports `targetscan_70.pl` 7.0: per species (only those in which the family is annotated, `$FIND_SITES_ALL_SPECIES = 0`) the UTR is T→U-converted, its first 14 nucleotide letters are masked to N (`$BEG_UTR_MASK_LENGTH`), and the 8mer / 7mer-m8 / 7mer-A1 / 6mer seed matches are found on the upper-cased row with any number of `-` between consecutive nucleotides, at every start; shorter matches inside longer ones are dropped; sites of different species are grouped all-vs-all in ordinal order of `species::MSA_start::MSA_end` when they share a start or end column, overlap by ≥ 2 columns (`$REQUIRED_OVERLAP`) or nest — an overlapping pair joins the group of whichever site already has one, and two existing groups are never merged (`groupThisPair`, reproduced because it defines TargetScan's groups). Each group reports the group type (`7mer-m8+8mer-1a` …), the group's species and, per site type, the species having that type or a longer type containing it (`Species_in_this_group_with_this_site_type`) — the list `targetscan_70_BL_PCT.pl` computes the branch length of. The ORF-overlap flag reproduces the script's window (`substr($alignment, $start, $len)` with the 1-based start, i.e. one column right of the site) because it defines TargetScan's `ORF_overlap` column and PCT "NA". `AssignTargetScanUtrBin` ports `targetscan_70_BL_bins.pl`: for every column where the reference (9606) has a non-gap character, the unrounded connecting-subtree length on `Tree.generic.txt` of the species whose character equals it (case-insensitive); median (Statistics::Lite = `StatisticsHelper.Median`); bin = largest k with median > `@BL_thresholds[k−1]` (0, 1.21207417, …, 5.78196252). The ten bin trees and the generic tree are bundled as single-line Newick (generated `MiRnaAnalyzer.TargetScanTrees.cs`, SHA-256 of each kit file in its header). `ScoreTargetScanConservation` chains groups → bin → Bls on `Tree.bin_NN` → PCT (bundled sigmoid; "NA" for seeds without parameters or ORF-overlapping sites) → Conserved. From the raw kit alignment (`UTR_Sequences_sample.txt`, 28 rows, 6 genes, 14 families) the library reproduces all 213 rows of `targetscan_70_output.txt` verbatim (incl. group numbers continued over genes and families), the 6 bins and medians of `UTRs_median_BLs_bins.txt` (`%.15g` text equal) and every Bls / PCT / Conserved column of `targetscan_70_output.BL_PCT.txt` (213 Bls, 161 PCT + 52 "NA", no conserved site in the sample). Perl quirks not copied (no effect on valid input): species id "0" read as missing; the bins script reuses the previous gene's reference row when a gene lacks one (here an error); the 7mer-m8 / 7mer-A1 hash-order tie that exists only for seed UUUUUUU (7mer-m8 kept). On the kit's sample (213 rows) the library, using the bundled table, reproduces every Bls (213/213) and every PCT (161/161, incl. the 28 non-zero values) of the shipped `targetscan_70_output.BL_PCT.txt`, and "NA" exactly for the 52 rows of seeds without parameters. **SA** (the unpaired probability of the 14 nt ending 7 nt 3' of the seed-match start; `_lunp` row `utrStart+7`, or the last row when the UTR ends sooner, column `l=14`) is taken from `ContextPlusPlusInputs.PlfoldUnpairedProbability` when the caller supplies the RNAplfold value — then it equals TargetScan's. Otherwise it is estimated with the Turner-2004 McCaskill partition function on one 80-nt window centred on the region. That is not RNAplfold: RNAplfold `-L 40 -W 80` averages the partition functions of every 80-nt window holding the region and limits pairs to a 40-nt span, and the canonical partition function has no span limit. On the reference UTRs the estimate differed from RNAplfold by up to 3 log10 units (e.g. log10 p −0.86 vs −2.01 for the 8mer at 31), which moves the SA contribution by up to its whole coefficient. A zero probability contributes 0 (perl `! $plfold`). All context++ coefficients are the verbatim values from `Agarwal_2015_parameters.txt` [6] and every feature/scaling formula is ported from `targetscan_70_context_scores.pl` [7]. The returned `ContextScorePartial` is therefore a partial context++ score, not the published headline CS.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Canonical 8mer, 7mer-m8, 7mer-A1, 6mer, and offset-6mer seed classes [1][2][3][5].
- Antiparallel seed reverse-complement matching between miRNA and target [1][2].
- Priority ordering that prefers stronger canonical seed classes over weaker offset-6mer interpretations at overlapping positions [1][3].
- (Opt-in `ScoreTargetSiteContextPlusPlus`) The context++ features computable from the miRNA and the supplied 3'UTR — site-type Intercept, Local_AU, sRNA position-1/8 and target site-position-8 nucleotide-identity indicators, **3' supplementary pairing (3P_score), minimum distance to the 3'UTR end (Min_dist), 3'UTR length (Len_3UTR), and offset-6mer count (Off6m)** — plus **SPS and TA_3UTR** (bundled TargetScan 7.0 seed-region table, or caller values) and **Len_ORF, ORF8m** (from the ORF sequence as `targetscan_count_8mers.pl`, or caller values); all with the **verbatim fitted coefficients** from `Agarwal_2015_parameters.txt` [6] and the per-site-type min-max scaling / DP / indicator logic ported exactly from `targetscan_70_context_scores.pl` [7] (3P_score raw values cross-checked against the reference perl).
- (Opt-in `ScoreTargetSiteContextPlusPlus`, `ComputeBranchLengthScore`, `PctFromBranchLength`) The **PCT** (probability of conserved targeting) feature: the Friedman et al. (2009) **branch-length score** — the total branch length of the minimal subtree connecting the species in which the site is conserved [8] — and its mapping to a PCT value via the published logistic `PCT(Bls) = b0 + b1/(1 + e^(−b2·Bls + b3))`, 4 decimals, truncated at 0, 0 for Bls = 0 and 6mer sites [9], then entering context++ with the verbatim Agarwal PCT coefficient + min-max scaling [6][7]. Computed when the caller supplies a `Conservation` input (tree + conserved-species set + the published per-site-type `b0..b3`).
- (Opt-in `ComputeTa3Utr` / `CountSeedSites3Utr`) The **`TA_3UTR`** (target-site abundance) feature, computed from a supplied 3'UTR set: `TA = log10(N)`, `N` = total non-overlapping 8mer + 7mer-m8 + 7mer-A1 sites of the seed across the 3'UTRs, exactly the count defined by Garcia et al. (2011) [10] and stored (as log10) in TargetScan's `TA_SPS_by_seed_region.txt` [7]. The computed value enters context++ through `ContextPlusPlusInputs.Ta` (verbatim Agarwal TA coefficient + per-site-type min-max scaling).

**Intentionally simplified:**

- `Score` is the fitted TargetScan 5.0 context score of Grimson et al. (2007), negated and floored at 0 (`max(0, −ContextScore)`): site-type contribution −0.310 (8mer) / −0.161 (7mer-m8) / −0.099 (7mer-A1) + 3' pairing (slope × 3P score + intercept) + local AU (slope × position-weighted A/U fraction of 2 × 30 nt, 3-dp rounded) + position (slope × min(distance to nearest 3'UTR end, 1500), 3-dp rounded), with the constants of `targetscan_50_context_scores.pl` ("Regression info from Grimson et al., 2007 supp data"). **Oracle:** the script was re-run on its kit sample (`targetscan_50_context_scores.zip`, `vert_50` data-download page; 26 aligned UTRs, 76 mature miRNAs): all 400 rows with UTR start ≥ 32 are reproduced to 1e-9 for the 3' pairing, local AU, position contributions and the context score. **Consequences:** the 6mer and offset 6mer have no TargetScan 5.0 coefficients (Score 0, `ContextScore` `null`); the sequence passed in is treated as the 3'UTR (sites with 1-based start < 15 are `too_close`, no score); the 3' pairing alignment is the context++ port (shared). The opt-in `ScoreTargetSiteContextPlusPlus` provides the newer Agarwal (2015) context++ score. TargetScan 5.0 passes a negative `substr` offset for sites within 31 nt of the UTR start (the single kit row at start 21 differs); TargetScan fixed that in 5.3 by clipping, and so does this port.
- (Opt-in context++) **SA** without a supplied RNAplfold value is an approximation (one 80-nt McCaskill window, no 40-nt span limit, no window averaging; see above). **Consequence:** SA and the score can differ from TargetScan; supply `PlfoldUnpairedProbability` to reproduce it. Blocked on a span-limited partition function in `RnaSecondaryStructure` (request to the RNA-structure owner).
- (Opt-in context++) The PCT **sigmoid parameters** `b0..b3` (TargetScan 7.0 kit, 111 seed families × 3 site types; caller rows override) and the ten per-UTR-bin **trees** + generic tree are bundled; the **species with the site** and the **UTR bin** are derived from a multi-species 3'UTR alignment (`FindTargetScanSiteGroups`, `AssignTargetScanUtrBin`, end-to-end `ScoreTargetScanConservation`; B13 F35/F36) or still accepted from the caller (`PctConservation`). **TA_3UTR / SPS** come from the bundled TargetScan table (or caller values; `ComputeTa3Utr` computes TA from your own 3'UTR set); **Len_ORF / ORF8m** need the transcript's ORF (`OrfSequence`, or the two numbers). The only context++ inputs the library cannot derive or bundle are the multi-species alignment itself and the exact RNAplfold SA value. **Consequence:** without an alignment (and without an ORF), `ContextScorePartial` is a partial context++ score, reported alongside the explicit `OmittedFeatures` list of whatever remains residual for that call.
- Duplex free energy is the Turner 2004 energy of the fixed-register ungapped duplex (MIRNA-PAIR-001); **consequence:** a site whose 3′ pairing needs a bulge is scored in the seed register only — use `FoldMiRnaDuplex` for the MFE gapped duplex.
- Context and accessibility are exposed only as separate helper methods; **consequence:** the main finder does not integrate AU context, conservation, or accessibility into site discovery or ranking.
  - `CalculateSiteAccessibility` = P(site entirely unpaired) = Z_open/Z from the Turner-2004 McCaskill partition function, delegated to the canonical `RnaSecondaryStructure.CalculateRegionUnpairedProbability` over an RNAplfold-style local context of W = max(80, site length) nt centred on the site (same helper as the context++ SA feature). RNAplfold's L = 40 span limit and averaging over all W-windows are not applied. Cross-check vs ViennaRNA 2.x (Z_c/Z, hard-constrained unpaired site): 0.9758 vs 0.9770 (dangles 2) / 0.9938 (dangles 0); 0.9674 vs 0.9701 / 0.9930; helix-locked site 1.28e-4 vs 1.28e-4 / 1.9e-3 (canonical d2 dangles + McCaskill grammar, B12 RNA-STRUCT-001 F4/F7; before: 0.9721, 0.9630, 3.2e-4). (Replaced 2026-09-28 a Watson-Crick pair-density heuristic that returned 0 for all three.)
  - `AnalyzeTargetContext` implements Grimson et al. (2007) context determinants [3]: `AuContent` = TargetScan 5.0 `getLocalAUcontribution` weighted A/U fraction of ≤ `contextWindow` (30) nt per flank, site excluded, with the weights of the `siteType` argument (default 7mer-m8: upstream 1/(i+1), downstream ½ for i = 0 then 1/(i+1)); `NearStart` = 1-based start < 15 (`$MIN_DIST_TO_CDS`); `NearEnd` = end > 0.85·length (descriptive); `ContextScore` = 0 when `NearStart` or when `siteType` has no TargetScan 5.0 coefficients, else the UTR-context part of the Grimson score, round3(slope_AU·AuContent + intercept_AU) + round3(slope_pos·min(d5′,d3′,1500) + intercept_pos) in log2 units (negative = repressed). The site-type and 3' pairing terms need the miRNA: `ScoreTargetSiteGrimson`. (Re-based 2026-11-25 on the TargetScan 5.0 coefficients; the earlier 0.5·AU + 0.5·end-proximity combination was an unfitted heuristic.)

**Not implemented:**

- Bulged / mismatched seed sites and 3'-compensatory sites; **users should rely on:** an external miRNA target-prediction tool for broader site classes. 3'-supplementary pairing is not a site class: Grimson et al. (2007) score it as the 3' pairing term of a seed site (`ScoreTargetSiteGrimson(...).ThreePrimePairingScore`).
- The full context++ headline score needs inputs that are not in a single 3'UTR: the transcript ORF (`OrfSequence`), the multi-species 3'UTR alignment for PCT (from which `FindTargetScanSiteGroups` / `ScoreTargetScanConservation` derive the species with the site and the bin tree) and, for TargetScan's exact SA, the RNAplfold value. TA_3UTR, SPS and the PCT sigmoid parameters are bundled. **Users should rely on:** supplying these via `ContextPlusPlusInputs`, the official TargetScan distribution, or treating `ContextScorePartial` as a partial score.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `TargetSiteType.Supplementary` is never emitted by any method (3'-supplementary pairing is the 3' pairing term of a seed site, Grimson 2007, not a site class); `TargetSiteType.Centered` is emitted only by `FindCenteredSites` (Shin et al. 2010), never by `FindTargetSites`. | Deviation | The public enum keeps both members for backward compatibility. | accepted | `Supplementary` documented as never emitted; centered sites: ≥ 11 contiguous Watson–Crick pairs at miRNA nt 4–14 / 5–15 / 4–15, nt 2–7 not all paired, G:U not counted (abstract + reprint wording; no reference implementation, the oracle is the hand construction in the tests); the "lacks 3'-compensatory pairing" criterion is not tested. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or `null` mRNA | Returns no sites. | The finder exits early for empty input. |
| Empty `miRna.Sequence` | Returns no sites. | The finder exits before seed matching if the miRNA sequence is empty. |
| Seed shorter than 7 nt | Returns no sites. | Canonical seed reverse-complement matching requires the 7-mer seed. |
| DNA target sequence | `T` is treated as `U`. | mRNA input is normalized with `Replace('T', 'U')`. |
| `minScore` above the largest context score | Returns no sites. | `Score = max(0, −ContextScore)` has no fixed ceiling; ≈ 0.65 for a strong 8mer in an AU-rich 3'UTR. |
| `minScore ≤ 0` | Every candidate site is emitted. | Scores are compared directly (`Score ≥ minScore`). |

### 6.2 Limitations

The finder reports only exact canonical seed-complement sites and offset-6mer sites. It does not model bulges, gaps, supplementary-only or 3'-compensatory sites, evolutionary conservation, transcript isoforms, or the 6mer / offset-6mer effect (no TargetScan 5.0 coefficients). `FreeEnergy` is an implementation-level ranking aid, not a biophysical confidence measure; `Score` is the Grimson (2007) fitted model.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
using Seqeron.Genomics.Annotation;

var let7a = MiRnaAnalyzer.CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU");
string mrna = "GGGGGCUACCUCAGGGGG";

var sites = MiRnaAnalyzer.FindTargetSites(mrna, let7a, minScore: 0.0).ToList();

// sites[0].Type == TargetSiteType.Seed8mer
// sites[0].SeedMatchLength == 8

// Opt-in TargetScan context++ score (Agarwal et al. 2015) for the discovered site:
var ctx = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(mrna, let7a, sites[0]);
// ctx.ContextScorePartial    -> partial CS (Intercept + Local_AU + sRNA1/8 + Site8
//                               + 3P_score + Min_dist + Len_3UTR + Off6m + TA + SPS from the bundled table)
//                               = -0.6031574254937664 here
// ctx.OmittedFeatures        -> still-residual features (here SA, PCT, Len_ORF, ORF8m)
// ctx.ContextScore           -> null here: the site starts at nt 6 (< 15, TargetScan "too_close");
//                               otherwise the 3-dp, ceiling-capped score TargetScan prints

// Compute TA_3UTR from a transcriptome (a 3'UTR set) per Garcia 2011: TA = log10(N),
// N = non-overlapping 8mer + 7mer-m8 + 7mer-A1 sites of the seed across the 3'UTRs.
double ta = MiRnaAnalyzer.ComputeTa3Utr(let7a, threePrimeUtrSet);  // IEnumerable<string>

// The ORF gives Len_ORF and ORF8m (as targetscan_count_8mers.pl); TA/SPS default to the bundled
// TargetScan table and may be overridden (here TA from your own 3'UTR set):
var withData = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(
    mrna, let7a, sites[0],
    new MiRnaAnalyzer.ContextPlusPlusInputs(Ta: ta, OrfSequence: cdsSequence,
        PlfoldUnpairedProbability: lunpRowUtrStartPlus7Col14)); // RNAplfold -L 40 -W 80 -u 20 value
// MiRnaAnalyzer.TryGetTargetScanTaSps("GAGGUAG", MiRnaAnalyzer.TargetSiteType.Seed8mer, out var ta7, out var sps7); // 3.393, -9.25

// PCT from multi-species conservation (Friedman 2009 Bls → published sigmoid → context++):
var tree = PhylogeneticAnalyzer.ParseNewick("((A:1.0,B:2.0):0.5,(C:1.5,D:3.0):4.0);");
var conservation = new MiRnaAnalyzer.PctConservation(
    tree,
    new[] { "A", "B" });                                       // species in which the site is conserved
// b0..b3 default to the bundled TargetScan 7.0 row for let-7a's seed GAGGUAG and the site type;
// pass a third argument new MiRnaAnalyzer.PctSigmoidParameters(b0, b1, b2, b3) to override.
// MiRnaAnalyzer.TryGetTargetScanPctParameters("GAGGUAG", MiRnaAnalyzer.TargetSiteType.Seed8mer, out var p8);
var withPct = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(
    mrna, let7a, sites[0], new MiRnaAnalyzer.ContextPlusPlusInputs(Conservation: conservation));
// withPct.BranchLengthScore == 3.0 ; withPct.Pct == sigmoid(3.0) to 4 dp (0 for a 6mer; PctIsNa for a seed outside the 111 PCT families) ; withPct.PctContribution folded into the score
// MiRnaAnalyzer.IsConservedByBranchLength(withPct.BranchLengthScore, site.Type) -> TargetScan "Conserved" column
```

**Numerical walk-through (8mer, let-7a, all-G flanks):** with `mRNA = GGGGG + CUACCUCA + GGGGG` (length 18), the 8mer occupies `Start=5,End=12`; both 30 nt flanks are all-G so the local-AU fraction is 0. The realised contributions are: `Intercept(8mer)=-0.589`; `Local_AU = -0.254 × ((0 − 0.308)/(0.814 − 0.308)) = +0.154608695652174`; `sRNA8G(8mer)=+0.015` (nt1=U ⇒ sRNA1 unscored; Site8 undefined for 8mer); `3P_score`: raw=0 ⇒ `-0.040 × ((0 − 1)/(3.5 − 1)) = +0.016`; `Min_dist`: nearest end = min(5, 18−13) = 5 ⇒ `0.118 × ((log10 5 − 1.415)/(3.113 − 1.415)) = -0.049759446106213065`; `Len_3UTR`: `0.310 × ((log10 18 − 2.392)/(3.637 − 2.392)) = -0.2830405810586145`; `Off6m`: pattern `CUACCU` occurs once ⇒ `-0.020 × 1 = -0.020`. These sum to -0.7561913315126536; the bundled TargetScan row of seed GAGGUAG (TA 3.393, 8mer SPS −9.250) adds `TA_3UTR = 0.222 × ((3.393 − 3.113)/(3.865 − 3.113)) = +0.08265957446808503` and `SPS = 0.210 × ((−9.250 + 11.13)/(−5.52 + 11.13)) = +0.07037433155080215`, so `ContextScorePartial` = **-0.6031574254937664** (SA omitted because the 14-nt window does not fit; no ORF ⇒ Len_ORF/ORF8m omitted; no conservation ⇒ PCT omitted). `ContextScore` is `null`: the site starts at nt 6, inside TargetScan's 15-nt `too_close` zone. All numbers are the verbatim Agarwal_2015_parameters.txt coefficients.

**Numerical walk-through (Grimson 2007 / TargetScan 5.0 `Score`, 8mer, let-7a, A flanks):** `mRNA = A×40 + CUACCUC + A + A×40` (length 88); the 8mer occupies `Start=40,End=47`. Site type −0.310; local AU: all flank nt are A/U ⇒ fraction 1 ⇒ round3(−0.64 × 1 + 0.365) = −0.275; position: 5.0 coordinates `utrStart=41`, `utrEnd=47` (the A1 lies outside) ⇒ distance min(40, 88−47=41) = 40 ⇒ round3(0.000172 × 40 − 0.07) = −0.063; 3' pairing raw score 2 ⇒ 0.011 − 0.0041 × 2 = +0.0028. Context score −0.6452, `Score` 0.6452. The same construction with a 7mer-m8 (`…CUACCUC + G + A×40`) gives −0.387, a 7mer-A1 −0.3082, a 6mer has no score.

**Oracle run (TargetScan 5.0):** `perl targetscan_50_context_scores.pl miR_for_context_scores.txt UTR_Sequences_sample.txt targetscan_50_output.txt out.txt` (kit SHA-256 72e95029…90c6, perl 5.38) gives 401 rows; the library equals 400 of them to 1e-9 (3' pairing, local AU, position, context score). The 401st (FNDC3A 9615 cfa-miR-101, 7mer-A1 at 21) hits perl's negative-`substr` defect (printed AU −0.095) and is clipped here (−0.147, independent derivation in the tests). The shipped output file differs from the re-run only in 65 7mer-m8 rows at the 1500-nt position cap (0.0995 prints 0.099 there, 0.100 on x86-64).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MiRnaAnalyzer_TargetPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/MiRnaAnalyzer_TargetPrediction_Tests.cs)
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
7. TargetScan distribution. `targetscan_70_context_scores.pl` (context++ feature computation and scaling: `getAgarwalContribution`, `getLocalAU_contribution`, `get_sRNA1_8_contributions`, `getSite8_contribution`, `getPCT_contribution`). https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_context_scores.pl
8. Friedman RC, Farh KKH, Burge CB, Bartel DP. 2009. Most mammalian mRNAs are conserved targets of microRNAs. Genome Research. 19(1):92-105. doi:10.1101/gr.082701.108. https://pmc.ncbi.nlm.nih.gov/articles/PMC2612969/ (branch-length score; PCT = E[(S−B)/S]). PCT documentation: https://www.targetscan.org/docs/pct.html
9. TargetScan distribution. `targetscan_70_BL_PCT.pl` (`getBranchLength`, `getPCT`, `calculatePCTthisBL`: `PCT = b0 + b1/(1 + e^(−b2·BL + b3))`, `sprintf("%.4f")`, truncated at 0; `%siteTypeToConsThreshold`; per-seed-family `b0..b3` from `PCT_parameters/*_PCT_parameters.txt`). https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_BL_PCT.pl ; kit with parameter files, trees and sample output: https://www.targetscan.org/vert_70/vert_70_data_download/targetscan_70_BL_PCT.zip
10. Garcia DM, Baek D, Shin C, Bell GW, Grimson A, Bartel DP. 2011. Weak seed-pairing stability and high target-site abundance decrease the proficiency of lsy-6 and other microRNAs. Nature Structural & Molecular Biology. 18(10):1139-1146. doi:10.1038/nsmb.2115. https://pmc.ncbi.nlm.nih.gov/articles/PMC3190056/ (TA = number of non-overlapping 3′UTR 8mer/7mer-m8/7mer-A1 sites; stored as log10 in TargetScan's `TA_SPS_by_seed_region.txt`).
11. TargetScan 7.0 context-score kit `TargetScan7_context_scores.zip` (www.targetscan.org/vert_70/vert_70_data_download/; SHA-256 bc79560d…118e): `TA_SPS_by_seed_region.txt` (bundled as `MiRnaAnalyzer.TargetScanTaSps.cs`), `targetscan_count_8mers.pl` (ORF length and ORF 8mer counts), `README_CONTEXT_SCORES_70.txt`, and the sample output `targetscan_70_context_scores_output.txt` used as the oracle.
12. TargetScan 7.0 target-prediction kit `targetscan_70.zip` (www.targetscan.org/vert_70/vert_70_data_download/; SHA-256 6e23c72a…d749; opened 2026-11-25): `targetscan_70.pl` 7.0 (site finding and cross-species grouping; same logic as the nsoranzo mirror up to an unused `%pairToDistance` and an equivalent split of the same-start/same-end test), sample `UTR_Sequences_sample.txt`, `miR_Family_info_sample.txt`, `targetscan_70_output.txt`; and `targetscan_70_BL_bins.pl`, `UTRs_median_BLs_bins.txt`, `PCT_parameters/Tree.bin_01..10.txt`, `Tree.generic.txt` from `targetscan_70_BL_PCT.zip` (SHA-256 5ecd0e90…5885).
