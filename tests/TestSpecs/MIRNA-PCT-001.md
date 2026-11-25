# Test Specification: MIRNA-PCT-001

**Test Unit ID:** MIRNA-PCT-001
**Area:** MiRNA
**Algorithm:** TargetScan PCT (Branch-Length Conservation)
**Status:** ☑ Complete — re-validated 2026-11-25 (B13 review: Stage A PASS-WITH-NOTES / Stage B FAIL→FIXED, F19–F22)
**Last Updated:** 2026-11-25 (WP-2: bundled PCT parameters; WP-3: species with the site, bin and trees from a UTR alignment)

> **Re-validated (B13, 2026-11-25).** TargetScan 7.0 `targetscan_70_BL_PCT.pl` was run (perl 5.38.2 + BioPerl) on
> the official BL/PCT kit (`targetscan_70_BL_PCT.zip`, targetscan.org vert_70); it reproduces the shipped
> `targetscan_70_output.BL_PCT.txt` exactly. The library now reproduces all 213 Bls and all 161 PCT values of that
> output (bin trees + kit parameter rows). Fixed: 4-dp rounding of Bls and PCT (F20), `getPCT` rules Bls = 0 ⇒ 0 and
> 6mer ⇒ 0 (F19), species absent from the tree ignored (F21), argument validation + negative Bls ⇒ NA, and the
> "Conserved" call (F22). The 2026-06-25 claims "PCT(0) is the untruncated logistic floor" and "4-dp rounding is a
> benign divergence" were wrong (they were not checked against `getPCT`).
>
> **Bundled parameters (B13 WP-2, F32, 2026-11-25).** The kit's three `PCT_parameters/*_PCT_parameters.txt` tables
> (111 seed families × 8mer / 7mer-m8 / 7mer-1a) are bundled verbatim in
> `MiRnaAnalyzer.TargetScanPctParameters.cs` (generated; SHA-256 of each source in the header) and are the default
> when no caller row is given. With the bundled table the unit test reproduces all 213 Bls, all 161 PCT and the 52
> "NA" rows of the kit output.
>
> **PCT evidence from a multi-species UTR alignment (B13 WP-3, F35/F36, 2026-11-25).** `FindTargetScanSiteGroups`
> ports `targetscan_70.pl` 7.0 (site finding + cross-species grouping), `AssignTargetScanUtrBin` ports
> `targetscan_70_BL_bins.pl`, the 10 bin trees + generic tree are bundled as single-line Newick
> (`MiRnaAnalyzer.TargetScanTrees.cs`), and `ScoreTargetScanConservation` chains alignment → groups → bin → Bls →
> PCT → Conserved. From the raw kit alignment the library reproduces `targetscan_70_output.txt` (213 rows verbatim),
> `UTRs_median_BLs_bins.txt` (6 genes) and every Bls / PCT / Conserved of `targetscan_70_output.BL_PCT.txt`. The
> caller-supplied path (`PctConservation`) is unchanged.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Friedman, Farh, Burge & Bartel (2009) Genome Res 19:92 — Bls = total branch length of the tree connecting the species with the site perfectly aligned; PCT |
| 2 | `targetscan_70_BL_PCT.pl` (raw.githubusercontent.com/nsoranzo/targetscan/main, and the kit copy — identical apart from paths/comments): `getBranchLength`, `lookAtBranchLength`, `getPCT`, `calculatePCTthisBL`, `%siteTypeToConsThreshold`, `selectConservationForThisSite` |
| 3 | `targetscan_70_BL_PCT.zip` (targetscan.org/vert_70/vert_70_data_download): `PCT_parameters/{8mer,7mer_m8,7mer_1a}_PCT_parameters.txt` (111 seed families each), `Tree.bin_01..10.txt`, `Tree.generic.txt`, sample input/output, `README_BL_PCT_70.txt` |
| 4 | `targetscan_70_context_scores.pl` `readTargets`/`getPCT_contribution` (PCT "NA" ⇒ contribution 0) |

## 2. Canonical Method(s)

- `MiRnaAnalyzer.ComputeBranchLengthScore(PhyloNode tree, IReadOnlyCollection<string> speciesWithSite)` — Bls, 4 dp
- `MiRnaAnalyzer.PctFromBranchLength(double bls, PctSigmoidParameters p)` and `(…, TargetSiteType siteType)` — PCT, 4 dp
- `MiRnaAnalyzer.IsConservedByBranchLength(double bls, TargetSiteType, BranchLengthConservationThresholds?)`
- `MiRnaAnalyzer.TryGetTargetScanPctParameters(string seedRegion, TargetSiteType, out PctSigmoidParameters)` /
  `IsTargetScanPctSeedFamily(string seedRegion)` — bundled TargetScan 7.0 rows (`MiRnaAnalyzer.TargetScanPctParameters.cs`)
- wired into `ScoreTargetSiteContextPlusPlus` via `ContextPlusPlusInputs.Conservation`
- `MiRnaAnalyzer.FindTargetScanSiteGroups(string seedOrMature, IReadOnlyDictionary<string,string> alignedUtrs, IReadOnlyCollection<string>? speciesWithMiRna)` — `targetscan_70.pl` sites + groups (`TargetScanSiteGroup`, `TargetScanAlignedSite`)
- `MiRnaAnalyzer.AssignTargetScanUtrBin(alignedUtrs, referenceSpecies = "9606", tree = null)` — `targetscan_70_BL_bins.pl`; `TargetScanUtrBinThresholds`
- `MiRnaAnalyzer.GetTargetScanBinTree(int)` / `GetTargetScanBinTreeNewick(int)` / `GetTargetScanGenericTree()` / `TargetScanGenericTreeNewick` — bundled trees (`MiRnaAnalyzer.TargetScanTrees.cs`)
- `MiRnaAnalyzer.ScoreTargetScanConservation(seedOrMature, alignedUtrs, speciesWithMiRna, referenceSpecies, utrBin, thresholds)` — end to end (`MiRnaAnalyzer.TargetScanConservation.cs`)
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
- **Test fixtures:** `Unit/Annotation/MiRnaAnalyzer_PctTargetScanParity_Tests.cs` (TargetScan parity; kit rows and
  bin trees in the generated `Unit/Annotation/MiRnaPctTargetScanKitData.cs`),
  `Unit/Annotation/MiRnaAnalyzer_TargetPrediction_Tests.cs` (region CTX-PCT), `Fuzzing/MiRnaPctFuzzTests.cs`,
  `Unit/Annotation/MiRnaAnalyzer_TargetScanAlignment_Tests.cs` (alignment pipeline; kit files verbatim in the generated
  `Unit/Annotation/MiRnaTargetScanAlignmentKitData.cs`)

## 3. Contract / Invariants

- Bls = total length of the minimal subtree connecting the species that are leaves of the tree; 0 for < 2 such
  species; independent of the root position; reported `sprintf("%.4f")`.
- PCT = `sprintf("%.4f", b0 + b1/(1+e^(−b2·Bls+b3)))`, then `< 0 ⇒ 0`; Bls = 0 ⇒ 0; 6mer ⇒ 0; Bls < 0 / NaN / ∞ ⇒
  `ArgumentOutOfRangeException` (perl "NA"); non-finite b ⇒ `ArgumentException`.
- R: PCT ∈ [0,1] for published-style parameters; M: non-decreasing in Bls for b1, b2 > 0; D: deterministic.
- Conserved ⇔ Bls ≥ threshold(site type); TargetScanHuman 7 = 1.8 / 2.8 / 3.6 (8mer / 7mer-m8 / 7mer-1a), mouse
  0.6 / 1.8 / 2.5; 6mer never.
- context++: PCT contribution = Agarwal min-max scaled 4-dp PCT; negative Bls ⇒ contribution 0 and PCT listed omitted.
- Parameter lookup (`readAllPCTdata` / `readMiRNAs`): key = 7-nt seed region (miRNA nt 2–8), `s/T/U/gi` (lower case
  also accepted), ordinal; rows for 8mer / 7mer-m8 / 7mer-1a only (6mer and other types ⇒ not found); 111 families.
  No caller row ⇒ bundled row for the miRNA's seed; seed outside the table ⇒ PCT "NA" (`PctIsNa`, contribution 0,
  not omitted); a caller row overrides the table.
- Sites (`targetscan_70.pl`): T/t→U; first 14 A/C/G/T/U letters (any case) → N; 8mer = revcomp(seed)+A, 7mer-m8 =
  revcomp(seed), 7mer-A1 = its last 6 + A, 6mer = its last 6; '-' allowed between consecutive letters; every start;
  only the species in which the family is annotated (null = all). 8mer drops the 7mer-m8/7mer-A1/6mer it contains,
  7mer the 6mer. MSA start/end 1-based columns; UTR start = non-gap characters in columns 1..start; UTR end = start +
  length − 1. ORF_overlap = lower case in columns start+1 … end+1 (the script's off-by-one window, copied).
- Groups: all-vs-all over ordinal "species::start::end", different species; same start/end, ≥ 2-column overlap or
  nesting; a pair joins an existing group of either site, existing groups never merged; singletons numbered last.
  Group type = distinct type names in id order (7mer-1a, 7mer-m8, 8mer-1a, 6mer) joined by "+"; per-type species =
  sites of that type or a longer type containing it; lists unique, perl-numeric order.
- Bin (`targetscan_70_BL_bins.pl`): reference (9606) non-gap columns; unrounded connecting-subtree length on
  Tree.generic of the species whose character equals the reference's (case-insensitive); median; bin = largest k with
  median > threshold[k−1]; no reference nucleotide ⇒ NaN / bin 1; reference row missing or not a tree leaf ⇒
  `ArgumentException`.
- End to end: Bls on Tree.bin_NN of the site type's species list; PCT "NA" (null) when the seed has no PCT row or
  the site overlaps the ORF; else `getPCT`; Conserved per thresholds.

## 4. Cross-check / Differential Oracle

| ID | Case | Oracle (perl output) |
|---|---|---|
| PCT-KIT-01..09 | kit bin-1 rows: let-7 6mer {10090,10116}; let-7 8mer 3/5 spp; let-7 7mer-m8/7mer-1a; miR-18 7mer-m8; miR-10 7mer-m8 (`-0.0000`); miR-101 7mer-1a; miR-1 single species | 0.1513/0.0, 1.1469/0.0, 1.6379/0.0, 1.6379/0.0058, 1.5614/0.0, 1.6379/0.0547, 0.9265/−0.0000, 0.9265/0.0, 0/0.0 |
| PCT-RND-01 | tree `((9606:1.60004,10090:1.6):0.5,(9031:1.5,9615:3.0):4.0)`, miR-30 8mer, {9606,10090} / {9606,9031} | 3.2000 / 0.1787 / x ; 7.6000 / 0.4605 (unrounded Bls would give 0.1788) |
| PCT-ABS-01 | same tree, "9606 10090 99999" / "9606 99999" / "9606 9606 10090" | 3.2000 / 0.0000 / 3.2000 |
| PCT-ZERO-01 | (0.1,0.8,2,4): single species; 6mer 3.2; 7mer-1a 12.2 | 0.0 ; 0.0 ; 0.9000 x |
| PCT-CONS-01 | thresholds | 7mer-1a 3.2000 not conserved; 8mer 3.2 conserved |
| CTX-PCT-003..010 | worked tree `((A:1,B:2):0.5,(C:1.5,D:3):4)` with (0,1,1,0) and miR-30 8mer row | 0.9526, 0.9991, 6mer 0.0, Bls sweep 0→0.0 … 8→0.9997, miR-30 3/4/5 → 0.0/0.4605/0.4605 |

| PCT-KIT-ALL | every row of `targetscan_70_output.BL_PCT.txt` (bin 1/2 trees, species lists from `targetscan_70_output.txt`), PCT from the **bundled** row | 213/213 Bls, 161/161 PCT (28 non-zero: 0.0058×12, 0.0547×10, 0.0239×5, 0.0071×1), 52 NA ⇔ seed not in table, no "x" |
| PCT-TAB-01 | spot rows (let-7 ×3, miR-18, miR-10, miR-101, miR-1, first 8mer row UCCCUUU) | equal to the kit file text |
| PCT-TAB-02 | key: GAGGTAG / gaggtag found; 6 / 8 nt, family name, 6mer, offset-6mer not found; miR-601…609 seeds have no PCT | perl `s/T/U/gi` + hash lookup |
| PCT-CTX-TAB | let-7a 7mer-m8, kit species "9031 9606 9615 10090 10116", bin-1 tree, no caller row | Bls 1.6379, PCT 0.0058 (kit), contribution −0.048×0.0058/0.364; caller row (0,1,1,0) ⇒ 0.8372 |
| PCT-CTX-NA | miR-601 (GGUCUAG) 7mer-m8, no caller row | PCT NA ⇒ 0, contribution 0, not omitted |
| PCT-ALN-KIT | kit `UTR_Sequences_sample.txt` × `miR_Family_info_sample.txt` (6 genes × 14 families, perl `sort keys` order) | all 213 rows of `targetscan_70_output.txt` verbatim (group numbers 1…121 continued over genes/families) |
| PCT-ALN-BIN | kit alignment, generic tree | `UTRs_median_BLs_bins.txt`: CDC2L6 0.64763119 1, FNDC3A 1.48570309 2, LIN28B 1.12468545 1, LPHN1 0.64763119 1, NLRP1 0 1, ZNF197 0.28661355 1 (`%.15g` text equal) |
| PCT-ALN-E2E | raw kit alignment → `ScoreTargetScanConservation` | every row of `targetscan_70_output.BL_PCT.txt`: UTR coords, group, Bls 213/213, PCT 161/161 (28 non-zero), 52 NA, Conserved 0 |
| PCT-ALN-01 | constructed 9606/10090 alignment (masked lower-case/t/N/gap prefix, gapped 8mer, lower-case first letter vs letter after the site, final `tUACCUCA`) | 12 perl rows verbatim (ORF_overlap 0 / 1 by the shifted window; unpaired sites groups 6, 7) |
| PCT-ALN-02 | species 1–5, gapped 8mer of 5 overlapping two existing groups | perl: groups {1,2,5} and {3,4} (not merged) |
| PCT-ALN-03 | constructed 11-row alignment (lower case, short row, N, gaps, 99999 not in tree) | perl bins: 2.55195311, bin 3 |
| PCT-ALN-04 | FNDC3A miR-1 8mer group, per-type lists → `ComputeBranchLengthScore(GetTargetScanBinTree(2), …)` | equals the end-to-end Bls; 8mer list 9606 9615 10090 10116, 7mer-m8 list adds 9031 |

## 5. Validation Checklist

- [x] Stage A: formula, Bls definition, `getPCT` rules, rounding, thresholds and parameter-column order re-read in the perl and the kit.
- [x] Stage B: perl run on the kit sample and on constructed inputs; library equal on every number above.
- [x] WP-3: `targetscan_70.pl` (official 7.0 and mirror both reproduce the kit output) and `targetscan_70_BL_bins.pl` (+ Statistics::Lite) re-run; library equal on PCT-ALN rows.
- [x] Targeted tier green (MiRna; the 4 MIRNA-PAIR duplex failures pre-date this unit — see B13.md), Mcp.Annotation green; heavy-tier PCT classes re-locked and green.
