# Validation Report: MIRNA-PCT-001 — TargetScan PCT (branch-length score → logistic)

- **Re-validated:** 2026-11-25 (B13 review; supersedes the 2026-06-25 report)   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.ComputeBranchLengthScore`, `MiRnaAnalyzer.PctFromBranchLength` (+ site-type overload),
  `MiRnaAnalyzer.IsConservedByBranchLength` / `BranchLengthConservationThresholds`; wired into `ScoreTargetSiteContextPlusPlus`
  through `ContextPlusPlusInputs.Conservation` (`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`).
- **Stage A verdict:** PASS-WITH-NOTES   **Stage B verdict:** FAIL → FIXED (B13 F19–F22)   **State:** FIXED

## Sources opened (2026-11-25)

1. Friedman et al. (2009) Genome Res 19:92 — Bls definition (Methods), PCT (via the 2026-06-25 quotation; PMC blocked for fetch).
2. `targetscan_70_BL_PCT.pl` — raw.githubusercontent.com/nsoranzo/targetscan/main and the copy inside the kit (same code; only
   the `PCT_parameters` path and two comments differ).
3. `targetscan_70_BL_PCT.zip` from targetscan.org (vert_70 data download): `README_BL_PCT_70.txt`, `PCT_parameters/8mer|7mer_m8|7mer_1a_PCT_parameters.txt`
   (111 seed families each), `Tree.bin_01..10.txt`, `Tree.generic.txt`, `miR_Family_info_sample.txt`, `targetscan_70_output.txt`,
   `UTRs_median_BLs_bins.txt`, `targetscan_70_output.BL_PCT.txt`.
4. `targetscan_70_context_scores.pl` (`readTargets`: PCT "NA" ⇒ contribution 0; `getPCT_contribution`).
5. PHYLIP Newick description (phylipweb.github.io/phylip/newicktree.html): a tree "can continue to subsequent lines"; blanks may appear anywhere except inside names / branch lengths.

## Reference run

perl 5.38.2 + BioPerl (`libbio-perl-perl`). `perl targetscan_70_BL_PCT.pl miR_Family_info_sample.txt targetscan_70_output.txt UTRs_median_BLs_bins.txt`
reproduces the shipped `targetscan_70_output.BL_PCT.txt` (diff empty). A C# harness fed the same species lists, bin trees and parameter rows to the
library: **213/213 Bls and 161/161 PCT equal** (exact double equality after the fixes). Constructed inputs (own `PCT_parameters` dir, trees with
NCBI taxonomy ids) supplied the edge cases below; every expected test value is a number the perl printed.

## Stage A — Description

- Bls: `getBranchLength` climbs from each species to the first already-included node; this is the minimal connecting-subtree (Steiner) length,
  independent of rooting. It returns `sprintf("%.4f")`, and that rounded string is what `getPCT`, the thresholds and the output use. **Note:** the
  description claimed full precision as "more precise, benign"; TargetScan's reported value is the rounded one.
- PCT: `getPCT` evaluates `calculatePCTthisBL` only for Bls > 0 and a non-6mer type; Bls = 0 ⇒ "0.0"; 6mer ⇒ "0.0"; non-numeric ⇒ "EMPTY";
  negative ⇒ NA. `calculatePCTthisBL`: `sprintf("%.4f")`, then `< 0 ⇒ "0.0"`. **Note:** the description claimed "PCT(0) = b0 + b1/(1+e^b3),
  used as-is"; that contradicts `getPCT`.
- Conserved call: `%siteTypeToConsThreshold` 1.8 / 2.8 / 3.6 (TargetScanHuman, TS7 trees; mouse block 0.6 / 1.8 / 2.5; 6mer 100). The README's
  "0.8 / 1.3 / 1.6" is stale; the code values produced the shipped output. Not described before.
- Parameter rows: columns 2–5 = b0..b3 (confirmed). The tables are downloadable (kit), not "unpublished"; doc corrected (no 6mer table:
  6mer PCT is 0). **Update (WP-2, F32):** now bundled — see "Bundled parameter tables" below.
- Effective grouping: `processPreviousGroup` tests `$siteType ne "8mer"` while the type is named "8mer-1a", so the multi-BL branch (8mer group with
  m8/1a sites, best of three PCTs) never runs; every site uses the species list of its own site type. The library's single species set matches this.

## Stage B — Implementation (findings, all fixed)

| F | Defect | Repro (library before → perl) | Fix |
|---|---|---|---|
| F19 | `getPCT` rules missing | (0.1,0.8,2,4), Bls 0: 0.1144 → 0.0; 6mer {A,B} (0,1,1,0): PCT 0.9526, contribution +0.0247 → 0.0 / 0 | Bls 0 ⇒ 0; overload with site type (6mer ⇒ 0) used by context++ |
| F20 | No 4-dp rounding of Bls / PCT | tree with 9606:1.60004, miR-30 8mer: Bls 3.20004, PCT 0.17881 (0.1788) → 3.2000 / 0.1787 | `ComputeBranchLengthScore` and `PctFromBranchLength` round as `sprintf("%.4f")` (IEEE-exact "F4"); `-0.0000` ⇒ 0 |
| F21 | A name that is not a leaf still counted | "9606 99999": 2.10004 → 0.0000; "9606 10090 99999": 3.70004 → 3.2000 | counts only species found in the tree; O(n) single walk (was O(n²)) |
| F22 | No argument validation / NA; no Conserved call | Bls −2.0 → finite PCT; NaN → NaN | negative / NaN / ∞ Bls ⇒ `ArgumentOutOfRangeException`, non-finite b ⇒ `ArgumentException`; context++ treats a negative Bls as NA (contribution 0, listed omitted); `IsConservedByBranchLength` + `BranchLengthConservationThresholds.TargetScanHuman7/TargetScanMouse7` |

Reference defect, not reproduced: when the FIRST species of the list is not in the tree, `getBranchLength` loses the path of the next species
("99999 9606 10090" prints 1.6000; the connecting length is 3.2000). The README says unknown species are ignored; the library ignores them.

## Tests

- New `MiRnaAnalyzer_PctTargetScanParity_Tests` (15): 9 kit rows (bin-1 tree embedded), rounding-before-sigmoid, unknown/duplicate species,
  root independence, `getPCT` zero rules, thresholds (human + mouse, boundaries), argument validation.
- Re-locked to perl values: CTX-PCT-003/004/005/006/007/008/009/010 (`MiRnaAnalyzer_TargetPrediction_Tests`), `MiRnaAlgebraicTests` PCT identity,
  `MiRnaCombinatorialTests` PCT grid (4-dp tolerance), `MiRnaPctFuzzTests` (INV-PCT-1/3/5: b = 0 ⇒ 0, 4-dp reference, negative Bls ⇒ throws).

## Boundaries (documented, not limitations)

- `PctConservation` still takes a caller tree + species-with-site list; since WP-3 (F35/F36) both can instead be derived from a multi-species
  UTR alignment (`FindTargetScanSiteGroups`, `AssignTargetScanUtrBin`, `GetTargetScanBinTree`, `ScoreTargetScanConservation`). The b0..b3 row
  defaults to the bundled TargetScan 7.0 table (WP-2, F32); a caller row overrides it.
- `PhylogeneticAnalyzer.ParseNewick` rejects the kit's multi-line tree files (whitespace in front of a `(`: `"(A:1, (B:2,C:3):1);"` throws
  `FormatException`, Biopython 1.88 reads 3 leaves); the bundled trees are single-line, so the library does not depend on it (cross-batch
  request to B16 stands for caller-supplied files).

## Bundled parameter tables (B13 WP-2, F32, 2026-11-25)

- Source: kit `targetscan_70_BL_PCT.zip` (https://www.targetscan.org/vert_70/vert_70_data_download/targetscan_70_BL_PCT.zip, archive dated
  2015-09-03, SHA-256 5ecd0e90…5885), `PCT_parameters/8mer|7mer_m8|7mer_1a_PCT_parameters.txt`, 111 data rows each (same 111 seeds in all three;
  header has an extra `#kmer` column, data rows have 5 fields). Terms: TargetScan FAQ (www.targetscan.org/faqs.Release_7.html, opened 2026-11-25)
  — "free to use TargetScan data or code for your research … as long as you cite at least one TargetScan reference" (BL/PCT: Friedman et al. 2009).
  The kit README has no licence text; the perl headers say "Copyright(c) 2009-2015 Whitehead Institute … All Rights Reserved" (the code is not
  copied; only the parameter values, with the citation, in `MiRnaAnalyzer.TargetScanPctParameters.cs`).
- Keying (`targetscan_70_BL_PCT.pl`): `readAllPCTdata` stores `$familyPlusType2coeff{$f[0]}{$siteType}` with `$f[0]` = the 7-nt seed; the targets'
  family name is mapped to its seed via the miRNA file (`readMiRNAs`, `s/T/U/gi`); `$thisMiRNAhasPCT{seed}` decides PCT vs "NA". The library
  key is the 7-nt seed region (nt 2–8), T→U, upper-cased (perl is case-sensitive; every key is upper case).
- Generated by a script from the kit files (exact decimal text; file names, row counts and SHA-256 in the header). Public API:
  `TryGetTargetScanPctParameters`, `IsTargetScanPctSeedFamily`, `TargetScanPctSeedFamilyCount`; `PctConservation.SigmoidParameters` is now
  optional (`null` ⇒ bundled row for the miRNA's seed); `ContextPlusPlusScore.PctIsNa`.
- Oracle: with the bundled table, the new unit test `KitSample_AllRows_BranchLengthAndPct_WithBundledParameters_EqualTargetScan` reproduces all
  213 Bls, all 161 PCT (28 of them non-zero, i.e. dependent on the row values) and the 52 "NA" rows of the shipped
  `targetscan_70_output.BL_PCT.txt` (exact double equality).

## Species with the site, UTR bin and trees from an alignment (B13 WP-3, F35/F36, 2026-11-25)

- Sources opened: `targetscan_70.zip` (www.targetscan.org/vert_70/vert_70_data_download/, SHA-256 6e23c72a…d749): `targetscan_70.pl` 7.0,
  `README_70.txt`, `UTR_Sequences_sample.txt`, `miR_Family_info_sample.txt`, `targetscan_70_output.txt`; mirror
  raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70.pl (same logic; the official copy has an unused `%pairToDistance`, an
  unused `%GROUP_TYPES_LIST_2_GROUP_TYPE` lookup and splits the same-start/same-end test into equivalent branches); BL/PCT kit
  `targetscan_70_BL_bins.pl`, `UTRs_median_BLs_bins.txt`, `PCT_parameters/Tree.bin_01..10.txt`, `Tree.generic.txt`; Statistics::Lite 3.62
  `median` (apt `libstatistics-lite-perl`).
- Reference run: perl 5.38.2 — both `targetscan_70.pl` copies reproduce `targetscan_70_output.txt` byte for byte from the sample;
  `targetscan_70_BL_bins.pl` reproduces `UTRs_median_BLs_bins.txt` byte for byte; `targetscan_70_BL_PCT.pl` the BL/PCT output (F19–F22).
- F35 (`targetscan_70.pl` port): `FindTargetScanSiteGroups` — T→U, 14-letter 5' mask, gap-tolerant 8mer / 7mer-m8 / 7mer-A1 / 6mer at every
  start, subset removal, MSA + ungapped coordinates, all-vs-all grouping in ordinal key order with REQUIRED_OVERLAP 2 and the script's
  non-merging `groupThisPair`, group type, group species and per-type species lists (`makeListNonRedundant` order), ORF_overlap with the script's
  one-column-right window (copied: it defines the ORF_overlap column and PCT "NA"; constructed perl run: lower-case first letter ⇒ 0,
  lower-case letter right after the site ⇒ 1).
- F36 (`targetscan_70_BL_bins.pl` + end to end): `AssignTargetScanUtrBin` (reference-anchored per-column connecting-subtree length on the generic
  tree, unrounded — shared private `ConnectingSubtreeLength` with `ComputeBranchLengthScore`; median via `StatisticsHelper.Median`; active
  `@BL_thresholds`); the 10 bin trees + generic tree bundled as single-line Newick in the generated `MiRnaAnalyzer.TargetScanTrees.cs` (SHA-256 per
  file); `ScoreTargetScanConservation` = groups → bin → Bls (bin tree, type list) → PCT (bundled row; NA for no family / ORF overlap) → Conserved.
- Oracle (tests `MiRnaAnalyzer_TargetScanAlignment_Tests`, kit files verbatim in the generated `MiRnaTargetScanAlignmentKitData.cs`): 213/213 rows of
  `targetscan_70_output.txt` verbatim; 6/6 bins with the median text (`%.15g`) equal; end to end 213 Bls, 161 PCT (28 non-zero), 52 NA, 0 conserved
  equal to `targetscan_70_output.BL_PCT.txt`; constructed perl runs: 12-row mask/gap/subset/ORF-window case, non-merged groups, an 11-row bin case
  (2.55195311, bin 3).
- Not copied (no effect on valid input; documented in the XML docs): species id "0" read as missing; the bins script reusing the previous gene's
  reference row when a gene has none (error here); the hash-order tie between 7mer-m8 and 7mer-A1 that exists only for seed UUUUUUU (7mer-m8 kept).
