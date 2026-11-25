# Test Specification: MIRNA-CLEAVAGE-001

**Test Unit ID:** MIRNA-CLEAVAGE-001
**Area:** MiRNA
**Algorithm:** Drosha/Dicer Cleavage-Site Prediction
**Status:** ☑ Validated — B13 re-review (2026-11-25): Stage A PASS-WITH-NOTES / Stage B FAIL→FIXED (F27–F30)
**Last Updated:** 2026-11-25

> The measuring rules are now applied to the folded stem: 11 base pairs from the basal junction (Han 2006),
> 22 nt from the 5' (default) or 3' end (Park 2011), and the opposite-arm product from the 2-nt 3'-overhang
> partner rule (miRDeep2 `find_star`). The former linear "star" (5p shifted by 2 nt on the same arm) and the
> nucleotide (not base-pair) Drosha count are gone. See `docs/Validation/reports/MIRNA-CLEAVAGE-001.md`.

---

## 1. Evidence Summary

| # | Source | Confirms |
|---|--------|----------|
| 1 | Han et al. (2006), *Cell* 125:887 | Drosha cleaves ~11 bp from the basal stem–ssRNA junction. |
| 2 | Auyeung et al. (2013), *Cell* 152:844 (PDF opened) | "an 11 bp basal stem was preferred"; positions numbered from the Drosha sites (− upstream of the 5p cut, + downstream of the 3p cut); basal U at −14 / G at −13; CNNC enriched at +16..+18; cleavage sites inferred from RNAfold (hairpin ± 20 nt) with a 2-nt 3' overhang. Apical UGU/GUG only as loop enrichment (no positional rule). |
| 3 | Park et al. (2011), *Nature* 475:201 | 5' counting (~22 nt from the 5' end); Dicer also anchors the 3' end (3' counting). |
| 4 | Lee et al. (2003) / Han 2006 | RNase III leaves a 2-nt 3' overhang. |
| 5 | miRDeep2 `find_star` (rajewsky-lab/mirdeep2 `src/miRDeep2_core_algorithm.pl`, opened) | Star from mature: start = bp(end − 2 − off) − off, end = bp(start + off) + off + 2. |
| 6 | Gu et al. (2012), *Cell* 151:900 (PMC3499986) | Loop counting = precision condition (cut 2 nt from a loop/bulge), not a placement rule. |
| 7 | Zeng, Yi & Cullen (2005), *EMBO J* 24:138 (PMC544904, opened) | "From the junction of the loop and the adjacent stem, Drosha then cleaves approximately two helical RNA turns" (~22 bp along the stem, "allowance for distortions"): the opt-in `DroshaCountingRule.ApicalJunction` ruler (test DD2b). |
| 8 | miRBase hairpins MI0000077, MI0000060, MI0000070, MI0000088, MI0000681 + UCSC hg38 flanks + ViennaRNA 2.8 RNA.fold | Real pri-miRNA oracle (see §4). |

## 2. Canonical Method(s)

`MiRnaAnalyzer.PredictDroshaDicerCleavage(string sequence, int basalJunction)` (canonical MFE fold, 5' counting) and
`PredictDroshaDicerCleavage(string sequence, int basalJunction, string? dotBracket, DicerCountingRule dicerRule = FivePrimeCounting)` → `DroshaDicerCleavage?`

- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
- **Tests:** `Unit/Annotation/MiRnaAnalyzer_PreMiRna_Tests.cs` (DD1–DD14, DD2b), `MiRnaAnalyzer_StructureCleavageKernels_Tests.cs`, `MiRnaAnalyzer_ContextScoreAssembly_Tests.cs` (CNNC bound), heavy tier `Fuzzing/MiRnaCleavageFuzzTests.cs`, `Properties/MiRnaProperties.cs`, `Metamorphic/MiRnaMetamorphicTests.cs`, `Algebraic/MiRnaAlgebraicTests.cs`, `Combinatorial/MiRnaCombinatorialTests.cs`.

## 3. Contract / Invariants

- Input normalised (upper, T→U); structure = caller dot-bracket (must be balanced and same length, else `ArgumentException`) or `RnaSecondaryStructure.CalculateMfeStructure`.
- `basalJunction` must be the 5' base of a pair; the unbranched stem from it is walked pair by pair (interior loops and bulges crossed, a multiloop ends it).
- `DroshaCountingRule.ApicalJunction` (opt-in 5-argument overload): `DroshaCut5'` = 5' base of the 22nd pair counted from the apical junction pair (innermost pair of the unbranched stem); needs ≥ 22 pairs. Apical UGU/GUG stays unflagged (BLOCKED, no reproducible window: see B13.md WP-8).
- `DroshaCut5' = 5' base of pair 12` (11 bp below it); `DroshaCut3' = partner(DroshaCut5') + 2`.
- 5' counting: 5p = [DroshaCut5', +21]; 3p = [find_star start, DroshaCut3']. 3' counting: 3p = [DroshaCut3' − 21, DroshaCut3']; 5p end = partner(3p start + off) + off + 2.
- `MatureEnd < StarStart` (opposite arms). `HasBasalUgMotif`: U at DroshaCut5' − 14, G at − 13. `HasCnncMotif`: C..C starting at DroshaCut3' + 16..18.
- `null`: null/empty, junction ∉ [0,len), junction not a 5' pair base, < 12 stem pairs, a product off the sequence or not on opposite arms. No LimitationPolicy guard (nothing approximate is returned).

## 4. Cross-check / Differential Oracle

Python re-implementation on RNAfold 2.8 structures of the miRBase hairpin ± 30 nt hg38 flank (mir-30a ± 20 nt), junction = pair 11 bp below the annotated 5p start. Predicted − annotated (hairpin 1-based), 5' counting, 5p start/end | 3p start/end:

| Hairpin | 5p | 3p | 3' counting 5p end / 3p start | UG | CNNC |
|---|---|---|---|---|---|
| hsa-mir-21 | 0 / −2 | +1 / −1 | −1 / 0 | – | – |
| hsa-let-7a-1 | +1 / 0 | −1 / 0 | 0 / −1 | – | – |
| hsa-mir-16-1 | 0 / 0 | 0 / +1 | −1 / +1 | ✓ | ✓ |
| hsa-mir-30a | 0 / 0 | +2 / 0 | +2 / 0 | ✓ | ✓ |
| hsa-mir-155 | 0 / −2 | +2 / 0 | 0 / 0 | – | – |

Canonical MFE equals RNAfold for 4/5 (mir-16-1 differs: B12 energy discrepancy, see report).

## 5. Validation Checklist

- [x] Stage A: sources opened (Auyeung PDF, miRDeep2 code, miRBase pages, PMC3499986); constants confirmed.
- [x] Stage B: fold-based implementation; oracle numbers locked in DD1–DD3; hand-derived synthetic geometry DD4–DD13.
