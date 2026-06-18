# Validation Report: ONCO-LOH-001 — Loss of Heterozygosity (HRD-LOH) detection

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.DetectLOH(IEnumerable<AlleleSpecificSegment>)`,
  `OncologyAnalyzer.CalculateHrdLohScore(...)`, `OncologyAnalyzer.CalculateLOHFraction(..., chromosome)`,
  internal `IsLohSegment`, `MergeAdjacentSameState`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN

---

## Stage A — Description

### Sources opened this session (retrieved live, not trusting the repo Evidence)

| Source | URL | What it confirms |
|--------|-----|------------------|
| Abkevich et al. 2012, full text (PMC open) | https://pmc.ncbi.nlm.nih.gov/articles/PMC3493866/ | HRD score definition; 15 Mb cutoff; three-feature LOH length distribution; whole-chromosome peak |
| Abkevich et al. 2012, PubMed abstract | https://pubmed.ncbi.nlm.nih.gov/23047548/ | "number of these regions"; intermediate-size LOH ↔ BRCA1/2 deficiency |
| scarHRD `calc.hrd.R` (full file) | https://raw.githubusercontent.com/sztup/scarHRD/master/R/calc.hrd.R | Exact counting algorithm + operation order |
| scarHRD `shrink.seg.ai.wrapper.R` | https://raw.githubusercontent.com/sztup/scarHRD/master/R/shrink.seg.ai.wrapper.R | Merge wrapper (per chr) |
| scarHRD `shrink.seg.ai.R` | https://raw.githubusercontent.com/sztup/scarHRD/master/R/shrink.seg.ai.R | Merge criterion: cols 7&8 (A,B) equal on consecutive rows |
| scarHRD `preprocess.hrd.R` | https://raw.githubusercontent.com/sztup/scarHRD/master/R/preprocess.hrd.R | Sex-chr removal; no gap-filling |
| oncoscanR `score_loh` | https://rdrr.io/github/yannchristinat/oncoscanR/man/score_loh.html | Independent corroboration; merge "overlapping or neighbor (1bp)" |

### Formula check (against the sources, verbatim)

1. **HRD-LOH score (Abkevich, PMC, verbatim):** *"the homologous recombination deficiency (HRD) score
   was defined as the number of these regions"* where *these regions* = *"LOH regions >15 Mb, but less
   than a whole chromosome in length."* → count of LOH regions, lower bound 15 Mb (strict), upper bound
   < whole chromosome.
2. **15 Mb cutoff (Abkevich, PMC, verbatim):** *"The value of 15 Mb was selected somewhat arbitrarily,
   but further analysis showed that the exact value of this cut-off does not have significant impact on
   the results."*
3. **LOH-segment criterion (scarHRD `calc.hrd.R`, verbatim):**
   `segLOH <- segSamp[segSamp[,nB] == 0 & segSamp[,nA] != 0,,drop=F]` → minor (B) allele CN = 0 AND major
   (A) allele CN ≠ 0. This is **allelic loss including copy-neutral LOH**, excluding homozygous deletion.
4. **Strict size filter (scarHRD, verbatim):** `segLOH <- segLOH[segLOH[,4]-segLOH[,3] > sizelimit1,,drop=F]`,
   `sizelimitLOH = 15e6`. Length = end − start. Strict `>`. oncoscanR independently: *"larger than 15Mb"*.
5. **Whole-chromosome exclusion (scarHRD, verbatim):**
   `for(j in unique(segSamp[,2])){ if(all(segSamp[segSamp[,2]==j,nB]==0)) chrDel <- c(chrDel,j) }`
   then `segLOH <- segLOH[!segLOH[,2] %in% chrDel,,drop=F]` → a chromosome where **every** segment has
   minor==0 is "global LOH" and is dropped. Matches Abkevich *"less than a whole chromosome"*.
6. **Exact scarHRD operation order (retrieved verbatim):** compute `chrDel` on raw segments → cap
   `nA>1 → 1` → `shrink.seg.ai.wrapper` merge → LOH filter → size filter → chrDel exclusion → `nrow`.

### Edge-case semantics (all sourced)

- Length exactly 15 Mb → **not** counted (strict `>`). [scarHRD]
- Homozygous deletion (minor==0 & major==0) → **not** LOH (`nA != 0`). [scarHRD]
- Heterozygous retained (minor≠0) → **not** LOH (`nB == 0`). [scarHRD]
- Whole-chromosome LOH → excluded (`chrDel`). [scarHRD / Abkevich]
- Copy-neutral LOH (minor==0, major==2, total==2) → **is** LOH (criterion is allele-specific, not total CN). [scarHRD]

### Independent cross-check (hand computation of the Evidence 7-segment dataset)

Applying the scarHRD rule by hand: chr1 (20 Mb LOH, partial chr) → **count**; chr2 (10 Mb LOH, single seg
⇒ whole-chr LOH and <15 Mb) → 0; chr3 (16 Mb, single seg ⇒ whole-chr LOH) → 0; chr4 (homdel) → 0; chr5
(15 Mb LOH not >15 Mb, plus het) → 0. **Total = 1.** Independently matches the Evidence and test M1.

### Findings / divergences (Stage A)

- **PASS.** The biology/maths of the description is correct and matches the primary paper plus two
  reference implementations. The Evidence's verbatim quotes were all reproduced from the live sources.
- **Description fix applied:** the by-area Registry row (ALGORITHMS_CHECKLIST_V2.md §ONCO-LOH-001) still
  listed the obsolete signature `DetectLOH(tumorVcf, normalVcf)`. Stage A establishes (Abkevich/scarHRD/
  oncoscanR all operate on allele-specific CN segments) that the correct canonical input is allele-specific
  segments; raw VCF/segmentation is upstream (ONCO-CNA-001). Updated the Registry row to the segment-based
  signatures + a validation note. (Not a code change — a description correction, per protocol.)

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:1864–2109`
(`AlleleSpecificSegment`, `IsLohSegment`, `DetectLOH`, `CalculateHrdLohScore`, `CalculateLOHFraction`,
`GroupValidatedByChromosome`, `MergeAdjacentSameState`, `ValidateSegment`).

### Formula realised correctly?

- LOH criterion `IsLohSegment` = `MinorCopyNumber == 0 && MajorCopyNumber != 0` — exactly scarHRD
  `nB==0 & nA!=0` (line 1916). ✓
- Size filter `merged.Length > HrdLohMinRegionLengthBp` with `HrdLohMinRegionLengthBp = 15_000_000`,
  `Length = End - Start`, strict `>` (lines 1874, 1894, 1955). ✓ (matches `> sizelimit1`, length = col4−col3).
- Whole-chromosome exclusion `group.All(s => s.MinorCopyNumber == 0)` computed on the **raw** per-chromosome
  segments before merging (lines 1946–1950) — exactly scarHRD's `chrDel` (which is also computed on
  `segSamp` before the cap/merge). ✓
- Score = number of surviving regions (`regions.Count`, line 1962) = scarHRD `nrow(segLOH)`. ✓
- Copy-neutral LOH handled (criterion is allele-specific, never `total<2`). ✓ (now locked by a test).

### Cross-verification table recomputed vs code (run + hand-trace)

| Case | Source-derived expected | Code result | Match |
|------|-------------------------|-------------|-------|
| Evidence 7-segment dataset | 1 (hand-computed via scarHRD rule) | 1 | ✓ |
| 18 Mb copy-neutral LOH (major=2) on partial chr | 1 | 1 | ✓ |
| Exactly 15 Mb LOH | 0 (strict `>`) | 0 | ✓ |
| 30 Mb homozygous deletion | 0 (`major!=0`) | 0 | ✓ |
| 40 Mb het retained + 20 Mb LOH | 1 | 1 | ✓ |
| Whole-chr LOH (all minor=0) | 0 (`chrDel`) | 0 | ✓ |
| Same LOH + het ⇒ partial chr | 1 | 1 | ✓ |
| Two adjacent 8 Mb LOH → 16 Mb | 1 (merge) | 1 | ✓ |
| LOH fraction 20M/60M | 0.3333… | 0.3333… | ✓ |

### Variant/delegate consistency

`CalculateHrdLohScore` = `DetectLOH(...).Score` (line 1976) — consistent by construction. ✓

### Numerical robustness

Coordinates `long`; fraction is `(double)lohLength/totalLength` with a `totalLength==0 → 0` guard (no
div-by-zero). No overflow on chromosome-scale bp (≤ ~2.5e8). ✓

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** M3 (exactly 15 Mb → 0) fails a `>=` impl; M6/M7 contrast pair pins the
  `chrDel` logic; S1 fails without the merge step; M4 fails a "total<2" impl; the new **M2b copy-neutral
  LOH** test fails any impl keying on total copy number. Each asserts the exact sourced value.
- **No green-washing:** all assertions use exact equality (`Is.EqualTo`) or exact fractions
  (`Within(1e-10)`); no Greater/AtLeast/Contains substituted for known exact values; no skips/ignores.
  (The two `Is.InRange(0,1)` assertions in `CalculateLOHFraction_MixedChromosome` are *accompanied* by the
  exact `Is.EqualTo(1/3)` / `Is.EqualTo(0.3)` assertions, so the invariant check does not weaken the exact check.)
- **Coverage:** all three public methods + both error paths (null, non-positive length, negative CN) +
  whole-chr / homdel / het / strict-boundary / merge / copy-neutral / empty / unknown-chr / order-invariance.
- **Gap closed this session:** added `DetectLOH_CopyNeutralLoh_IsCounted` (M2b) to explicitly lock the
  source distinction "LOH = minor==0 & major!=0, independent of total CN."
- **Honest green:** full unfiltered suite **6633 passed / 0 failed** (run 5×; one transient non-LOH failure
  on an early run was not reproducible in 4 subsequent full runs and the LOH suite passes in isolation,
  19/19). Build 0 errors; the only warnings are pre-existing and unrelated to this unit.

### Findings / defects (Stage B)

- **No code defect.** The implementation faithfully realises the validated Abkevich/scarHRD algorithm.
- **NOTE 1 (merge semantics, documented divergence — not a defect):** the repo's `MergeAdjacentSameState`
  merges on **LOH state** (`IsLohSegment` boolean) and requires **gap ≤ 1 bp**, whereas scarHRD's
  `shrink.seg.ai` merges consecutive rows whose **(capped-major, minor)** pair is equal, with **no gap
  check**. For LOH counting these are equivalent on contiguous segmentation (scarHRD's intended input),
  because (a) the LOH/non-LOH partition is what drives the count, and (b) real CN segmentations tile each
  chromosome contiguously so no gaps arise. They diverge only on the artificial case of two same-state LOH
  segments separated by a >1 bp gap on a partial chromosome (scarHRD would merge them; the repo would not).
  The repo's behaviour matches **oncoscanR** ("merge overlapping or neighbor LOH segments at 1 bp"), which
  is an authoritative source explicitly cited in the code, so this is a defensible, documented design
  choice — recorded in the TestSpec Assumption register.
- **NOTE 2 (`CalculateLOHFraction` is a definitional API choice):** Abkevich/scarHRD define only the
  *count*, not a per-chromosome fraction. The fraction (length-weighted LOH burden, no 15 Mb filter, no
  whole-chr exclusion) is a Registry-mandated quantity satisfying `0 ≤ f ≤ 1`; its source-backed core is
  the LOH-segment criterion. Documented as Assumption #1 in the TestSpec/Evidence. Acceptable.

---

## Verdict & follow-ups

- **Stage A: PASS** — description independently confirmed against Abkevich 2012 (PMC full text) and two
  reference implementations (scarHRD, oncoscanR). One stale Registry signature corrected.
- **Stage B: PASS-WITH-NOTES** — code matches the validated algorithm exactly; two documented,
  source-defensible divergences (oncoscanR-style merge; the definitional LOH-fraction API), neither a defect.
- **End-state: ✅ CLEAN** — no defect; one test added to lock copy-neutral LOH; description corrected;
  `dotnet build` 0 errors and the FULL unfiltered suite Failed: 0.
- **Test-quality gate: PASS** — sourced expectations, no green-washing, full branch coverage, honest green.
- No findings logged in FINDINGS_REGISTER (no defect); a docs-correction note is recorded there for traceability.
