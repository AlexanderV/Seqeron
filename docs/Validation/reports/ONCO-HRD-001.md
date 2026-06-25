# Validation Report: ONCO-HRD-001 — Homologous Recombination Deficiency (HRD) composite genomic-scar score

- **Validated:** 2026-06-24   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CalculateHRDScore(int,int,int)`, `ClassifyHRDStatus(int)`, `DetectHRD(HrdComponents)`, `DetectHRD(IEnumerable<AlleleSpecificSegment>, int tai, int lst)`, `CalculateHrdTaiScore(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)`, `CalculateHrdLstScore(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)`, `DetectHRD(IEnumerable<AlleleSpecificSegment>, ReferenceGenome)` (all-derived). Plus embedded `GRCh38Centromeres`/`GRCh37Centromeres` tables and helpers.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** CLEAN (one by-design residual documented: odd-ploidy TAI branch)

## Context

This re-validation covers the limitations-campaign extension (commits `601cfbdf`, `7bb6dfb2`) that turned ONCO-HRD-001 from a composite-sum-only unit into an **end-to-end** one: `DetectHRD(segments, genome)` now derives all three scarHRD components — HRD-LOH (Abkevich 2012), HRD-TAI (Birkbak 2012 / scarHRD `calc.ai_new`), HRD-LST (Popova 2012 / scarHRD `calc.lst`) — directly from allele-specific copy-number segments, using embedded GRCh38/GRCh37 centromere coordinates from UCSC cytoBand `acen`.

## Stage A — Description

### Sources opened & what they confirm (retrieved this session via web tools)

1. **scarHRD `calc.ai_new.R`** (raw GitHub) — confirmed verbatim the even-ploidy TAI path:
   - `seg <- seg[seg[,4]-seg[,3] >= min.size,]` (min.size = 1e6 drop).
   - even-ploidy AI: `c(0,2)[match(seg[,7]==seg[,8], ...)]` — AI present ⟺ major (col 7) ≠ minor (col 8).
   - p-telomeric: `if(seg[1,'AI']==2 & nrow != 1 & seg[1,4] < chrominfo[i,2]) → AI<-1` (first segment END < centromere start).
   - q-telomeric (quoted verbatim this session): `if(seg[nrow,'AI']==2 & nrow != 1 & seg[nrow,3] > chrominfo[i,3]) → AI<-1` (last segment START > centromere end).
   - whole-chromosome: `if(nrow==1 & AI!=0) → AI<-3` (single imbalanced segment is NOT telomeric).
   - returned TAI = `nrow(seg[seg[,'AI']==1,])` — count of AI==1 segments.
2. **scarHRD `calc.lst.R`** (raw GitHub) — confirmed the `chr.arm='no'` path: p-arm = segments with end ≤ `chrominfo[i,2]`, q-arm = segments with start ≥ `chrominfo[i,3]`; clamp `p.arm[last,4]<-chrominfo[i,2]`, `q.arm[1,3]<-chrominfo[i,3]`; `while(length(n.3mb)>0){ remove <3e6; shrink.seg.ai() }` smoothing; count adjacent pairs both ≥ 10e6 with gap < 3e6.
3. **Popova 2012 (Cancer Res 72(21):5454)** — LST = chromosomal break between adjacent regions each ≥ 10 Mb, obtained after smoothing/filtering < 3 Mb variation. Matches scarHRD.
4. **Birkbak 2012 (Cancer Discov 2(4):366)** — NtAI = number of subchromosomal AI regions extending to the telomere (not crossing the centromere). Matches scarHRD telomeric definition.
5. **Telli 2016 / Stewart 2022** — HRD = unweighted sum LOH+TAI+LST; HR-deficient at ≥ 42 inclusive. (Confirmed in the prior pass; unchanged.)
6. **UCSC cytoBand `acen` (api.genome.ucsc.edu)** — spot-checked 3 chromosomes against the embedded tables (see cross-check below).

### Formula check
- TAI: telomeric-AI count over autosomes, even-ploidy path, with min.size 1 Mb, same-state merge, first/last-segment-vs-centromere geometry, single-segment = whole-chr exclusion — matches `calc.ai_new` exactly.
- LST: per-arm split at centromere with boundary clamp, iterative < 3 Mb smoothing + re-merge, adjacent ≥ 10 Mb pair with < 3 Mb gap — matches `calc.lst` (`chr.arm='no'`) exactly.
- Sum + ≥ 42 cutoff — matches Telli 2016 (unchanged from prior pass).

### Centromere coordinate cross-check (embedded vs UCSC cytoBand acen, retrieved this session)

| Build | Chr | UCSC acen (p11 chromStart / q11 chromEnd) | Embedded table | Match |
|-------|-----|-------------------------------------------|----------------|-------|
| hg38  | 1   | 121,700,000 / 125,100,000 | `new(121_700_000, 125_100_000)` | ✅ exact |
| hg38  | 21  | 10,900,000 / 13,000,000   | `new(10_900_000, 13_000_000)`   | ✅ exact |
| hg19  | 1   | 121,500,000 / 128,900,000 | `new(121_500_000, 128_900_000)` | ✅ exact |

(p-arm end = centromere start = p11 acen chromStart; q-arm start = centromere end = q11 acen chromEnd — exactly the code's convention.)

### Edge-case semantics
- Sub-1 Mb terminal fragment dropped before AI assignment (min.size). Single imbalanced segment = whole-chr (AI=3), not telomeric. Interstitial imbalance not counted. First segment crossing centromere (end ≥ cen start) not telomeric. Sex chromosomes excluded from both TAI and LST (autosome-only centromere table; scarHRD removes chr23/24/X/Y). All sourced and well-defined.

### Findings / divergences
- One **by-design** residual: the odd-ploidy TAI branch of `calc.ai_new` (`seg[,7]+seg[,8]==ploidy & seg[,8]!=0`) is not reproduced because `AlleleSpecificSegment` lacks the ASCAT per-sample ploidy / aberrant-cell-fraction columns it needs. The implemented even-ploidy path (AI ⟺ major≠minor) is the dominant path and matches Birkbak's plain "regions of allelic imbalance". Documented in the TestSpec Assumption Register; not a defect.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `CalculateHrdTaiScore` (2320–2373): drop < 1 Mb, `MergeAdjacentSameAlleleState`, single-segment skip, first `IsAllelicImbalance && End < cen.Start`, last `IsAllelicImbalance && Start > cen.End`. Matches `calc.ai_new` verbatim.
- `CalculateHrdLstScore` (2397–2442) + `CountArmTransitions` (2449–2493) + `SmoothShortSegments` (2500–2542): p/q split (`Start <= cen.Start` / `End >= cen.End`), clamp inner edge, iterative < 3 Mb removal with re-merge + re-clamp, adjacent both-≥10 Mb gap-<3 Mb count. Matches `calc.lst`.
- Centromere tables (2196–2253), `GetCentromeres`/`TryGetCentromere` (2263–2285) autosome-only via `TryGetAutosomeNumber`.
- All-derived `DetectHRD(segments, genome)` (2597–2611): materializes once, `loh=DetectLOH.Score`, `tai=CalculateHrdTaiScore`, `lst=CalculateHrdLstScore`, then `DetectHRD(HrdComponents)` → sum + classify. Realizes scarHRD `sum_HRD0`.
- Constants: `HrdTaiMinSegmentLengthBp=1e6`, `HrdLstSmoothingLengthBp=3e6`, `HrdLstLargeSegmentLengthBp=1e7`, `HrdHighScoreThreshold=42` — all match sources.

### Hand cross-check recomputed vs code (GRCh38 chr1 cen 121.7M/125.1M)

| Case | Segments | Code path result | Expected | Match |
|------|----------|-------------------|----------|-------|
| TAI both arms (M13) | [0,50M]2:1; [50M,121.7M]1:1; [130M,248.9M]2:1 | first end 50M<121.7M → +1; last start 130M>125.1M → +1 | 2 | ✅ |
| TAI first crosses cen (M15) | [0,130M]2:1; [130M,end]1:1 | first end 130M ≥ 121.7M → not telomeric | 0 | ✅ |
| TAI whole-chr (M16) | [0,end]2:1 | single merged segment → skip | 0 | ✅ |
| TAI sub-1Mb (M17) | [0,0.5M]2:1; [0.5M,end]1:1 | 0.5M dropped → one balanced segment | 0 | ✅ |
| LST adjacent large (M18) | [0,40M]2:1; [40M,80M]1:1 | both ≥10M, gap 0<3M | 1 | ✅ |
| LST 5Mb breaks pair (M19) | 40M/5M/65M | no adjacent both-large pair | 0 | ✅ |
| LST 2Mb smoothed (M20) | [0,40M]2:1; [40M,42M]4:0; [42M,90M]1:1 | 2M<3M removed, two ≥10M, gap 2M<3M | 1 | ✅ |
| LST q-arm (M22) | [125.1M,180M]2:1; [180M,248.9M]1:1 | both ≥10M q-arm, gap 0<3M | 1 | ✅ |
| TAI/LST sex chr (S8/S10) | chrX | TryGetCentromere false → skip | 0 | ✅ |
| TAI GRCh37 (S9) | [0,50M]1:1; [130M,end]2:1, hg19 | last start 130M>128.9M → +1 | 1 | ✅ |

### Note (faithful, not a defect)
A segment spanning the centromere (start ≤ cenStart AND end ≥ cenEnd) is selected onto BOTH the p- and q-arm lists in `CalculateHrdLstScore` and clamped on each side. This mirrors scarHRD's independent p/q selection + clamp of a spanning segment, so it is faithful to the reference.

### Variant/delegate consistency
- `DetectHRD(segments, genome)` reuses `DetectLOH` + `CalculateHrdTaiScore` + `CalculateHrdLstScore` and feeds `DetectHRD(HrdComponents)` (INV-8 holds by construction; locked by M24).
- `DetectHRD(segments, tai, lst)` derives LOH only, INV-6 (M11).
- TAI/LST order-independence (P1) via per-chromosome sort.

### Test quality audit
`OncologyAnalyzer_CalculateHRDScore_Tests.cs` — 85 tests, **all pass** (filtered run: Failed 0, Passed 85). Every assertion uses `Is.EqualTo` with an exact sourced value or exact enum; null→throw, empty→0, negative→throw guards covered. No range/tolerance/`no-throw`-only assertions on the maths. M13–M24, S7–S11, P1 lock the TAI/LST/all-derived paths and the GRCh37 table.

### Findings / defects
None. No code or test change required.

## Verdict & follow-ups
- **Stage A: PASS** — TAI (`calc.ai_new`), LST (`calc.lst`), sum/cutoff all match the primary papers and the verbatim scarHRD source; 3/3 centromere coordinates match UCSC cytoBand acen exactly.
- **Stage B: PASS** — code realizes the validated formulas; 10 hand-computed TAI/LST cross-checks reproduce the code's results; full unfiltered suite green **18213 passed / 0 failed**.
- **End-state: CLEAN.** The only residual (odd-ploidy TAI branch needing ASCAT columns absent from `AlleleSpecificSegment`) is by-design and documented in the Assumption Register; the implemented even-ploidy path is the dominant scarHRD path.
- **Follow-ups:** none.
