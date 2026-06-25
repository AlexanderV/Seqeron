# Test Specification: META-CHECKM-001

**Test Unit ID:** META-CHECKM-001
**Area:** Metagenomics
**Algorithm:** CheckM Marker-Gene Completeness/Contamination
**Status:** ☑ Validated — Stage A PASS / Stage B PASS / CLEAN (2026-06-25)
**Last Updated:** 2026-06-25

> Independently re-validated under the two-stage protocol. The CheckM completeness/contamination
> formula is implemented verbatim from the reference `markerSets.py` (`MarkerSet.genomeCheck`),
> reconfirmed by hand (83.33% / 11.11% on a synthetic bin), and the real-marker detection path
> was cross-checked against the reference HMMER engine via pyhmmer 0.12.1. See
> `docs/Validation/reports/META-CHECKM-001.md`.

---

## 1. Evidence Summary

| # | Source | What it establishes |
|---|--------|---------------------|
| 1 | Parks et al. (2015) *Genome Res* 25:1043 — CheckM | completeness = unique SCGs present / expected; contamination = SCGs in multiple copies; collocated marker SETS grouping + per-set averaging |
| 2 | Ecogenomics/CheckM `checkm/markerSets.py` `MarkerSet.genomeCheck` | exact arithmetic: `comp += present/|s|`, `cont += multiCopy/|s|`, `percComp = 100·comp/|M|`, `percCont = 100·cont/|M|`; `multiCopy += N−1` for `N>1` |
| 3 | Parks et al. (2018) *Nat Biotechnol* 36:996 (GTDB) + GTDB r214 Pfam HMM listing | bac120 = 120 markers (6 Pfam + 114 TIGRFAM), ar122 = 122 (35 Pfam + 87 TIGRFAM); bundled Pfam accessions verified |
| 4 | Pfam = CC0 (InterPro docs) | bundled HMMs are public-domain redistributable |
| 5 | pyhmmer 0.12.1 hmmsearch (this session) | PF00410→uS8 (166.58 bits, GA 24.0), PF01025→GrpE (163.42 bits, GA 25.8); specificity confirmed |

## 2. Canonical Method(s)

`EstimateBinQualityFromMarkerCounts`, `EstimateBinQualityFromMarkers`, `DetectMarkers`,
`LoadBundledRibosomalMarkerHmms`, `LoadBundledBacterialMarkerHmms`, `LoadBundledArchaealMarkerHmms`,
`Bundled{Ribosomal,Bacterial,Archaeal}MarkerSets`, `LoadMarkerHmms`

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_MarkerGeneQuality_Tests.cs`

## 3. Contract / Invariants

- Completeness = 100·(1/|M|)·Σ_{s∈M} present_s/|s|; Contamination = 100·(1/|M|)·Σ_{s∈M} (Σ_{g∈s} N_g−1)/|s|.
- 0 ≤ completeness ≤ 100; contamination ≥ 0 (may exceed 100 under heavy duplication).
- A multi-copy marker counts once toward completeness and N−1 toward contamination.
- Empty marker set excluded from |M|; |M| = 0 ⇒ 0/0 (no div-by-zero). Deterministic.
- A marker is "present" when its protein Plan7 Viterbi bit score ≥ the Pfam GA1 gathering threshold.

## 4. Cross-check / Differential Oracle

- **Formula oracle:** CheckM `markerSets.py` `MarkerSet.genomeCheck`; hand computation gives
  83.3333% / 11.1111% on `M={{A,B},{C,D,E},{F}}`, copies A1 B0 C2 D1 E1 F1.
- **Detection oracle:** pyhmmer 0.12.1 hmmsearch on the bundled HMMs vs E. coli uS8 (P0A7W7) and
  GrpE (P09372).

## 5. Validation Checklist (restored to ☑)

- [x] Stage A: sources retrieved; formula matches markerSets.py + the 2015 paper; hand-derived 83.33%/11.11%.
- [x] Stage B: code matches the formula; pyhmmer confirms detection; cross-check table recomputed vs code.
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [x] `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the `docs/checklists/*.md`.
