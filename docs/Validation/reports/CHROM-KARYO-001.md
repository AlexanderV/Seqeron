# Validation Report: CHROM-KARYO-001 — Karyotype Analysis

- **Validated:** 2026-06-24   **Area:** Chromosome Analysis
- **Canonical method(s):** `ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes, expectedPloidyLevel)`, `ChromosomeAnalyzer.DetectPloidy(normalizedDepths, expectedDiploidDepth)` (`src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

> **Scope clarification.** The session prompt framed this unit around the Levan, Fredga & Sandberg (1964) arm-ratio / centromeric-index chromosome classification (metacentric / submetacentric / subtelocentric / telocentric). That model is **not** CHROM-KARYO-001. Per `ALGORITHMS_CHECKLIST_V2.md` lines 73–74 the Levan classification is a separate unit, **CHROM-CENT-001** (`ChromosomeAnalyzer_Centromere_Tests.cs`, `CalculateArmRatio`/`ClassifyChromosomeByArmRatio`/`AnalyzeCentromere`). CHROM-KARYO-001 is the **karyotype / ploidy / aneuploidy** unit: `AnalyzeKaryotype` + `DetectPloidy`, sourced to Wikipedia Karyotype / Ploidy / Aneuploidy. This report validates the actual CHROM-KARYO-001 unit as defined by its TestSpec, Evidence doc, checklist block, and test file.

## Stage A — Description

### Sources opened
- **Wikipedia — Aneuploidy** (fetched 2026-06-24). Confirms the nomenclature table by **absolute copy count**: Nullisomy=0, Monosomy=1, Disomy=2, Trisomy=3, Tetrasomy=4, Pentasomy=5. Explicitly states disomy is the *normal* condition for diploids but is aneuploid "in organisms that normally have three or more copies" — i.e. terminology is absolute, not relative to expected ploidy.
- **Wikipedia — Ploidy** (fetched 2026-06-24). Confirms haploid=1 set, diploid=2 sets, triploid=3, tetraploid=4; human 2n=46, n=23.

### Formula / model check
1. **Aneuploidy nomenclature** (`GetAneuploidyTerm`, lines 186–195): maps 0→Nullisomy, 1→Monosomy, 2→Disomy, 3→Trisomy, 4→Tetrasomy, 5→Pentasomy, ≥6→"Polysomy (N copies)". Matches the Wikipedia Aneuploidy table exactly. Absolute-count basis confirmed by source.
2. **Aneuploidy trigger** (`AnalyzeKaryotype`, lines 161–168): a chromosome group is abnormal iff `count != expectedPloidyLevel`. Correct generalisation of euploidy/aneuploidy (a group deviating from the organism's set count is aneuploid). In tetraploid context (expected=4), 2 copies → "Disomy", 3 copies → "Trisomy" — absolute terms, matching the source.
3. **Summary fields**: TotalGenomeSize = Σ lengths; MeanChromosomeLength = total/count; sex/autosome split by `IsSexChromosome`. All standard karyotype summary quantities (Wikipedia Karyotype).
4. **DetectPloidy** (lines 217–241): ratio = median(depths)/expectedDiploidDepth; ploidy = round(ratio×2) clamped to [1,8]; confidence = max(0, 1 − |ratio×2 − ploidy|×2). This is a depth-ratio heuristic (not a literature formula but a reasonable engineering model); diploid ratio≈1, tetraploid ratio≈2 agree with the ploidy definition. Median is a *true* median (averages two middle elements for even counts, lines 227–229).

### Edge-case semantics
- Empty chromosomes → empty karyotype, no aneuploidy (DD1) — defined and tested.
- Empty depths → (2, 0): default diploid, zero confidence (DD2) — defined and tested.
- Nullisomy (0 copies) is unreachable via `GroupBy` (DD4) — absent chromosomes form no group; term mapped for completeness. Documented, correct.
- Ploidy clamp [1,8] (DD3) — documented engineering bound.

### Independent cross-check (hand)
DetectPloidy, recomputed by hand and matching the documented examples:

| ratio | round(ratio×2) | clamp | conf = max(0, 1−\|ratio×2−ploidy\|×2) |
|-------|----------------|-------|----------------------------------------|
| 1.0 | 2 | 2 | 1.0 |
| 2.0 | 4 | 4 | 1.0 |
| 0.5 | 1 | 1 | 1.0 |
| 1.5 | 3 | 3 | 1.0 |
| 1.2 | 2 | 2 | 1 − 0.4×2 = 0.2 |
| 10.0 | 20 | 8 | max(0, 1 − 12×2) = 0 |
| 0.1 | 0 | 1 | max(0, 1 − 0.8×2) = 0 |

### Findings / divergences
None affecting CHROM-KARYO-001. Note: confidence uses the **clamped** ploidy as the reference, so for out-of-range depths the |·| term is the distance to the clamp (giving 0 confidence) rather than to the raw round — internally consistent and explicitly tested (`ExtremeHighDepth`/`ExtremeLowDepth`).

## Stage B — Implementation

### Code path reviewed
- `AnalyzeKaryotype` — `ChromosomeAnalyzer.cs:136–180`
- `GetAneuploidyTerm` — `ChromosomeAnalyzer.cs:186–195`
- `GetChromosomeBaseName` — `ChromosomeAnalyzer.cs:200–212` (strips trailing numeric `_N` copy suffix only)
- `DetectPloidy` — `ChromosomeAnalyzer.cs:217–241`
- `Karyotype` record — `ChromosomeAnalyzer.cs:45–53`

### Formula realised correctly?
Yes. The aneuploidy term switch matches the sourced table verbatim; the trigger is `count != expectedPloidyLevel`; DetectPloidy implements the documented ratio/round/clamp/confidence steps with a true median. Verified by the seven hand computations above against the test assertions.

### Cross-verification vs code
All 36 unit tests in `ChromosomeAnalyzer_Karyotype_Tests.cs` assert exact values (e.g. `Abnormalities[0] == "Disomy chr1"`, confidence `0.2 .Within(1e-10)`, ploidy/confidence pairs) and pass. The hand-computed DetectPloidy table reproduces the in-test comments exactly.

### Variant/delegate consistency
TestSpec lists no wrappers/delegates. A sibling method `IdentifyWholeChromosomeAneuploidy` (lines 881–917) uses the same term vocabulary but **hardcodes a diploid baseline** (`CopyNumber != 2`, no `2 => "Disomy"` case) — this is correct *for that method's diploid-only contract* and belongs to a different unit (DetectAneuploidy / CHROM-ANEU-001), so it is not a divergence within CHROM-KARYO-001.

### Test quality audit
Real, deterministic, value-exact assertions covering all M/S/C cases and both invariant sets (Total = Autosomes + Sex; HasAneuploidy ↔ non-empty Abnormalities; Ploidy∈[1,8]; Confidence∈[0,1]). Edge cases (empty inputs, single value, clamps, between-ploidy, absolute-terminology in tetraploid context) are all covered.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN.** No code changed.
- Cross-check categories independently confirmed: aneuploidy nomenclature by absolute copy count (Disomy=2 … Pentasomy=5) per Wikipedia Aneuploidy; ploidy/depth-ratio mapping per Wikipedia Ploidy; seven DetectPloidy ratio/ploidy/confidence values hand-recomputed.
- Process note: the prompt's Levan arm-ratio framing applies to **CHROM-CENT-001**, not this unit; validate that separately.
