# Validation Report: CHROM-KARYO-001 — Karyotype Analysis

- **Validated:** 2026-06-12   **Area:** Chromosome Analysis
- **Canonical method(s):** `ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes, expectedPloidyLevel=2)`, `ChromosomeAnalyzer.DetectPloidy(normalizedDepths, expectedDiploidDepth=1.0)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

Source file: `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs` (lines 42–53, 131–242)
Test file: `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Karyotype_Tests.cs` (37 tests)

---

## Stage A — Description

### Sources opened & what they confirm

1. **Wikipedia — Ploidy** (https://en.wikipedia.org/wiki/Ploidy)
   - "Ploidy is the number of complete sets of chromosomes in a cell." Haploid/monoploid = 1 set, diploid = 2, triploid = 3, tetraploid = 4.
   - Human: **2n = 46**, **n = 23** (diploid reference confirmed).
   - Euploidy = chromosome number is an exact multiple of the haploid number; aneuploidy = one or more individual chromosomes absent or present in more than usual copies.

2. **Wikipedia — Aneuploidy** (https://en.wikipedia.org/wiki/Aneuploidy)
   - Nomenclature is by **absolute copy count**: Nullisomy 0, Monosomy 1, Disomy 2, Trisomy 3, Tetrasomy 4, Pentasomy 5.
   - "Sex chromosome tetrasomy and pentasomy have been reported in humans, including XXXX, XXXY, XXYY, XXXXX, XXXXY, and XYYYY." → terms for 4 and 5 copies are real and standard.

3. **Wikipedia — Karyotype** (https://en.wikipedia.org/wiki/Karyotype)
   - Formula notation: female 46,XX; male 46,XY. Down syndrome (trisomy 21) = 47,XX,+21 / 47,XY,+21.
   - Normal human = 22 autosome pairs + 1 sex-chromosome pair = 46.

### Definition check

- **Ploidy = count of complete chromosome sets** — matches spec/evidence exactly.
- **Aneuploidy term = function of absolute per-chromosome copy count** — matches the implementation's `GetAneuploidyTerm` switch (0→Nullisomy … 5→Pentasomy) exactly. The map is verbatim correct against the Aneuploidy article.
- **Disomy (2 copies)** is normal for a diploid (never aneuploid when expected=2) but IS an abnormality when the expected ploidy is higher — consistent with the absolute-count convention (DD5 in the evidence doc).

### Worked example (hand computation)

- Diploid human-style set {chr1×2, chr2×2, chrX, chrY}, expected ploidy 2: each autosome group has count 2 = expected → no aneuploidy; total 6, autosomes 4, sex 2. ✓
- Trisomy 21 {chr21×3}, expected 2: group count 3 ≠ 2 → absolute term for 3 = "Trisomy" → "Trisomy chr21". ✓
- Tetraploid context {chr1×3}, expected 4: count 3 ≠ 4 → absolute term for 3 = "Trisomy" (NOT "Monosomy"). Confirms absolute, not relative, terminology. ✓
- `DetectPloidy`: depth ratio r = median/expectedDiploidDepth; ploidy = round(2r); confidence = max(0, 1 − 2·|2r − ploidy|). For r=1.0 → ploidy 2, conf 1.0; r=2.0 → 4, conf 1.0; r=0.5 → 1, conf 1.0; r=1.2 → round(2.4)=2, conf=1−2·0.4=0.2. All match.

### Findings / divergences

None at the description level. The "ploidy from read depth" model (ploidy ≈ round(2 × depth-ratio relative to the diploid baseline) with a triangular confidence around integer ploidy) is a reasonable, standard depth-based estimator; it is not from a single canonical paper but is mathematically sound and matched by the tests.

---

## Stage B — Implementation

### Code path reviewed

`AnalyzeKaryotype` (lines 136–180), `GetAneuploidyTerm` (186–195), `GetChromosomeBaseName` (200–211), `DetectPloidy` (216–240).

### Formula realised correctly?

- **Set counting → ploidy/aneuploidy:** autosomes are grouped by base name (`GroupBy(GetChromosomeBaseName)`); each group's `Count()` is compared to `expectedPloidyLevel`; any mismatch sets `HasAneuploidy` and adds `"{term} {group}"` where `term = GetAneuploidyTerm(count)`. This faithfully implements set-counting with absolute-copy-count terminology. ✓
- **Chromosome number total:** `TotalChromosomes = chromList.Count`; partitioned into `AutosomeCount` and `SexChromosomes` by the `IsSexChromosome` flag. ✓
- **Genome size / mean:** `TotalGenomeSize = Σ length`; `MeanChromosomeLength = total / count` (double). ✓
- **DetectPloidy:** true median (averages the two middle elements for even counts), ratio, `round(2r)`, clamp [1,8], `confidence = max(0, 1 − 2·|2r − ploidy|)`. ✓

### Cross-verification table (recomputed vs code, all PASS)

| Input | Expected (sourced) | Code result |
|-------|--------------------|-------------|
| {chr1×2, chr2×2, chrX, chrY}, exp 2 | no aneuploidy; 6/4/2 | ✓ |
| {chr21×3}, exp 2 | "Trisomy chr21" | ✓ |
| {chr1×1}, exp 2 | "Monosomy chr1" | ✓ |
| {chr18×4}, exp 2 | "Tetrasomy" (not Trisomy) | ✓ |
| {chr21×5}, exp 2 | "Pentasomy" (not Trisomy) | ✓ |
| {chr1×3}, exp 4 | "Trisomy" (absolute, not Monosomy) | ✓ |
| {chr1×2}, exp 4 | "Disomy chr1" | ✓ |
| depth 1.0×100, exp 1.0 | ploidy 2, conf 1.0 | ✓ |
| depth 2.0×100 | ploidy 4, conf 1.0 | ✓ |
| depth 0.5×100 | ploidy 1, conf 1.0 | ✓ |
| depth 1.5×100 | ploidy 3, conf 1.0 | ✓ |
| depth 1.2×100 | ploidy 2, conf 0.2 | ✓ |
| depth 10×100 | ploidy 8 (clamped), conf 0 | ✓ |
| depth 0.1×100 | ploidy 1 (clamped), conf 0 | ✓ |
| empty depths | (2, 0) | ✓ |

### Invariants verified

- `TotalChromosomes == AutosomeCount + SexChromosomes.Count` — holds (simple partition). ✓
- `TotalGenomeSize == Σ length`; `MeanChromosomeLength == total/count`. ✓
- `HasAneuploidy ⇔ Abnormalities.Count > 0`. ✓
- `PloidyLevel ∈ [1,8]` (explicit clamp); `Confidence ∈ [0,1]` (frac ≥ 0 ⇒ conf ≤ 1; `Math.Max(0, …)` ⇒ conf ≥ 0). ✓

### Test quality audit

37 tests assert exact sourced values (exact term strings, exact counts, exact confidence 0.2 within 1e-10), cover every Stage-A edge case (empty, single chromosome/monosomy, ratios 0.5/1.0/2.0, clamp high/low, between-ploidy, custom expected depth, tetra/pentasomy, tetraploid absolute terminology, disomy in non-diploid context), and check the invariants with parameterised cases. Tests are deterministic and not tautological.

### Findings / notes (non-defects)

1. **No ISCN formula string emitted.** `Karyotype` reports counts + ploidy + a list of textual abnormality descriptions (e.g. "Trisomy chr21"); it does not assemble the compact ISCN string such as "47,XX,+21". This is within the documented scope of the spec/evidence (output = counts and per-group terms) and is internally consistent — a reporting choice, not a correctness defect.
2. **`GetChromosomeBaseName` strips only numeric `_N` suffixes.** Its comment mentions "chr1a, chr1b" but the code only handles `_`-delimited integer suffixes. No spec/test relies on letter suffixes; the documented copy-naming convention is `chrN_k`. Cosmetic comment/behaviour mismatch only.

These are noted as PASS-WITH-NOTES and warrant no code change under this protocol (no defect, no failing assertion, no divergence from the validated description).

---

## Verdict & follow-ups

- **Stage A: PASS** — ploidy/aneuploidy/karyotype definitions and the diploid human reference (2n=46) confirmed against Wikipedia Ploidy/Aneuploidy/Karyotype; aneuploidy nomenclature is absolute-copy-count and matches the implementation map exactly.
- **Stage B: PASS-WITH-NOTES** — implementation realises set-counting ploidy, total chromosome number, genome size/mean, and absolute-copy-count aneuploidy terms correctly; every worked example and cross-check value matches the code and the 37 tests. Notes: no compact ISCN formula string (by design) and a cosmetic comment mismatch in `GetChromosomeBaseName`.
- **State: CLEAN** — no defect found; no code change required. Karyotype filter: 37/37 passed. Full suite: 4484/4484 passed.
