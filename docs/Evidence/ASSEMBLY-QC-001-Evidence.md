# Evidence Artifact: ASSEMBLY-QC-001

**Test Unit ID:** ASSEMBLY-QC-001
**Algorithm:** `GenomeAssemblyAnalyzer.CalculateLocalQuality`, `FindSuspiciousRegions` (+ additive `GcSkewCalculator.CalculateWindowedGcContent(..., includePartialWindow)`); MCP `local_quality`, `find_suspicious_regions`
**Date Collected:** 2026-11-25 (B10 sweep U5)

---

## Online Sources

### lh3/sdust — `sdust.c` (symmetric DUST, Morgulis et al. 2006, Bioinformatics 22:134)

**URL:** https://raw.githubusercontent.com/lh3/sdust/master/sdust.c (+ `sdust.h`, `kdq.h`, `kvec.h`, `kalloc.h`, `ketopt.h`)
**Accessed:** 2026-11-25 (curl, neutral User-Agent; the C file read in full, **not compiled or run** — the sandbox refused to build downloaded code)
**Authority rank:** 3 (reference implementation)

1. `sdust_core`: word = triplet (`SD_WLEN` 3), window W = 64 bases, T = 20 (score threshold ×10), non-ACGT letters break the input into independent pieces.
2. Output intervals are `[start, finish)` (printed `start`, `finish`); overlapping or adjacent perfect intervals are merged in `save_masked_regions` (`p->start <= f`).
3. Perfect-interval test `r·10 > T·l`, r = Σ c_t(c_t−1)/2 over the interval's triplets, l = triplets − 1 … (`find_perfect`).

### Biopython 1.88 — `Bio.SeqUtils.gc_fraction`, `GC_skew`

Installed package (no network). `gc_fraction(seq[i:i+w])`, default `ambiguous="remove"`; `GC_skew` slices `seq[i:i+window]` for every i = 0, window, … < len, i.e. the last window is truncated — the rule of the new `includePartialWindow`.

### QUAST `basic_stats.GC_content` / misassembly definitions

Via docs/Evidence/ASSEMBLY-STATS-001-Evidence.md (F10): GC = (G+C)/(A+C+G+T). QUAST / GAGE misassemblies (relocation, inversion, translocation) are defined from alignments to a reference; no sequence-only GC / N / complexity rule with thresholds is published, so `FindSuspiciousRegions` is documented as a screen.

### Linguistic complexity

Canonical `SequenceComplexity.CalculateLinguisticComplexity` (Orlov & Potapov 2004; Rosalind LING; docs/Evidence/ASSEMBLY-REPEAT-001-Evidence.md): LC = Σ V_i / Σ min(a^i, N − i + 1), m = 10, a = 4 extended by every other symbol (N → 5).

---

## Reference numbers (`tools/review/oracles/B10/qc_oracle.py`, `python3 -I`)

Independent Python: Biopython `gc_fraction` windows, LC from the formula (a = 5 with N), a line-by-line port of `sdust_core`, region union/merge. Fixtures use a deterministic LCG mirrored in the C# tests.

| Fixture | Expected (oracle = C#) |
|---------|------------------------|
| `LocalQuality` random 1500, w 1000 | (0,1000) GC .529 LC .913347087763847; (500,1000) .516 / .9154102523408982; (1000,500) .49 / .9509239624356256 |
| `LocalQuality` 400 rnd + 50 N + 300 rnd, w 300, step 100 | 8 windows; window 200: GC .48, N 50, LC .7844554902846814 |
| Suspicious: random 2000 + A×2000 + random 2000 | [0,749] GC dev .16883; [1000,4000] GC dev .33317 + DUST 999.5; [4750,5999] GC dev .21483 |
| Suspicious: 3000 + A×300 + 3000 | [1711,1717] DUST 2.5; [2750,3499] GC .25686 + DUST 149; [5372,5379] 3.0; [6237,6244] 3.0 |
| Suspicious: GC island 0.95 | [2473,2479], [2750,4249] GC .476 + DUST 48.5453, [4851,4858] |
| Suspicious: random + N×400 + random | [1250,2249] N .8 |
| poly(A)600 | [0,599] DUST 299 = 598·597/2/597 |

All C# values agree with the oracle to 1e-9. (The oracle's first LC draft used the Troyanskaya product form and disagreed with the canonical unit, which is the Orlov-Potapov/Rosalind sum form; the oracle was corrected to the sum form, not the code.)

## Findings

F58–F61 (see B10 summary): infinite loop for windowSize 1 / invalid sizes; inconsistent GC denominators and private GC loop; unsourced `FindSuspiciousRegions` thresholds, score and region ends; additive Core `includePartialWindow`.
