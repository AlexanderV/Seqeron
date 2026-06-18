# Validation Report: ASSEMBLY-COVER-001 — Coverage (Depth) Calculation

- **Validated:** 2026-06-15   **Area:** Assembly
- **Canonical method(s):** `SequenceAssembler.CalculateCoverage(string reference, IReadOnlyList<string> reads, int minOverlap = 20)` → `int[]`; internal helper `FindBestAlignment` (tested indirectly).
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent of repo artifacts)

| Source | Retrieved | What it confirms |
|--------|-----------|------------------|
| Metagenomics Wiki — SAMtools breadth of coverage (https://www.metagenomics.wiki/tools/samtools/breadth-of-coverage) | WebFetch 2026-06-15 | Breadth = covered bases / reference length. Worked example: **32876 covered / 45678 total = 0.719** (71.9%). |
| Daniel Cook — Depth/Breadth from a BAM (https://www.danielecook.com/calculate-depth-and-breadth-of-coverage-from-a-bam-file/) | WebFetch 2026-06-15 | Depth of coverage = **`Sum of Depths / genome size`**; breadth = **`Bases Mapped / genome size`**; depth = "average number of reads aligned to an individual base". |
| Illumina — Sequencing Coverage for NGS (https://sapac.illumina.com/.../coverage.html) | WebFetch 2026-06-15 | NGS coverage = "average number of reads that align to / 'cover' known reference bases"; Lander/Waterman **C = LN/G** with C=coverage, L=read length, N=#reads, G=haploid genome length. |

### Formula check
- Per-base depth[i] = #{reads spanning position i} — matches Metagenomics Wiki ("number of reads mapping to a specific reference position").
- Average depth = Σ depth / G — matches Cook verbatim (`Sum of Depths / genome size`).
- Breadth = #{i : depth[i] ≥ 1} / G — matches Cook (`Bases Mapped / genome size`) and Metagenomics Wiki.
- C = LN/G — matches Illumina verbatim. Cross-checked against the per-base average on the worked dataset: L=5, N=3, G=10 → C = 15/10 = 1.5, **identical** to Σdepth/G = 1.5. Two independent methods agree.

### Edge-case semantics
- Empty reads / no aligned reads → all-zero depth, average 0, breadth 0 (Cook, Metagenomics Wiki). Defined & sourced.
- Empty reference → empty array (one value per reference base; zero bases → length-0). Consistent with the per-base definition.
- Unaligned read (below minOverlap) → contributes 0. Sourced (unaligned read adds no depth).

### Independent cross-check (hand computation of M1)
Reference `ACGTTGCAAT` (len 10); reads `ACGTT`/`TTGCA`/`GCAAT`. Tracing the implementation's `FindBestAlignment` by hand:
`ACGTT`→pos 0 (5 matches, covers [0,5)); `TTGCA`→pos 3 (ref[3..7]=`TTGCA`, 5 matches, [3,8)); `GCAAT`→pos 5 (ref[5..9]=`GCAAT`, 5 matches, [5,10)). Depth = **[1,1,1,2,2,2,2,2,1,1]**, Σ=15, average 1.5, breadth 10/10=1.0. Matches the spec's M1/M7/M8 exactly and the source definitions.

### Findings / divergences (NOTE)
- **Boundary-clipping convention not realised by the placement model.** Sources describe a read overhanging the reference end contributing only its overlapping portion (clip at the boundary). The repository's `FindBestAlignment` only scans positions `0 … reference.Length − read.Length`, i.e. it places a read **only where it fits entirely**. Consequently the `min(p + L, reference.Length)` clip in the depth loop is genuinely **dead/defensive** code and the "overhang → partial contribution" case is unreachable; an over-long read simply fails to place (depth 0, M4). This is an honest, documented placement-model limitation (description §4.2, §5.3 "intentionally simplified") — **not** an arithmetic error. The depth-counting arithmetic itself is exactly source-faithful. Hence PASS-WITH-NOTES rather than PASS.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs:772-819`.
- Null `reference`/`reads` → `ArgumentNullException` (ThrowIfNull).
- Allocates `int[reference.Length]` zeroed; for each read, `FindBestAlignment` returns best position (≥ minOverlap matches, leftmost tie) or −1; increments `[pos, min(pos+L, refLen))`.
- `FindBestAlignment` scans only positions where the read fits → over-long reads return −1.

### Formula realised correctly?
Yes. depth[i] = count of placed reads spanning i, matching the validated per-base definition. Average/breadth are correctly left to the caller as Σ/len and covered-count/len (the tests derive them exactly).

### Cross-verification table recomputed vs code (test run)
| Case | Expected (sourced/hand) | Code result | Match |
|------|-------------------------|-------------|-------|
| M1 depth | [1,1,1,2,2,2,2,2,1,1] | same | ✅ |
| M7 Σ=15, avg=1.5 | 15 / 1.5 | same | ✅ |
| M8 breadth=1.0 | 10/10 | same | ✅ |
| M3 single read | [0,0,0,0,0,1,1,1,1,1,0,0,0,0,0] | same | ✅ |
| M4 over-long → all 0 | [0,0,0,0] | same | ✅ |
| M5 unmatched → all 0 | all 0 | same | ✅ |
| M6 empty reads | all 0 len 10 | same | ✅ |
| S1 partial breadth/avg=0.5 | [1,1,1,1,1,0,0,0,0,0] | same | ✅ |
| Empty reference → [] | len 0 | same | ✅ (new test) |

### Variant/delegate consistency
No public overload of this method in `SequenceAssembler`; the MCP wrapper (`AlignmentTools.CalculateCoverage`) delegates to it; `BedParser.CalculateCoverage` is a separate, unrelated unit. No divergence.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M1/M7/M8/S1 assert exact arrays/values traced to the per-base definition + Cook/Metagenomics formulas and hand-computation, not to code output. A deliberately-wrong impl (e.g. inclusive intervals, off-by-one, wrong sum) would fail them.
- **No green-washing:** all assertions are exact (`Is.EqualTo` on arrays / `.Within(1e-10)` on derived doubles), no Greater/AtLeast/Contains, no widened tolerances, no skips.
- **Coverage:** all 8 MUST + 2 SHOULD + 1 COULD cases, plus null-reference and null-reads validation. **Gap found & fixed:** the Stage-A-documented "empty reference → empty array" edge case (description §3.3) was untested → added `CalculateCoverage_EmptyReference_ReturnsEmptyArray`.
- **Honest green:** full unfiltered suite **6532 → 6533 passed, Failed: 0** (1 pre-existing skipped benchmark), build 0 errors (4 pre-existing warnings in unrelated files).

### Findings / defects
No algorithm defect. One test-coverage gap (empty reference) fixed in session.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — formulas and edge semantics all confirmed against three external sources + hand computation; one documented placement-model note (boundary-clip path unreachable; intentional, recorded in the description).
- **Stage B: PASS** — code realises the validated per-base depth exactly; tests are exact and sourced; added the missing empty-reference test.
- **End-state: ✅ CLEAN** — no defect; the single test gap was fully closed; build + full suite green.
- **Test-quality gate: PASS** (after adding the empty-reference test).
