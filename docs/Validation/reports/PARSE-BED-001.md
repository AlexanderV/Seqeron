# Validation Report: PARSE-BED-001 — BED File Parsing

- **Validated:** 2026-06-12   **Area:** FileIO
- **Canonical method(s):** `BedParser.Parse(content/reader/file, BedFormat)`, `ParseFile`, `FilterByChrom`, `FilterByRegion`, `MergeOverlapping`, `Intersect` (+ supporting filters, interval ops, BED12 block ops, stats, writer)
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.IO/BedParser.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/BedParserTests.cs` (88 tests)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **UCSC Genome Browser FAQ — BED format** (https://genome.ucsc.edu/FAQ/FAQformat.html):
  - chromStart: "The first base in a chromosome is numbered 0." → 0-based start.
  - chromEnd: "The chromEnd base is not included in the display of the feature" → exclusive (half-open) end.
  - 3 required fields: chrom, chromStart, chromEnd. 9 optional in order: name(4), score(5), strand(6), thickStart(7), thickEnd(8), itemRgb(9), blockCount(10), blockSizes(11), blockStarts(12).
  - score: "between 0 and 1000". strand: "." (no strand), "+" or "-".
  - BED12: "The first blockStart value must be 0"; "the final blockStart position plus the final blockSize value must equal chromEnd"; blockStarts are "calculated relative to chromStart".
- **Wikipedia — BED (file format)** (https://en.wikipedia.org/wiki/BED_(file_format)):
  - "zero-based for the coordinate start and one-based for the coordinate end" = zero-based, half-open interval. Length = end − start (no +1). Contrasts with GFF (1-based inclusive).
  - Field order must be respected; intermediate columns must be filled.
  - "separated by spaces or tabs, the latter being recommended."
  - Header/comment lines: `browser`, `track`, `#`.
- **Quinlan & Hall (2010), BEDTools** (Bioinformatics 26(6):841–842): standard merge/intersect/subtract interval semantics (same-chromosome merge; intersection = overlap region).

### Coordinate convention (the error-prone part)
BED is **0-based, half-open**: chromStart 0-based, chromEnd exclusive. A feature covering the first 100 bases is chromStart=0, chromEnd=100; **length = chromEnd − chromStart = 100 (NO +1)**. Confirmed against both UCSC and Wikipedia.

### Worked example
Line `chr1\t0\t100\tfeature1\t960\t+`:
- chrom = `chr1`, chromStart = 0, chromEnd = 100, **length = 100 − 0 = 100**
- name = `feature1`, score = 960, strand = `+`

All values match the spec; length uses `end − start` with no off-by-one.

### Findings / divergences
None. Description (TestSpec + Evidence doc) is faithful to authoritative sources. Minor note: the Evidence doc phrases chromEnd as "1-based, non-inclusive" — this is the standard Wikipedia framing of the half-open end and is equivalent to "exclusive 0-based end"; not a defect.

## Stage B — Implementation

### Code path reviewed
`BedParser.cs:24-43` (BedRecord, `Length => ChromEnd - ChromStart`), `Parse` (102-149), `CountFields` (151-157), `ParseLine` (159-252), `ParseIntList` (254-259).

### Coordinate convention realised correctly
- `BedRecord.Length => ChromEnd - ChromStart` (line 39) — half-open, no +1. ✔
- `GenomicInterval.Length => End - Start` (line 59) — consistent. ✔
- Worked example recomputed in code: parses to ChromStart=0, ChromEnd=100, Length=100, Name="feature1", Score=960, Strand='+'. ✔

### Field order / optional-column handling
- Tab split first (line 161), space fallback if <3 fields (165). ✔
- Required: chrom (0), chromStart (1), chromEnd (2); non-numeric coords → null/skip. ✔
- Optional read by index in correct order: name[3], score[4], strand[5], thickStart[6], thickEnd[7], itemRgb[8], block fields[9..11]. ✔
- score clamped to [0,1000] via `Math.Clamp` (187). ✔
- strand accepts only `+`/`-`/`.` (193); others leave Strand null (lenient, spec-compatible). ✔

### Edge cases (Stage-A) traced/tested
- 3-column minimal (BED3): parsed; tested `Coordinate_*`. ✔
- 12-column full (BED12): blockCount/array-length match, first blockStart==0, last block reaches feature length, no overlap (lines 211-245). ✔
- track/browser/comment/empty lines skipped (125-132). ✔
- start>end rejected (179-180); zero-length start==end allowed (insertions). ✔
- negative coords: `int.TryParse` accepts negatives; a negative start with start≤end would parse, but is out of biological range — not asserted by spec as an error and not a stated defect. Negative score is clamped to 0. ✔
- blockStarts relative to chromStart: `ExpandBlocks`/`GetIntrons` add `ChromStart + BlockStarts[i]` (504, 538). ✔
- Column consistency in Auto mode: first data line sets expected field count; mismatched lines skipped (135-143). ✔

### Cross-verification (recomputed vs code, via passing tests)
- `chr1 0 100 first_100_bases` → Length 100 (`Coordinate_EndIsNonInclusive_LengthCorrect`). ✔
- `chr1 0 1 first_base` → start 0, length 1. ✔
- `chr1 5 5 insertion_point` → length 0. ✔
- UCSC BED12 `chr22 1000 5000 cloneA 960 + ... 2 567,488, 0,3512` parses (first start 0, 3512+488==4000==5000-1000). ✔
- Merge same-chrom overlapping, intersect overlap region, diff-chrom not merged — match BEDTools semantics. ✔

### Test quality audit
Tests assert exact sourced values (coordinates, lengths, scores, strands, block coords) with UCSC-cited rationale strings, not just no-throw. Coordinate-system block (730-795) directly locks 0-based half-open. Deterministic.

### Findings / defects
None.

## Verdict & follow-ups
- Stage A: **PASS** — 0-based half-open, field order, score/strand, block constraints all confirmed against UCSC + Wikipedia + BEDTools.
- Stage B: **PASS** — implementation realises the validated spec; length = end − start with no off-by-one; field order and optional-column handling correct; worked example recomputed against code.
- **State: CLEAN.** No defect found. Build succeeds; Bed-filtered tests 88/88 pass; full suite 4486/4486 pass. No code changed.
