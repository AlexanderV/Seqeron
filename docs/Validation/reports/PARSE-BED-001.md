# Validation Report: PARSE-BED-001 — BED File Parsing

- **Validated:** 2026-06-24   **Area:** FileIO
- **Canonical method(s):** `BedParser.Parse(content/reader/file, BedFormat)`, `ParseFile`, `FilterByChrom`, `FilterByRegion`, `MergeOverlapping`, `Intersect` (+ supporting filters, interval ops, BED12 block ops, stats, writer)
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.IO/BedParser.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/BedParserTests.cs` (103 BedParser-matched tests)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (fetched fresh, not from memory)
- **UCSC Genome Browser FAQ — BED format** (https://genome.ucsc.edu/FAQ/FAQformat.html), verbatim:
  - chromStart: *"The first base in a chromosome is numbered 0"* → 0-based start.
  - chromEnd: *"The chromEnd base is not included in the display of the feature"* → exclusive (half-open) end.
  - Three required fields: chrom, chromStart, chromEnd.
  - score: *"A score between 0 and 1000."*
  - strand: *"Either '.' (=no strand) or '+' or '-'."*
  - BED12: *"the first blockStart value must be 0... the final blockStart position plus the final blockSize value must equal chromEnd. Blocks may not overlap."*
  - zero-length: *"chromStart and chromEnd can be identical, creating a feature of length 0, commonly used for insertions."*
- **Wikipedia — BED (file format)** (https://en.wikipedia.org/wiki/BED_(file_format)), verbatim:
  - *"the system used by the BED format is zero-based for the coordinate start and one-based for the coordinate end."*
  - Length rationale: *"this calculation being based on the simple subtraction of the end coordinates (column 3) by those of the start (column 2)"* → length = end − start.
  - Contrast: *"Unlike the coordinate system used by other standards such as GFF..."* (GFF is 1-based inclusive).
  - Minimum 3 columns; *"These columns must be separated by spaces or tabs, the latter being recommended..."*; *"Each row of a file must have the same number of columns."*; header lines begin with `browser` / `track` / `#`.
- **Quinlan & Hall (2010), BEDTools** (Bioinformatics 26(6):841–842): standard merge/intersect/subtract interval semantics (same-chromosome merge; intersection = overlap region).

### Coordinate convention (the error-prone part)
BED is **0-based, half-open**: chromStart 0-based, chromEnd exclusive. First 100 bases → chromStart=0, chromEnd=100; **length = chromEnd − chromStart = 100 (NO +1)**. The distinction vs GFF (1-based inclusive) is explicit in both sources.

### Hand cross-checks
- `chr1 0 100` → start 0, end 100 (exclusive), length 100. ✔
- `chr1 0 1` → single base 0, length 1. ✔
- `chr1 5 5` → zero-length insertion at 5. ✔
- UCSC BED12 `chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512`: blockStarts relative to chromStart, first=0 ✔; final blockStart+blockSize = 3512+488 = 4000 = chromEnd−chromStart = 5000−1000 ✔; blocks don't overlap (0+567=567 ≤ 3512) ✔.

### Findings / divergences
None. Minor note: the Evidence doc phrases chromEnd as *"1-based, non-inclusive"* — this is exactly Wikipedia's *"one-based for the coordinate end"* framing of the half-open end, equivalent to "exclusive 0-based end". Not a defect.

## Stage B — Implementation

### Code path reviewed
`BedParser.cs:24-43` (BedRecord, `Length => ChromEnd - ChromStart`), `Parse` (102-149), `CountFields` (151-157), `ParseLine` (159-252), `ParseIntList` (254-259), `GenomicInterval` (57-78), block ops `ExpandBlocks`/`GetIntrons` (494-552).

### Coordinate convention realised correctly
- `BedRecord.Length => ChromEnd - ChromStart` (line 39) — half-open, no +1. ✔
- `GenomicInterval.Length => End - Start` (line 59) — consistent. ✔
- `chromStart > chromEnd` rejected (179-180); `start == end` allowed (zero-length insertions). ✔

### Field order / optional-column handling
- Tab split first (161), space fallback if <3 fields (165); non-numeric coords → skip. ✔
- Optional fields read by correct index: name[3], score[4], strand[5], thickStart[6], thickEnd[7], itemRgb[8], blockCount[9], blockSizes[10], blockStarts[11]. ✔
- score clamped to [0,1000] via `Math.Clamp` (187). ✔
- strand accepts only `+`/`-`/`.` (193); others leave Strand null (lenient, spec-compatible). ✔

### BED12 block constraints realised correctly
- blockCount must equal both array lengths (219). ✔
- first blockStart == 0 (225). ✔
- `starts[bc-1] + sizes[bc-1] == chromEnd - chromStart` (231) — the correct relative form of UCSC's "final blockStart + final blockSize must equal chromEnd", because blockStarts are stored relative to chromStart. ✔
- no overlap: `starts[i] < starts[i-1] + sizes[i-1]` rejected (235-239). ✔
- `ExpandBlocks`/`GetIntrons` add `ChromStart + BlockStarts[i]` to map to absolute coords (504, 538). ✔

### Header / consistency handling
- `track `, `browser ` (case-insensitive), `#`, and whitespace lines skipped (125-132). ✔
- Auto mode: first data line sets expected field count; mismatched lines skipped (135-143) — realises "same number of columns". ✔

### Cross-verification (recomputed vs code, via passing tests)
- `Coordinate_EndIsNonInclusive_LengthCorrect`: `chr1 0 100` → Length 100. ✔
- `Coordinate_ZeroLength_ValidForInsertions`: `chr1 5 5` → Length 0. ✔
- `Coordinate_BrowserToZeroBased_ConversionCorrect`: browser `chr7:127471196-...` → chromStart 127471195 (−1 on start, end unchanged) — correct 1-based-inclusive → 0-based-half-open conversion. ✔
- BED12 first-blockStart / final-block / overlap validation tests (806-922) all pass. ✔

### Test quality audit
Tests assert exact sourced coordinates, lengths, scores, strands, and block coords, each with a UCSC/Wikipedia evidence rationale string — not no-throw tautologies. Deterministic. Coordinate-system block (730-796) directly locks 0-based half-open length = end − start.

### Findings / defects
None.

## Verdict & follow-ups
- Stage A: **PASS** — 0-based half-open, field order, score/strand ranges, and BED12 block constraints confirmed against fresh UCSC + Wikipedia fetches + BEDTools; hand cross-checks reproduced.
- Stage B: **PASS** — implementation realises the validated spec; `Length = end − start` with no off-by-one; field order, optional-column, and block validation correct.
- **State: CLEAN.** No defect found; no code changed. Build succeeds; BedParser-filtered tests 103/103 pass.
