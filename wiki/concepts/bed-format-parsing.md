---
type: concept
title: "BED format parsing (0-based half-open intervals)"
tags: [file-io, algorithm]
sources:
  - docs/algorithms/FileIO/BED_Parsing.md
  - docs/Evidence/PARSE-BED-001-Evidence.md
  - docs/Validation/reports/ANNOT-GFF-001.md
source_commit: 1edb4d1fcb74936366d343cafa919645c6555e45
created: 2026-07-10
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: parse-bed-001-evidence
      evidence: "Test Unit ID: PARSE-BED-001 ... Algorithm Group: FileIO ... BED File Parsing"
      confidence: high
      status: current
---

# BED format parsing (0-based half-open intervals)

**BED** (Browser Extensible Data) is the UCSC tab-delimited format for genomic **intervals** —
each line is one feature on a chromosome. Parsing it correctly hinges on one non-obvious fact:
BED coordinates are **0-based, half-open**. This is the anchor page for Seqeron's file-parsing
**FileIO** (`PARSE-*`) family; the first ingested member is test unit **PARSE-BED-001**, whose
literature-traced record is [[parse-bed-001-evidence]]. [[test-unit-registry]] tracks the unit,
[[algorithm-validation-evidence]] describes the evidence-artifact pattern, and [[fuzzing]]
explains why parsers like this one are the campaign's highest-priority malformed-input target.

## The coordinate system (the thing that bites)

- **`chromStart` is 0-based** — the first base of a chromosome is position 0.
- **`chromEnd` is 1-based and exclusive** — it is the position *after* the last included base.
- A feature therefore spans the **half-open interval `[chromStart, chromEnd)`** and has length
  `chromEnd − chromStart`.

| Record | `chromStart` | `chromEnd` | Bases covered |
|--------|-------------|-----------|---------------|
| First 100 bases | 0 | 100 | 0–99 (100 bases) |
| Single base at position 0 | 0 | 1 | base 0 only |
| Zero-length insertion | 5 | 5 | none — an insertion point *at* position 5 |

The browser display `chr1:1-100` (1-based inclusive) corresponds to BED `chromStart=0,
chromEnd=100`. Mixing the two conventions is the classic off-by-one source when converting
between BED and 1-based formats (GFF/GTF/VCF are 1-based inclusive), so coordinate/coordinate
conversions must always account for the half-open, 0-based origin. The sibling GFF/GTF
tab-delimited annotation format (9-column, 1-based inclusive) is traced in
[[parse-gff-001-evidence]], and the tab-delimited **VCF** variant format (8 fixed columns +
optional FORMAT/genotype samples, 1-based `POS`) in [[parse-vcf-001-evidence]]. GFF3 has a
**second, annotation-layer** code path (`GenomeAnnotator.ToGff3` / `ParseGff3`, distinct from the
FileIO parser above) — its own concept [[gff3-io]], whose validation, including the per-transcript
**cumulative CDS phase** formula `(3 − Σ preceding lengths mod 3) mod 3` on both strands, is
recorded in [[annot-gff-001-report]].

## The column ladder: BED3 → BED12

Three required fields plus up to nine optional ones (12 total), and **every row in a file must
carry the same number of columns**:

1. **BED3 (required):** `chrom`, `chromStart`, `chromEnd`.
2. **Optional (BED4–BED12), in fixed order:** `name`, `score` (0–1000, grayscale),
   `strand` (`+` / `−` / `.`), `thickStart`, `thickEnd`, `itemRgb` (`R,G,B`), `blockCount`,
   `blockSizes` (comma-separated), `blockStarts` (comma-separated, **relative to `chromStart`**).

BED12 encodes multi-block features (e.g. exons of a transcript). Tab separation is recommended;
spaces are tolerated but break names that contain spaces.

## Validation rules (all from the UCSC spec)

- **`chromStart ≤ chromEnd`** — equal coordinates are legal (zero-length insertion); a line with
  `chromStart > chromEnd` is rejected (Seqeron's `ParseLine` returns null).
- **`score` clamped to `[0, 1000]`**.
- **`strand` ∈ {`+`, `−`, `.`}** (`.` = no/unknown strand).
- **Column consistency** — Auto mode locks the column count to the first data line.
- **BED12 block constraints:** first `blockStart` must be **0**; `blockCount` must equal the
  `blockSizes` / `blockStarts` array lengths; blocks may not overlap
  (`blockStarts[i] ≥ blockStarts[i−1] + blockSizes[i−1]`); the last block must reach the feature
  end (`lastBlockStart + lastBlockSize == chromEnd − chromStart`).

Header/metadata lines (`track `, `browser `, and `#` comments) are skipped, not parsed as
features. Known failure modes the parser guards against: non-numeric coordinates, fewer than 3
columns, block-count/array mismatch, block overflow past the feature length, and invalid strand
characters — the malformed inputs [[fuzzing]] exercises.

## Downstream use

Parsed BED intervals feed interval arithmetic (BEDTools-style merge / intersect / subtract /
complement). The evidence artifact cites BEDTools (Quinlan & Hall 2010) for that usage context
and for the round-trip parse→write→parse and zero-length-feature edge-case test practice.

## Parser surface & behavioral notes (`BedParser`)

The primary spec (`docs/algorithms/FileIO/BED_Parsing.md`) documents the concrete
`BedParser` (`src/Seqeron/Algorithms/Seqeron.Genomics.IO/BedParser.cs`) beyond the coordinate
rules above. `Parse` / `ParseFile` scan the input in **`O(n)` time, `O(1)` auxiliary** space,
splitting each line on **tabs first and falling back to whitespace only when fewer than 3 tab
fields are present**. Beyond `ParseLine`, the type ships an interval toolkit and BED12 helpers:

| Operation | Time | Behavioral note (distinct from the pure spec) |
|-----------|------|-----------------------------------------------|
| `FilterByChrom` / `FilterByRegion` / `FilterByStrand` / `FilterByLength` / `FilterByScore` | `O(r)` | Record-level predicates over `r` parsed records. |
| `MergeOverlapping` | `O(r log r)` | **Merges *touching* intervals** — treats `next.ChromStart <= current.ChromEnd` as mergeable (adjacent, not just overlapping). Dominated by the sort. |
| `Intersect` | `O(a·b)` worst case | Grouped by chromosome but still quadratic in the worst case. |
| `Subtract` | `O(a·b)` worst case | Per-record interval subtraction. |
| `ExpandIntervals` | — | **Swaps upstream/downstream on negative-strand records** (strand-aware flanking). |
| `CalculateCoverage` | — | Emits **depth change-points**, not one row per base position. |
| `ExpandBlocks` / `GetTotalBlockLength` / `GetIntrons` | — | BED12 block helpers. Each exon expands as `exonStart = chromStart + blockStarts[i]`, `exonEnd = exonStart + blockSizes[i]`; `GetIntrons` derives the gaps between blocks. |

## `Auto` vs. explicit format modes

`BedFormat` defaults to **`Auto`**, which reads the first non-header **data** line to fix an
expected field count and then **skips later data lines whose width differs** — so a mixed-width
BED file is *partially* parsed, not reinterpreted per line. Explicit non-`Auto` values
(`Bed3`/`Bed6`/`Bed12`/…) do **not** currently force their nominal field count: the parser still
follows whatever field count each line actually presents. Optional display fields
(`thickStart`, `thickEnd`, `itemRgb`) are parsed syntactically but not semantically validated,
and **bigBed conversion / full UCSC-toolchain compatibility is out of scope** (implementation
status: *Simplified*) — those workflows need external UCSC tooling.
