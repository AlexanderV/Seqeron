---
type: source
title: "Evidence: PARSE-BED-001 (BED file parsing — 0-based half-open intervals)"
tags: [validation, file-io]
doc_path: docs/Evidence/PARSE-BED-001-Evidence.md
sources:
  - docs/Evidence/PARSE-BED-001-Evidence.md
source_commit: 44edb4ad55f6742c06654061979db03670e8dcfa
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PARSE-BED-001

The validation-evidence artifact for test unit **PARSE-BED-001** — **BED file parsing**:
read UCSC BED interval records (`chrom` / `chromStart` / `chromEnd` + up to nine optional
columns) into typed features, enforcing the format's coordinate and block constraints. This
is a **FileIO** (file-parsing `PARSE-*`) family Evidence file and the **first** ingested
unit of that family; it is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The format itself — the
0-based half-open coordinate model, the BED3→BED12 column ladder, and the validation
rules — is summarized in [[bed-format-parsing]]. See [[test-unit-registry]] for how units
are tracked and [[fuzzing]] for why parsers are the family's hottest malformed-input target.

## What this file records

- **Online sources:**
  - **UCSC Genome Browser FAQ — BED format** (authoritative, `FAQformat.html#format1`) — the
    3 required + 9 optional fields (12 total); the **critical coordinate system**:
    `chromStart` **0-based**, `chromEnd` **1-based non-inclusive** (so `chr1:1-100` in the
    browser = `chromStart=0, chromEnd=100`, spanning bases 0–99); **zero-length features**
    where `chromStart == chromEnd` (insertion points); the BED12 **block constraints** (first
    `blockStart` must be 0, `lastBlockStart + lastBlockSize == chromEnd − chromStart`, blocks
    may not overlap); header lines (`track `, `browser `, `#` comment).
  - **Wikipedia — BED (file format)** — tab separation recommended (space allowed), **all
    rows must have the same column count**, and the formal **GA4GH BEDv1** spec (2021).
  - **BEDTools manual** (Quinlan & Hall 2010) — interval operations (merge / intersect /
    subtract / complement) as the usage context, and the round-trip / edge-case test practice.
- **Coordinate & column facts (from the artifact):**
  - Required BED3 = `chrom`, `chromStart` (0-based), `chromEnd` (exclusive). Optional BED4–BED12
    add `name`, `score` (0–1000, grayscale), `strand` (`+`/`−`/`.`), `thickStart`, `thickEnd`,
    `itemRgb` (`R,G,B`), `blockCount`, `blockSizes`, `blockStarts` (relative to `chromStart`).
  - `chromStart=0, chromEnd=100` → 100 bases (0–99); `0,1` → single base 0; `5,5` → zero-length
    insertion at position 5.
- **Validation rules enforced (per UCSC FAQ, §6 of the artifact):**
  - `chromStart ≤ chromEnd` — `ParseLine` returns null when `chromStart > chromEnd` (zero-length
    equal-coordinate features are allowed).
  - `score` clamped to `[0, 1000]` (`Math.Clamp`).
  - `strand` restricted to `+` / `−` / `.`.
  - **Column consistency** — Auto mode locks the column count to the first data line.
  - BED12 block rules — first `blockStart == 0`; `blockCount` equals the `blockSizes` /
    `blockStarts` array lengths; no block overlap (`blockStarts[i] ≥ blockStarts[i−1] +
    blockSizes[i−1]`); last block reaches the feature end.
- **Documented failure modes:** `chromStart > chromEnd` (non-zero-length), non-numeric
  coordinates, fewer than 3 columns, `blockCount` mismatch, block overflow past feature length,
  invalid strand character.

## Deviations and assumptions

**Deviations: none** (§7 of the artifact: *"All parsing rules are directly derived from UCSC
FAQ and Wikipedia."*). The half-open, 0-based `chromStart` / 1-based-exclusive `chromEnd`
convention and the block/score/strand/column-consistency rules are taken verbatim from the
UCSC specification. The only implementation-shaped behaviours are the **null-on-invalid-line**
contract (`ParseLine` returns null rather than throwing on a rule violation) and header-line
skipping (`track ` / `browser ` / `#`), both API-contract details outside the format spec.
No source contradictions.
