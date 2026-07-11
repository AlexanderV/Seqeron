---
type: concept
title: "GFF3 I/O — annotation-layer lightweight parser/exporter (per-transcript CDS phase)"
tags: [annotation, file-io, algorithm]
mcp_tools:
  - parse_gff3
  - to_gff3
sources:
  - docs/algorithms/Annotation/GFF3_IO.md
source_commit: d4bfbd6d7b77287fda1b9ce616720a53312b9779
created: 2026-07-11
updated: 2026-07-11
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: annot-gff-001-report
      evidence: "Test Unit ID: ANNOT-GFF-001 (GFF3_IO.md); canonical methods GenomeAnnotator.ParseGff3 / ToGff3"
      confidence: high
      status: current
---

# GFF3 I/O — annotation-layer lightweight parser/exporter

**GFF3** (Generic Feature Format v3, Sequence Ontology, Stein 2020) is the tab-delimited genome
**annotation** interchange format: one feature per line, nine ordered columns, `#`/`##` comment
and directive lines. Seqeron has **two** GFF3 code paths. This page is the **annotation layer's**
lightweight helper surface — `GenomeAnnotator.ParseGff3(...)` / `ToGff3(...)` /
`EncodeGff3Value(...)`, test unit **ANNOT-GFF-001** (Implementation Status: *Simplified*). The
**FileIO** full parser (`GffParser.Parse` / `ParseFile` / `WriteToStream` / `BuildGeneModels`,
GFF3 **and** GTF/GFF2 dialects, hierarchical gene models) is the separate **PARSE-GFF-001** unit,
traced on [[parse-gff-001-evidence]]. The shared 9-column schema, attribute-dialect facts, and the
`Parent`/child hierarchy live there; this page carries only what is **distinct** about the
annotation-layer helper — a deliberately reduced record shape, a `GeneAnnotation`-only exporter,
and the load-bearing **per-transcript cumulative CDS phase** algorithm.

The two-stage validation verdict for this unit is [[annot-gff-001-report]] (**Stage A/B PASS,
State CLEAN**, no defect). See [[test-unit-registry]] for how the unit is tracked and
[[algorithm-validation-evidence]] for the evidence-artifact pattern.

## The nine columns and the coordinate model

`seqid`, `source`, `type` (a Sequence Ontology term), `start`, `end`, `score`, `strand`
(`+`/`-`/`.`/`?`), `phase` (`0`/`1`/`2` for CDS, else `.`), and `attributes` (`;`-separated
`key=value`). GFF3 coordinates are **1-based, inclusive** in the external text — the opposite of
[[bed-format-parsing]]'s 0-based half-open convention (that anchor spells out the off-by-one
contrast; GFF/GTF/VCF are the 1-based-inclusive side), and the same 1-based-inclusive convention
as the INSDC flat-file location grammar on [[insdc-feature-location]]. `.` is the undefined
sentinel for `score` and `phase`.

**Coordinate handling is the sharp edge of the two-path split.** Internally `GeneAnnotation.Start`
is **0-based**; `ToGff3` emits column 4 as `Start + 1` and leaves `End` unchanged (the internal
half-open `[Start, End)` has 1-based-inclusive length `End − Start`). `ParseGff3` preserves the
file's 1-based `start`/`end` **verbatim** into the returned `GenomicFeature` (no conversion on the
read path). So the exporter and parser use *different* record types and coordinate spaces.

## The reduced record shape (why round-trips are lossy)

`ParseGff3(IEnumerable<string>)` returns a **flat `GenomicFeature`** carrying only `FeatureId`,
`Type`, `Start`, `End`, `Strand`, `Score`, `Phase`, and `Attributes`. It **drops `seqid` and
`source`** (columns 1–2) entirely — so a `GenomeAnnotator` parse→export is **not lossless** even
on valid GFF3. Parse behaviour: skips blank and every `#`-prefixed line (both comments **and**
`##gff-version`/directive lines are ignored), skips any line with **fewer than 9** tab-separated
fields, `.`→`null` for score/phase, percent-decodes attribute values with
`Uri.UnescapeDataString`, and fills a **missing `ID`** with a running `feature_n` counter. Numeric
parsing is tolerant (`TryParse`; malformed numbers skip the line, no throw).

`ToGff3(IEnumerable<GeneAnnotation>, seqId = "seq1")` writes `##gff-version 3` then one row per
annotation: `{seqId}\t{Source|"."}\t{Type}\t{Start+1}\t{End}\t{Score|"."}\t{Strand}\t{phase}\t{attrs}`.
It exports **only** from `GeneAnnotation` (not from the parsed `GenomicFeature`), so it is **not** a
general-purpose GFF3 writer. `source` and `score` carry the feature's **real** values when present
and fall back to `.` when absent (`Source` defaults to `.`, `Score` to `null`); `seqId` is written
**verbatim with no escaping** (callers must supply a compliant identifier).

## Load-bearing algorithm: per-transcript cumulative CDS phase

Column 8 **phase** is *required* for CDS features — the number of bases to skip from the segment's
5′ end to reach the next codon boundary. `ToGff3` precomputes it per transcript:

1. Group CDS rows by `GeneId`.
2. Order each group **5′→3′**: ascending genomic `start` on `+`, **descending on `−`** (the 5′ end
   of a minus-strand feature is its `end`).
3. The **5′-most** segment → phase **0**.
4. Each later segment → **`(3 − (Σ preceding-segment lengths) mod 3) mod 3`**, accumulating
   `End − Start` (the 1-based-inclusive segment length) over all preceding segments.

Phase depends only on *preceding* segment lengths, and each transcript is independent. Non-CDS
features emit `.`. Worked oracles hand-checked against the SO v1.26 canonical EDEN gene
(`annot-gff-001-report`): plus-strand `cds00003` (lengths 602/501/601) → **0, 1, 1**; a run of
`≡0 mod 3` segments → **0, 0, 0, 0**; a single-segment CDS → **0**.

**Assumption — no per-feature phase field.** Because `GeneAnnotation` has no explicit phase field,
the 5′-most CDS segment is *always* treated as phase 0. A **5′-partial CDS** starting mid-codon
(non-zero opening phase) cannot be represented — matching the SO canonical example (all `cdsNNNNN`
segments start at phase 0) but a real simplification. Documented, not a defect.

## Column-9 percent-encoding (export)

`EncodeGff3Value` follows the GFF3 rule that **only required characters may be encoded**. It
encodes exactly tab, newline, carriage return, `%`, control characters (`<0x20`, `0x7F`), and the
column-9 reserved set `;` `=` `&` `,` — and **leaves spaces literal** (spaces are legal inside a
field; parsers split on tabs, not spaces, so `ID=gene 1` is emitted un-encoded). Attribute
**order** on export: `ID` and `product` first, then the rest, **omitting `translation`** (large
translated payloads are not preserved in lightweight output). Parse **decodes** percent-escapes.

## What this surface intentionally does not do

- **Lossy parse:** drops `seqid`/`source`; not a round-trip-faithful reader.
- **`GeneAnnotation`-only export:** cannot serialize arbitrary feature records.
- **No hierarchy:** no `Parent`/child gene-model reconstruction, no multi-parent handling — for
  those, and for **GTF/GFF2** dialects, use the FileIO `GffParser` ([[parse-gff-001-evidence]]).
- **No FASTA-section** handling; **no** faithful classic quoted-GFF2-attribute parsing (no repo
  surface covers the latter).
- **No validation** of Sequence Ontology terms, `seqId` encoding, or parent-child consistency; the
  helpers assume non-null enumerables (no explicit null guards).

These are Framework/Simplified [[research-grade-limitations|limitations]] of a lightweight
interoperability helper, not the full toolkit.
