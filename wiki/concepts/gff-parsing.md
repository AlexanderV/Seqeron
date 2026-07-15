---
type: concept
title: "GFF/GTF parsing (FileIO GffParser — 9-column parser, GFF3+GTF+GFF2 dialects, gene models)"
tags: [file-io, algorithm]
sources:
  - docs/algorithms/FileIO/GFF_Parsing.md
  - docs/Evidence/PARSE-GFF-001-Evidence.md
source_commit: c2007120f9b30878e0a7a27ea4b6ccd9302979bf
created: 2026-07-13
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: parse-gff-001-evidence
      evidence: "Test Unit ID: PARSE-GFF-001 | Algorithm Group: FileIO (GffParser)"
      confidence: high
      status: current
---

# GFF/GTF parsing (FileIO GffParser — 9-column parser, GFF3+GTF+GFF2 dialects, gene models)

**GFF/GTF** are the nine-column tab-delimited genome-**annotation** formats: one feature per line,
`seqid`/`source`/`type`/`start`/`end`/`score`/`strand`/`phase`/`attributes`, with coordinates
**1-based and fully closed**. Seqeron's FileIO **`GffParser`** (test unit **PARSE-GFF-001**, status
*Simplified*) is the **fuller** of the repository's two GFF3 code paths: it parses **GFF3, GTF, and
legacy GFF2-style** rows from strings, files, and readers into typed `GffRecord` values, then layers
filtering, hierarchical gene-model building, statistics, writing, subsequence extraction, and
interval merging on top.

This concept is the **parser-surface** synthesis of the primary algorithm spec
`docs/algorithms/FileIO/GFF_Parsing.md` — it owns the `GffParser` contract, the line-oriented parse,
and the gene-model builder. It deliberately **does not re-derive** the literature-traced format
facts (the 9-column schema tables, the attribute-dialect / RFC 3986 percent-escape tables, the SO
GFF3 v1.26 EDEN oracle) — those live in the Evidence source page [[parse-gff-001-evidence]] and are
**cross-linked, not duplicated** here. This is a member of the FileIO file-parsing **`PARSE-*`**
family anchored by [[bed-format-parsing]] (where the **coordinate contrast** — GFF is 1-based
inclusive, BED is 0-based half-open — is spelled out); its format siblings are [[fasta-parsing]] and
[[fastq-parsing]]. [[test-unit-registry]] tracks the unit, [[algorithm-validation-evidence]]
describes the evidence-artifact pattern, and [[fuzzing]] explains why parsers are the campaign's
highest-priority malformed-input target.

## Two GFF3 code paths — this is the FileIO one

Seqeron has **two** GFF3 surfaces, and they are **not** the same unit:

- **This one — FileIO `GffParser`** (**PARSE-GFF-001**): the fuller parser. GFF3 **and** GTF/GFF2
  dialects, hierarchical `Parent`/child gene models, `< 8`-field rejection, merge/extract/statistics
  utilities.
- **The annotation-layer `GenomeAnnotator.ParseGff3` / `ToGff3`** (**ANNOT-GFF-001**), a
  deliberately reduced lightweight helper — see [[gff3-io]] (reduced `GenomicFeature` record,
  `GeneAnnotation`-only exporter, load-bearing per-transcript cumulative CDS-phase algorithm, GFF3
  only, `< 9`-field skip, drops `seqid`/`source`).

The shared 9-column schema, attribute-dialect facts, and `Parent` hierarchy concepts belong to the
format and are documented on [[parse-gff-001-evidence]]; the two concept pages each carry only what
is **distinct about their own implementation**.

## The format (what the parser consumes)

```text
seqid  source  type  start  end  score  strand  phase  attributes
```

Coordinates are **1-based, inclusive** (a feature spans both `start` and `end`). Column 9 encodes
attributes, and its **encoding differs by dialect** — GFF3 uses `ID=gene00001;Name=EDEN` while GTF
uses `gene_id "ENSG00001"; gene_name "TestGene";`. The full column/encoding/percent-escape tables
are on [[parse-gff-001-evidence]].

### Invariants (from the spec)

| ID | Invariant |
|----|-----------|
| INV-01 | GFF/GTF coordinates are 1-based and fully closed. |
| INV-02 | `start <= end` for valid features (closed intervals). |
| INV-03 | `phase` is meaningful for CDS and is one of `0`, `1`, `2`, or `.`. |
| INV-04 | Reserved characters in GFF3 attributes are percent-encoded (URL-escape model). |

## The parser (line-oriented, format-resolved)

`Parse(content, format)` / `ParseFile(filePath, format)` / `Parse(TextReader, format)` read the
input **line by line**:

1. **Skip** blank lines, `##` directives, and `#` comment lines.
2. In **`Auto`** mode, inspect `##gff-version` directives when present (see the gotcha below).
3. **Split each data row on tabs** and parse columns 1–8. A line is **rejected when it has fewer
   than 8 tab-separated fields**, or when `start`/`end` do not parse as integers.
4. Normalize sentinels: a `.` **score → `null`**, a `.` **phase → `null`**.
5. **Parse attributes** (column 9) per the resolved format; when column 9 is **absent**, the record
   receives an **empty attribute dictionary**.
6. Yield `GffRecord` values and expose the filter / gene-model / write / extract / merge helpers.

### Attribute parsing by dialect

| Format | Parsing rule |
|--------|--------------|
| **GFF3** | Split on `;`, split each part on the **first `=`**, then **URL-unescape** keys and values. |
| **GTF** | Split on `;`, trim each part, split on the **first space**, and **strip surrounding quotes** from the value. |
| **GFF2** | Reuses the **non-GTF** branch in the implementation. |

### Gotcha — `Auto` only consults `##gff-version`

`Auto` format detection is **intentionally conservative**: it uses **only** the `##gff-version`
directive to distinguish GFF3 from legacy GFF, and **otherwise defaults to GFF3**. It does **not**
infer GTF from the `key "value"` attribute syntax. **Consequence:** GTF inputs that lack a
`##gff-version` directive must be selected **explicitly** (`format = GffFormat.Gtf`) or their
`gene_id "…"` attributes will be mis-parsed as GFF3.

## Contract & surface (`GffParser`, `GffRecord`, `GeneModel`)

Implementation: `src/Seqeron/Algorithms/Seqeron.Genomics.IO/GffParser.cs`.

| Group | Entry points | Notes |
|-------|--------------|-------|
| **Parse** | `Parse(string, GffFormat)`, `ParseFile(string, …)`, `Parse(TextReader, …)` | String / file / reader; `format` defaults to **`Auto`**. |
| **Filter** | `FilterByType`, `FilterBySeqid`, `FilterByRegion`, `GetGenes`, `GetExons`, `GetCDS` | `O(n)` scan, `O(m)` matches. |
| **Hierarchy / attrs** | `BuildGeneModels`, `GetAttribute`, `GetGeneName` | Gene-model tree from `ID`/`Parent`. |
| **Statistics / write** | `CalculateStatistics`, `WriteToStream`, `WriteToFile` | `WriteToStream` emits `##gff-version 3` **only** for GFF3 output; formats GTF attributes with trailing `;`. |
| **Sequence / interval** | `ExtractSequence`, `MergeOverlapping` | Feature-level subsequence (interval + strand); interval merge. |

`GffRecord` exposes `Seqid`, `Source`, `Type`, `Start`/`End` (1-based inclusive `int`),
`Score` (`double?`), `Strand` (`char`), `Phase` (`int?`), and `Attributes`
(`IReadOnlyDictionary<string,string>`). `BuildGeneModels` returns a `GeneModel` aggregating gene,
transcript, exon, CDS, and UTR collections.

### Complexity

| Operation | Time | Space |
|-----------|------|-------|
| `Parse` / `ParseFile` | `O(n)` | `O(n)` |
| `FilterByType` / `FilterByRegion` | `O(n)` | `O(m)` (m = matches) |
| `BuildGeneModels` | `O(n)` | `O(n)` (index by `ID`/`Parent`) |
| `MergeOverlapping` | `O(n log n)` | `O(n)` (dominated by sorting) |

## Gene-model construction (`BuildGeneModels`)

The builder reconstructs the `ID`/`Parent` hierarchy over a **restricted feature vocabulary**:

1. Group top-level **`gene`** records.
2. Collect direct **transcript-bearing children** of type **`mRNA`, `transcript`, or `ncRNA`**.
3. Under each transcript, collect children whose type is **`exon`, `CDS`, or contains `utr`** (any
   type containing `utr` is treated as a UTR).

When a child record lists **multiple `Parent` IDs** (comma-separated), the **same child is attached
under each referenced parent**. Richer Sequence-Ontology hierarchies beyond this vocabulary remain
available only as raw `GffRecord` rows.

## Edge cases & intentional simplifications

- **Null / empty input** → **no records** (explicit early-return guards).
- **Directive (`##`) / comment (`#`) lines** → skipped; blank lines → skipped.
- **Row with fewer than 8 columns** → **skipped** (minimum structural validation). *Note the
  difference from the annotation-layer [[gff3-io]], which requires `< 9` fields to skip.*
- **`.` score / `.` phase** → stored as **`null`**; unknown strand `?` is valid.
- **Column 9 missing** → record gets an **empty attribute dictionary** (8-column rows are accepted).
- **Multiple `Parent` IDs** → child attached to **each** parent.
- **`Auto` detection is conservative** → GTF without `##gff-version` must be selected explicitly (see
  the gotcha above); the attribute syntax is **not** used to infer the dialect.
- **Round-trip is lossy for metadata** → `WriteToStream`/`WriteToFile` do **not** preserve original
  directives, comments, or full metadata beyond the generated `##gff-version 3` header.
- **No conformance validator** → the surface focuses on core row parsing plus practical annotation
  utilities, not on validating Sequence-Ontology terms or full parent-child consistency.

For validation of the sibling annotation-layer path see [[annot-gff-001-report]]; the format-facts,
dialect tables, and testing-methodology categories for **this** unit live on
[[parse-gff-001-evidence]].
