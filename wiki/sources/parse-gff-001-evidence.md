---
type: source
title: "Evidence: PARSE-GFF-001 (GFF/GTF parsing — 9-column tab schema + GFF3/GTF attribute dialects)"
tags: [validation, file-io]
doc_path: docs/Evidence/PARSE-GFF-001-Evidence.md
sources:
  - docs/Evidence/PARSE-GFF-001-Evidence.md
source_commit: 0dc9d22d7d1e7c57ca0bd01652ff34963a0b08a7
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-13
---

# Evidence: PARSE-GFF-001

The validation-evidence artifact for test unit **PARSE-GFF-001** — **GFF/GTF parsing**
(`Parse(content)` / `ParseFile(filePath)` / `ToGff3(features)`): read General Feature Format
annotation records — **nine tab-delimited columns** (`seqid` / `source` / `type` / `start` /
`end` / `score` / `strand` / `phase` / `attributes`) — into typed features with parsed
attributes, plus a round-trip GFF3 writer. This is a **FileIO** (file-parsing `PARSE-*`) family
Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. Like [[parse-bed-001-evidence]] it
is a **tab-delimited interval** format (not an INSDC flat file such as
[[parse-genbank-001-evidence]] / [[parse-embl-001-evidence]]), so it shares the family anchor
[[bed-format-parsing]]; the **coordinate contrast** that matters here — GFF is **1-based,
inclusive** whereas BED is 0-based half-open — is already spelled out on that anchor page. See
[[test-unit-registry]] for how units are tracked and [[fuzzing]] for why parsers are the
family's hottest malformed-input target. The `GffParser` **implementation surface** (the parse
state machine, `GffRecord`/`GeneModel` contract, the `Auto`-only-`##gff-version` gotcha,
`BuildGeneModels`, and the merge/extract utilities) is synthesized in the concept [[gff-parsing]];
this Evidence page records the **literature-traced format facts** below and the coordinate model on
the shared anchor. Note the **separate annotation-layer** GFF3 path — `GenomeAnnotator.ParseGff3` / `ToGff3`,
test unit ANNOT-GFF-001 — is its own concept [[gff3-io]] (reduced record shape,
`GeneAnnotation`-only export, per-transcript CDS phase); this evidence file covers the fuller
`GffParser` (GFF3 **and** GTF/GFF2, hierarchical gene models).

## What this file records

- **Authoritative / online sources (all accessed 2026-02-05, confidence HIGH):**
  - **Wikipedia — General Feature Format** — general overview; the **9-column** tab-delimited
    structure; all dialects (GFF2, GFF3, GTF) share the same first 8 fields; **GFF2/GTF is
    deprecated in favour of GFF3**.
  - **UCSC Genome Browser — GFF/GTF Format** (`FAQformat.html`) — the GFF2 field descriptions
    and the **GTF extension**: attributes end in a semicolon, separated by exactly one space,
    with **mandatory `gene_id "value"` and `transcript_id "value"`**; score `0–1000` drives gray
    shading; fields must be tab-separated (spaces within a field are allowed).
  - **Sequence Ontology — GFF3 Specification v1.26** (Lincoln Stein, 18 Aug 2020,
    `The-Sequence-Ontology/Specifications`) — the authoritative GFF3 grammar: UTF-8 recommended;
    `##` directives, `#` comments, blank lines ignored; RFC 3986 escaping; predefined attribute
    tags and their semantics; the phase-field rule.

- **The nine columns (1-based coordinates):** `seqid`, `source` (generating program/procedure),
  `type` (e.g. gene/exon/CDS), `start` (1-based), `end` (1-based **inclusive**), `score`
  (numeric or `.`), `strand` (`+` / `-` / `.` unstranded / `?` unknown), `phase` (`0`/`1`/`2`
  for CDS, else `.`), `attributes` (column 9, tag-value pairs).

- **Phase field (column 8, GFF3):** required for CDS features; number of bases to remove from
  the 5′ end to reach the next codon — **phase 0** = codon begins at the first nucleotide,
  **1** = second, **2** = third. `.` for non-CDS.

- **Attribute dialects (the distinguishing feature of the family):**
  - **GFF3:** `tag=value` pairs, `;`-separated (e.g. `ID=gene00001;Name=EDEN`).
  - **GTF (GFF2 extension):** `tag "value";` pairs — quoted values, each terminated by `;` and
    separated by one space (e.g. `gene_id "ENSG00001"; transcript_id "ENST00001";`).
  - **GFF2:** a free-form `group` field (column 9) linking lines that share a group.
  - Predefined GFF3 tags: `ID` (required for features with children), `Name`, `Alias`,
    **`Parent`** (part-of relationship — the hierarchy link), `Target`, `Gap` (CIGAR-like),
    `Derives_from`, `Note`, `Dbxref`, `Ontology_term`, `Is_circular`.

- **Hierarchical feature model (GFF3):** `Parent` establishes gene → mRNA → exon/CDS part-of
  trees. **Multiple parents** are allowed via comma-separated `Parent` values
  (e.g. `Parent=mRNA00001,mRNA00002`); **discontinuous features** repeat the same `ID` across
  lines (canonical for multi-segment CDS).

- **RFC 3986 percent-escaping (GFF3):** tab `%09`, newline `%0A`, CR `%0D`, percent `%25`, and
  control chars `%00`–`%1F` / `%7F` must be escaped everywhere; column 9 additionally reserves
  `;` `%3B`, `=` `%3D`, `&` `%26`, `,` `%2C`. Parsers must **unescape** on read.

- **Directives / structure:** `##gff-version 3.x.x` must be the first line;
  `##sequence-region seqid start end` bounds a sequence; `##FASTA` marks embedded FASTA;
  `###` signals all forward references resolved. `#` comment and blank lines are ignored.

- **Canonical oracle (GFF3 spec):** the `ctg123` EDEN gene record — gene `1000..9000` →
  mRNA → exons / CDS with `ID`/`Parent`/`Name` attributes — and an Ensembl/UCSC **GTF** example
  (`gene`/`transcript`/`exon` with `gene_id`/`transcript_id`/`exon_number`).

## Documented corner cases and edge cases

- **From the GFF3 spec:** zero-length features (`start == end`, site to the right of the base);
  circular genomes (`Is_circular=true`); discontinuous features (same `ID` on multiple lines);
  multiple parents (comma-separated); polycistronic and trans-spliced transcripts.
- **Implementation edge cases:** empty **or** null content → empty collection; **malformed lines
  with < 8 fields → skip**; missing score (`.`) → null; missing phase (`.`) → null; unknown
  strand `?` → valid; URL-encoded characters → unescaped; `#` comment lines → skip; `##`
  directive lines → processed as metadata; blank lines → skip.

## Testing methodology

The artifact's required test categories mirror the format: column parsing (all 9 columns);
attribute parsing across both dialects (**GFF3 `key=value;` vs GTF `key "value";`**); RFC 3986
escape/unescape; hierarchical `Parent`/child linking; **format auto-detection** (GFF3 vs GTF vs
GFF2); directive processing (version, sequence-region); and **round-trip** (write via `ToGff3`
then re-parse preserves data).

## Deviations and assumptions

The artifact is source-first — all parsing rules trace directly to the Wikipedia overview, the
UCSC GFF/GTF FAQ, and the Sequence Ontology GFF3 v1.26 spec — and records **no source
contradictions** (stated confidence: HIGH). The implementation-shaped points are API-contract
behaviours outside the format spec: **null/empty content → empty collection** (defensive),
**lines with fewer than 8 fields are skipped rather than throwing**, and missing `score`/`phase`
(`.`) map to null. The 1-based-inclusive coordinate model is taken verbatim from the specs; the
BED-vs-GFF off-by-one contrast is documented on [[bed-format-parsing]].
