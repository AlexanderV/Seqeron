---
type: source
title: "Evidence: PARSE-GENBANK-001 (GenBank flat-file parsing — LOCUS/FEATURES/ORIGIN + INSDC locations)"
tags: [validation, file-io]
doc_path: docs/Evidence/PARSE-GENBANK-001-Evidence.md
sources:
  - docs/Evidence/PARSE-GENBANK-001-Evidence.md
source_commit: 6dbe2cbf99c5f8ee4eb8eb90b438d7b10e3d5021
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PARSE-GENBANK-001

The validation-evidence artifact for test unit **PARSE-GENBANK-001** — **GenBank flat-file
parsing**: read NCBI GenBank sequence records — keyword-delimited **sections** (`LOCUS` /
`DEFINITION` / `ACCESSION` / `VERSION` / `KEYWORDS` / `SOURCE`+`ORGANISM` / `REFERENCE` /
`FEATURES` / `ORIGIN` / `//` terminator) — into typed entries with features, qualifiers, and a
normalized sequence. This is a **FileIO** (file-parsing `PARSE-*`) family Evidence file and one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern. It cross-links to the family anchor [[bed-format-parsing]] (the first ingested `PARSE-*`
unit) and its format siblings [[parse-embl-001-evidence]] (its INSDC close cousin),
[[parse-fasta-001-evidence]], and [[parse-fastq-001-evidence]]; see [[test-unit-registry]] for how
units are tracked and [[fuzzing]] for why parsers are the family's hottest malformed-input target.
GenBank and EMBL are the two INSDC flat-file dialects — same feature-table + location grammar,
different line syntax — so the shared descriptor grammar (complement/join/order/partial/remote) is
factored into its own concept, [[insdc-feature-location]]; only the GenBank-specific format and
testing facts live here.

## What this file records

- **Authoritative / online sources:**
  - **NCBI GenBank Sample Record** (`ncbi.nlm.nih.gov/Sitemap/samplerecord.html`, accessed
    2026-02-05) — the official field-by-field explanation of the flat-file format and the
    canonical example record (accession **U49845**).
  - **Wikipedia — GenBank** (`en.wikipedia.org/wiki/GenBank`, accessed 2026-02-05) — GenBank as
    NCBI's open-access nucleotide database (started 1982; doubles ~every 18 months; 34+ trillion
    bases as of Oct 2024).
  - **INSDC Feature Table Definition** (`insdc.org/documents/feature_table.html`) — the DDBJ/ENA/
    GenBank shared feature keys, qualifiers, and location-descriptor grammar (referenced by NCBI).
- **Record structure (from the artifact):** keyword-per-section layout — `LOCUS`, `DEFINITION`,
  `ACCESSION`, `VERSION` (`accession.version`, e.g. `U49845.1`), `KEYWORDS` (often just `.`),
  `SOURCE` with nested `ORGANISM` (scientific name + taxonomy lineage), optional `REFERENCE`
  blocks, `FEATURES`, `ORIGIN` (sequence marker), and the `//` record terminator.
- **`LOCUS` line** (whitespace-delimited): `LOCUS` keyword, locus name (≤16 chars), sequence
  length + `bp`/`aa`, molecule type (DNA/RNA/mRNA/…), optional topology (`linear`/`circular`),
  **GenBank division** (3-letter), and modification date `DD-MMM-YYYY`. Example:
  `LOCUS  SCU49845  5028 bp  DNA  PLN  21-JUN-1999`.
- **18 GenBank division codes:** PRI, ROD, MAM, VRT, INV, PLN, BCT, VRL, PHG, SYN, UNA, EST, PAT,
  STS, GSS, HTG, HTC, ENV. (Compare EMBL's taxonomic-division vocabulary — overlapping but not
  identical.)
- **`ORIGIN` sequence section:** each line = a position number + 60 nucleotides in 6 groups of
  10; sequence stored **lowercase** in standard format; numbers/spaces/newlines stripped on parse;
  uppercase normalization is common practice.
- **Feature locations** — the INSDC descriptor grammar (`n..m`, single `n`, `<n..m`/`n..>m`
  partials, `complement(…)`, `join(…)`, and `complement(join(…))`) is shared verbatim with EMBL;
  the full rule set and worked assembly oracles live in [[insdc-feature-location]].

## Documented edge cases and failure modes

- **From NCBI:** locus names no longer guaranteed to follow historical naming conventions; dates
  in `DD-MMM-YYYY` **or** `DD-MMM-YY`; empty `KEYWORDS` (`.`); abbreviated taxonomy lineage for
  long lineages; older records may lack some fields.
- **From implementation analysis:** empty content → empty enumerable; **null content → empty
  enumerable (defensive)**; missing `LOCUS` line → skip record; malformed location strings handled
  gracefully; missing `ORIGIN` → empty sequence; non-alphabetic characters in the sequence
  filtered out.
- **Multi-line fields** (`DEFINITION`, `ORGANISM`) span lines with continuation indentation;
  **multiple records** per file split on `//`.

## Invariants and testing methodology

Stated invariants: declared `LOCUS` length should match the actual sequence length; every record
ends with `//`; location `Start ≤ End` (equal for single positions); a feature's extracted
subsequence matches its location bounds; the complement flag is detected for every `complement(…)`
location. The artifact's recommended test categories mirror the sections — header/`LOCUS` parsing,
metadata extraction, organism/taxonomy, reference parsing, feature + location parsing, sequence
extraction/normalization, multi-record handling, and defensive edge cases (empty/null/missing/
malformed). Recommendations: use NCBI field formats as ground truth, test all 18 divisions and all
location variants, and verify the length invariant.

## Deviations and assumptions

The artifact is source-first (NCBI sample record + INSDC feature table as ground truth) and records
no source contradictions. The only implementation-shaped points are the defensive API-contract
behaviours above (null/empty → empty enumerable; missing `LOCUS` → skip; missing `ORIGIN` → empty
sequence; non-alphabetic sequence characters filtered) — contract choices outside the format spec,
not departures from it.
