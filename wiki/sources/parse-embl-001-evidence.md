---
type: source
title: "Evidence: PARSE-EMBL-001 (EMBL flat-file parsing — ID/AC/FT/SQ records + INSDC locations)"
tags: [validation, file-io]
doc_path: docs/Evidence/PARSE-EMBL-001-Evidence.md
sources:
  - docs/Evidence/PARSE-EMBL-001-Evidence.md
source_commit: 6fe820d709820582a9f2476ebbfae0034dc0d5e1
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PARSE-EMBL-001

The validation-evidence artifact for test unit **PARSE-EMBL-001** — **EMBL flat-file
parsing** (`EmblParser.Parse` / `EmblParser.ParseFile`): read EMBL/ENA sequence entries —
two-character **line-type records** (`ID` / `AC` / `DT` / `DE` / `FH`/`FT` feature table /
`SQ` sequence / `//` terminator) — into typed entries with features, qualifiers, and
normalized sequence. This is a **FileIO** (file-parsing `PARSE-*`) family Evidence file and
one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. It cross-links to the family
anchor [[bed-format-parsing]] (the first ingested `PARSE-*` unit); see [[test-unit-registry]]
for how units are tracked and [[fuzzing]] for why parsers are the family's hottest
malformed-input target. EMBL is a **flat-file record** format — categorically different from
BED's tab-delimited intervals — but shares the INSDC feature-table grammar with GenBank (a
close cousin; the parser even exposes a `ToGenBank` conversion), so no separate concept page
is warranted yet.

## What this file records

- **Online / authoritative sources:**
  - **EBI EMBL User Manual, Release 143 (March 2020)** —
    `ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt`; ENA/EBI's official flat-file spec.
    Supplies the line-type table, the per-record grammar (§3.4.1–3.4.21), the data-class and
    taxonomic-division vocabularies (§3.1/§3.2), and the IUPAC base codes (Appendix A).
  - **INSDC Feature Table Definition v11.3 (Oct 2024)** — `insdc.org/files/feature_table.html`;
    the DDBJ/ENA/GenBank shared feature-table + **location-descriptor** grammar. Re-retrieved
    this session (2026-06-24 and 2026-06-26) from the EBI/INSDC/DDBJ mirrors for the operator
    and remote-reference rules.
- **Record grammar (from the artifact):**
  - **Line structure** — every line starts with a 2-char type code + 3 blanks; content begins
    at column 6. `XX` spacer lines and `//` terminator delimit blocks/records.
  - **`ID` line** (§3.4.1): `ID <accession>; SV <version>; <topology>; <mol_type>; <data_class>;
    <division>; <length> BP.` — e.g. `X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.` Topology
    is `circular` or `linear`.
  - **`AC`** semicolon-separated accessions (primary first; ranges like `X00001-X00005` allowed);
    **`DT`** two dates (created / last-updated + version); **`DE`/`KW`/`OS`/`OC`/`OG`** metadata;
    reference block ordered `RN, RC, RP, RX, RG, RA, RT, RL`; `DR` cross-refs; `CC` comments.
  - **`FH`/`FT` feature table** — fixed-format header, then features: `FT <key> <location>` with
    `/qualifier="value"` lines; qualifier continuation starts at column 22.
  - **`SQ` line** carries the composition summary (`Sequence <len> BP; <A> A; <C> C; …`); sequence
    data lines are 60 bases/line in groups of 10, 5'→3', stored **lowercase** (normalized on read).
- **INSDC feature-location descriptors** (shared with GenBank):
  - Simple `100..200`, single base `467`, between-bases site `123^124`.
  - Partial `<1..200` (5′), `100..>500` (3′), `<1..>500` (both).
  - Operators: `complement(loc)`, `join(loc,…)` (end-to-end into one contiguous span),
    `order(loc,…)`. **Nesting rule (verbatim):** `complement` may combine with `join`/`order`,
    but **nested `join`/`order` within the same location is illegal**.
  - **Remote references** `accession[.version]:span` (e.g. `J00194.1:100..202`) — 1-based
    inclusive on the *remote* entry — and may appear nested inside an operator, e.g.
    `join(1..100,J00194.1:100..202)`.
- **Vocabularies:** data classes CON/PAT/EST/GSS/HTC/HTG/WGS/TSA/STS/STD; taxonomic divisions
  PHG/ENV/FUN/HUM/INV/MAM/VRT/MUS/PLN/PRO/ROD/SYN/TGN/UNC/VRL; IUPAC ambiguity codes (R/Y/M/K/
  S/W/H/B/V/D/N).
- **Edge cases:** null → throw/empty, empty/whitespace → empty collection; minimal valid record
  (ID+AC+DT×2+DE+KW+OS+OC+FT source+SQ+data+`//`); multiple records per string split on `//`;
  circular vs linear; lowercase-sequence normalization; multi-line continuation of
  DE/KW/OS/OC/RA/RT/RL and FT qualifier values.

## Deviations and assumptions

The parser follows the EBI manual and INSDC feature-table spec directly; the notable
implementation-shaped points documented in the artifact are:

- **Remote-reference parsing fix:** a nested `accession[.version]:` prefix is captured
  per-segment and stripped **before** the local numeric span parse — otherwise the version
  digit (`.1`) is mis-read by the shared range regex as a spurious single-base part. Captured
  into `Location.RemoteParts` (accession, version, span); top-level remote refs stay in
  `RemoteAccession`/`RemoteVersion`.
- **Remote-aware location→sequence assembly (2026-06-26 enhancement):** the library is
  **offline-first and does no network I/O**; the **caller supplies a resolver delegate** for a
  remote entry's sequence, and `FeatureLocationHelper.ResolveLocationSequence` /
  `EmblParser.ResolveLocationSequence` perform the full assembly. Rules taken verbatim from
  INSDC §3.4/§3.5: `join`/`order` concatenate in listed order; `complement(...)`
  reverse-complements the assembled string; `complement(join(a,b)) == join(complement(b),
  complement(a))` (outer complement reverses element order); remote spans are 1-based inclusive.
  Worked oracles are given (e.g. `join(1..3,7..9)` over `ACGTACGTAC` → `ACGGTA`;
  `complement(Y.1:1..4)` with `Y→AAAC` → `GTTT`).
- **ASSUMPTION — `<`/`>` partials slice the stated number verbatim** (the spec gives no other
  coordinate; the partial flag does not move the slice bounds).
- **ASSUMPTION — missing/null resolver contributes an empty segment** (the grammar does not
  define an unavailable remote entry; local segments are still assembled in place, matching the
  clamp-not-throw behaviour of local out-of-range spans).

No source contradictions are recorded.
