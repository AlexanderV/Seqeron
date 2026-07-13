---
type: source
title: "Evidence: PARSE-FASTA-001 (FASTA parsing — >defline + sequence lines, multi-record, opt-in alphabets)"
tags: [validation, file-io]
doc_path: docs/Evidence/PARSE-FASTA-001-Evidence.md
sources:
  - docs/Evidence/PARSE-FASTA-001-Evidence.md
source_commit: 5dfac5a498e70d7d66967d830f9f18fbedda0224
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PARSE-FASTA-001

The validation-evidence artifact for test unit **PARSE-FASTA-001** — **FASTA parsing**
(`FastaParser.Parse` / `ParseFile` / `ParseFileAsync`, plus `ToFasta` / `WriteFile` output):
read one or more FASTA records — a `>` **description line (defline)** followed by one or more
**sequence lines** — into typed `FastaEntry` values (`Id`, `Description`, sequence), and write
them back with configurable line wrapping. This is a **FileIO** (file-parsing `PARSE-*`) family
Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. It cross-links to the family
anchor [[bed-format-parsing]] (the first ingested `PARSE-*` unit) and its EMBL sibling
[[parse-embl-001-evidence]]; see [[test-unit-registry]] for how units are tracked and
[[fuzzing]] for why parsers are the family's hottest malformed-input target. FASTA is the
simplest and most ubiquitous of the sequence formats — a `>`-header + free-form residue lines,
with no coordinate model (BED) or record grammar (EMBL) to summarize. This Evidence page records the
literature-traced **format facts and character-set tables**; the parser's **implementation surface**
(state machine, invariants, contract, complexity, opt-in `SequenceAlphabet`/`FastaRecord` overloads),
synthesized from the primary algorithm spec, lives in the concept [[fasta-parsing]].

## What this file records

- **Authoritative / online sources:**
  - **Wikipedia — FASTA format** (`en.wikipedia.org/wiki/FASTA_format`, accessed 2026-02-05) —
    the defline + sequence-line structure, the ≤80-characters-per-line recommendation,
    whitespace-in-sequence ignoring, and the "lower-case letters are accepted and mapped into
    upper-case" rule.
  - **NCBI BLAST Help — Query Input** (`blast.ncbi.nlm.nih.gov/doc/blast-topics/`, accessed
    2026-02-05) and the **NCBI FASTA specification** (`ncbi.nlm.nih.gov/blast/fasta.shtml`) —
    "blank lines are not allowed in the middle of FASTA input", the NCBI identifier formats
    (`gi|`, `gb|`, pipe-delimited), and the amino-acid code set for BLASTP/TBLASTN
    (incl. `U` selenocysteine, `X` any, `*` stop).
  - **Original papers:** Lipman DJ & Pearson WR (1985), *Science* 227:1435–41; Pearson WR &
    Lipman DJ (1988), *PNAS* 85:2444–8 — the FASTA search tools that named the format.
  - **Alphabet code sets** (re-retrieved 2026-06-24): **NC-IUB (1985)** "Nomenclature for
    Incompletely Specified Bases…", *Nucleic Acids Research* 13(9):3021–3030 (via Wikipedia
    "Nucleic acid notation") for the IUPAC nucleotide table; **bioinformatics.org** IUPAC
    nucleotide (`/sms/iupac.html`) and amino-acid (`/sms2/iupac.html`) tables, cross-confirmed
    against NCBI BLAST topics.
- **Format facts (from the artifact):**
  - **Description line (defline):** begins with `>`; the **first word** (up to the first
    space/tab) is the sequence **identifier**, the remainder is the optional **description**;
    all on a single line.
  - **Sequence lines:** follow immediately after the header; **≤80 characters recommended**
    (original format ≤120, usually ≤80); may be **multi-line (interleaved)** or **single-line
    (sequential)**; **whitespace within the sequence is ignored**.
  - **Multi-FASTA:** multiple records in one file, each starting with `>` — a concatenation of
    single-FASTA entries.
- **Character sets (opt-in `SequenceAlphabet` modes):** the default parser is **strict DNA-only**
  (`A`/`C`/`G`/`T`). Opt-in modes accept, verbatim from the sources: **IUPAC nucleotide**
  (A C G T U R Y S W K M B D H V N + gap `-`; NC-IUB codes W=A,T · S=C,G · M=A,C · K=G,T ·
  R=A,G · Y=C,T · B=C,G,T · D=A,G,T · H=A,C,T · V=A,C,G · N=any); **RNA** (A C G U); **Protein**
  (20 standard + ambiguity B/Z/J/X, rare U=Sec, O=Pyl, and stop `*`).
- **Edge cases (documented):** empty content → no sequences; **blank lines skipped** (NCBI:
  "not allowed in the middle of FASTA input" — defensive skip); whitespace in sequence ignored;
  **lower-case mapped to upper-case**; header-with-no-sequence is implementation-specific;
  special header characters (NCBI `gi|`/`gb|` pipes, colons); output line-width historically 80,
  modern configurable 60–80.
- **Output / round-trip:** `ToFasta` and `WriteFile` re-emit records with configurable line
  wrapping (default width **80**); the artifact's must-test list includes single/multi-record
  parsing, multi-line sequence concatenation, ID/description split, line-width wrapping on
  output, and **round-trip integrity** (parse → write → parse).

## Deviations and assumptions

The artifact records no source contradictions. The implementation-shaped choices in
`FastaParser.cs` are:

- **Default alphabet is strict DNA (`A/C/G/T`)** — broader IUPAC-nucleotide / RNA / protein
  code sets are **opt-in** via `SequenceAlphabet`; the sources define the codes but not which is
  the parser default, so DNA-only is a Seqeron design choice.
- **A header with no following sequence is not yielded** — grounded in NCBI ("the line after the
  FASTA definition line begins the nucleotide sequence"); an API-contract detail, not a format
  rule.
- **Blank lines are skipped** defensively even though NCBI states they are "not allowed" mid-input
  (tolerant-read behaviour), and **case mapping to upper-case** is delegated to the
  `DnaSequence` value type, matching the Wikipedia/NCBI upper-casing rule.
