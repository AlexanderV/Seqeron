---
type: source
title: "Evidence: PARSE-FASTQ-001 (FASTQ parsing — 4-line records + Phred+33/+64 quality encoding)"
tags: [validation, file-io]
doc_path: docs/Evidence/PARSE-FASTQ-001-Evidence.md
sources:
  - docs/Evidence/PARSE-FASTQ-001-Evidence.md
source_commit: d977e955ad5bf9f2eea32e0bd6c12987ab01edbc
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PARSE-FASTQ-001

The validation-evidence artifact for test unit **PARSE-FASTQ-001** — **FASTQ parsing**
(`FastqParser.cs`): read FASTQ records — the **4-line unit** `@header` / sequence / `+` /
quality — into typed reads, decode the per-base **Phred** quality string, and re-emit
(round-trip) records. This is a **FileIO** (file-parsing `PARSE-*`) family Evidence file and one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern. It cross-links to the family anchor [[bed-format-parsing]] (the first ingested `PARSE-*`
unit) and its format siblings [[parse-fasta-001-evidence]] and [[parse-embl-001-evidence]]; see
[[test-unit-registry]] for how units are tracked and [[fuzzing]] for why parsers are the
family's hottest malformed-input target. Unlike FASTA/EMBL, FASTQ carries a **quality-score
encoding** with real algorithmic content, so that nuance is factored into its own concept,
[[phred-quality-encoding]]; only the format-and-testing facts specific to this unit live here.

## What this file records

- **Authoritative / online sources:**
  - **Wikipedia — FASTQ format** (`en.wikipedia.org/wiki/FASTQ_format`, accessed 2026-02-05) —
    the 4-line record structure, the Phred+33 / Phred+64 encodings, the `Q = −10·log₁₀(p)`
    formula, the ASCII ranges, the offset auto-detection heuristic, and the documented parsing
    edge cases.
  - **Cock et al. (2009), "The Sanger FASTQ file format…"**, *Nucleic Acids Research*
    38(6):1767–1771 (DOI 10.1093/nar/gkp1137) — the authoritative Sanger FASTQ definition and
    the Solexa/Illumina encoding variants and their history.
  - **NCBI Sequence Read Archive — File Format Guide** (`ncbi.nlm.nih.gov/sra/docs/submitformats/`)
    — SRA submission requirements, paired-end conventions (`/1`, `/2` suffixes), and encoding
    requirements.
- **Format facts (from the artifact):**
  - **Record = 4 lines:** `@<identifier> [description]`, the sequence, a `+` separator line
    (optionally repeating the identifier), and the quality string.
  - **Core invariant:** **sequence length == quality-string length** per record.
  - **Quality encoding** — Phred+33 (Sanger / Illumina 1.8+, Q 0–93) and Phred+64
    (Illumina 1.3–1.7, Q 0–62); the score math, ASCII ranges, boundary characters, and the
    `< '@'` → +33 / `> 'I'` → +64 auto-detection heuristic are detailed in
    [[phred-quality-encoding]].
  - **Illumina conventions:** header
    `@INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y READ:FILTER:CONTROL:INDEX`; paired-end indicators
    `/1`,`/2` (or `1:`,`2:`); interleaved R1/R2 records.
- **Documented edge cases:** legacy **multi-line** sequence/quality (Sanger files may wrap a
  record across lines); a literal **`@` inside a quality string** (makes multi-line parsing
  ambiguous); a **`+` inside sequence** data (unusual but allowed); **empty records / blank
  lines** skipped gracefully; **encoding-detection failure defaults to Phred+33**.
- **Boundary quality values** (test oracles): Phred+33 `!`=Q0, `I`=Q40, `~`=Q93 (max);
  Phred+64 `@`=Q0, `h`=Q40, `~`=Q62 (max).
- **Must-test surface:** single/multi-record and multi-line parsing; ID/description header split
  (first space separates `Id` and `Description`, specials + unicode preserved); per-record
  encoding auto-detection; **statistics** (total reads/bases, GC%, Q20/Q30 %, per-position mean
  quality — 1-based, population stddev); **filtering** by quality threshold and by length range;
  **trimming** (quality-based end trim, adapter removal = 3′ prefix + internal match); and
  **round-trip integrity** (parse → write → parse yields equivalent data).

## Deviations and assumptions

The artifact records **no source contradictions**. The implementation-shaped choices in
`FastqParser.cs`:

- **`ErrorProbabilityToPhred` caps at Q93** for `p ≤ 0` — the maximum score representable in a
  Sanger Phred+33 string (a representability bound, not a spec rule).
- **Encoding auto-detection is per-record and deterministic**, and the ambiguous `@`–`I`
  character window **defaults to Phred+33** (the current standard), per the Wikipedia heuristic.
- **Header parsing** splits `Id` from `Description` on the first space and **preserves special
  characters and unicode** — an API-contract detail, not a format constraint.
