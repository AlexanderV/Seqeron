---
type: source
title: "Validation report: ANNOT-GFF-001 (annotation-layer GFF3 I/O — export source/score + per-transcript CDS phase)"
tags: [validation, annotation, file-io, governance]
doc_path: docs/Validation/reports/ANNOT-GFF-001.md
sources:
  - docs/Validation/reports/ANNOT-GFF-001.md
source_commit: ffe651dc64c6b2478550efc53671d6f5e2c1ccbe
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-11
---

# Validation report: ANNOT-GFF-001

The two-stage **validation write-up** for test unit **ANNOT-GFF-001** — GFF3 **input/output
inside the annotation layer** (`GenomeAnnotator.ToGff3` / `ParseGff3`), validated 2026-06-25.
This is the *report* artifact that feeds one row of the [[validation-ledger]]; it records the
validator's verdict on both the format description and the shipped code, following the
[[validation-protocol]] Stage A / Stage B method that underpins the [[validation-and-testing]]
campaign. The annotation-layer format/algorithm itself — the reduced `GenomicFeature` record, the
`GeneAnnotation`-only exporter, and the per-transcript CDS-phase rule below — is synthesized on the
concept page [[gff3-io]] (which cites the primary spec `docs/algorithms/Annotation/GFF3_IO.md`).

**Scope note — distinct from PARSE-GFF-001.** This unit is the **annotation layer's own** GFF3
emitter/parser on `GenomeAnnotator` (record `GeneAnnotation`, 0-based half-open internally,
carrying `Source`/`Score`), *not* the FileIO `GffParser` traced in [[parse-gff-001-evidence]]
and anchored on [[bed-format-parsing]]. Same 9-column [Sequence Ontology GFF3 v1.26] format,
different code path and a different record type (`GeneAnnotation` for export vs `GenomicFeature`
for parse, file 1-based). There is no separate `annot-gff-001-evidence` artifact; this report is
the sole wiki record for the unit.

## Verdict

**Stage A: PASS · Stage B: PASS · State: ✅ CLEAN.** A **fresh re-validation** (the unit was
reset to pending after a campaign export-fidelity fix: real `source`/`score` columns + per-
transcript cumulative CDS phase on both strands). The authoritative **SO GFF3 Specification
v1.26** was retrieved live this session and quoted directly rather than trusting the repo's own
tests. **No defect found; no code or test change required.** GFF3 fixture **46/46 green**; full
`dotnet test Seqeron.sln` green (Seqeron.Genomics.Tests 18783 passed / 0 failed).

## Canonical methods validated

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs`:

- `ToGff3(IEnumerable<GeneAnnotation>, string)` — writes `##gff-version 3` then per row
  `{seqId}\t{source}\t{Type}\t{Start+1}\t{End}\t{score}\t{Strand}\t{phase}\t{attrs}`. Emits the
  internal 0-based `Start` as **1-based** (`Start+1`); `source = IsNullOrEmpty ? "." : Source`;
  `score = HasValue ? invariant(value) : "."`; `phase = CDS ? computed : "."`.
- `ParseGff3(IEnumerable<string>)` — skips blank / `#` / `<9`-column lines; tolerant `TryParse`
  (malformed numerics → skip, no throw); `.` → null score/phase; empty strand → `.`;
  `Uri.UnescapeDataString` decode; auto-ID `feature_{n}`; preserves the file's 1-based start/end
  verbatim into `GenomicFeature`.
- `ComputeCdsPhases` — the load-bearing algorithm (see below).
- Helpers `FormatGff3Attributes`, `EncodeGff3Value`, `ParseGff3Attributes`; record fields
  `GeneAnnotation.Source` / `GeneAnnotation.Score`.

## The load-bearing algorithm — per-transcript cumulative CDS phase

Column 8 phase is REQUIRED for CDS features: the number of bases to skip from the 5′ end to reach
the next codon boundary. `ComputeCdsPhases` groups CDS rows by `GeneId` (Ordinal), orders each
group **5′→3′** — ascending genomic start on `+`, **descending on `−`** (the 5′ end of a minus-
strand feature is its *end*) — then assigns:

- first (5′-most) segment → **phase 0**;
- each later segment → **`(3 − (Σ preceding segment lengths) mod 3) mod 3`**, accumulating
  `End − Start` (= 1-based inclusive length for the half-open internal `[Start, End)`).

Phase depends only on *preceding* segment lengths, and each transcript is independent.

**Cross-checks (hand-recomputed against the SO v1.26 canonical EDEN gene, ctg123):**

| Case | Phase sequence | Test |
|------|----------------|------|
| Plus-strand `cds00003` (602/501/601) | **0, 1, 1** — exact spec match | M28 |
| Plus-strand `cds00001` (all ≡0 mod 3) | **0, 0, 0, 0** — exact spec match | M30 |
| Minus-strand (hand-derived; spec gives no `−` CDS example), emitted in input order | **2, 2, 0** | M29 |
| Per-transcript independence, single CDS | 0, 0 / 0 | M30, M21 |

## Other confirmed column rules (SO GFF3 v1.26)

- **1-based coordinates**, length = end − start + 1; internal 99 → emitted 100 (M15/M17). The
  BED-vs-GFF off-by-one contrast lives on [[bed-format-parsing]].
- **Real `source` (col 2)** and **`score` (col 6, InvariantCulture)** emitted; absent → `.`
  (M23–M26). Non-CDS phase → `.` (M20).
- **Col-9 escaping** (`EncodeGff3Value`): encodes exactly `\t \n \r % ; = & ,` + control chars
  (<0x20, 0x7F); everything else passes through — matching the spec's "no other characters may be
  encoded" and "unescaped **spaces are allowed** within fields; parsers split on tabs, not spaces"
  (reserved `a%3Bb%3Dc%26d%2Ce` M19; control `a%09b%0Ac%0D%25d` M22; `ID=gene 1` un-encoded M16).
- **Round-trip** source/score survive parse (M27, S4); parse decodes percent-escapes (M11);
  all strands `+ − . ?` preserved (M18); malformed/short lines skipped (M12).

## Test quality

46 tests, all green, all asserting **exact spec-sourced values** (column values, encoded byte
sequences, coordinates, hand-derived phase sequences) — no code-echoes, no "no-throw"
tautologies, deterministic. One **cosmetic** note (not a defect): the M28 doc-comment labels its
third segment "len 602" whereas the internal `(6999,7600)` is length 601; the third segment's own
length does not enter its phase, so the asserted `0,1,1` is unaffected. Numerics are locale-safe
(InvariantCulture) and the parser is exception-free on malformed input.

## Findings

- **No findings.** Both stages PASS, State ✅ CLEAN. The export-fidelity fix (real source/score +
  both-strand cumulative CDS phase) is faithfully realised against SO GFF3 v1.26. Zero code
  change; no follow-up.
