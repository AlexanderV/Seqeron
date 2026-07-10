---
type: source
title: "Evidence: QUALITY-PHRED-001 (Phred score handling — parse/encode + Phred+33 ↔ Phred+64 conversion)"
tags: [validation, file-io]
doc_path: docs/Evidence/QUALITY-PHRED-001-Evidence.md
sources:
  - docs/Evidence/QUALITY-PHRED-001-Evidence.md
source_commit: 25fe7f865ba8c3ce681652d15fc633919907a6e5
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: QUALITY-PHRED-001

The validation-evidence artifact for test unit **QUALITY-PHRED-001** — **Phred score handling**:
decode a FASTQ quality string to per-base Phred scores (`ParseQualityString`), re-encode scores to
characters (`ToQualityString`), and convert a quality string between the two live ASCII offsets
(`ConvertEncoding`, Phred+33 ↔ Phred+64). This is a **QUALITY**-family Evidence file and one
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern.

Its subject is entirely the **Phred quality-score encoding** already factored into the shared
concept [[phred-quality-encoding]] (first traced by the FASTQ parsing unit
[[parse-fastq-001-evidence]]). This unit adds a **primary-literature anchor** for that concept —
Cock et al. (2010), the de-facto FASTQ specification — and the **cross-variant conversion**
representability rules; both are enriched into [[phred-quality-encoding]] rather than duplicated
into a new concept. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Authoritative source:** **Cock et al. (2010)**, "The Sanger FASTQ file format for sequences
  with quality scores, and the Solexa/Illumina FASTQ variants", *Nucleic Acids Research*
  38(6):1767–1771 (DOI 10.1093/nar/gkp1137, retrieved via PMC2847217, accessed 2026-06-13) —
  authority rank 1, the peer-reviewed and de-facto FASTQ format specification.
- **Encoding facts (verbatim from the paper):**
  - **Sanger / Phred+33:** ASCII 33–126 encodes Phred Q 0–93 (offset 33). Char = `chr(Q + 33)`;
    `Q = ord(char) − 33`.
  - **Illumina 1.3+ / Phred+64:** ASCII 64–126 encodes Phred Q 0–62 (offset 64).
  - **Phred score definition:** `Q = −10·log₁₀(P)`, `P` = base-call error probability.
  - **Cross-variant conversion is a pure re-offset:** the Phred score is *invariant* across the
    two variants — converting fastq-sanger ↔ fastq-illumina only changes the ASCII offset (shift
    each byte by ±31). No numeric rescaling.
  - **Solexa is out of scope:** Solexa uses the odds-ratio score `Q = −10·log₁₀(P/(1−P))`
    (offset 64, ASCII 59–126, scores −5..62); Solexa→Phred conversion is lossy. This unit handles
    only the two straight-Phred variants.
- **Conversion representability (corner cases):**
  - **Phred+64 → Phred+33 always representable:** Phred+64 holds Q ∈ [0, 62] ⊆ Phred+33's [0, 93].
  - **Phred+33 → Phred+64 may overflow:** a Phred+33 score in (62, 93] exceeds the Phred+64
    maximum (62) and cannot be encoded — `ArgumentOutOfRangeException`.
  - **Below-offset byte → negative Q:** a character below the variant's offset decodes to a
    negative Phred score, which is invalid for either variant (Phred Q ≥ 0) — malformed input.
- **Worked oracles (test fixtures):**
  - Phred+33: `!`→Q0, `5`→Q20, `?`→Q30, `I`→Q40, `~`→Q93.
  - Phred+64: `@`→Q0, `h`→Q40, `~`→Q62.
  - Phred+64→Phred+33 (score preserved): `@h~` → `!I_` (Q 0/40/62).
  - Phred+33→Phred+64: `!I` → `@h` (Q 0/40).
- **Must-test surface:** decode both variants at boundary + interior chars; round-trip
  parse∘encode = identity; `ConvertEncoding` both directions preserving Phred score; out-of-range
  decode throws; Phred+33→Phred+64 overflow (Q>62, e.g. `~`=93) throws; empty/null handling;
  seeded round-trip property check.

## Deviations and assumptions

The artifact records **no source contradictions**. Two API-shape assumptions (range bounds are
source-backed; only the exception *type* is chosen):

- **Malformed-byte decode** (char decoding outside a variant's valid Phred range) raises
  `ArgumentOutOfRangeException` — surfaces malformed input rather than emitting a negative/invalid
  score.
- **Phred+33 → Phred+64 overflow** (Q > 62, per the representability corner case) raises
  `ArgumentOutOfRangeException`.
