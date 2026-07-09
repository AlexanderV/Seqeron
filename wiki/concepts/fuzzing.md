---
type: concept
title: "Fuzzing — random/malformed/boundary input for crashes and unhandled exceptions"
tags: [testing, validation, methodology]
sources:
  - docs/checklists/03_FUZZING.md
source_commit: 08ebf05f070b0cf9bc90d7ef1b1083b07a391606
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:validation-and-testing
      source: fuzzing-checklist
      evidence: "Priority P2. Feeds random/invalid/boundary data to surface crashes, hangs, and unhandled exceptions; 'критично важливо для парсерів файлових форматів і точок валідації' — critical for file-format parsers and input-validation points."
      confidence: high
      status: current
---

# Fuzzing

Fuzzing feeds **random, invalid, or boundary inputs** to surface crashes, hangs, and unhandled
exceptions — failures that correctness-oriented tests never provoke because they only feed
well-formed data. The property under test is weak but essential: the code must **fail cleanly**
(a typed exception or a defined empty result), never crash, hang, or corrupt. In genomics the
hot targets are **file-format parsers** and **input-validation points** — garbage sequences,
truncated files, invalid characters, null bytes. This is a **P2** member of the
[[validation-and-testing]] program; the checklist record is [[fuzzing-checklist]].

## The strategy taxonomy

- **RB — Random Bytes**: arbitrary byte streams into parsers and validators.
- **TF — Truncated Fields**: files cut off mid-record (a FASTQ with quality shorter than the
  sequence, a VCF row missing columns).
- **MC — Malformed Content**: structurally invalid but plausible (pattern longer than the
  subject, negative coordinates, out-of-vocabulary symbols).
- **BE — Boundary Exploitation**: the extremes — `0`, `-1`, `MaxInt`, empty string, single
  character, extremely long input.
- **INJ — Injection**: special characters, null bytes, Unicode, control characters.
- **OVF — Overflow**: extreme lengths and nesting depths.

## Coverage — a gap that was closed

The checklist's starting point was stark: **"Zero for Seqeron.Genomics"** — the only fuzz
coverage was `SuffixTreeFuzzTests` (corruption headers) in the SuffixTree library. The campaign
brought the genomics surface to **258 / 258**, prioritised as **12 high** (parsers + validation
points), **45 medium** (boundary inputs), **29 lower** (algorithm-specific edge cases). This is
the methodology most directly aligned with the [[mutation-testing]] finding that **null-handling
mutations survived in the parsers** — fuzzing and null-argument boundary tests attack the same
weak surface from opposite directions. The clean-failure contract it enforces is part of what the
[[research-grade-limitations|research-grade]] envelope depends on: the library degrades
predictably at its edges rather than crashing.
