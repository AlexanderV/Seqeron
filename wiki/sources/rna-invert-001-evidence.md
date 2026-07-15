---
type: source
title: "Evidence: RNA-INVERT-001 (RNA inverted repeats / potential stem regions — antiparallel complementary arms)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-INVERT-001-Evidence.md
sources:
  - docs/Evidence/RNA-INVERT-001-Evidence.md
source_commit: b36d042573850cc6e5d5f099aa7f01a4be79cdb7
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-INVERT-001

The validation-evidence artifact for test unit **RNA-INVERT-001** — **RNA inverted
repeats** (the stem-forming, antiparallel **reverse-complement arms** that make a
potential stem-loop), on `RnaSecondaryStructure` in the RNA secondary-structure family.
It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; [[test-unit-registry]]
tracks the unit.

**No new concept is warranted — this is a reuse.** An inverted repeat is *exactly* the
`W G W̄ᴿ` model already synthesized on the repeats-family anchor
[[repetitive-element-detection]] (its *Inverted repeats* section, from the same IUPACpal
source); the RNA twist is only that the arm complement is the RNA base-pairing rule
{A-U, G-C} of [[rna-base-pairing]] and that a non-zero loop makes the object a **stem-loop**
— the same object detected structurally by [[pre-mirna-hairpin-detection]] and written in
[[rna-dot-bracket-notation]]. This page records what the artifact adds; the model lives on
those concept pages.

## What this file records

- **Online sources:**
  - **IUPACpal** (Alamro et al., *BMC Bioinformatics* 2021; PMC7866733; rank 1) — the
    operative formal definitions: an inverted repeat is *"a single stranded sequence of
    nucleotides with a subsequent downstream sequence consisting of its reverse
    complement"*; the **gapped model `W G W̄ᴿ`** with left arm `W`, spacer/loop `G`
    (`|G| ≥ 0`), right arm = reverse complement of `W`; the **mismatch model** *k*-Hamming
    (`δ_H(W, W̄ᴿ) ≤ k`, perfect = `k = 0`); Watson-Crick complement basis (RNA A⟷U, C⟷G);
    tool parameters minimum/maximum arm length, maximum gap size, maximum mismatches.
  - **Wikipedia "Inverted repeat"** (citing Ussery, Wassenaar & Borini 2008; rank 4) — the
    same definition and the worked example `5'---TTACGnnnnnnCGTAA---3'` (right arm `CGTAA`
    = reverse complement of left arm `TTACG`); **zero-gap ⇒ palindrome** (`5'TTACGCGTAA3'`
    equals its own reverse complement).
  - **EMBOSS `einverted` manual** (rank 3, reference tool) — confirms **inverted repeat =
    stem-loop**: it finds "regions of local alignment of the input sequence and its reverse
    complement that exceed a threshold score"; two arms complementary, intervening region =
    loop/bulge, mismatches/gaps = bulges. Its scored (match +3 / mismatch −4 / gap 12 /
    min-score 50) DP variant is noted but **not** implemented here.
- **Datasets (documented oracles):**
  - *Wikipedia/Ussery example, adapted to RNA* — `UUACGAAAAAACGUAA` (len 16): left arm `UUACG`
    (0–4), 6-nt loop `AAAAAA` (5–10), right arm `CGUAA` (11–15); right arm = revcomp(left)
    (UUACG → complement AAUGC → reverse CGUAA), arm length 5.
  - *Palindromic IR with minimal RNA loop* — `GGCCAAAGGCC` (len 11): left arm `GGCC` (0–3),
    3-nt loop `AAA` (4–6), right arm `GGCC` (7–10); GGCC → CCGG → GGCC.
- **Corner cases / failure modes:**
  - **Antiparallel, not parallel** — the right arm is the *reverse* complement; a sequence
    whose arms match in parallel (complement-only, not reversed) is not an IR.
  - Zero-gap palindrome (`|G| = 0`) is out of scope of the RNA default (a loop is required —
    biologically a hairpin loop needs ≥ 3 unpaired nt).
  - Arm shorter than the minimum length is not reported; a perfect arm is extended to its
    maximal complementary run.
  - Loop length below `minSpacing` / above `maxSpacing` excludes the IR; no `W G W̄ᴿ`
    decomposition ⇒ empty result.

## Deviations and assumptions

**Deviations from the reference tools: a documented scope restriction, not invented
behaviour.** The repository reports **perfect (zero-mismatch), ungapped arms only** — the
`k = 0` special case of the IUPACpal model (`W G W̄ᴿ` with `δ_H = 0`), via strict
Watson-Crick + IUPAC `GetComplement`/`GetRnaComplementBase`. EMBOSS `einverted`'s
scored mismatch/gap DP variant is **Not Implemented** (out of scope). Three assumptions:

1. **Perfect ungapped arms (`k = 0`).** A documented special case of the source model, so a
   scope restriction rather than a new rule.
2. **Loop bounds via `minSpacing`/`maxSpacing`.** These public parameters bound `|G|`,
   matching IUPACpal's "maximum gap size" / einverted loop concept; defaults
   (minSpacing 3, maxSpacing 100) are repository API conventions (3-nt minimal hairpin
   loop) and do not affect formal correctness of a reported IR.
3. **Maximal-arm + non-overlapping greedy reporting.** For each accepted left start the
   longest perfect antiparallel arm is reported and reported repeats do not overlap —
   mirrors einverted's maximal-local-alignment reporting; does not change which positions
   are valid complementary pairs.

Recommended coverage (MUST): a known IR (`UUACGAAAAAACGUAA`) returns arms at the exact
positions with right arm = reverse complement of left (antiparallel); a palindromic IR
`GGCCAAAGGCC` returns left 0–3 / right 7–10 / length 4; a parallel-only sequence is **not**
reported; no IR present → empty. SHOULD: arm below `minLength` excluded, exactly `minLength`
reported; loop outside `[minSpacing, maxSpacing]` excluded. COULD: a perfect arm longer than
`minLength` reported at full length. No source contradictions — IUPACpal, Ussery/Wikipedia,
and EMBOSS einverted agree on the reverse-complement-arms definition.
