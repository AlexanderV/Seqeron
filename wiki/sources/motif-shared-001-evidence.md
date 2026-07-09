---
type: source
title: "Evidence: MOTIF-SHARED-001 (Shared motifs via fixed-k word enumeration + matching-sequence quorum)"
tags: [validation, motif]
doc_path: docs/Evidence/MOTIF-SHARED-001-Evidence.md
sources:
  - docs/Evidence/MOTIF-SHARED-001-Evidence.md
source_commit: 9f3180f840fb594bb106edef7ac44083d6d57c8a
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: MOTIF-SHARED-001

The validation-evidence artifact for test unit **MOTIF-SHARED-001** — **shared motifs**
(`FindSharedMotifs`): enumerate every fixed-length (`k`) exact word across a *set* of
sequences and report each word present in at least `minSequences` of them, keyed by the
oligo-analysis **"matching sequences"** quorum. It is one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the model,
contract, invariants, worked oracle, and the not-LCSM contrast are synthesized in
[[shared-motifs]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **RSAT `oligo-analysis` manual** (rank 3, reference implementation of van Helden's
    method) — the operative verbatim definitions: **occurrences** = *"a simple count ...
    Overlapping matches are detected and summed"*; **matching sequences** = *"the number of
    sequences from the input set which contain at least one occurrence of the
    oligonucleotide"* (the quorum decision variable, each sequence counted at most once);
    fixed oligo length enumerated across the whole set; and *"only the first occurrence of
    each sequence is taken into consideration"* for the matching-sequence statistic.
  - **Das & Dai 2007**, *A survey of DNA motif finding algorithms* (rank 1, BMC
    Bioinformatics 8(Suppl 7):S21) — classifies van Helden oligo-analysis as a
    word-based/enumerative method and states the exact-word limitation: *"there are no
    variations allowed within an oligonucleotide"* (justifies exact, non-degenerate matching).
    Companion WebSearch: one approach *"records the number of sequences containing occurrences
    of each k-mer"* and reports those over a threshold.
  - **van Helden, André & Collado-Vides 1998** (rank 1, J Mol Biol 281(5):827–842) — the
    named primary (over-represented oligonucleotides via non-coding frequency tables); direct
    article HTTP 403, so RSAT + the survey carry the verbatim definitions.
  - **ROSALIND "Finding a Shared Motif" (LCSM)** (rank 4) — cited **only** to delineate the
    *alternative* longest-common-substring framing this unit does **not** implement (common
    substring = present in all; LCSM = the maximum-length such; non-unique).
- **Datasets (documented oracles):**
  - *Hand-traced k=3 quorum* — S0=`ATGATG`, S1=`ATGCCC`, S2=`CCCGGG`, `minSequences=2`:
    `ATG` → matching sequences {0,1} (count 2, Prevalence 2/3, note `ATG` twice in S0 still
    counts 1), `CCC` → {1,2}. Shared motifs = {`ATG`,`CCC`}.
  - *Rosalind LCSM sample* (contrast only) — GATTACA/TAGACCA/ATACA: LCSM answer `AC` (len 2);
    this unit at `k=2,minSeq=3` instead reports **all** 2-mers present in all three (`AC`,
    `TA`, …) — fixed-k quorum, not a single longest substring.
- **Corner cases / failure modes:**
  - Within-sequence repeat contributes **1** to the matching-sequence count (presence/absence).
  - Overlapping occurrences are summed only for the raw occurrence count, not the quorum.
  - Exact-word only — a one-substitution near-miss is not matched (Das & Dai).
  - `k` longer than the shortest sequence → that sequence yields no words.
  - Empty collection → no motifs; `k < 1` → throws.

## Deviations and assumptions

**Deviations: none** — the matching-sequence quorum, fixed-k exact enumeration, and
per-sequence presence/absence counting follow the RSAT/van Helden definitions. Two
**assumptions**, both presentation/API conveniences the sources do not prescribe:

1. **Default `k = 6`, `minSequences = 2`** — in RSAT's permitted range (any oligo length,
   any quorum) but not source-mandated; treated as caller-supplied (the algorithm is correct
   for any `k ≥ 1`, `minSequences ≥ 1`).
2. **`Prevalence = matchingSequences / totalSequences`** — RSAT reports raw matching-sequence
   counts; the fraction (a value in (0, 1]) is a presentation convenience consistent with the
   definition, not itself a source formula.

Recommended coverage (MUST): a word in exactly the quorum number of distinct sequences reports
the correct `SequenceIndices` + count; a word repeated within one sequence contributes 1;
below-quorum words excluded; `Prevalence` = matchingSequences/totalSequences exactly. SHOULD:
exact-word semantics (one substitution not matched); `k` > shortest sequence yields no words
from it. COULD: empty input → no motifs, `k < 1` throws. No contradictions among sources.
