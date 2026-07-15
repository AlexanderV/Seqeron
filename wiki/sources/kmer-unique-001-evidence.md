---
type: source
title: "Evidence: KMER-UNIQUE-001 (unique / singleton k-mers + min-count frequency filtering)"
tags: [validation, analysis]
doc_path: docs/Evidence/KMER-UNIQUE-001-Evidence.md
sources:
  - docs/Evidence/KMER-UNIQUE-001-Evidence.md
source_commit: 8107afde5814b7194214921227878bd429f3fc04
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: KMER-UNIQUE-001

The validation-evidence artifact for test unit **KMER-UNIQUE-001** — **unique k-mers / k-mers with
minimum count** (`KmerAnalyzer.FindUniqueKmers` + `FindKmersWithMinCount`): the two k-mer
frequency-filtering operations that select, respectively, the singletons (`Count == 1`) and the
recurrent k-mers (`Count ≥ minCount`, ordered by count descending). This is one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the operations
themselves are synthesized in [[unique-and-mincount-kmers]]. See [[test-unit-registry]] for how units
are tracked.

Both methods filter the shared `CountKmers` sliding-window multiset by the per-k-mer **Count** — the
number of overlapping occurrences of a k-mer as a substring (Compeau & Pevzner `Count(Text, Pattern)`).

## What this file records

- **Online sources:**
  - **Wikipedia — K-mer** (rank 4) — the foundational definition (k-mers = length-k substrings), the
    **total** count `L − k + 1` (overlapping, step +1) and universe size `n^k`, and the `AGAT` worked
    example (k=1→4, k=2→3 {AG,GA,AT}, k=3→2, k=4→1). Also the `k > L ⇒ L−k+1 ≤ 0 ⇒ 0 k-mers` corner.
  - **BioInfoLogics — k-mer counting, part I (Clavijo et al. 2018)** (rank 3) — the precise
    **total / distinct / unique** terminology: *distinct* = counted once even if repeated; ***unique*
    = "appear only once"** (frequency exactly 1) — the canonical definition `FindUniqueKmers` must
    implement. Worked `ATCGATCAC` k=3 table: 7 total, 6 distinct, 5 unique {TCG,CGA,GAT,TCA,CAC};
    `ATC` occurs twice (positions 0,4) so is distinct but NOT unique.
  - **Compeau & Pevzner — Bioinformatics Algorithms (2nd ed., 2015)** (rank 1) — `Count(Text, Pattern)`
    (overlapping occurrences), *most-frequent* k-mers maximise Count, and thresholding `Count ≥ t` to
    isolate recurrent/over-represented k-mers — the basis for `FindKmersWithMinCount`.

- **Documented corner cases / failure modes:** `k > L` ⇒ no k-mers ⇒ empty for both methods; a k-mer
  appearing ≥2 times (e.g. `ATC` in `ATCGATCAC`) is excluded from the unique set (unique ≡ the
  frequency-1 subset strictly).

- **Datasets (documented oracles):**
  - `ATCGATCAC` (L=9), k=3 → 7 total / 6 distinct / 5 unique {TCG,CGA,GAT,TCA,CAC}; ATC Count 2.
  - `AGAT`, k=2 → `FindUniqueKmers` = {AG, GA, AT} (all three unique).
  - `ACGTACGT`, k=4 → total 5, counts ACGT=2/CGTA=1/GTAC=1/TACG=1;
    `FindKmersWithMinCount(…,2)` = {(ACGT,2)}; `FindKmersWithMinCount(…,1)` = all 4 distinct, Count-desc
    (ACGT first); `FindUniqueKmers(…,4)` = {CGTA, GTAC, TACG}.
  - `AAAAA`, k=3 → zero unique (single distinct 3-mer `AAA`, Count 3).

- **Test-coverage recommendations:** MUST — the four oracle calls above + homopolymer→∅ unique.
  SHOULD — empty / `k > L` → empty for both methods; `k ≤ 0` → `ArgumentOutOfRangeException`.
  COULD — case-insensitivity (lower-case input yields the same unique set).

## Deviations and assumptions

Two **assumptions**, both source-consistent and non-correctness-affecting:

- **`minCount ≤ 1`** — `Count ≥ minCount` is then satisfied by every observed k-mer, so
  `FindKmersWithMinCount` returns **all distinct** k-mers ordered by count. This is the mathematically
  consistent extension of the `Count ≥ t` predicate, treated as defined behaviour rather than an
  unknown (authoritative sources define min-count filtering only for `t ≥ 1`).
- **Case normalisation** — input is upper-cased (consistent with the sibling `KmerAnalyzer` methods)
  so case variants count as the same k-mer; no source contradicts this.

No source contradictions — Wikipedia's `L − k + 1`, BioInfoLogics' distinct-vs-unique distinction, and
Compeau & Pevzner's `Count(Text, Pattern)` are mutually consistent.
