---
type: source
title: "Evidence: KMER-BOTH-001 (Both-strand / additive k-mer counting)"
tags: [validation, analysis]
doc_path: docs/Evidence/KMER-BOTH-001-Evidence.md
sources:
  - docs/Evidence/KMER-BOTH-001-Evidence.md
source_commit: c85acb02f865639aee5a5c5cc2c2257a1295fe7c
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: KMER-BOTH-001

The validation-evidence artifact for test unit **KMER-BOTH-001** — **both-strand k-mer counting**
(`KmerAnalyzer.CountKmersBothStrands`): for every k-mer, the sum of its forward-strand and
reverse-complement-strand occurrences. This is the second **K-mer** family Evidence file (after
[[kmer-async-001-evidence|KMER-ASYNC-001]]) and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the operation itself is synthesized in
[[both-strand-kmer-counting]]. See [[test-unit-registry]] for how units are tracked.

It validates the **additive ("balance")** strand convention — `count[w] = forward[w] +
forward[RC(w)]`, keeping a key per observed k-mer — and is careful to distinguish it from the
**canonical-collapsing** convention (Jellyfish `-C` / Mash), which it does **not** implement.

## What this file records

- **Online sources:**
  - **kPAL — Methodology** (authority rank 3) — the additive both-strand operation verbatim: kPAL
    balances a profile "by adding the values of each k-mer to its reverse complement", which
    "enforce[s] balance between sequence information from the minus or plus strand".
  - **Anvar et al. 2014 — kPAL paper** (Genome Biology 15:555, authority rank 1) — grounds the
    additive semantics: "The balance function uses a sum of k-mers and their reverse complements".
  - **Shporer et al. 2016 — Inversion symmetry of DNA k-mer counts** (BMC Genomics, rank 1) —
    the identity behind INV-01: counts of a length-`k` string on one strand equal the counts of its
    reverse-complement k-mer, so occurrences of `w` on the RC strand (read 5'→3') = occurrences of
    `RC(w)` on the forward strand ⇒ `count[w] = forward[w] + forward[RC(w)]`.
  - **Marçais & Kingsford 2011 — Jellyfish** (Bioinformatics 27(6):764–770, rank 1) — the
    single-strand counting primitive ("number of occurrences of every k-mer in a long string") on
    which both-strand counting builds, and the **canonical `-C`** contrast (collapses `w`/`RC(w)`
    onto one key) — explicitly **not** what this unit does.
  - **Mash issue #45 / Ondov et al.** (rank 3) — canonical definition ("only the lexicographically
    smaller of the forward and reverse complement representations of a k-mer is hashed"), cited to
    contrast the collapsing approach.
  - **BioInfoLogics — k-mer counting part I (Clavijo 2018)** (rank 4) — biological motivation for
    strand-aware counting only.

- **Documented corner cases / failure modes:** grand total `2·(L − k + 1)` (each strand contributes
  `L − k + 1` windows); reverse-complement palindromes (`RC(w)=w`, e.g. `AT`, `GC`, `ACGT`) get both
  contributions on one key ⇒ `count[w] = 2·forward[w]`; `k > L` ⇒ `L − k + 1 ≤ 0` ⇒ empty result on
  each strand ⇒ empty combined result.

- **Datasets (documented oracles):**
  - `ATGGC`, k=2 → `{ AT:2, TG:1, GG:1, GC:2, CC:1, CA:1 }` (forward `{AT,TG,GG,GC}` + `RC=GCCAT`
    `{GC,CC,CA,AT}`), grand total `8 = 2·(5−2+1)`.
  - `ACGT`, k=2 → `{ AC:2, CG:2, GT:2 }` (every 2-mer is an RC palindrome ⇒ each key doubled).
  - `AAA`, k=2 → `{ AA:2, TT:2 }` (`RC(AAA)=TTT`; equal counts on complementary keys).

- **Test-coverage recommendations:** MUST — the three datasets exactly; grand total `2·(L − k + 1)`;
  the per-key identity `count[w] = forward[w] + forward[RC(w)]`. SHOULD — case-insensitivity
  (lowercase == uppercase), `k = L` (single window per strand), `DnaSequence` overload delegates to
  the `string` overload. COULD — empty/null → empty, `k > L` → empty, `k ≤ 0` → throws.

## Deviations and assumptions

**Two ASSUMPTIONS**, both API-shape and non-correctness-affecting on valid input:

- **Empty/short input returns an empty dictionary** — no source explicitly defines the both-strand
  result for empty input or `k > L`; resolved consistently with the `L − k + 1 ≤ 0 ⇒ no windows`
  formula and sibling `CountKmers` behaviour.
- **`k ≤ 0` throws `ArgumentOutOfRangeException`** — sources define `k` as a positive substring
  length but prescribe no exception type; resolved to match the sibling
  `CountKmers` / `GenerateAllKmers` contract in this repository.

No source contradictions — kPAL "balance" and inversion symmetry give the identical additive
semantics, and the Jellyfish/Mash canonical wording is cited only to contrast the not-implemented
collapsing mode.
