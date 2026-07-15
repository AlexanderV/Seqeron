---
type: source
title: "Evidence: ASSEMBLY-MERGE-001 (Contig Merging)"
tags: [validation, assembly]
doc_path: docs/Evidence/ASSEMBLY-MERGE-001-Evidence.md
sources:
  - docs/Evidence/ASSEMBLY-MERGE-001-Evidence.md
source_commit: 35a37c03842e8267eaa63463955c27f25f2c9ca5
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ASSEMBLY-MERGE-001

The validation-evidence artifact for test unit **ASSEMBLY-MERGE-001** — contig merging: the
suffix–prefix overlap collapse primitive `MergeContigs(contig1, contig2, overlapLength)` that
removes one copy of the shared region to form a superstring. One instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the overlap definition,
merge arithmetic, oracles and assumptions are summarized in [[contig-merge-overlap-collapse]], the
anchor for the assembly MERGE family. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (all rank 1, accessed 2026-06-13; PDFs fetched and text-extracted with
  `pdftotext -layout`):
  - **Langmead (JHU) — "Assembly & shortest common superstring" notes** — the overlap definition
    (length-`l` suffix of X = length-`l` prefix of Y), the `suffixPrefixMatch(x, y, k)` primitive
    (longest suffix ≥ k matching a prefix, else 0; guards `len < k`), superstring-by-collapse ("just
    concatenate them" without the shortest requirement; collapse keeps one copy of the overlap), and
    the greedy worked traces (`BAA`+`AAB` → `BAAB`; the `{AAA,AAB,ABB,BBB,BBA}` chain → `AAABBBA`).
  - **Langmead (JHU) — "Overlap Layout Consensus assembly" notes** — contig terminology, the OLC
    stages (Overlap graph / Layout bundling / Consensus), "an overlap is a suffix/prefix match",
    "report only the longest suffix/prefix match", and the generalized-suffix-tree overlap
    *discovery* step (O(N + a)) — noted as discovery, not the merge step.
  - **MIT 7.91J Foundations of Computational and Systems Biology, Lecture 6 (Spring 2014, MIT OCW)**
    — independent corroboration of the same overlap definition and the OLC consensus stage.
  - **Reference-only:** Compeau, Pevzner & Tesler (2011), *Nature Biotechnology* 29:987–991
    (DOI 10.1038/nbt.2023) — cited as background on overlap/de Bruijn assembly.
- **Datasets (published oracles):**
  - `BAA` + `AAB`, overlap 2 → `BAAB` (length 4 = 3 + 3 − 2) — greedy-SCS trace.
  - `{AAA, AAB, ABB, BBB, BBA}`, each overlap 2, chained → SCS `AAABBBA` (length 7).
  - `BAA` + `AAB`, overlap 0 → `BAAAAB` (plain concatenation, length 6).
- **Corner cases / failure modes** — no overlap (`suffixPrefixMatch` = 0) ⇒ merge equals `X + Y`;
  overlap bounded by `min(|X|, |Y|)` (cannot take a suffix/prefix longer than the string); when
  multiple overlaps exist only the longest is used, removing exactly one copy of the region.
- **Recommended coverage** — MUST: exact collapse `BAA`+`AAB` (ov 2) → `BAAB`; chained SCS →
  `AAABBBA`; overlap 0 → `X + Y`; overlap = `min(|X|,|Y|)` full-containment boundary. SHOULD:
  overlapLength > min ⇒ concatenation; negative overlap ⇒ concatenation; null contig ⇒
  `ArgumentNullException` (mirrors `AssembleOLC` / `BuildDeBruijnGraph`). COULD: empty + non-empty
  (ov 0) returns the non-empty contig (identity); length invariant `|merge| = |c1|+|c2|−usedOverlap`.

## Assumptions (from the artifact)

Two assumption records, both API-contract rather than algorithm-correctness: (1) **the caller-supplied
overlap length is trusted, not re-verified** — `MergeContigs` is a low-level collapse primitive;
verifying the region truly matched is the separate job of `FindOverlap` / `FindAllOverlaps`, and an
unverified length changes only whether the collapsed region matched, not the merge arithmetic.
(2) **out-of-range overlap (≤ 0 or > min length) → plain concatenation** — follows directly from the
source facts that overlap 0 means "no overlap → concatenate" and a valid overlap is bounded by
`min(|X|, |Y|)`.

No contradictions among the sources — Langmead's SCS notes, Langmead's OLC notes and MIT 7.91J
Lecture 6 all give the identical suffix-of-X / prefix-of-Y overlap definition, so the three
corroborate one another. The Compeau 2011 paper is cited for background only.
