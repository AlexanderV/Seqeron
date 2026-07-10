---
type: concept
title: "Shared motifs across a sequence set (fixed-k word enumeration + matching-sequence quorum)"
tags: [motif, algorithm]
mcp_tools:
  - find_shared_motifs
sources:
  - docs/Evidence/MOTIF-SHARED-001-Evidence.md
  - docs/algorithms/Motif_Discovery/Shared_Motifs.md
source_commit: 9f3180f840fb594bb106edef7ac44083d6d57c8a
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: motif-shared-001-evidence
      evidence: "Test Unit ID: MOTIF-SHARED-001 ... Algorithm: Shared Motifs via fixed-length word enumeration with matching-sequence quorum (oligo-analysis matching sequences)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:longest-common-substring
      source: motif-shared-001-evidence
      evidence: "MOTIF-SHARED-001-Evidence.md contrasts this unit with the ROSALIND LCSM framing: 'LCSM requires a substring present in all sequences with variable length. The unit under test instead fixes the length (k) and uses a quorum (>= minSequences), so LCSM is documented as a related-but-distinct algorithm, not the contract here.' Both find substrings common to multiple sequences but differ: fixed-k + quorum returning ALL qualifying words vs the single longest substring present in all, via a generalized suffix tree."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:overrepresented-kmer-discovery
      source: motif-shared-001-evidence
      evidence: "Both are the van Helden oligo-analysis word-enumeration family (Das & Dai 2007): enumerate every fixed-length exact word and count. MOTIF-DISCOVER-001 counts occurrences in ONE sequence and ranks by observed/expected enrichment; MOTIF-SHARED-001 counts how many of a SET of sequences each word appears in (matching-sequence quorum). overrepresented-kmer-discovery names FindSharedMotifs as the cross-sequence sibling."
      confidence: high
      status: current
---

# Shared motifs across a sequence set (fixed-k word enumeration + matching-sequence quorum)

**Shared-motif finding** surfaces the words that recur **across a *set* of sequences**:
it enumerates every fixed-length (`k`) exact word in the input collection and reports each
word present in **at least `minSequences`** of them. This is the **oligo-analysis
"matching sequences"** statistic of van Helden, André & Collado-Vides (1998) — as
implemented by RSAT `oligo-analysis` — where the decision variable is *"the number of
sequences from the input set which contain at least one occurrence of the
oligonucleotide."* Seqeron exposes it as `FindSharedMotifs`. Validated under test unit
**MOTIF-SHARED-001**; the validation record is [[motif-shared-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
artifact pattern.

## The matching-sequence quorum (THE defining rule)

Each input sequence is counted **at most once** per word — presence/absence, not
multiplicity:

- A word's **matching-sequence count** = the number of *distinct* input sequences containing
  **≥ 1** occurrence of it (RSAT: *"only the first occurrence of each sequence is taken into
  consideration"*). A word repeated many times *within a single sequence* still contributes
  exactly **1** to that count.
- A word is a **shared motif** iff its matching-sequence count **≥ `minSequences`** (the
  quorum). Words below the quorum are excluded.
- Each shared motif carries the **set of sequence indices** it was found in and a
  **`Prevalence` = matchingSequences / totalSequences** ∈ (0, 1].

This is distinct from the raw **occurrence** count (RSAT: *"overlapping matches are detected
and summed"*) — the occurrence count sums all hits including within-sequence overlaps, but it
does **not** drive the shared-motif decision.

## Fixed-k, exact words

The method enumerates a **fixed** word length `k` throughout the whole input set and matches
**exactly** — Das & Dai (2007) note the *"greatest shortcoming"* of van Helden oligo-analysis
is that *"there are no variations allowed within an oligonucleotide."* No degenerate/IUPAC
positions, no substitutions: a one-mismatch near-miss is **not** a match. Defaults are
`k = 6`, `minSequences = 2`, both **API ergonomics** (any `k ≥ 1` and `minSequences ≥ 1` is
valid) — RSAT permits any oligo length and any quorum, so the defaults are documented but
treated as caller-supplied parameters, not biological constants.

## Worked oracle

S0 = `ATGATG`, S1 = `ATGCCC`, S2 = `CCCGGG`, with `k = 3`, `minSequences = 2`:

| Word | Found in | Matching-sequence count | Shared? |
|------|----------|-------------------------|---------|
| `ATG` | S0 (pos 0, 3), S1 (pos 0) | 2 → indices {0, 1} | ✔ (Prevalence 2/3) |
| `CCC` | S1 (pos 3), S2 (pos 0) | 2 → indices {1, 2} | ✔ |

Result = `{ATG` → [0,1]`, CCC` → [1,2]`}`. Note `ATG` occurs **twice** in S0 but still
counts **1** toward its matching-sequence total.

## Contract, invariants, and corner cases

| Aspect | Behaviour |
|--------|-----------|
| Inputs | a collection of sequences, `k` (default 6, ≥ 1), `minSequences` (default 2, ≥ 1) |
| Per-motif output | the word, its **`SequenceIndices`** (distinct sequences containing it), matching-sequence count, `Prevalence` = count/total |
| Within-sequence repeats | contribute **1**, not their occurrence multiplicity (matching-sequence semantics) |
| Below quorum | a word in fewer than `minSequences` sequences is **excluded** |
| `k` > shortest sequence | that sequence yields **no** windows — contributes to no word's count |
| Empty collection | no motifs |
| `k < 1` | throws (input-validation contract) |

## Not LCSM — the contrast to nail

The nearby **longest common substring / LCSM** framing ([[longest-common-substring]]) is
deliberately **not** what this unit computes:

| | Shared motifs (this unit, MOTIF-SHARED-001) | LCSM / LCS ([[longest-common-substring]]) |
|---|---|---|
| Word length | **fixed** `k` | **variable** — the maximum length achievable |
| Membership | **quorum**: present in ≥ `minSequences` | present in **all** sequences |
| Output | **all** qualifying words + their sequence sets | a **single** longest common substring |
| Engine | enumerate-count over fixed-`k` windows | generalized suffix tree (deepest common node) |

ROSALIND's own example makes the divergence concrete: on `GATTACA` / `TAGACCA` / `ATACA`,
LCSM returns the single string `AC` (length 2), whereas this unit at `k = 2, minSequences = 3`
reports **every** 2-mer present in all three (e.g. `AC`, `TA`). LCS is the "single longest,
present in all" specialization; shared motifs is the "all words, present in a quorum" family.

## Siblings in the motif / word-enumeration family

- [[overrepresented-kmer-discovery]] — the **single-sequence** member of the same
  van Helden oligo-analysis word-enumeration family: it counts each word's occurrences in
  *one* sequence and ranks by the observed/expected (O/E) enrichment ratio. Shared-motif
  finding instead counts *how many sequences* a word appears in across a **set** (the
  cross-sequence "matching sequences" quorum). Same enumerate-count machinery, different
  count and different question.
- [[known-motif-search]] — matches a caller-supplied set of **known** query motifs by exact
  substring in one subject (motif is an *input*). Shared-motif finding *derives* the motifs
  as *output* from the sequence set.

## Sources and deviations

**RSAT `oligo-analysis` manual** (rank 3, reference implementation of the van Helden method —
the verbatim "matching sequences", "at least one occurrence", "only the first occurrence of
each sequence", and "overlapping matches ... summed" definitions), **Das & Dai 2007** (rank 1,
BMC Bioinformatics survey — places oligo-analysis in the word-enumeration family and states the
exact-word "no variations allowed" limitation), **van Helden, André & Collado-Vides 1998**
(rank 1, J Mol Biol 281(5):827–842 — the named primary; direct article HTTP 403, so the RSAT
manual + survey carry the operative definitions), and **ROSALIND LCSM** (rank 4, cited only to
delineate the alternative longest-common-substring framing this unit does **not** implement).
**No source contradictions.** Two flagged **assumptions**, both presentation/API conveniences
outside the source formula: the default `k = 6` / `minSequences = 2` (in-range but not
prescribed by any source — treated as caller-supplied), and expressing the RSAT raw
matching-sequence count as a `Prevalence` fraction of total sequences (a value in (0, 1],
consistent with the definition but not itself a source formula).
