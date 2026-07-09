---
type: concept
title: "De novo motif discovery via over-represented k-mers (observed/expected enrichment)"
tags: [motif, algorithm]
sources:
  - docs/Evidence/MOTIF-DISCOVER-001-Evidence.md
  - docs/algorithms/Motif_Discovery/Overrepresented_Kmer_Discovery.md
source_commit: ca2ab6940b31b0df620d2b47f5d84f4cc49652d5
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: motif-discover-001-evidence
      evidence: "Test Unit ID: MOTIF-DISCOVER-001 ... Algorithm: Motif Discovery via Overrepresented k-mers (observed/expected enrichment)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:known-motif-search
      source: motif-discover-001-evidence
      evidence: "Discovery finds UNKNOWN over-represented k-mers by enumerate-count-rank (no query supplied); known-motif search matches a caller-supplied set of KNOWN motifs by exact substring. Both share overlapping all-occurrences window enumeration but answer opposite questions (which words recur? vs where is this word?)."
      confidence: high
      status: current
---

# De novo motif discovery via over-represented k-mers

**De novo motif discovery** surfaces *unknown* candidate motifs in a single DNA
sequence by **enumerating every length-`k` k-mer, counting its occurrences, and ranking
by how far the observed count exceeds the count expected by chance** — the observed/expected
(O/E) enrichment ratio. It is the "you do **not** know what you are looking for" end of
motif analysis: no query motif is supplied. Seqeron exposes it as
`MotifFinder.DiscoverMotifs(DnaSequence, k, minCount)`. Validated under test unit
**MOTIF-DISCOVER-001**; the validation record is [[motif-discover-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes
the artifact pattern.

## The model: observed vs expected under a uniform background

For a sequence of length `N` there are `N − k + 1` length-`k` windows. Under a
**zero-order i.i.d. uniform** background (each of A/C/G/T drawn independently with
probability `1/4`, Compeau & Pevzner), the **expected** number of occurrences of any
specific k-mer is

```
E = (N − k + 1) / 4^k
```

— the window count divided by the `4^k` distinct DNA k-mers. The **overrepresentation**
(enrichment) of a k-mer observed `c` times is the **observed/expected ratio**

```
enrichment = c / E = c / ((N − k + 1) / 4^k)
```

A value **> 1** means the k-mer recurs more often than random DNA predicts — the signal
that regulatory motifs (TF binding sites, etc.) tend to be locally over-represented. The
method is **deterministic and exact** (no sampling): it does one linear pass with a
hash-map from k-mer → start-position list, so each returned record carries its `Count`,
its **0-based** `Positions`, and its `Enrichment`.

## Worked oracles

- **Tandem repeat** `ATGCATGCATGC` (N=12, k=4): windows = 9, `E = 9/4⁴ = 9/256 =
  0.03515625`; `ATGC` occurs at positions **0, 4, 8** (Count 3) → enrichment
  `3 / (9/256) = 768/9 ≈ 85.333`.
- **Homopolymer** `AAAAAAAAAA` (N=10, k=3): windows = 8, `E = 8/4³ = 0.125`; `AAA`
  occurs at positions **0..7** (Count 8) → enrichment `8 / 0.125 = 64.0`.

## Contract, invariants, and corner cases

| Aspect | Behaviour |
|--------|-----------|
| Inputs | `sequence` (non-null), `k` (default 6, ≥1), `minCount` (default 2) |
| Record | `Sequence` (the k-mer), `Count`, `Positions` (0-based window starts, all overlaps), `Enrichment` = Count/E |
| Overlaps | counted at every window (INV-01) — `AAA` in `AAAAAAAAAA` is 8, not fewer |
| `minCount` filter | **presentation-only**: excludes k-mers with Count < minCount; it never alters a returned record's Count/Positions/Enrichment (the formal O/E is defined for every k-mer) |
| Enrichment > 0 | always, since `E > 0` whenever any k-mer was counted (INV-04) |
| `k > N` | empty result — no length-`k` windows exist |
| Null sequence / `k < 1` | `ArgumentNullException` / `ArgumentOutOfRangeException` |
| Order | hash-map enumeration order; caller sorts if needed |

## Scope, deviations, and siblings

This is the **enumerative, exact-word, single-sequence** discovery method — the counting
family of de novo discovery. It is a genuinely **distinct** operation from the two sibling
motif concepts and is modelled `alternative_to` the known-motif matcher:

- [[known-motif-search]] — matches a caller-supplied set of **known** motifs by exact
  substring (the "where is this word?" question). Discovery answers the opposite: "which
  words recur more than chance?" — no query. Both share the overlapping all-occurrences
  window enumeration; they differ in whether the motif is an input or an output.
- [[consensus-from-alignment]] — collapses an already-aligned motif instance set to one
  plurality consensus. Discovery operates on a *single unaligned* sequence and produces
  ranked candidate words, not a consensus.

**Intentionally simplified / not implemented** (per the evidence + algorithm doc):

- **Background is zero-order uniform only.** Higher-order Markov backgrounds (as in
  monaLisa `getKmerFreq`) are not modelled — on compositionally biased (e.g. GC-rich)
  sequences the O/E ratio over/under-states true enrichment (ASM-01).
- **No statistical p-value / E-value.** The closed-form probability `Pr(N,4,k,t)` for ≥ t
  occurrences (Compeau & Pevzner) is *not* computed — it is an approximation that ignores
  self-overlap, and only affects a *probability* statistic, not the exact Count or the O/E
  denominator. Rank by the deterministic Count / Enrichment, or use an external
  significance tool.
- **Single-sequence, DNA-only.** Cross-sequence shared-motif discovery is a separate unit,
  [[shared-motifs]] (`FindSharedMotifs`, MOTIF-SHARED-001 — the same fixed-`k` exact
  word-enumeration family, but counting *how many sequences of a set* contain each word via
  the van Helden / RSAT "matching sequences" quorum instead of one sequence's O/E enrichment);
  other de novo families (greedy / median-string / Gibbs-sampling motif search) are separate
  methods not implemented here.

**No source contradictions.** The only assumption is that `minCount` is a presentation
threshold, not part of the published statistic — it changes which rows appear, never a
returned record's correctness-affecting fields. Compeau & Pevzner (expected-count formula,
worked `ATGC` example) and the monaLisa / PeerJ O/E-ratio corroboration agree that
`E = length/4^k` and enrichment = O/E under the uniform model.
