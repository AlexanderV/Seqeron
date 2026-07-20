---
type: concept
title: "Linguistic complexity (vocabulary-usage ratio)"
tags: [analysis, algorithm]
mcp_tools:
  - complexity_linguistic
  - linguistic_complexity
sources:
  - docs/algorithms/Sequence_Composition/Linguistic_Complexity.md
source_commit: e961cd9f4ee3fb5796e5ae9593ec6d74c99df8b2
created: 2026-07-17
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: linguistic-complexity
      evidence: "Test Unit ID: SEQ-COMPLEX-001; Algorithm: Linguistic Complexity (summation-variant vocabulary-usage ratio)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:windowed-sequence-complexity-profile
      source: linguistic-complexity
      evidence: "CalculateWindowedComplexity delegates per-window LC to the same CalculateLinguisticComplexityCore; the windowed profile emits this LC per window"
      confidence: high
      status: current
---

# Linguistic complexity (vocabulary-usage ratio)

The **whole-sequence scalar** measure of how much of a sequence's *possible* vocabulary it
actually uses: the ratio of **distinct observed subwords** to the **maximum number that could
occur**. Low values flag simple/tandem repeats and other low-complexity tracts; high values track
diverse, coding-like content (Trifonov 1990; Gabrielian & Bolshoy 1999; Troyanskaya et al. 2002).
Validated as test unit **SEQ-COMPLEX-001** (CLEAN, PASS/PASS 2026-06-24); [[test-unit-registry]]
tracks the unit and [[algorithm-validation-evidence]] describes the artifact pattern. The
implementation is [[research-grade-limitations|research-grade]] and **DNA-oriented** (hard-coded
alphabet size 4).

## Where it sits in the complexity family

This is the **scalar linguistic-complexity member** of the `SEQ-COMPLEX-*` sequence
complexity / entropy family — a *vocabulary-usage ratio*, distinct from every sibling measure:

- vs. **compression complexity** — [[sequence-complexity-compression-lempel-ziv]] (LZ76) counts
  variable-length phrases in a left-to-right exhaustive-history parse and normalizes by
  `n/log_b n`. LC instead counts **distinct fixed-length subwords for each length `1..m`** and
  divides by the *combinatorial* maximum vocabulary. Both are low on repetitive input, but LC's
  denominator is `min(4^i, N−i+1)`, not a randomness bound.
- vs. **DUST triplet score** — [[dust-low-complexity-score]] fixes `k = 3` and sums
  `∑ c(c−1)/2 / (L−2)` where a *high* score means *low* complexity (opposite direction). LC
  aggregates over word lengths `1..m` and a *high* LC means *high* complexity.
- vs. **k-mer k-entropy** — [[k-mer-statistics]] (`CalculateKmerEntropy`, SEQ-COMPLEX-KMER-001)
  reduces a fixed-`k` count profile to Shannon entropy `−Σ pᵢ log₂ pᵢ`. LC uses **presence/absence
  of distinct subwords** (vocabulary count `Vᵢ`), not their frequency distribution, and sums across
  several lengths.
- vs. **windowed profile** — [[windowed-sequence-complexity-profile]] (SEQ-COMPLEX-WINDOW-001) is
  the *profiling* member: it runs **this same LC per sliding window** (delegating to the shared
  `CalculateLinguisticComplexityCore`) alongside a per-window Shannon entropy, emitting a
  `ComplexityPoint` profile. This unit is the **standalone whole-sequence scalar** version.
- vs. **protein SEG** — [[protein-low-complexity-seg]] is the amino-acid-side low-complexity
  detector (Shannon entropy of composition); LC is the DNA vocabulary-richness analogue.

## Definition — the summation variant

Seqeron implements the **summation form** over word lengths `1 … m`:

```
LC = ( Σᵢ Vᵢ ) / ( Σᵢ Vmax,i )        i = 1 … m,   m = min(maxWordLength, N)
Vmax,i = min(K^i, N − i + 1),          K = 4 (DNA alphabet, hard-coded)
```

- `Vᵢ` = number of **distinct** length-`i` subwords observed (overlapping windows), counted with a
  `HashSet<string>` per length.
- `Vmax,i` = the maximum distinct subwords of length `i` that *could* appear, bounded by both the
  **alphabet** (`K^i`) and **positional availability** (`N − i + 1` windows).
- `m` = the effective maximum word length, capped at the sequence length.

There is **one aggregate LC scalar per call** (summed across all lengths), not a per-`k` vector; the
`maxWordLength` parameter (default `10`) controls the upper length included in the summation. This is
the direct definition with hash-based subword enumeration — **not** the suffix-tree optimization of
Troyanskaya et al. (2002), which is deliberately not implemented.

## Method surface and contract

Primary spec `docs/algorithms/Sequence_Composition/Linguistic_Complexity.md`; implementation
`SequenceComplexity.CalculateLinguisticComplexity` on `Seqeron.Genomics.Analysis`
(`SequenceComplexity.cs`), two overloads sharing `CalculateLinguisticComplexityCore`:

| Overload | Returns | Validation |
|----------|---------|------------|
| `CalculateLinguisticComplexity(DnaSequence, int maxWordLength = 10)` | `double` LC | null ⇒ `ArgumentNullException`; `maxWordLength < 1` ⇒ `ArgumentOutOfRangeException` |
| `CalculateLinguisticComplexity(string, int maxWordLength = 10)` | `double` LC | null/empty ⇒ `0`; uppercases input; **no alphabet validation** |

The MCP tools `complexity_linguistic` / `linguistic_complexity` are the thin wrappers over the
scalar method (see [[mcp-tool-catalog]]). The per-window sibling `windowed_complexity` belongs to
[[windowed-sequence-complexity-profile]].

## Invariants

- **INV-01** `0 ≤ LC ≤ 1` for DNA-alphabet inputs — observed distinct counts cannot exceed the
  `K = 4` theoretical maximum when the input matches the hard-coded denominator.
- **INV-02** empty sequence ⇒ `0` (short-circuits before accumulating).
- **INV-03** word lengths capped at `min(maxWordLength, N)` (loop bound `m`).

## Edge cases

| Case | Behavior | Rationale |
|------|----------|-----------|
| Empty string | `0` | no subword vocabulary |
| Null `DnaSequence` | `ArgumentNullException` | explicit guard |
| `maxWordLength < 1` (typed) | `ArgumentOutOfRangeException` | explicit validation |
| Single nucleotide `A` | positive value | one distinct 1-mer exists |
| Homopolymer (`AAAA…`) | low value | only one word per length observed |
| Random-like sequence | high value (→ 1) | observed vocabulary approaches the maximum |
| Raw-string, non-ACGT symbols | **may exceed** the `[0,1]` DNA interpretation | denominator stays `4^i`; observed words can include out-of-alphabet symbols |

**Worked oracle** (from the windowed sibling, same core): `ACGTACGT` (N = 8, all four bases)
→ `Σ Vᵢ = 23`, `Σ Vmax,i = 29` ⇒ **LC = 23/29 ≈ 0.7931**; `AAAAAAAA` poly-A → `Σ Vᵢ = 6`,
`Σ Vmax,i = 29` ⇒ **LC = 6/29 ≈ 0.2069**. A summary cross-check reports `ATGCATGC` →
`CalculateLinguisticComplexity ≈ 0.83968…`.

## Complexity and limitations

Time `O(n · k²)` effective, space `O(u)` (distinct-subword count): the direct path allocates and
hashes fresh substrings of lengths `1..k` for each window rather than indexing a suffix tree. The
metric is **DNA-specific** — the denominator is fixed to alphabet size 4, so the raw-string overload
accepts arbitrary uppercase symbols without reconciling the denominator, and callers should treat
non-ACGT results as DNA-oriented (potentially outside `[0,1]`). Suffix-tree acceleration
(Troyanskaya et al. 2002) is not implemented; direct enumeration is authoritative here.

## Applications

Low-complexity region detection in repetitive DNA; characterizing microsatellites, tandem repeats,
and palindrome/hairpin-rich segments; contrasting low-complexity vs coding-like regions. Consumed
per-window by [[windowed-sequence-complexity-profile]] and complementary to the explicit-repeat
[[repetitive-element-detection]].

## References

Trifonov E.N. (1990) *Making sense of the human genome*; Troyanskaya O.G. et al. (2002)
*Bioinformatics* 18(5):679–688 (PMID 12050064); Orlov Y.L. & Potapov V.N. (2004) *Nucleic Acids
Research* 32(W):W628–W633; Gabrielian A. & Bolshoy A. (1999) *Computers & Chemistry* 23(3–4):263–274.
