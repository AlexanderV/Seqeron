---
type: concept
title: "Shannon entropy (per-symbol sequence information content)"
tags: [analysis, algorithm]
mcp_tools:
  - complexity_shannon
  - shannon_entropy
sources:
  - docs/algorithms/Sequence_Composition/Shannon_Entropy.md
source_commit: a444b72a7bdf4b704592af25b5fc00838a9f36cc
created: 2026-07-17
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: shannon-entropy-spec
      evidence: "Test Unit ID: SEQ-ENTROPY-001; Algorithm Group: Sequence Composition; Algorithm: Shannon Entropy for Biological Sequences"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:k-mer-statistics
      source: shannon-entropy-spec
      evidence: "K-mer entropy H_k = −Σ p(kmer) log₂ p(kmer) is the higher-order extension of per-base Shannon entropy; the per-base version is effectively k=1 over the 4-nucleotide alphabet. SequenceComplexity.CalculateKmerEntropy is the k-mer member (SEQ-COMPLEX-KMER-001, home k-mer-statistics)."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:base-composition
      source: shannon-entropy-spec
      evidence: "Shannon entropy is computed over the base-frequency (composition) distribution; base composition is its input distribution."
      confidence: medium
      status: current
    - predicate: relates_to
      object: concept:windowed-sequence-complexity-profile
      source: shannon-entropy-spec
      evidence: "The windowed complexity profile emits a per-window Shannon entropy of base composition (uniform ⇒ 2.0, homopolymer ⇒ 0); this scalar unit is the whole-distribution counterpart of that per-window channel."
      confidence: medium
      status: current
---

# Shannon entropy (per-symbol sequence information content)

**Shannon entropy** measures the information content / uncertainty of a symbol distribution:
`H = −Σ p·log₂ p` in **bits**. Applied to a biological sequence it quantifies how evenly the symbols
are used — a diversity/complexity signal used to flag low-complexity and repetitive regions, to build
sequence logos, and to compare compositionally constrained versus unconstrained regions. It is the
**base per-symbol entropy** member of the sequence complexity/entropy family (siblings
[[sequence-complexity-compression-lempel-ziv]], [[dust-low-complexity-score]], [[linguistic-complexity]],
and the fixed-`k` k-entropy of [[k-mer-statistics]]). Validated as test unit **SEQ-ENTROPY-001**
(CLEAN, PASS/PASS 2026-06-24); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern. The implementation is
[[research-grade-limitations|research-grade]].

## Core model

For a symbol distribution `X` with probabilities `p(xᵢ)`:

```
H(X) = −Σᵢ p(xᵢ) · log₂ p(xᵢ)          (bits, convention 0·log 0 = 0)
```

Probabilities are the observed **relative symbol frequencies** `p(xᵢ) = countᵢ / total`. Base-2
logarithms make the unit **bits**. For DNA over the 4-base alphabet the **maximum entropy is
`log₂ 4 = 2` bits** (reached when all four bases are equiprobable); the minimum is `0` (one symbol).

The **k-mer extension** applies the identical formula to the frequency distribution of overlapping
length-`k` substrings, `H_k = −Σ p(kmer) log₂ p(kmer)`, whose maximum over the full DNA k-mer alphabet
is `log₂(4ᵏ) = 2k` bits. That higher-order form is the entropy member of the complexity family and is
owned by [[k-mer-statistics]] (`SequenceComplexity.CalculateKmerEntropy`, SEQ-COMPLEX-KMER-001); per-base
Shannon entropy here is effectively its **`k = 1`** special case over just the four nucleotides —
composition-only, blind to order.

## Two implementations: DNA-canonical vs general-alphabet

The spec documents **two distinct entry points** with deliberately different alphabet handling:

| Method (location) | Alphabet counted | Range (DNA) | Use when |
|-------------------|------------------|-------------|----------|
| `SequenceComplexity.CalculateShannonEntropy(DnaSequence \| string)` | **only `A/T/G/C`** (canonical DNA; `N` and other ambiguity codes ignored) | `[0, 2]` bits | DNA composition entropy |
| `SequenceStatistics.CalculateShannonEntropy(string)` | **all letters** | alphabet-dependent | non-DNA alphabets or general text |

The canonical `SequenceComplexity` path is DNA-focused: it upper-cases, counts the four standard bases,
converts to probabilities over the **counted total** (so ambiguity codes never contribute), uses base-2
logs, and returns `0.0` when no counted DNA base remains after filtering. The `SequenceStatistics`
alternative counts every letter and is therefore the one to reach for when all symbols must participate
(general-alphabet entropy is deliberately **not implemented** in the canonical DNA path). Both live in
`Seqeron.Genomics.Analysis` (`SequenceComplexity.cs`, `SequenceStatistics.cs`). The k-mer overload
`SequenceComplexity.CalculateKmerEntropy(DnaSequence, int k = 2)` is documented on [[k-mer-statistics]].
The MCP tools `complexity_shannon` (canonical DNA) and `shannon_entropy` (general-purpose) are the thin
wrappers over the two methods (see [[mcp-tool-catalog]]).

## Contract

| Name | Type | Default | Constraints |
|------|------|---------|-------------|
| `sequence` | `DnaSequence` or `string` | required | null `DnaSequence` → `ArgumentNullException`; null/empty string → `0.0` (raw-string path uppercases before analysis) |
| `k` | `int` | `2` | `CalculateKmerEntropy` overload only; typed overload throws `ArgumentOutOfRangeException` when `k < 1`; sequence shorter than `k` → `0.0` |

**Output:** `entropy` (`double`) — Shannon entropy in bits.

## Invariants

- **INV-01** DNA entropy from the canonical implementation lies in **`[0, 2]` bits** (4-symbol
  alphabet counted).
- **INV-02** Empty sequences return **`0.0`** (short-circuits before computing frequencies).
- **INV-03** Homopolymers return **`0.0`** (one symbol has probability `1.0`, all others `0`).

## Reference values (DNA)

| Distribution | Entropy (bits) | Meaning |
|--------------|----------------|---------|
| Uniform (25% each base) | `2.0` | maximum DNA entropy |
| Two bases 50/50 | `1.0` | intermediate |
| Single base 100% | `0.0` | minimum entropy |

## Edge cases

| Input | Result | Rationale |
|-------|--------|-----------|
| Empty sequence | `0.0` | no information content |
| Single repeated base / homopolymer `AAAA` | `0.0` | no uncertainty (INV-03) |
| Equal use of all four bases | `2.0` | maximum DNA uncertainty |
| Sequence length `< k` (k-mer entropy) | `0.0` | no k-mers extractable |

## Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateShannonEntropy` | `O(n)` | `O(1)` | counts four fixed DNA bases |
| `CalculateKmerEntropy` | `O(n)` | `O(u)` | builds a k-mer count dictionary (`u` distinct k-mers) before the entropy sum |

## Relationship to neighbouring measures

- **k-mer k-entropy** — [[k-mer-statistics]] holds the higher-order `H_k` (fixed-`k` frequency
  distribution, max `log₂(4ᵏ)`); per-base Shannon entropy here is its `k = 1` composition-only case.
- **Windowed profile** — [[windowed-sequence-complexity-profile]] emits a per-window Shannon entropy of
  base composition along the sequence (a *profile*); this unit is the whole-distribution scalar. Keep it
  distinct from the Statistics-domain windowed **entropy profile** consumer (`entropy-profile`, still
  pending) that scans this same measure across sliding windows.
- **Other complexity-family scalars** — the compression [[sequence-complexity-compression-lempel-ziv]]
  (variable-length LZ76 phrases, order-sensitive), the fixed-`k=3` [[dust-low-complexity-score]] masker
  (high score ⇒ low complexity — opposite numeric direction), and the vocabulary-usage
  [[linguistic-complexity]] are *different* scalars over the same sequence; a perfectly periodic string
  has flat composition entropy per window yet non-trivial LZ/linguistic structure. All are low exactly
  where low-complexity/repeat tracts are.
- **Composition** — [[base-composition]] is the input frequency distribution this entropy reduces;
  entropy is its evenness summary. The general-purpose `SequenceStatistics.CalculateShannonEntropy` is
  also one of the fields bundled by `SummarizeNucleotideSequence` (SEQ-SUMMARY-001).

## Applications

- Low-complexity region detection (entropy-like masking, DUST/SEG-style workflows).
- Sequence-logo information content: `2 − H` bits per DNA position (Schneider & Stephens 1990).
- Comparing coding vs non-coding composition patterns.

## References

Shannon, C.E. (1948) "A Mathematical Theory of Communication", *Bell System Technical Journal*
27(3):379–423. Cover, T.M. & Thomas, J.A. (1991) *Elements of Information Theory*, Wiley.
Schneider, T.D. & Stephens, R.M. (1990) "Sequence Logos: A New Way to Display Consensus Sequences",
*Nucleic Acids Res.* 18(20):6097–6100.
