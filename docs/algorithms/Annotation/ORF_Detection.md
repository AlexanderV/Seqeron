# ORF Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Annotation |
| Test Unit ID | ANNOT-ORF-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Open reading frame (ORF) detection identifies contiguous coding candidates in DNA by finding in-frame start-to-stop intervals. In this repository, the canonical implementation is `GenomeAnnotator.FindOrfs(...)`, which scans the forward strand and, by default, the reverse complement as well. The method is deterministic and exact with respect to its fixed start-codon and stop-codon sets, but it is still a simplified coding-region heuristic rather than a trained gene finder. Related wrappers exist elsewhere in the repository, but they do not all share the same parameter semantics.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An ORF is a continuous stretch of codons that begins with a start codon and ends with an in-frame stop codon, with no intervening in-frame stop between those boundaries. Because DNA can be read in three frames on each strand, double-stranded DNA exposes six possible reading frames: `+1`, `+2`, `+3`, `-1`, `-2`, and `-3`. The current document and repository evidence use the common prokaryotic start codons `ATG`, `GTG`, and `TTG`, and the standard stop codons `TAA`, `TAG`, and `TGA`. The background references also note that minimum ORF-length thresholds are heuristic and context-dependent, with common practical cutoffs around 100 codons for gene-finding workflows.

### 2.2 Core Model

For a start codon beginning at nucleotide index $s$ and an in-frame stop codon beginning at index $t$, an ORF spans the half-open interval $[s, t + 3)$. The translated amino-acid length before the terminal stop is:

$$
L_{aa} = \frac{t - s}{3}
$$

under the condition that $t - s$ is divisible by $3$. ORF enumeration therefore consists of scanning each frame in codon steps, tracking candidate start positions, and emitting every candidate that reaches a valid in-frame stop while satisfying the minimum-length threshold.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every returned ORF begins with `ATG`, `GTG`, or `TTG` when `requireStartCodon = true` | The canonical implementation only adds pending starts from that fixed start-codon set |
| INV-02 | Every returned ORF ends with `TAA`, `TAG`, or `TGA` when `requireStartCodon = true` | ORFs are emitted only when the scan reaches a stop codon from the fixed stop-codon set |
| INV-03 | Every returned ORF span is divisible by `3` nucleotides | The scan advances in codon-sized steps within a fixed frame |
| INV-04 | Returned coordinates satisfy `0 <= Start < End <= dnaSequence.Length` | Forward coordinates are emitted directly and reverse-strand coordinates are remapped back to the original sequence bounds |
| INV-05 | When both strands are searched, reported frame labels are `1`, `2`, `3` on the forward strand and `-1`, `-2`, `-3` in the longest-per-frame view | The repository groups reverse-complement ORFs under negative frame keys in `FindLongestOrfsPerFrame(...)` |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `dnaSequence` | `string` | required | DNA sequence to scan for ORFs | Case-insensitive; null or empty input yields no ORFs in the canonical API |
| `minLength` | `int` | `100` | Minimum ORF length in amino acids for `FindOrfs(...)` | Compared against amino-acid length excluding the terminal stop |
| `searchBothStrands` | `bool` | `true` | Whether to scan the reverse complement in addition to the forward strand | When `false`, only frames `1`, `2`, and `3` are searched |
| `requireStartCodon` | `bool` | `true` | Whether a start codon is required before an ORF can be emitted | When `false`, trailing frame segments without a terminating stop may also be emitted if they satisfy `minLength` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | 0-based inclusive ORF start coordinate relative to the original input sequence |
| `End` | `int` | 0-based exclusive ORF end coordinate relative to the original input sequence |
| `Frame` | `int` | Reading frame number `1`, `2`, or `3` in `OpenReadingFrame`; `FindLongestOrfsPerFrame(...)` uses negative keys for reverse-complement frames |
| `IsReverseComplement` | `bool` | Indicates whether the ORF was found on the reverse-complement strand |
| `Sequence` | `string` | DNA sequence of the ORF, including the terminal stop codon when present |
| `ProteinSequence` | `string` | Translation of `Sequence`; tests confirm the canonical implementation includes the terminal `*` stop symbol |

### 3.3 Preconditions and Validation

`GenomeAnnotator.FindOrfs(...)` returns an empty sequence for null or empty input. It uppercases each scanned codon internally, so lowercase and mixed-case DNA are accepted. The canonical ORF finder does not pre-validate the alphabet; codons containing characters such as `N` simply fail start/stop matching and are treated as ordinary non-boundary codons. When `requireStartCodon = true`, an ORF without a terminating stop codon is not emitted. Reverse-strand results are remapped to the coordinate system of the original forward input before they are returned.

## 4. Algorithm

### 4.1 High-Level Steps

1. If `dnaSequence` is null or empty, return no ORFs.
2. Scan each requested reading frame in codon steps and record every encountered start codon as a pending ORF start.
3. When an in-frame stop codon is found, emit every pending start in that frame whose stop-delimited amino-acid length is at least `minLength`, then clear the pending-start list.
4. If `searchBothStrands = true`, repeat the scan on the reverse complement and remap the resulting coordinates back to the original sequence.
5. For `FindLongestOrfsPerFrame(...)`, group emitted ORFs by frame and keep the longest translated result in each frame bucket.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The canonical implementation uses two fixed codon sets:

| Category | Values |
|----------|--------|
| Start codons | `ATG`, `GTG`, `TTG` |
| Stop codons | `TAA`, `TAG`, `TGA` |

Within each frame, `GenomeAnnotator.FindOrfsInFrame(...)` keeps a list of pending start positions rather than a single active start. This allows nested or overlapping ORFs that share the same terminating stop codon to be emitted independently when they all meet the length threshold.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindOrfs(...)` | `O(n)` | `O(p)` plus output | `n` = sequence length; `p` = pending start positions in the current frame; forward and reverse scans only change the constant factor |
| `FindLongestOrfsPerFrame(...)` | `O(n + m log m)` | `O(m)` | `m` = number of emitted ORFs materialized before longest-per-frame selection |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomeAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs), [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs), [Translator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs)

- `GenomeAnnotator.FindOrfs(string, int, bool, bool)`: Canonical ORF finder used by the annotation layer.
- `GenomeAnnotator.FindLongestOrfsPerFrame(string, bool)`: Returns the longest ORF found in each requested frame.
- `GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int)`: Alternate helper with different semantics from the canonical API.
- `Translator.FindOrfs(DnaSequence, GeneticCode?, int, bool)`: Wrapper surface that accepts an explicit genetic code.

### 5.2 Current Behavior

Repository-specific behavior confirmed by source and tests:

- `GenomeAnnotator.FindOrfs(...)` always translates with `GeneticCode.Standard`.
- Returned `ProteinSequence` values include the terminal `*` stop symbol for stop-terminated ORFs.
- Multiple start codons before the same stop codon produce multiple ORFs, so nested ORFs can be returned from a single frame.
- `FindLongestOrfsPerFrame(...)` internally calls `FindOrfs(...)` with `minLength: 1` and `requireStartCodon: true` before selecting the longest translated ORF in each frame.
- `GenomicAnalyzer.FindOpenReadingFrames(...)` is not equivalent to the canonical API: it recognizes only `ATG` as a start codon, always searches both strands, and compares its `minLength` argument to nucleotide span in the source/tests rather than to amino-acid length.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Start-to-stop ORF detection within fixed reading frames.
- Six-frame search when both strands are requested.
- Recognition of the standard stop codons `TAA`, `TAG`, and `TGA`.
- Recognition of the common prokaryotic start codons `ATG`, `GTG`, and `TTG` in the canonical annotation API.

**Intentionally simplified:**

- The canonical annotation API hard-codes the standard genetic code; **consequence:** callers cannot change codon-translation behavior through `GenomeAnnotator.FindOrfs(...)`.
- ORF detection is based only on start/stop structure and minimum length; **consequence:** coding potential, codon bias, and organism-specific evidence do not influence which ORFs are returned.
- Alternate APIs in `GenomicAnalyzer` and `Translator` are not normalized to the same contract; **consequence:** wrapper results may differ from the canonical annotation surface for the same sequence and threshold values.

**Not implemented:**

- Configurable genetic-code selection in the canonical annotation API; **users should rely on:** [Translator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs).
- Organism-specific start-codon sets or trained coding-potential models; **users should rely on:** no current alternative in this repository.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Alternate ORF surfaces are not contract-equivalent | Deviation | `GenomeAnnotator`, `GenomicAnalyzer`, and `Translator` can return different results for the same input and threshold values | accepted | Most visible differences are start-codon handling and `minLength` semantics |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null or empty sequence | Returns no ORFs | `FindOrfs(...)` checks `string.IsNullOrEmpty(...)` and `yield break`s |
| Very short sequence (`Length < 3`) | Returns no ORFs | No complete codon can be scanned |
| No valid start codon with `requireStartCodon = true` | Returns no ORFs | No pending ORF start is ever recorded |
| Valid start but no terminating stop with `requireStartCodon = true` | Returns no ORFs | Stop codon is required before emission |
| Lowercase or mixed-case DNA | Handled identically to uppercase DNA | The implementation uppercases codons before set membership checks |
| `N` inside a candidate start or stop codon | That codon is not treated as a start or stop | The fixed codon sets contain only exact `A/C/G/T` triplets |
| `N` inside an internal coding codon | ORF can continue until a recognized stop codon is found | Only start/stop boundaries are matched against fixed codon sets |
| Nested start codons before the same stop | Multiple ORFs can be emitted | Pending starts are accumulated in a list until the stop codon clears them |

### 6.2 Limitations

This implementation is an ORF enumerator, not a full gene-prediction model. The canonical API does not estimate coding likelihood, does not incorporate codon-usage bias, and does not expose a configurable genetic code. The repository also contains alternate ORF-related entry points whose semantics differ from the canonical `GenomeAnnotator` contract, so callers should choose the entry point deliberately.

## 8. References

1. Wikipedia contributors. Open reading frame. https://en.wikipedia.org/wiki/Open_reading_frame
2. Rosalind. Open Reading Frames. https://rosalind.info/problems/orf/
3. NCBI ORF Finder. https://www.ncbi.nlm.nih.gov/orffinder/
4. Deonier R, Tavare S, Waterman M. Computational Genome Analysis: An Introduction. Springer-Verlag, 2005.
5. Claverie JM. Computational methods for the identification of genes in vertebrate genomic sequences. Human Molecular Genetics. 1997;6(10):1735-1744.
6. Sieber P, Platzer M, Schuster S. The Definition of Open Reading Frame Revisited. Trends in Genetics. 2018;34(3):167-170.
