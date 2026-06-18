# Average Nucleotide Identity (ANI)

| Field | Value |
|-------|-------|
| Algorithm Group | Comparative Genomics |
| Test Unit ID | COMPGEN-ANI-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Average Nucleotide Identity (ANI) measures whole-genome relatedness between two prokaryotic
genomes as the mean nucleotide identity of their conserved (alignable) genomic fragments [1]. It
is the in-silico replacement for laboratory DNA-DNA hybridization (DDH): an ANI of ≈95 %
corresponds to the 70 % DDH species boundary [1], with ≈94 % reported by an independent study [2].
The query genome is cut into fixed-length fragments, each fragment is aligned to the reference,
and qualifying fragments (passing identity and alignable-length cut-offs) are averaged [1]. The
method is heuristic with respect to alignment but the averaging formula is exact.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Bacterial/archaeal species have traditionally been delineated by DDH ≥ 70 %. DDH fragments
genomic DNA to ~1 kb pieces and measures cross-hybridization. ANI reproduces this in silico by
fragmenting one genome into ~1 kb pieces and aligning them to the other, then averaging the
identities of the conserved fragments [1][2].

### 2.2 Core Model

Let the query genome `Q` be cut into consecutive non-overlapping fragments `f_1, …, f_n` of
length `L` (Goris et al. use `L = 1020` nt) [1]. For each fragment `f_i` let `a_i` be its best
alignment against the reference `R`, with identity recalculated over the whole fragment length:

```
id(f_i) = (matching bases of f_i in a_i) / L
cov(f_i) = (alignable length of f_i) / L
```

A fragment qualifies iff `id(f_i) > 0.30` and `cov(f_i) ≥ 0.70` [1]. Then

```
ANI(Q, R) = mean over qualifying f_i of id(f_i)
```

— "the mean identity of all BLASTN matches that showed more than 30 % overall sequence identity
(recalculated to an identity along the entire sequence) over an alignable region of at least 70 %
of their length" [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 ≤ ANI ≤ 1` | each `id(f_i) ∈ [0,1]`; the mean preserves the bound [1] |
| INV-02 | Identical genomes → ANI = 1.0 | every fragment is a perfect substring, `id = 1` [1] |
| INV-03 | Only fragments with `id > 0.30` and `cov ≥ 0.70` contribute | Goris cut-off clause [1] |
| INV-04 | Fragmentation is consecutive, non-overlapping; trailing partial fragment (< L) is dropped | "consecutive 1020 nt fragments" [1] |
| INV-05 | No qualifying fragment / empty / null input → 0 | mean over the empty set is reported as 0 |

### 2.5 Comparison with Related Methods

| Aspect | ANIb (this) | OrthoANI [3] |
|--------|-------------|--------------|
| Fragment matching | best one-way BLASTN hit | reciprocal orthologous fragment pairs |
| Both genomes fragmented | no (query only) | yes |
| Symmetry | direction-dependent | symmetric |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `querySequence` | `string` | required | Query genome (the one fragmented) | DNA; case-insensitive |
| `referenceSequence` | `string` | required | Reference genome searched against | DNA; case-insensitive |
| `fragmentLength` | `int` | 1020 | Fragment length in nt [1] | > 0 |
| `minIdentity` | `double` | 0.30 | Min per-fragment identity (over fragment length) [1] | [0,1] |
| `minAlignableFraction` | `double` | 0.70 | Min fraction of fragment that must align [1] | [0,1] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `double` | ANI as a fraction in [0, 1]: mean identity of qualifying fragments, else 0 |

### 3.3 Preconditions and Validation

Null or empty `querySequence`/`referenceSequence` → 0. `fragmentLength ≤ 0` →
`ArgumentOutOfRangeException`. Input is uppercased (`ToUpperInvariant`) before comparison; the
DNA alphabet is assumed (no IUPAC degeneracy handling). Coordinates are 0-based internally.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; uppercase both sequences.
2. Cut the query into consecutive non-overlapping fragments of `fragmentLength` (drop trailing partial).
3. For each fragment, find the best ungapped full-length placement in the reference; compute
   identity = matching bases / fragmentLength, and alignable fraction.
4. Keep fragments with `identity > minIdentity` and `alignableFraction ≥ minAlignableFraction`.
5. Return the mean identity of the kept fragments (0 if none).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Constant | Value | Source |
|----------|-------|--------|
| Fragment length | 1020 nt | Goris et al. 2007 [1] |
| Min identity | 0.30 | Goris et al. 2007 [1] |
| Min alignable fraction | 0.70 | Goris et al. 2007 [1] |
| Species boundary (interpretation only) | ≈95 % ANI | Goris et al. 2007 [1]; ≈94 % [2] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateANI | O((n/L) × m × L) = O(n × m) | O(1) | n = query length, m = reference length, L = fragmentLength; ungapped sliding placement per fragment |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ComparativeGenomics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs)

- `ComparativeGenomics.CalculateANI(query, reference, fragmentLength, minIdentity, minAlignableFraction)`: public ANI entry point.
- `ComparativeGenomics.BestUngappedFragmentMatch(...)` (private): best ungapped full-length placement and identity for one fragment.

### 5.2 Current Behavior

Fragmentation is consecutive and non-overlapping (a trailing fragment shorter than
`fragmentLength` is ignored). Per fragment, the implementation slides the fragment across every
full-length offset of the reference and counts matching bases at the best offset; identity is that
count divided by the fragment length. Because the placement is full-length, the alignable fraction
is 1.0 whenever the reference is at least as long as the fragment, and 0 otherwise (reference too
short → fragment excluded by the 70 % cut-off).

**Search-reuse decision:** the repository suffix tree (`SuffixTree`) supports *exact* occurrence
enumeration only; ANI requires *identity over mismatches* (approximate, scoring-based) which the
suffix tree does not provide. The previous implementation misused `LongestCommonSubstring` as an
identity proxy, which is not nucleotide identity and was a defect. A direct best-placement scan is
the correct primitive for per-fragment identity, so the suffix tree is not used here.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Consecutive non-overlapping query fragmentation at a fixed length (default 1020 nt) [1].
- Per-fragment identity recalculated over the whole fragment length [1].
- ANI = mean of fragment identities passing `> 30 %` identity and `≥ 70 %` alignable cut-offs [1].
- Default fragment length and cut-offs taken directly from Goris et al. 2007 [1].

**Intentionally simplified:**

- Per-fragment alignment uses an **ungapped** best full-length placement instead of gapped BLASTN
  (X=150, q=-1) [1]; **consequence:** indel-containing homologous regions score lower identity
  than a gapped aligner would report. For substitution-only divergence the recalculated-over-fragment
  identity matches the definition.

**Not implemented:**

- Reciprocal/averaged two-direction ANI and the OrthoANI orthologous-pair refinement [3];
  **users should rely on:** dedicated tools (pyani, OrthoANI, FastANI) for taxonomy-grade values.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Ungapped placement vs gapped BLASTN | Assumption | Lower identity for indel regions; exact for substitutions | accepted | §5.3 simplified; Evidence ASSUMPTION 1 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null / empty input | return 0 | Validation contract |
| `fragmentLength ≤ 0` | `ArgumentOutOfRangeException` | Validation contract |
| Query shorter than `fragmentLength` | return 0 | No full fragment fits |
| Reference shorter than `fragmentLength` | return 0 | No fragment can align ≥ 70 % (INV-03) |
| No fragment passes cut-offs | return 0 | Mean over empty set (INV-05) |

### 6.2 Limitations

Does not model gaps/indels (ungapped placement), IUPAC degenerate bases, reverse-complement
matching, or reciprocal averaging. Not suitable for taxonomy-grade ANI on real genomes; intended
as an exact realization of the ANI averaging formula over conserved fragments.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Identical genomes → ANI = 1.0
double ani = ComparativeGenomics.CalculateANI(
    "AAAACCCCGGGGTTTT", "AAAACCCCGGGGTTTT", fragmentLength: 4);
// One substitution in the last fragment ("TTTA"): (1+1+1+0.75)/4 = 0.9375
double ani2 = ComparativeGenomics.CalculateANI(
    "AAAACCCCGGGGTTTA", "AAAACCCCGGGGTTTT", fragmentLength: 4);
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ComparativeGenomics_CalculateANI_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_CalculateANI_Tests.cs) — covers INV-01..INV-05
- Evidence: [COMPGEN-ANI-001-Evidence.md](../../../docs/Evidence/COMPGEN-ANI-001-Evidence.md)
- TestSpec: [COMPGEN-ANI-001.md](../../../tests/TestSpecs/COMPGEN-ANI-001.md)

## 8. References

1. Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, Tiedje JM. 2007. DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. Int J Syst Evol Microbiol 57(1):81-91. https://doi.org/10.1099/ijs.0.64483-0
2. Konstantinidis KT, Tiedje JM. 2005. Genomic insights that advance the species definition for prokaryotes. Proc Natl Acad Sci USA 102(7):2567-2572. https://doi.org/10.1073/pnas.0409727102
3. Lee I, Ouk Kim Y, Park S-C, Chun J. 2016. OrthoANI: An improved algorithm and software for calculating average nucleotide identity. Int J Syst Evol Microbiol 66(2):1100-1103. https://doi.org/10.1099/ijsem.0.000760
4. Pritchard L, et al. pyani ANIb documentation. https://pyani.readthedocs.io/en/latest/api/pyani.anib.html
