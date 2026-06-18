# Contig Merging (Suffix–Prefix Overlap Collapse)

| Field | Value |
|-------|-------|
| Algorithm Group | Extended Assembly |
| Test Unit ID | ASSEMBLY-MERGE-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Contig merging joins two sequences whose suffix/prefix already overlap into a single longer
sequence (a superstring), keeping exactly one copy of the shared region. It is the elementary
"collapse" step used in the consensus stage of overlap-layout-consensus (OLC) assembly and in
shortest-common-superstring (greedy SCS) assembly, where each round merges the pair with the
longest overlap. The operation is exact and deterministic: given two strings and a known overlap
length `l`, the result is `contig1` concatenated with `contig2` minus its length-`l` prefix [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In de novo assembly, reads (and later contigs) that originate from overlapping stretches of the
genome share a suffix/prefix. The overlap graph places a directed edge from X to Y whose weight is
the length of the longest suffix of X that exactly matches a prefix of Y; assembling the genome
means collapsing chains of such overlaps into longer contiguous sequences ("fragments are contigs,
short for contiguous") [2]. This unit covers only the per-pair collapse arithmetic, not overlap
discovery or layout.

### 2.2 Core Model

An overlap is "a length-`l` suffix of X [that] matches a length-`l` prefix of Y, where `l` is
given" [1]. To merge two overlapping strings into a superstring, keep one copy of the overlapping
region:

```
merge(X, Y, l) = X + Y[l:]        (the length-l prefix of Y is dropped)
|merge(X, Y, l)| = |X| + |Y| − l
```

When there is no qualifying overlap the strings are simply concatenated: "without requirement of
'shortest', it's easy: just concatenate them" — i.e. `merge(X, Y, 0) = X + Y` [1]. A valid overlap
length is bounded by the shorter string, because it is simultaneously a suffix of X and a prefix of
Y (`suffixPrefixMatch` guards `if len(x) < k or len(y) < k: return 0`) [1]. When multiple overlaps
exist, only the longest suffix/prefix match is used [2].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | For a valid overlap `0 < l ≤ min(|c1|,|c2|)`: result = `c1 + c2[l..]`, length = `|c1|+|c2|−l` | Collapse keeps one copy of the length-`l` overlap [1] |
| INV-02 | `l = 0` ⇒ result = `c1 + c2` (plain concatenation) | Source defines the no-overlap case as concatenation [1] |
| INV-03 | `l ≤ 0` or `l > min(|c1|,|c2|)` ⇒ result = `c1 + c2` | Such `l` is not a valid suffix/prefix overlap; bound is `min(|x|,|y|)` [1] |
| INV-04 | Result always starts with `c1` and ends with the non-overlapped tail of `c2` | Concatenation never deletes content from `c1` or from `c2[l..]` [1] |

### 2.5 Comparison with Related Methods

| Aspect | Contig Merging (this unit) | Full OLC assembly |
|--------|----------------------------|-------------------|
| Scope | Collapse one ordered pair at a known overlap | Build overlap graph, layout, consensus over many reads [2] |
| Overlap discovery | Caller-supplied length (via `FindOverlap`) | Computed for all read pairs |
| Determinism | Exact, O(n) | Greedy heuristic (NP-hard SCS/Hamiltonian path) [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `contig1` | `string` | required | Left contig; its suffix overlaps the prefix of `contig2`. | Non-null; any alphabet; case preserved verbatim |
| `contig2` | `string` | required | Right contig; its length-`overlapLength` prefix is the shared region. | Non-null |
| `overlapLength` | `int` | required | Length of the suffix(contig1)/prefix(contig2) overlap to collapse. | Treated as valid only when `0 < l ≤ min(|c1|,|c2|)` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `string` | The merged superstring `c1 + c2[l..]` for a valid `l`, otherwise `c1 + c2`. |

### 3.3 Preconditions and Validation

Null `contig1` or `contig2` throws `ArgumentNullException`. The method is case-sensitive and
alphabet-agnostic (it manipulates characters, not nucleotides) — no T↔U normalization is performed.
A non-positive `overlapLength`, or one exceeding `min(|c1|,|c2|)`, is treated as "no usable overlap"
and yields plain concatenation rather than an error. The supplied overlap length is trusted and not
re-verified against the actual characters; verifying the overlap is the role of `FindOverlap`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate that `contig1` and `contig2` are non-null.
2. If `overlapLength ≤ 0` or `overlapLength > min(|contig1|, |contig2|)`, return `contig1 + contig2`.
3. Otherwise return `contig1 + contig2.Substring(overlapLength)`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `MergeContigs` | O(n) | O(n) | n = `|c1| + |c2|`; a single substring + concatenation, matching the checklist O(n). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.MergeContigs(string contig1, string contig2, int overlapLength)`: collapses the
  two contigs at the given overlap, or concatenates when the overlap is not usable.

### 5.2 Current Behavior

The overlap region is removed from the *front of `contig2`* (the prefix), not the back of `contig1`,
so `contig1` is preserved verbatim as a prefix of the result. The single named constant `NoOverlap`
(= 0) documents the no-overlap sentinel. No suffix tree is used here (see §5.3): this is a single
fixed-length substring/concatenation on a known overlap, not an occurrence-search, so the repository
`SuffixTree` would add construction cost without changing the O(n) result; the suffix tree is the
right tool for *overlap discovery* (`FindOverlap`/`FindAllOverlaps`), not for this collapse step [2].

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Merge as `X + Y[l:]`, keeping one copy of the length-`l` overlap, length `|X|+|Y|−l` [1].
- No-overlap (l = 0) and out-of-range-overlap cases collapse to plain concatenation `X + Y` [1].
- Overlap bounded by `min(|X|,|Y|)` [1].

**Intentionally simplified:**

- (none) — the collapse arithmetic is exact for any supplied overlap length.

**Not implemented:**

- Overlap *verification*: the method does not re-check that `contig1`'s suffix equals `contig2`'s
  prefix; users should rely on `FindOverlap` / `FindAllOverlaps` to compute a verified overlap [2].
- Multi-contig layout/consensus over an overlap graph; **users should rely on**
  `SequenceAssembler.AssembleOLC` / `AssembleDeBruijn`.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Caller-supplied overlap length is trusted | Assumption | A wrong `l` collapses a region that may not match | accepted | API-contract boundary; `FindOverlap` computes/verifies `l` [2] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `overlapLength = 0` | `c1 + c2` | No overlap ⇒ concatenate [1] |
| `overlapLength < 0` | `c1 + c2` | Non-positive is not a valid overlap [1] |
| `overlapLength > min(|c1|,|c2|)` | `c1 + c2` | Overlap bounded by the shorter string [1] |
| `overlapLength = min(|c1|,|c2|)` | full collapse of the shorter prefix | Boundary of a valid overlap [1] |
| empty `c1` or `c2`, `l = 0` | the other contig | Concatenation with the empty string |
| null `c1` or `c2` | `ArgumentNullException` | Input validation (sibling-method convention) |

### 6.2 Limitations

This is a single-pair primitive: it neither discovers overlaps nor resolves repeats, mismatches, or
indels in the overlap (it assumes an exact, already-known overlap). For approximate overlaps or
full assembly use the OLC / de Bruijn entry points.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// BAA + AAB, suffix "AA" of BAA == prefix "AA" of AAB, overlap length 2.
string merged = SequenceAssembler.MergeContigs("BAA", "AAB", 2); // "BAAB"
```

**Numerical / biological walk-through:**

From the Langmead greedy-SCS trace [1], chaining {AAA, AAB, ABB, BBB, BBA} at overlap 2 each:
`AAA`+`AAB`→`AAAB`, +`ABB`→`AAABB`, +`BBB`→`AAABBB`, +`BBA`→`AAABBBA` (the SCS, length 7).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_MergeContigs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_MergeContigs_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ASSEMBLY-MERGE-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-MERGE-001-Evidence.md)
- Related algorithms: [Overlap_Layout_Consensus](../Assembly/Overlap_Layout_Consensus.md)

## 8. References

1. Langmead, B. (n.d.). *Assembly & shortest common superstring* (JHU computational genomics lecture
   notes). https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_scs.pdf (accessed 2026-06-13).
2. Langmead, B. (n.d.). *Overlap Layout Consensus assembly* (JHU computational genomics lecture
   notes). https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf (accessed 2026-06-13).
3. MIT 7.91J / Burge C., Fraenkel E. (2014). *Foundations of Computational and Systems Biology,
   Lecture 6: Genome Assembly*. MIT OpenCourseWare.
   https://ocw.mit.edu/courses/7-91j-foundations-of-computational-and-systems-biology-spring-2014/e885f0eb376ea6c2045eb9d8847f106f_MIT7_91JS14_Lecture6.pdf (accessed 2026-06-13).
4. Compeau, P.E.C., Pevzner, P.A., Tesler, G. (2011). *How to apply de Bruijn graphs to genome
   assembly*. Nature Biotechnology 29:987–991. https://doi.org/10.1038/nbt.2023
