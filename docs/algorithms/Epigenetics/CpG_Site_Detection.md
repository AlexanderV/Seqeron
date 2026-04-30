# CpG Site Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Epigenetics |
| Test Unit ID | EPIGEN-CPG-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

CpG site detection identifies every position in a DNA sequence where a cytosine is immediately followed by a guanine in the 5'→3' direction [1][4]. The repository also computes the canonical CpG observed/expected ratio and detects CpG islands using sequence length, GC content, and CpG-density criteria [1][2]. Site scanning is exact and deterministic, while island detection is a window-based classification procedure whose reported boundaries depend on the configured window length and thresholds [1]. The implementation exposes separate entry points for site enumeration, ratio calculation, and island detection in the same analyzer class.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In single-stranded 5'→3' notation, a CpG site is a cytosine followed by guanine and must not be confused with GpC [1][4]. CpG islands are regions with elevated CpG density that are commonly formalized by minimum length, GC-content, and observed/expected CpG-ratio criteria [1][2].

### 2.2 Core Model

CpG site scanning evaluates each adjacent dinucleotide window and reports the 0-based position of the cytosine when the window is `CG` [4].

The CpG observed/expected ratio is defined as [1]:

$$
\mathrm{O/E} = \frac{\text{CpG count}}{\left(\frac{C_{\text{count}} \times G_{\text{count}}}{L}\right)}
$$

where CpG count is the number of `CG` dinucleotides, $C_{\text{count}}$ is the number of cytosines, $G_{\text{count}}$ is the number of guanines, and $L$ is the sequence length [1].

The cited default CpG-island criteria are [1]:

1. Length at least 200 bp.
2. GC content at least 50%.
3. CpG observed/expected ratio at least 0.6.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported CpG site marks a `C` whose next nucleotide is `G`. | This is the defining predicate of a CpG dinucleotide [1][4]. |
| INV-02 | The observed/expected ratio is uniquely determined by CpG count, C count, G count, and sequence length. | The formula is fixed by the Gardiner-Garden and Frommer definition [1]. |
| INV-03 | Site detection is deterministic for a given input sequence. | The algorithm is a literal scan over adjacent dinucleotides with no heuristic branching. |
| INV-04 | Any reported CpG island satisfies the configured minimum length, GC-content threshold, and CpG-ratio threshold. | Island classification is defined by those three criteria [1][2]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | DNA sequence to scan. | `null` and empty input are accepted; methods normalize to uppercase before scanning. |
| `minLength` | `int` | `200` | Minimum island length for `FindCpGIslands`. | Effective only for island detection; sequences shorter than this yield no islands. |
| `minGc` | `double` | `0.5` | Minimum GC fraction for `FindCpGIslands`. | Compared against the GC fraction computed for each candidate island. |
| `minCpGRatio` | `double` | `0.6` | Minimum observed/expected CpG ratio for `FindCpGIslands`. | Compared against the CpG O/E value computed for each candidate island. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| CpG positions | `IEnumerable<int>` | 0-based positions of the cytosine in each detected `CG` dinucleotide. |
| Observed/expected ratio | `double` | CpG O/E ratio for the full input sequence; `0.0` when the sequence is null, empty, shorter than 2, or has zero expected CpG count. |
| CpG islands | `IEnumerable<(int Start, int End, double GcContent, double CpGRatio)>` | Candidate islands returned by `FindCpGIslands`, where `Start` is 0-based inclusive and `End` is the exclusive end position used by the implementation. |

### 3.3 Preconditions and Validation

`FindCpGSites` and `FindCpGIslands` yield no results for `null` or empty input. `CalculateCpGObservedExpected` returns `0.0` for `null`, empty, or length-1 input, and also returns `0.0` when the expected CpG count is zero. The implementation uppercases input before scanning but does not enforce a DNA alphabet; characters other than `C` and `G` simply do not contribute to CpG counts unless they appear in a literal `CG` pair after normalization. Positions are 0-based, and the island tuple uses an exclusive `End` coordinate because the implementation stores `i + windowSize` as the region boundary before slicing the final substring.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the sequence to uppercase.
2. For CpG site detection, scan adjacent nucleotide pairs and yield the current index whenever the pair is `CG`.
3. For the observed/expected ratio, count cytosines, guanines, and CpG dinucleotides, then apply the Gardiner-Garden and Frommer formula.
4. For CpG islands, slide a window of length `minLength` one base at a time, mark contiguous windows whose GC fraction and CpG ratio satisfy the configured thresholds, then recompute GC content and CpG ratio on each merged candidate island before yielding it.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindCpGSites` | `O(n)` | `O(1)` plus yielded output | Single pass over adjacent dinucleotides. |
| `CalculateCpGObservedExpected` | `O(n)` | `O(1)` | Counts `C`, `G`, and `CG` across the sequence. |
| `FindCpGIslands` | `O(n × w)` | `O(w)` per active window plus yielded output | `w` is the effective window length because each window is rescanned to compute GC content and CpG ratio. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [EpigeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs)

- `EpigeneticsAnalyzer.FindCpGSites(string)`: Scans the normalized sequence and yields 0-based CpG positions.
- `EpigeneticsAnalyzer.CalculateCpGObservedExpected(string)`: Counts `C`, `G`, and `CG` and returns the CpG O/E ratio.
- `EpigeneticsAnalyzer.FindCpGIslands(string, int, double, double)`: Performs sliding-window island detection and returns merged island tuples.

### 5.2 Current Behavior

`FindCpGSites` and `CalculateCpGObservedExpected` both uppercase the input string before scanning. `FindCpGIslands` also uppercases the input, advances the window by one nucleotide at a time, merges consecutive qualifying windows into one candidate island, and recomputes GC content and CpG ratio on the merged substring before yielding it. In code, island thresholds are compared inclusively against the configured parameters (`gc >= minGc` and `cpgRatio >= minCpGRatio`). The returned island tuple uses an exclusive `End` coordinate because the stored boundary is `i + windowSize`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Exact CpG-site scanning as adjacent `C` then `G` in the 5'→3' direction [1][4].
- CpG observed/expected ratio calculation using $\text{CpG count} / ((C_{\text{count}} \times G_{\text{count}}) / L)$ [1].
- Default island parameters of length `200`, GC fraction `0.5`, and CpG ratio `0.6`, matching the cited Gardiner-Garden and Frommer thresholds numerically [1].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Dedicated presets for alternative CpG-island definitions such as the stricter Takai and Jones criteria; **users should rely on:** `FindCpGIslands(sequence, minLength, minGc, minCpGRatio)` with custom threshold arguments.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` or empty input | `FindCpGSites` and `FindCpGIslands` yield no results; `CalculateCpGObservedExpected` returns `0.0`. | Guard behavior implemented directly in the entry points. |
| Sequence length `< 2` | No CpG sites and ratio `0.0`. | A dinucleotide cannot be formed from fewer than two bases. |
| Sequence shorter than `minLength` | `FindCpGIslands` yields no islands. | The algorithm does not evaluate windows smaller than the configured minimum length. |
| Sequence with no `C` or no `G` | Observed/expected ratio is `0.0`. | Expected CpG count becomes zero, so the implementation returns `0.0`. |
| Adjacent CpGs such as `CGCG` | Both CpG positions are reported. | The scan evaluates every adjacent dinucleotide window independently. |
| `GpC` without `CpG` | No CpG is reported for the `GC` window. | The definition requires `C` followed by `G`, not the reverse order [1][4]. |

### 6.2 Limitations

The island detector is a sequence-only classifier: it does not infer methylation state, chromatin context, or promoter status from the returned regions. The implementation is also not optimized with rolling counts, so island detection rescans each window and has `O(n × w)` complexity rather than a fully incremental `O(n)` implementation. Alternative CpG-island criteria from later literature must be supplied explicitly through the method parameters rather than selected through named presets.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [EpigeneticsAnalyzer_CpGDetection_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_CpGDetection_Tests.cs)
- Test spec: [EPIGEN-CPG-001.md](../../../tests/TestSpecs/EPIGEN-CPG-001.md)
- Evidence: [EPIGEN-CPG-001-Evidence.md](../../../docs/Evidence/EPIGEN-CPG-001-Evidence.md)

## 8. References

1. Gardiner-Garden M, Frommer M. 1987. CpG islands in vertebrate genomes. Journal of Molecular Biology. https://doi.org/10.1016/0022-2836(87)90689-9
2. Takai D, Jones PA. 2002. Comprehensive analysis of CpG islands in human chromosomes 21 and 22. Proceedings of the National Academy of Sciences of the United States of America. https://doi.org/10.1073/pnas.052410099
3. Saxonov S, Berg P, Brutlag DL. 2006. A genome-wide analysis of CpG dinucleotides in the human genome. Proceedings of the National Academy of Sciences of the United States of America. https://doi.org/10.1073/pnas.0510310103
4. Wikipedia. 2026. CpG site. Wikipedia. https://en.wikipedia.org/wiki/CpG_site
