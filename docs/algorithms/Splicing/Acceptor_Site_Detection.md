# Acceptor Site Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Splicing |
| Test Unit ID | SPLICE-ACCEPTOR-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Acceptor splice-site detection identifies candidate 3' splice sites at the ends of introns. In this repository, `SpliceSitePredictor.FindAcceptorSites` scans a normalized RNA sequence for canonical `AG` acceptors and, when requested, U12-style `AC` acceptors. Canonical scoring combines a polypyrimidine-tract contribution with a sparse position-weight matrix, while U12 acceptors are scored against a `YCCAC`-style consensus plus upstream pyrimidine content. Candidates meeting the configured threshold are returned as `SpliceSite` records.[1][2]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The original document describes the 3' splice-site consensus as `(Y)nNCAG|G`, with an almost invariant terminal `AG` dinucleotide preceded by a polypyrimidine tract (PPT).[1][2]

| Feature | Repository description |
|---------|------------------------|
| Canonical acceptor | `AG` |
| Upstream context | PPT enriched in `C` and `U` |
| Key conserved positions | `-3` favors `C`, `-2` is `A`, `-1` is `G`, position `0` often favors `G` |
| U12-type acceptor | Optional `AC` with `YCCAC`-style upstream context |

### 2.2 Core Model

Canonical acceptor scoring combines two components:

1. PPT contribution: count pyrimidines in positions `[position - 15, position - 3)`, divide by `12`, and scale by `2`.
2. Sparse PWM contribution: add log-like weights from the repository's `AcceptorPwm` at offsets `-15`, `-10`, `-5`, `-4`, `-3`, `-2`, `-1`, and `0` relative to the splice site.

The combined value is normalized with:

```text
(score / (count + 1) + 2) / 4
```

U12 acceptors are scored by matching `C` at `-1` and `-2`, a pyrimidine at `-3`, and an upstream PPT fraction, then normalizing by the maximum score of `3.5`.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | A 15-base upstream PPT window captures enough local acceptor context. | Real acceptors with weaker or differently positioned upstream context can be underscored. |
| ASM-02 | The sparse acceptor PWM is sufficient for candidate prioritization. | Fine-grained sequence preferences outside the represented offsets are ignored. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Returned `Score` values are in `[0, 1]`. | Both canonical and U12 scorers clamp the normalized result. |
| INV-02 | Returned `Confidence` values are in `[0, 1]`. | Confidence is computed through a clamped linear interpolation. |
| INV-03 | Canonical `AG` sites use `Type = Acceptor`; optional `AC` sites use `Type = U12Acceptor`. | The public method assigns those enum values explicitly. |
| INV-04 | Returned `Position` for a canonical acceptor is `i + 1`, the index of the `G` in `AG`. | The method yields `Position: i + 1`. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | DNA or RNA sequence to scan. | Empty, `null`, or length `< 20` yields no results. |
| `minScore` | `double` | `0.5` | Minimum normalized score required to emit a site. | Applied after canonical or U12 scoring. |
| `includeNonCanonical` | `bool` | `false` | Whether to include U12-style `AC` acceptors. | Canonical `AG` scanning is always active. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | For canonical `AG`, the index of the `G` in the dinucleotide; for `AC`, the index of the `C`. |
| `Type` | `SpliceSiteType` | `Acceptor` or `U12Acceptor`. |
| `Motif` | `string` | Local context extracted around the first base of the candidate dinucleotide via `GetMotifContext(sequence, position - 1, 15, 2)`; `Position` itself reports the terminal intronic base. |
| `Score` | `double` | Normalized acceptor score. |
| `Confidence` | `double` | Clamped confidence derived from the score. |

### 3.3 Preconditions and Validation

The method uppercases the sequence and converts `T` to `U`. Canonical scanning begins at index `15` so that the upstream PPT window is available. Canonical `AG` and optional `AC` candidates are scored independently, and sites are emitted only when their normalized score meets `minScore`.[8]

## 4. Algorithm

### 4.1 High-Level Steps

1. Return no results for `null`, empty, or too-short input.
2. Normalize the sequence to uppercase RNA notation.
3. Starting at index `15`, scan for canonical `AG` candidates and, optionally, `AC` candidates.
4. For `AG`, compute PPT contribution and sparse-PWM contribution, then normalize the combined score.
5. For `AC`, compute the U12 `YCCAC`-style score plus upstream PPT fraction, then normalize.
6. Emit a `SpliceSite` when the candidate score is at least `minScore`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The current scoring rules are:

| Candidate class | Recognition rule | Scoring rule |
|-----------------|------------------|--------------|
| Canonical acceptor | `AG` at positions `i`, `i+1` | PPT contribution + sparse `AcceptorPwm`, normalized to `[0, 1]` |
| U12 acceptor | `AC` at positions `i`, `i+1` when non-canonical enabled | `YCCAC` component + PPT fraction, normalized by `3.5` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindAcceptorSites` | `O(n)` | `O(1)` auxiliary | Single scan over the normalized sequence. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SpliceSitePredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs)

- `SpliceSitePredictor.FindAcceptorSites(string, double, bool)`
- `SpliceSitePredictor.ScoreAcceptorSite(string, int)` (private helper)
- `SpliceSitePredictor.ScoreU12AcceptorSite(string, int)` (private helper)

### 5.2 Current Behavior

Canonical acceptor scoring uses a separate PPT component and applies the `AcceptorPwm` with offsets relative to the splice site, implemented as `position + 2 + offset` inside `ScoreAcceptorSite`. `FindAcceptorSites` reports the index of the terminal intronic nucleotide (`G` in `AG`) rather than the first exonic nucleotide. U12 acceptors are returned only when `includeNonCanonical` is enabled and are scored against `YCCAC`-style upstream context plus PPT quality.[8]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Canonical acceptor recognition is based on the `AG` terminus and upstream PPT-rich context of `(Y)nNCAG|G`.[1][2]
- Optional U12-style acceptor handling is based on `AC` and `YCCAC`-style context.[4][5][6]

**Intentionally simplified:**

- The canonical scorer uses a sparse eight-position PWM and a separate PPT term rather than a full statistical acceptor model; **consequence:** scoring is interpretable and fast, but less expressive than richer trained predictors.
- U12 acceptor scoring uses a short consensus and PPT fraction only; **consequence:** minor-spliceosome candidates are heuristic approximations.

**Not implemented:**

- Full-length acceptor models such as maximum-entropy or species-trained acceptor predictors; **users should rely on:** no current alternative.
- Automatic disambiguation between cryptic and functional `AG` sites beyond local scoring; **users should rely on:** caller-side thresholding and downstream intron pairing.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | The acceptor PWM is intentionally sparse, covering eight positions rather than a full continuous window. | Design choice | Reduces model complexity while preserving the strongest documented positions. | accepted | Confirmed in [SpliceSitePredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or `null` input | Returns no acceptor sites. | Guard clause. |
| Sequence shorter than 20 nt | Returns no acceptor sites. | The method requires upstream context. |
| DNA input | Same behavior as RNA input after normalization. | `T` is converted to `U`. |
| Lowercase input | Same behavior as uppercase input. | The sequence is uppercased before scanning. |
| Higher `minScore` | Returns a subset of lower-threshold results. | Filtering uses `score >= minScore`. |

### 6.2 Limitations

Acceptor detection in this repository is a local motif-and-PPT scorer. It does not use a full statistical acceptor model and does not decide exon or intron structure by itself.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SpliceSitePredictor_AcceptorSite_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_AcceptorSite_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [SPLICE-ACCEPTOR-001.md](../../../tests/TestSpecs/SPLICE-ACCEPTOR-001.md)
- Related algorithms: [Donor_Site_Detection.md](Donor_Site_Detection.md), [Gene_Structure_Prediction.md](Gene_Structure_Prediction.md)

## 8. References

1. Shapiro MB, Senapathy P. 1987. RNA splice junctions of different classes of eukaryotes. Nucleic Acids Research. N/A
2. Burge CB, Tuschl T, Sharp PA. 1999. The RNA World. N/A
3. Yeo G, Burge CB. 2004. Journal of Computational Biology. N/A
4. Patel AA, Steitz JA. 2003. Nature Reviews Molecular Cell Biology. N/A
5. Hall SL, Padgett RA. 1994. Journal of Molecular Biology. N/A
6. Jackson IJ. 1991. Nucleic Acids Research. N/A
7. Dietrich RC, Incorvaia R, Padgett RA. 1997. Molecular Cell. N/A
8. Test specification: [SPLICE-ACCEPTOR-001.md](../../../tests/TestSpecs/SPLICE-ACCEPTOR-001.md)
