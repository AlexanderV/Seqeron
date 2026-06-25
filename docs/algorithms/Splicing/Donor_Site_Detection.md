# Donor (5') Splice Site Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Splicing |
| Test Unit ID | SPLICE-DONOR-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Donor splice-site detection identifies candidate 5' splice sites at the start of introns. In this repository, `SpliceSitePredictor.FindDonorSites` scans a sequence for canonical `GU` donor dinucleotides and, when requested, optional `GC` and U12-style `AU` donors. Candidates are scored with a binary consensus model over the donor context and filtered by `minScore`. The method is deterministic, linear in sequence length, and returns a `SpliceSite` record for each retained candidate.[1][2]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The donor splice site marks the 5' end of an intron. The original document records the extended consensus motif `MAG|GURAGU`, where the canonical intronic start is `GU` and non-canonical `GC` and U12-type `AU` donors are recognized in smaller subsets of introns.[1][2]

| Position | Consensus |
|----------|-----------|
| `-3` | `M` (`A` or `C`) |
| `-2` | `A` |
| `-1` | `G` |
| `0` | `G` |
| `+1` | `U` |
| `+2` | `R` (`A` or `G`) |
| `+3` | `A` |
| `+4` | `G` |
| `+5` | `U` |

### 2.2 Core Model

The repository uses a binary consensus weight table rather than a log-odds model: each position in the donor motif contributes `1.0` for a consensus match and `0.0` otherwise, and the score is the fraction of matched positions among the positions available in the sequence window. U12 donors are scored separately against the fixed consensus `AUAUCC` using the fraction of matching positions.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | A short consensus-based donor motif is sufficient to prioritize candidate donor sites. | Scores can miss broader sequence determinants used by more sophisticated splice-site models. |
| ASM-02 | `GC` and U12-style `AU` donors should be optional rather than always included. | Canonical-only scans can miss these alternative donor classes by design. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Returned `Score` values are in `[0, 1]`. | Both donor-scoring helpers normalize by the number of matched positions. |
| INV-02 | Returned `Confidence` values are in `[0, 1]`. | `CalculateConfidence` clamps the result. |
| INV-03 | Canonical and `GC` donors are returned with `Type = Donor`; U12 `AU` donors use `Type = U12Donor`. | The public method assigns those enum values explicitly. |
| INV-04 | Higher `minScore` produces a subset of lower-threshold results. | Filtering is applied as `if (score >= minScore)`. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | DNA or RNA sequence to scan. | Empty, `null`, or length `< 6` yields no results. |
| `minScore` | `double` | `0.5` | Minimum donor score required to emit a site. | Applied after score normalization. |
| `includeNonCanonical` | `bool` | `false` | Whether to include `GC` donors and U12-style `AU` donors. | Canonical `GU` scanning is always active. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Index of the donor dinucleotide start in the normalized sequence. |
| `Type` | `SpliceSiteType` | `Donor` for canonical `GU` and `GC`, `U12Donor` for U12-style `AU`. |
| `Motif` | `string` | Context extracted by `GetMotifContext(sequence, position, 3, 6)`. |
| `Score` | `double` | Normalized donor score. |
| `Confidence` | `double` | Clamped confidence derived from the score. |

### 3.3 Preconditions and Validation

The method uppercases the input and converts `T` to `U` before scanning. No results are produced for sequences shorter than 6 bases because a donor context cannot be scored meaningfully within the method's windowing rules. Canonical scanning checks for `GU`, optional `GC` scanning reuses the same donor scorer, and optional U12 scanning looks for `AU` at the candidate position.[3]

## 4. Algorithm

### 4.1 High-Level Steps

1. Return no results for `null`, empty, or very short input.
2. Normalize the sequence to uppercase RNA notation.
3. Iterate across candidate dinucleotide starts.
4. Score canonical `GU` donors using the donor consensus table.
5. If `includeNonCanonical` is enabled, also score `GC` donors with the same scorer and `AU` donors with the U12 consensus scorer.
6. Emit each site whose score is at least `minScore`, together with context and confidence.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The donor scoring tables are implemented directly in `SpliceSitePredictor`.

| Candidate class | Recognition rule | Scoring rule |
|-----------------|------------------|--------------|
| Canonical donor | `GU` at positions `i`, `i+1` | Fraction of consensus matches over donor offsets `-3..+5` |
| `GC` donor | `GC` at positions `i`, `i+1` when non-canonical enabled | Same donor scorer as canonical sites |
| U12 donor | `AU` at positions `i`, `i+1` when non-canonical enabled | Fraction of matches to `AUAUCC` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindDonorSites` | `O(n)` | `O(1)` auxiliary | Single pass over the normalized sequence with streaming output. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SpliceSitePredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs)

- `SpliceSitePredictor.FindDonorSites(string, double, bool)`
- `SpliceSitePredictor.ScoreDonorSite(string, int)` (private helper)
- `SpliceSitePredictor.ScoreU12DonorSite(string, int)` (private helper)
- `SpliceSitePredictor.ScoreDonorMaxEnt(string)` (opt-in Yeo & Burge 2004 MaxEntScan `score5ss` maximum-entropy 5' donor score, in bits; embedded probability table)

### 5.2 Current Behavior

Canonical `GU` and optional `GC` donors are both emitted with `Type = Donor`; the distinction between `U2` and `GC-AG` introns is made later during intron classification rather than at donor-site emission time. `ScoreDonorSite` averages binary consensus matches across the 9 donor positions. `ScoreU12DonorSite` scores the six-base `AUAUCC` consensus by match fraction. Confidence is computed with `CalculateConfidence(score, 0.5, 1.0)`, and there is no extra `GC` penalty factor in the current implementation.[3]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Canonical donor detection is anchored on `GU` and the extended donor consensus around the exon-intron junction.[1][2]
- Optional support exists for `GC` donors and U12-style donor recognition.[2]
- The opt-in `ScoreDonorMaxEnt` implements the Yeo & Burge (2004) **MaxEntScan `score5ss`** maximum-entropy 5' donor model on a 9-nt window (3 exon + 6 intron, conserved `GT` at 0-based positions 3–4): `log2(P_maxent/P_background)`, removing the `GT` (scored by a consensus/background model) and looking up the maximum-entropy probability of the 7 remaining positions (`window[0:3] + window[5:9]`) directly in a single embedded table. Unlike score3, score5 is single-matrix (no overlapping sub-windows). The factorisation and table are the MIT-licensed maxentpy port; the canonical documented example `ScoreDonorMaxEnt("cagGTAAGT")` reproduces `10.86` bits exactly (provenance + licence: `Data/maxent_score5.LICENSE.md`).[4]

**Intentionally simplified:**

- The default donor scorer (`FindDonorSites` / `ScoreDonorSite`) uses binary consensus weights instead of a trained log-odds or maximum-entropy model; **consequence:** scores rank motif matches directly but do not capture richer statistical dependencies. (The opt-in `ScoreDonorMaxEnt` provides the true maximum-entropy model for callers who need it; the default scorer is unchanged.)
- Each donor site is scored independently from local motif context only; **consequence:** exon structure and longer-range splicing signals are not considered at this stage.

**Not implemented:**

- Species-trained donor predictors beyond the human/mammalian MaxEntScan score5 model; **users should rely on:** `ScoreDonorMaxEnt` (human-trained) or the MaxEntScan reference tool for other species.
- A distinct `SpliceSiteType` for `GC` donors; **users should rely on:** later intron classification in `DetermineIntronType(...)`.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No `GU` in canonical mode | Returns no donor sites. | No canonical donor candidate is available. |
| Sequence shorter than 6 nt | Returns no donor sites. | Guard clause. |
| DNA input with `T` | Same behavior as RNA input with `U`. | `T` is converted to `U`. |
| Lowercase input | Same behavior as uppercase input. | The sequence is uppercased before scanning. |
| Multiple donor motifs | Each candidate is scored independently. | The scan iterates over every candidate start position. |

### 6.2 Limitations

The default donor-site detection in this repository is a lightweight motif scorer. It does not use statistical training, does not model longer dependencies, and does not decide intron structure on its own. As an opt-in companion it offers the Yeo & Burge (2004) MaxEntScan `score5ss` maximum-entropy 5' donor score (`ScoreDonorMaxEnt`, on a 9-nt window, embedded MIT-licensed table — see §5.3). The score5 model is human/mammalian-trained; other species may have different donor statistics.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SpliceSitePredictor_DonorSite_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [SPLICE-DONOR-001.md](../../../tests/TestSpecs/SPLICE-DONOR-001.md)
- Related algorithms: [Acceptor_Site_Detection.md](Acceptor_Site_Detection.md), [Gene_Structure_Prediction.md](Gene_Structure_Prediction.md)

## 8. References

1. Shapiro MB, Senapathy P. 1987. RNA splice junctions of different classes of eukaryotes. Nucleic Acids Research. doi:10.1093/nar/15.17.7155
2. Burge CB, Tuschl T, Sharp PA. 1999. Splicing precursors to mRNAs by the spliceosomes. The RNA World. N/A
3. Test specification: [SPLICE-DONOR-001.md](../../../tests/TestSpecs/SPLICE-DONOR-001.md)
4. Yeo G, Burge CB. 2004. Maximum entropy modeling of short sequence motifs. Journal of Computational Biology. doi:10.1089/106652704773135290
5. Wikipedia contributors. 2026. RNA splicing. Wikipedia. https://en.wikipedia.org/wiki/RNA_splicing
6. Wikipedia contributors. 2026. Spliceosome. Wikipedia. https://en.wikipedia.org/wiki/Spliceosome
