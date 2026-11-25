# Gene Structure Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | Splicing |
| Test Unit ID | SPLICE-PREDICT-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Gene structure prediction in this repository infers exon and intron organization from splice-site candidates. `SpliceSitePredictor.PredictGeneStructure` orchestrates intron prediction, greedy non-overlapping intron selection, exon derivation, and spliced-sequence generation. The supporting public method `PredictIntrons` pairs donor and acceptor sites, optionally adds a branch-point contribution, classifies intron types, and exposes the intron candidates that feed the higher-level structure result. The implementation is deterministic and intentionally heuristic rather than a full probabilistic gene finder.[1][2][3][4]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Gene structure prediction partitions a pre-mRNA-like sequence into exons retained in the mature transcript and introns removed by splicing.[1][2] The original document records these exon and intron concepts:

| Concept | Repository description |
|---------|------------------------|
| Exon types | `Initial`, `Internal`, `Terminal`, `Single` |
| Intron types | `U2`, `U12`, `GcAg`, `Unknown` |
| Exon phase | Cumulative length of previous exons modulo `3` |

### 2.2 Core Model

The repository model is:

1. Predict donor and acceptor sites.
2. Pair donors and acceptors whose distance satisfies intron-length constraints.
3. Add the Gao et al. (2008) `yUnAy` branch-point term (canonical `FindAcceptorBranchPoint`; 0 when none is found).
4. Keep introns whose combined score meets the requested threshold.
5. Select a non-overlapping subset greedily by descending intron score.
6. Reject an intron that would leave an exon (between introns or at either end) shorter than `minExonLength`; derive exons from the gaps between selected introns, so exons + introns tile the sequence.
7. Build the spliced sequence by removing intron intervals from the normalized sequence.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Greedy highest-score non-overlap selection is acceptable for the use case. | The selected intron set can differ from the globally optimal combination under other objectives. |
| ASM-02 | A simple average of donor, acceptor, and optional branch-point scores is sufficient as an intron score. | Score ordering can differ from more formal probabilistic or HMM-based models. |
| ASM-03 | Exon phase can be measured from cumulative exon length starting at the beginning of the provided sequence. | Sequences containing untranslated regions or partial CDS context can produce phases that are locally consistent but biologically incomplete. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Selected introns do not overlap. | `SelectNonOverlappingIntrons` rejects any intron using a previously occupied position. |
| INV-02 | `OverallScore` is the mean of selected intron scores, or `0` when no introns are selected. | `PredictGeneStructure` computes it directly from `selectedIntrons`. |
| INV-03 | The first exon phase is `0` and later phases are cumulative exon length modulo `3`. | `CalculatePhase` sums lengths of previously added exons. |
| INV-04 | A sequence with no selected introns returns a single exon covering the full sequence. | `DeriveExons` explicitly creates a `Single` exon in that case. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[PredictGeneStructure] sequence` | `string` | required | DNA or RNA sequence to analyze. | Empty or `null` returns an empty `GeneStructure`. |
| `[PredictGeneStructure] minExonLength` | `int` | `30` | Minimum exon length retained in the `Exons` list. | Applied when deriving exons from selected introns. |
| `[PredictGeneStructure] minIntronLength` | `int` | `60` | Minimum intron length passed to intron prediction. | Used by `PredictIntrons`. |
| `[PredictGeneStructure] minScore` | `double` | `0.5` | Minimum combined intron score required for retained introns. | Also drives site thresholds indirectly. |
| `[PredictIntrons] maxIntronLength` | `int` | `100000` | Maximum intron length for donor-acceptor pairing. | Public on `PredictIntrons`; fixed to `100000` inside `PredictGeneStructure`. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `GeneStructure.Exons` | `IReadOnlyList<Exon>` | Exons derived from the selected non-overlapping introns. |
| `GeneStructure.Introns` | `IReadOnlyList<Intron>` | Selected introns sorted by start position. |
| `GeneStructure.SplicedSequence` | `string` | Sequence produced by removing selected intron intervals from the normalized input. |
| `GeneStructure.OverallScore` | `double` | Mean selected-intron score, or `0` if no introns are selected. |
| `Intron.Start` / `Intron.End` | `int` | Inclusive intron boundaries defined by donor and acceptor positions. |
| `Intron.Type` | `IntronType` | Intron classification (`U2`, `U12`, `GcAg`, or `Unknown`). |
| `Exon.Type` | `ExonType` | Exon classification (`Initial`, `Internal`, `Terminal`, or `Single`). |
| `Exon.Phase` | `int?` | Phase computed from cumulative preceding exon length. |

### 3.3 Preconditions and Validation

`PredictGeneStructure` and `PredictIntrons` normalize the sequence to uppercase RNA notation. `PredictIntrons` calls `FindDonorSites` and `FindAcceptorSites` with individual site thresholds of `minScore * 0.8` and enables non-canonical sites for both. Branch points are searched in the interval `[acceptor.Position - 50, acceptor.Position - 18]` with a minimum branch-point score of `0.4`. When no branch point is found, the combined intron score is the average of donor and acceptor scores; otherwise it is the mean of donor, acceptor, and branch scores.[6]

## 4. Algorithm

### 4.1 High-Level Steps

1. Return an empty `GeneStructure` for `null` or empty input.
2. Normalize the sequence to uppercase RNA notation.
3. Predict donor and acceptor sites and pair them into intron candidates that satisfy the length limits.
4. For each donor-acceptor pair, find branch points upstream of the acceptor and compute a combined intron score.
5. Keep intron candidates whose combined score is at least `minScore`.
6. Sort candidate introns by descending score and greedily select a non-overlapping subset.
7. Derive exons from the gaps between selected introns, subject to `minExonLength`.
8. Remove selected intron intervals from the sequence to build the spliced sequence.
9. Return exons, selected introns, spliced sequence, and the mean selected-intron score.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The current intron-pairing and classification rules are:

| Rule | Repository behavior |
|------|---------------------|
| Intron length | `acceptor.Position - donor.Position + 1` |
| Candidate acceptance | Length in range and combined score `>= minScore` |
| U2 / GC-AG score | `(donor + acceptor + branch) / 3`, branch term 0 when no branch point (fixed term set: a branch point can only raise the score) |
| U12 score | `(U12 donor + U12 branch-site acceptor) / 2` |
| U12 intron type | Only when BOTH the U12 donor (/RUAUCCUY) and the U12 acceptor (UUCCUUAAC branch site 8-21 nt upstream) score >= max(`minScore`, 0.75) with AU-AC or GU-AG termini; takes precedence over U2 for the same coordinates (overload `PredictIntrons(..., u12MinSiteScore)`) |
| `GcAg` intron type | `GC` donor + `AG` acceptor |
| `U2` intron type | `GU` donor + `AG` acceptor |
| Coordinates | 0-based inclusive; `Start` = first intron nt (donor G/A), `End` = last intron nt (acceptor G/C) |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictIntrons` | `O(D * A)` | `O(D + A + k)` | Pairs donor and acceptor candidates, where `k` is the number of introns yielded. |
| `PredictGeneStructure` | `O(D * A)` plus selection/exon derivation | `O(n)` | Pairwise intron generation dominates. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SpliceSitePredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs)

- `SpliceSitePredictor.PredictGeneStructure(string, int, int, double)`
- `SpliceSitePredictor.PredictIntrons(string, int, int, double)`
- `SpliceSitePredictor.SelectNonOverlappingIntrons(...)` (private helper)
- `SpliceSitePredictor.DeriveExons(...)` (private helper)
- `SpliceSitePredictor.GenerateSplicedSequence(...)` (private helper)

### 5.2 Current Behavior

`PredictGeneStructure` always uses a fixed `maxIntronLength` of `100000`. Candidate introns are selected greedily by descending score (ties: lower start); an intron is rejected when it overlaps a selected one or would leave an exon shorter than `minExonLength` (>= 1) between introns or at either sequence end. Exons are never dropped, so exons + introns tile the sequence and the spliced sequence is the input with exactly the selected introns removed. Exon phase = cumulative preceding exon length mod 3 (Ensembl start phase).

`DetectAlternativeSplicing` reports SUPPA/rMATS events relative to the selected (reference) introns: `ExonSkipping` (SE), `Alt5SS` (A5SS), `Alt3SS` (A3SS), `MutuallyExclusiveExons` (MXE); `IsWithinCodingRegion` walks upstream in-frame to an AUG without crossing a stop codon.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Gene structure is represented as alternating exonic and intronic regions derived from splice-site boundaries.[1][2]
- The repository distinguishes exon types (`Initial`, `Internal`, `Terminal`, `Single`) and intron classes (`U2`, `U12`, `GcAg`, `Unknown`).[1][5]
- Exon phase is computed from cumulative exon length modulo `3`.[6]

**Intentionally simplified:**

- Intron selection is greedy rather than globally optimized; **consequence:** a high-scoring intron can block an alternative combination with better global structure.
- Intron scoring is a simple arithmetic mean of site scores and optional branch-point score; **consequence:** scores are heuristic rather than probabilistic likelihoods.
- Exon derivation uses only intron gaps and minimum-length filtering; **consequence:** untranslated regions and richer transcript models are not represented.

**Not implemented:**

- HMM or generalized-HMM gene finders such as GenScan/Augustus-style dynamic programming; **users should rely on:** no current alternative.
- Transcript-evidence based intron retention (RI): `FindRetainedIntronCandidates` is a feature heuristic (short, weaker sites; Braunschweig 2014) only.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | The intron score has no published form; a fixed term set (branch term 0 when absent) is used. | Behavior | A branch point can only raise the score. | accepted | B11 F33. |
| 2 | The U12 site-score cut-off 0.75 and the 500 nt / 0.8 retained-intron cut-offs are heuristics. | Behavior | Tunable (`u12MinSiteScore` overload). | accepted | No published cut-off for fraction-of-consensus scores. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or `null` input | Returns empty exons, empty introns, empty spliced sequence, overall score `0`. | Explicit early return. |
| No introns pass threshold | Returns a single exon covering the whole sequence. | `DeriveExons` handles the no-intron case directly. |
| Introns shorter than `minIntronLength` | Not yielded. | Length filter in `PredictIntrons`. |
| Introns longer than `maxIntronLength` | Not yielded. | Length filter in `PredictIntrons`. |
| Higher `minScore` | Produces a subset of introns and potentially fewer exons. | Candidate and selection filtering are threshold-driven. |

### 6.2 Limitations

This implementation is a splice-site-driven heuristic structure predictor. It does not use a full gene-finding model, does not optimize globally over transcript structures, and, being splice-motif only, ranks many false GT-AG pairs: on Ensembl HBB / TP53 / CFTR (exons 2-5) all 15 annotated introns are among the candidates (typed U2) but the greedy selection recovers none of them (see the Evidence file).

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SpliceSitePredictor_GeneStructure_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/SpliceSitePredictor_GeneStructure_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [SPLICE-PREDICT-001.md](../../../tests/TestSpecs/SPLICE-PREDICT-001.md)
- Related algorithms: [Donor_Site_Detection.md](Donor_Site_Detection.md), [Acceptor_Site_Detection.md](Acceptor_Site_Detection.md)

## 8. References

1. Gilbert W. 1978. Why genes in pieces? Nature. doi:10.1038/271501a0
2. Breathnach R, Chambon P. 1981. Organization and expression of eucaryotic split genes. Annual Review of Biochemistry. doi:10.1146/annurev.bi.50.070181.002025
3. Shapiro MB, Senapathy P. 1987. RNA splice junctions of different classes of eukaryotes. Nucleic Acids Research. doi:10.1093/nar/15.17.7155
4. Burge CB, Tuschl T, Sharp PA. 1999. Splicing of precursors to mRNAs by the spliceosomes. The RNA World. N/A
5. Sakharkar MK et al. 2002. ExInt: an Exon Intron Database. Nucleic Acids Research. doi:10.1093/nar/30.1.191
6. Alberts B et al. 2002. Molecular Biology of the Cell. N/A
7. Test specification: [SPLICE-PREDICT-001.md](../../../tests/TestSpecs/SPLICE-PREDICT-001.md)
