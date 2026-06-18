# Guide RNA Design

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | CRISPR-GUIDE-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-18 |

## 1. Overview

Guide RNA (gRNA/sgRNA) design in this repository evaluates and generates CRISPR guide candidates using sequence-composition metrics that are intended to correlate with editing efficiency and specificity. The design surface combines PAM-based candidate extraction with a score-based evaluation over GC content, poly-T content, seed-region GC, self-complementarity, and common restriction-site presence. The resulting score is a heuristic quality measure rather than a learned or experimentally calibrated predictor.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An sgRNA consists of a spacer sequence that determines targeting specificity plus a scaffold sequence required for Cas binding. The original document cites Addgene for the canonical spacer-plus-scaffold structure and for the requirement that the target lie adjacent to a valid PAM. The same document also cites GC content, seed-region composition, and poly-T avoidance as practical guide-quality factors. Sources: Addgene CRISPR Guide, Wikipedia (Guide RNA, CRISPR gene editing, Protospacer adjacent motif).

### 2.2 Core Model

Guide candidates are extracted from PAM sites in a target region and scored from a base score of `100` with the following deductions preserved from the original document and confirmed in source:

| Rule | Deduction |
|------|-----------|
| GC below `MinGcContent` | `(MinGcContent - actual) × 2` |
| GC above `MaxGcContent` | `(actual - MaxGcContent) × 2` |
| Contains `TTTT` | `20` |
| Self-complementarity above `0.3` | `selfComplementarity × 30` |
| Seed-region GC outside `30-80%` | `5` |
| Common restriction site detected | `5` |

The seed region is defined as the last 10 nucleotides for PAM-after-target systems and the first 10 nucleotides for PAM-before-target systems. The final score is clamped to `>= 0`.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Guide scores are bounded below by `0` | `EvaluateGuideRna(...)` applies `Math.Max(0, score)` |
| INV-02 | `DesignGuideRnas(...)` yields only candidates with `Score >= MinScore` | The source filters candidates against `effectiveParams.MinScore` |
| INV-03 | The seed-region GC uses 10 nt at the PAM-proximal end chosen by system orientation | The source uses the last 10 bases for PAM-after-target systems and the first 10 for PAM-before-target systems |
| INV-04 | `FullGuideRna` appends a fixed scaffold to `Sequence` | `GuideRnaCandidate.FullGuideRna` is a derived property in source |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[DesignGuideRnas] sequence` | `DnaSequence` | required | Target DNA sequence containing the region of interest | Null input throws `ArgumentNullException` |
| `[DesignGuideRnas] regionStart` | `int` | required | Start of the target region | Must be `>= 0` and `< sequence.Length` |
| `[DesignGuideRnas] regionEnd` | `int` | required | End of the target region | Must be `>= regionStart` and `< sequence.Length` |
| `[DesignGuideRnas/EvaluateGuideRna] systemType` | `CrisprSystemType` | `SpCas9` | CRISPR system definition used for guide length and seed orientation | Must resolve through `CrisprDesigner.GetSystem(...)` |
| `[DesignGuideRnas/EvaluateGuideRna] parameters` | `GuideRnaParameters?` | `GuideRnaParameters.Default` | Optional design thresholds | Defaults are `MinGcContent = 40`, `MaxGcContent = 70`, `MinScore = 50`, `AvoidPolyT = true`, `CheckSelfComplementarity = true`; the two boolean flags are stored in the current parameter record but not honored by the scoring path |
| `[EvaluateGuideRna] guideSequence` | `string` | required | Guide sequence to score | Null or empty throws `ArgumentNullException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Sequence` | `string` | Guide sequence being evaluated |
| `Position` | `int` | `pamSite.TargetStart` for designed guides and `-1` for standalone evaluation; reverse-strand designed guides currently keep that reverse-complement coordinate rather than a forward-source remap |
| `IsForwardStrand` | `bool` | Strand orientation for designed candidates |
| `GcContent` | `double` | Overall GC percentage |
| `SeedGcContent` | `double` | GC percentage of the seed region |
| `HasPolyT` | `bool` | Whether the guide contains `TTTT` |
| `SelfComplementarityScore` | `double` | Internal self-complementarity measure |
| `Score` | `double` | Composite guide-quality score |
| `Issues` | `IReadOnlyList<string>` | Human-readable issues recorded during evaluation |
| `System` | `CrisprSystem` | CRISPR system metadata returned by the current evaluation path; designed guides for systems not explicitly remapped in the helper currently reuse `SpCas9` metadata |
| `FullGuideRna` | `string` | Guide plus the fixed scaffold sequence |

### 3.3 Preconditions and Validation

`DesignGuideRnas(...)` throws `ArgumentNullException` for a null sequence and `ArgumentOutOfRangeException` when the requested region falls outside the sequence. `EvaluateGuideRna(...)` throws `ArgumentNullException` for null or empty guide strings, uppercases the input, and scores it as given without enforcing a specific guide length. In contrast, `DesignGuideRnas(...)` derives the guide length from the selected CRISPR system's PAM definition and extracted target sequence.

## 4. Algorithm

### 4.1 High-Level Steps

1. Resolve the CRISPR system into its PAM pattern, guide length, and orientation.
2. Find PAM sites in the requested region.
3. Extract each candidate guide sequence from the PAM-adjacent target interval.
4. Compute GC content, seed-region GC, poly-T presence, self-complementarity, and restriction-site presence.
5. Apply the heuristic score deductions from a base score of `100`.
6. Yield only candidates whose score meets the configured minimum threshold.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Guide-design defaults confirmed in source:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MinGcContent` | `40` | Minimum acceptable GC percentage |
| `MaxGcContent` | `70` | Maximum acceptable GC percentage |
| `MinScore` | `50` | Minimum score for `DesignGuideRnas(...)` output |
| `AvoidPolyT` | `true` | Stored in the current parameter record, but the `TTTT` penalty is still applied even when the flag is `false` |
| `CheckSelfComplementarity` | `true` | Stored in the current parameter record, but self-complementarity penalties are still applied even when the flag is `false` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DesignGuideRnas` | `O(n)` | `O(k)` | The original document describes linear-time design over the target region |
| `EvaluateGuideRna` | `O(k)` | `O(1)` | The current source performs linear scans over guide length `k` for GC, seed, and self-complementarity checks |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CrisprDesigner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs)

- `CrisprDesigner.DesignGuideRnas(DnaSequence, int, int, CrisprSystemType, GuideRnaParameters?)`: Generates and filters guide candidates from PAM sites in a target region.
- `CrisprDesigner.EvaluateGuideRna(string, CrisprSystemType, GuideRnaParameters?)`: Scores a guide sequence with composition-based heuristics.
- `GuideRnaParameters.Default`: Provides the default GC and score thresholds.

### 5.2 Current Behavior

The current implementation matches the documented `40-70%` GC window, penalizes `TTTT` poly-T sequences, uses a 10-base seed region, and reports issues when the seed GC falls outside `30-80%`. `DesignGuideRnas(...)` first filters PAM sites by a cut-site-in-region check and then evaluates the extracted guide. Standalone evaluation accepts non-standard guide lengths and scores them naturally from their composition, while designed guides follow the guide length of the selected CRISPR system for PAM-adjacent extraction. The internal `EvaluateGuideRna(PamSite, ...)` helper remaps only `SpCas9`, `SaCas9`, and `Cas12a/Cpf1` by name, so designed guides for systems such as `SpCas9NAG`, `AsCas12a`, `LbCas12a`, and `CasX` currently reuse `SpCas9` scoring metadata in the returned candidate. The reported `Position` is copied from `pamSite.TargetStart`, which is not remapped back to the original forward-sequence coordinate for reverse-strand designs. The `AvoidPolyT` and `CheckSelfComplementarity` parameter flags are currently passive record fields: the scoring path still applies poly-T and self-complementarity penalties even when those booleans are set to `false`. The derived `FullGuideRna` appends a fixed scaffold sequence to the spacer.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- PAM-adjacent guide extraction for supported CRISPR systems.
- GC-range screening, poly-T detection, seed-region GC evaluation, and self-complementarity checks.
- A fixed scaffold appended to the designed spacer for `FullGuideRna`.

**Intentionally simplified:**

- The `EvaluateGuideRna` quality score is composition-based rather than position-weighted or experimentally fitted; **consequence:** it ranks candidates heuristically and does not itself reproduce learned activity models — for experimentally-calibrated on-target efficacy call the dedicated `CalculateOnTargetDoench2014` (Rule Set 1) or `CalculateOnTargetRuleSet2` (Rule Set 2 / Azimuth) over the 30-nt context.
- The seed region is fixed to 10 nt using the upper bound of the cited `8-10` range; **consequence:** guides are evaluated with one fixed PAM-proximal window rather than a variable seed definition.
- Designed-candidate metadata is only partially system-specific in the current helper path; **consequence:** some named CRISPR systems share `SpCas9` metadata in returned candidates even though guide extraction still follows the selected PAM geometry.
- `GuideRnaParameters.AvoidPolyT` and `GuideRnaParameters.CheckSelfComplementarity` are not honored by the current scoring path; **consequence:** callers can store those booleans in the parameter record, but they do not disable the corresponding penalties.
- Restriction-site screening is limited to the common-site check embedded in source; **consequence:** cloning constraints beyond that built-in list are not modeled.

**Not implemented:**

- Position-weighted scoring *within the `EvaluateGuideRna` heuristic* and chromatin-accessibility modeling; **users should rely on:** the dedicated `CalculateOnTargetDoench2014` / `CalculateOnTargetRuleSet2` methods for learned on-target efficacy prediction (both now implemented), and downstream experimental validation or external tools for chromatin accessibility.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty guide sequence | Throws `ArgumentNullException` | Explicit source guard in `EvaluateGuideRna(...)` |
| Null target sequence | Throws `ArgumentNullException` | Explicit source guard in `DesignGuideRnas(...)` |
| `regionStart < 0` | Throws `ArgumentOutOfRangeException` | Explicit source guard |
| `regionEnd` outside sequence bounds | Throws `ArgumentOutOfRangeException` | Explicit source guard |
| Guide length not equal to 20 bp in standalone evaluation | Accepted and scored naturally | The standalone evaluator does not enforce system length |
| All-A guide | Valid but low-scoring | GC penalties and issue reporting apply |
| All-G/C guide | Valid but low-scoring | High-GC penalties and issue reporting apply |
| Guide containing `TTTT` | Penalized and flagged | Poly-T detection is part of the score |

### 6.2 Limitations

The current guide-design surface does not implement position-weighted efficacy models, machine learning, chromatin context, or genome-wide uniqueness checks. It is a sequence-level heuristic evaluator whose main purpose is to rank or filter candidates rather than to predict editing outcomes directly.

## 7. Examples and Related Material

### 7.2 Applications and Use Cases (Optional)

Related material called out in the original document:

- `CRISPR-PAM-001`: PAM Site Detection (prerequisite).
- `CRISPR-OFF-001`: Off-Target Analysis (related).

## 8. References

1. Addgene CRISPR Guide. https://www.addgene.org/guides/crispr/
2. Wikipedia: Guide RNA. https://en.wikipedia.org/wiki/Guide_RNA
3. Wikipedia: CRISPR gene editing. https://en.wikipedia.org/wiki/CRISPR_gene_editing
4. Wikipedia: Protospacer adjacent motif. https://en.wikipedia.org/wiki/Protospacer_adjacent_motif
5. Doench et al. (2014). "Rational design of highly active sgRNAs". Nature Biotechnology.
6. Hsu et al. (2013). "DNA targeting specificity of RNA-guided Cas9 nucleases". Nature Biotechnology.
