# Guide RNA Design

| Field | Value |
|-------|-------|
| Algorithm Group | MolTools |
| Test Unit ID | CRISPR-GUIDE-001 |
| Related Projects | N/A |
| Implementation Status | Implemented: composite triage heuristic (declared) + published Rule Set 2 on-target model, Graf 2019 motifs and CHOPCHOP stem count |
| Last Reviewed | 2026-10-21 |

## 1. Overview

Guide RNA (gRNA/sgRNA) design in this repository evaluates and generates CRISPR guide candidates using sequence-composition metrics that are intended to correlate with editing efficiency and specificity. The design surface combines PAM-based candidate extraction with a score-based evaluation over GC content, poly-T content, seed-region GC, self-complementarity, and common restriction-site presence. The resulting `Score` is a heuristic quality measure rather than a learned or experimentally calibrated predictor; every designed candidate on a 20-nt NGG system whose 30-nt context fits in the sequence additionally carries the published Doench 2016 Rule Set 2 / Azimuth on-target score (`OnTargetScore`, from the 30-nt `Context30Mer`) and the Graf 2019 inefficiency motif (`GrafMotif`), and candidates can be ranked by that score (`GuideRnaParameters.Ranking = GuideRnaRanking.OnTargetRuleSet2`).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

An sgRNA consists of a spacer sequence that determines targeting specificity plus a scaffold sequence required for Cas binding. The original document cites Addgene for the canonical spacer-plus-scaffold structure and for the requirement that the target lie adjacent to a valid PAM. The same document also cites GC content, seed-region composition, and poly-T avoidance as practical guide-quality factors. Sources: Addgene CRISPR Guide, Wikipedia (Guide RNA, CRISPR gene editing, Protospacer adjacent motif).

### 2.2 Core Model

Guide candidates are extracted from PAM sites in a target region and scored as `Score = clamp(100 − penalty, 0, 100)`, where `penalty` is CHOPCHOP's own ascending penalty sum (B08 round 3 G5; before it the composite carried unsourced weights):

| Term | Effect on `Score` | Source |
|------|-------------------|--------|
| GC below `MinGcContent` or above `MaxGcContent` (default 40 / 70, edges inclusive) and `ScoreGcContent` | `+300` penalty, i.e. `Score = 0` (`CrisprDesigner.BadGcContentPenalty`) | CHOPCHOP `chopchop.py` `calcGCContent`: `if self.GCcontent > GC_HIGH or self.GCcontent < GC_LOW: self.score += SCORE['CRISPR_BAD_GC']` (`SCORE['CRISPR_BAD_GC'] = 300`, `GC_LOW = 40`, `GC_HIGH = 70`) |
| Each self-complementary 4-bp stem (`SelfComplementaryStems`), when `CheckSelfComplementarity` | `+1` per stem (`CrisprDesigner.SelfComplementarityStemPenalty`) | CHOPCHOP `self.score += self.folding * SCORE['FOLDING']`, `SCORE['FOLDING'] = 1` |
| Contains `TTTT` | none (issue only) | CRISPOR `crispor.py` L3026: warning "cannot be transcribed with a U6 or U3 promoter", no score; CHOPCHOP has no such term |
| Seed-region GC outside `30-80%` | none (issue only) | this library's own flag, no published threshold |
| Common restriction site detected | none (issue only) | this library's own flag |

**Provenance and the two declared deviations (B08 round 3 G5).** Fetched `chopchop.py` (sha256 `3c7a4344ccf2a77e70b638146989d7b9569665216bc15b476177db37630d5960`) ranks guides ascending by `score` (`sort_CRISPR_guides`), so every term is a penalty; the efficiency models enter the same sum as `score -= CoefficientsScore * SCORE['COEFFICIENTS']` (`COEFFICIENTS = 100`). Seqeron shows the penalty as a quality number `100 − penalty`, clamped to [0, 100] (the "100 minus" and the clamp are presentation, the weights are CHOPCHOP's). The 300 exceeds the whole range, so a guide outside the GC window scores 0; CHOPCHOP's flat amount is not proportional to the excess (the old Seqeron rule `2 × excess` was not from any tool). (1) CHOPCHOP ships the GC term **off** in every mode (`CRISPR_DEFAULT["SCORE_GC"] = False`, comment "this is already scored in many models!"; the only flag `-c` is `store_false`, so the command line cannot enable it): `GuideRnaParameters.ScoreGcContent` is `true` by default here because this composite has no efficiency model that covers GC; `false` reproduces CHOPCHOP's default. (2) `MinScore` (default 50) is not from any tool: CHOPCHOP and CRISPOR have no score threshold; their filters are the GC window (`filterGCmin/filterGCmax`, default 0..100 = off) and `filterSelfCompMax` (= `MaxSelfComplementaryStems`). With the GC term on, a 20-nt guide has `Score = 100 − stems` (stems ≤ 16) inside the window and 0 outside, so every `MinScore` from 1 to 84 selects exactly the guides inside the GC window. Oracle: the scoring lines copied verbatim and run in Python 3 + Biopython 1.88 on a guide table (GC 20–100 %, TTTT, seed GC, EcoRI), locked by `EvaluateGuideRna_Score_MatchesChopchopPenaltyOracle`. The raw-string `EvaluateGuideRna` and `CountSelfComplementaryStems` accept only A/C/G/T (`U` is read as `T`, case folded ASCII-only) and throw `ArgumentException` otherwise.

The seed region is defined as the last 10 nucleotides for PAM-after-target systems and the first 10 nucleotides for PAM-before-target systems. The 10-nt length is the upper bound of the Addgene CRISPR Guide's seed range, verified 2026-10-14 (https://www.addgene.org/guides/crispr/, "Basic Mechanisms of CRISPR"): *"the seed sequence (8–10 bases at the 3′ end of the gRNA targeting sequence)"*. That sentence describes Cas9; the PAM-before-target branch applies the same length at the PAM-proximal 5' end. The final score is clamped to `>= 0`.

### 2.3 Published on-target efficacy models (Doench 2014 Rule Set 1, Doench 2016 Rule Set 2)

The composition heuristic above is **not** the on-target activity model. Two published models are
implemented separately and are the ones to use for predicted cutting efficiency.

**Doench 2014 "Rule Set 1"** — `CrisprDesigner.CalculateOnTargetDoench2014(string context30Mer)`.

- **Source.** Doench, Hartenian, Graham, et al., *Nat Biotechnol* **32**:1262–1267 (2014), PMID
  25184501, doi:10.1038/nbt.3026. Coefficients are taken verbatim from the reference implementation
  `doenchScore.py` shipped at the root of the CRISPOR repository (Haeussler et al. 2016, *Genome
  Biol* **17**:148); the same table appears in CRISPOR's `crisporEffScores.py` as `doenchParams`.
- **Input convention — 30 nt, exactly.** `[4 nt upstream] + [20 nt protospacer] + [3 nt PAM] +
  [3 nt downstream]`, read on the protospacer strand, A/C/G/T only, case-insensitive. This is the
  window the reference itself cuts (`trimSeqs(seqs, -24, 6)` around the PAM start). The PAM must be
  **NGG**: offsets 25 and 26 of the 30-mer must both be `G`.
- **The model.** `score = 0.59763615` (intercept) `+ abs(10 − gcCount) × gcWeight` over the 20-nt
  protospacer `seq[4:24)`, with the paper's **two** GC coefficients — `gcWeight = −0.2026259` when
  `gcCount <= 10` and `−0.1665878` when `gcCount > 10` — `+` the 70 position-specific
  single-nucleotide and dinucleotide weights whose positions are 0-based offsets **into the 30-mer**
  (`seq[pos : pos+len(modelSeq)] == modelSeq`), all passed through the logistic sigmoid
  `1 / (1 + e^−score)`.
- **Output scale.** The sigmoid gives `(0, 1)`; this method returns it **× 100** and **unrounded**.
  CRISPOR's own table reports `int(100 × sigmoid)`, i.e. `Math.Floor` of what this returns.
- **The NGG guard is this library's explicit contract.** The bare reference function performs no PAM
  check (CRISPOR only ever calls it for SpCas9/NGG guides), but the published table itself presumes
  NGG: every dinucleotide feature at offset 24 ends in `G` (`(24,"AG"/"CG"/"TG")`) and the offset-26
  feature begins with `G` (`(26,"GT")`), so a non-NGG context would leave the model's domain.
- **Verification.** Coefficient table re-parsed out of the downloaded reference text and diffed
  entry by entry: 70/70 identical in order, intercept and both GC coefficients identical. Numerically
  cross-checked over 477 contexts (400 pseudo-random NGG 30-mers, fixed seed, plus the 77
  `Context30Mer` values `DesignGuideRnas` produces on three targets) — **0 mismatches, bit-equal**.
  Both GC branches and the `gcCount == 10` boundary are covered.

**Doench 2016 "Rule Set 2" / Azimuth** — `CalculateOnTargetRuleSet2`, same 30-nt context convention,
returned conventionally in `[0, 1]`. It is a trained gradient-boosted-tree model, not a coefficient
table, reproduced sklearn-free from Microsoft Research's Azimuth `V3_model_nopos` / `V3_model_full`;
`percentPeptide == -1` is Azimuth's "no gene position" sentinel and is answered with the
sequence-only model, as the reference does. `DesignGuideRnas` reports this score (and the 30-nt
context it was computed from) on every qualifying candidate and can rank by it
(`GuideRnaRanking.OnTargetRuleSet2`).

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
| `[DesignGuideRnas/EvaluateGuideRna] parameters` | `GuideRnaParameters?` | `GuideRnaParameters.Default` | Optional design thresholds | Defaults are `MinGcContent = 40`, `MaxGcContent = 70`, `MinScore = 50`, `AvoidPolyT = true`, `CheckSelfComplementarity = true`, `Ranking = Score`, `MaxSelfComplementaryStems = null` (no filter), `SelfComplementarityBackboneRegions = null`; both boolean flags are honoured by the scoring path (setting one to `false` removes its deduction — B08 F8) |
| `[EvaluateGuideRna] guideSequence` | `string` | required | Guide sequence to score | Null or empty throws `ArgumentNullException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Sequence` | `string` | Guide sequence being evaluated |
| `Position` | `int` | `pamSite.TargetStart` for designed guides and `-1` for standalone evaluation; since the 2026-09 review `PamSite.TargetStart` is a forward-strand, 0-based protospacer start on both strands (CRISPOR convention), so designed guides report forward-source coordinates for reverse-strand designs too |
| `IsForwardStrand` | `bool` | Strand orientation for designed candidates |
| `GcContent` | `double` | Overall GC percentage |
| `SeedGcContent` | `double` | GC percentage of the seed region |
| `HasPolyT` | `bool` | Whether the guide contains `TTTT` |
| `SelfComplementarityScore` | `double` | `SelfComplementaryStems / (length − 4)` in `[0, 1]`: CHOPCHOP's stem count over its own 4-bp window count (`range(0, len(fwd) − STEM_LEN)`); `0` when the guide is shorter than 5 nt. Replaced the former unsourced normalised complementary-pair fraction (B08 round 3 G13); not part of `Score` (the stems are) and not thermodynamic |
| `Score` | `double` | Composite quality number `clamp(100 − penalty, 0, 100)` with CHOPCHOP's penalty weights (1 per stem; flat 300 for a GC outside the window when `ScoreGcContent`); TTTT / seed GC / cloning site are issues only (B08 round 3 G5) — use `OnTargetScore` (Rule Set 2) for a published efficacy estimate |
| `Issues` | `IReadOnlyList<string>` | Human-readable issues recorded during evaluation |
| `System` | `CrisprSystem` | The CRISPR system the guide was evaluated for; designed guides carry the system of their own `PamSite` (SpCas9-NAG, SpCas9-NGA, SaCas9, SaCas9-NNGRRN, AsCas12a, LbCas12a, CasX keep their own PAM, guide length and seed orientation — B08 F4) |
| `SelfComplementaryStems` | `int` | CHOPCHOP's 4-bp self-complementary stem count (`CountSelfComplementaryStems`); drives the `Score` deduction and the optional `MaxSelfComplementaryStems` filter |
| `Context30Mer` | `string?` | Rule Set 2 context (4 nt + 20-nt protospacer + NGG PAM + 3 nt) read on the protospacer strand; `null` for standalone evaluation, non-20-nt/non-NGG systems, flanks outside the sequence, or a non-A/C/G/T base |
| `OnTargetScore` | `double?` | Doench 2016 Rule Set 2 / Azimuth score in `[0, 1]` from `Context30Mer`; `null` when the context is unavailable. Published model, independent of `Score` |
| `GrafMotif` | `GrafMotifType` | Graf 2019 inefficiency motif (`TtMotif` / `GccMotif` / `None`), evaluated for NGG systems only; reported as an issue, never deducted (as in CRISPOR) |
| `FullGuideRna` | `string` | Guide plus the fixed scaffold sequence |

### 3.3 Preconditions and Validation

`DesignGuideRnas(...)` throws `ArgumentNullException` for a null sequence and `ArgumentOutOfRangeException` when the requested region falls outside the sequence. `EvaluateGuideRna(...)` throws `ArgumentNullException` for null or empty guide strings, folds the input to upper case (ASCII only, `U` read as `T`), throws `ArgumentException` for any character other than A/C/G/T (B08 F70), and scores it as given without enforcing a specific guide length. In contrast, `DesignGuideRnas(...)` derives the guide length from the selected CRISPR system's PAM definition and extracted target sequence.

## 4. Algorithm

### 4.1 High-Level Steps

1. Resolve the CRISPR system into its PAM pattern, guide length, and orientation.
2. Find PAM sites in the requested region.
3. Extract each candidate guide sequence from the PAM-adjacent target interval.
4. Compute GC content, seed-region GC, poly-T presence, self-complementarity, and restriction-site presence.
5. Convert CHOPCHOP's penalty sum (stem count, flat GC penalty) into `Score = clamp(100 − penalty, 0, 100)`; list TTTT, seed-GC and cloning-site observations in `Issues`.
6. Yield only candidates whose score meets the configured minimum threshold.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Guide-design defaults confirmed in source:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MinGcContent` | `40` | Minimum acceptable GC percentage |
| `MaxGcContent` | `70` | Maximum acceptable GC percentage |
| `MinScore` | `50` | Minimum score for `DesignGuideRnas(...)` output. **Not from any tool** (CHOPCHOP/CRISPOR have no score threshold); for a 20-nt guide any value 1-84 selects exactly the guides inside the GC window |
| `AvoidPolyT` | `true` | When `true` a `TTTT` run deducts `20` and is flagged; when `false` neither happens (`HasPolyT` is still reported) — B08 F8 |
| `CheckSelfComplementarity` | `true` | When `true` each CHOPCHOP stem deducts `1` and is flagged; when `false` no deduction is made (CHOPCHOP `-SC/--noScoreSelfComp`); `SelfComplementaryStems` is still reported — B08 F8/F69 |
| `Ranking` | `Score` | Order of `DesignGuideRnas` output: descending `Score`, or descending `OnTargetScore` with `GuideRnaRanking.OnTargetRuleSet2` (candidates without a context rank last); ties by position, then forward strand first |
| `MaxSelfComplementaryStems` | `null` | Optional filter (CHOPCHOP `filterSelfCompMax`); `null` = report, do not filter |
| `SelfComplementarityBackboneRegions` | `null` | Optional scaffold regions the stem count may pair against (CHOPCHOP `--backbone`) |

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
- `CrisprDesigner.CalculateOnTargetDoench2014(string)`: the published Doench 2014 Rule Set 1 on-target
  score over the 30-nt context (see §2.3), returned on a 0–100 scale.
- `CrisprDesigner.CalculateOnTargetRuleSet2(string)` / `(string, int, double)`: the published Doench
  2016 Rule Set 2 / Azimuth on-target score (see §2.3).

### 5.2 Current Behavior

The current implementation matches the documented `40-70%` GC window, reports `TTTT` poly-T sequences (issue only, as CRISPOR does), uses a 10-base seed region, and reports an issue when the seed GC falls outside `30-80%` (no score effect). `DesignGuideRnas(...)` first filters PAM sites by a cut-site-in-region check and then evaluates the extracted guide. Standalone evaluation accepts non-standard guide lengths and scores them naturally from their composition, while designed guides follow the guide length of the selected CRISPR system for PAM-adjacent extraction. The guide-design pipeline passes `pamSite.System` straight into the shared evaluation core, so each designed candidate is scored with its own system's guide length and seed orientation (the earlier name-based remap to `SpCas9` is gone — B08 F4). The reported `Position` is copied from `pamSite.TargetStart`, which is a forward-strand coordinate on both strands (see `PAM_Site_Detection.md`). The `AvoidPolyT` and `CheckSelfComplementarity` parameter flags are honoured: each penalty is applied only when its flag is set (B08 F8). The derived `FullGuideRna` appends a fixed scaffold sequence to the spacer.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- PAM-adjacent guide extraction for supported CRISPR systems.
- GC-range screening, poly-T detection, seed-region GC evaluation, and self-complementarity checks.
- A fixed scaffold appended to the designed spacer for `FullGuideRna`.

**Intentionally simplified:**

- The `EvaluateGuideRna` quality score is CHOPCHOP's composition penalty sum (stems, GC window) rather than position-weighted or experimentally fitted; **consequence:** it ranks candidates heuristically and does not itself reproduce learned activity models — for experimentally-calibrated on-target efficacy call the dedicated `CalculateOnTargetDoench2014` (Rule Set 1) or `CalculateOnTargetRuleSet2` (Rule Set 2 / Azimuth) over the 30-nt context.
- The seed region is fixed to 10 nt using the upper bound of the cited `8-10` range; **consequence:** guides are evaluated with one fixed PAM-proximal window rather than a variable seed definition.
- ~~Designed-candidate metadata is only partially system-specific / some named systems share `SpCas9` metadata~~ — **withdrawn (B08 F4):** `PamSite` carries its real `System` and the evaluation uses it, so seed orientation and guide length follow the selected system.
- ~~`AvoidPolyT` and `CheckSelfComplementarity` are not honored~~ — **withdrawn (B08 F8):** both opt-outs are honoured.
- The Rule Set 2 context and `OnTargetScore` exist only for 20-nt NGG systems; **consequence:** for other systems (SaCas9, Cas12a, CasX, SpCas9-NAG/NGA) `Context30Mer` / `OnTargetScore` are `null` and only the composition `Score` (plus, for SaCas9, `CalculateSaCas9HitScore` off-target scoring) applies.
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
| All-A guide | Valid; Score 0 (GC outside the window: CHOPCHOP's flat 300) | issue reporting applies |
| All-G/C guide | Valid; Score 0 (GC 100 %: flat 300; with `ScoreGcContent = false`, 100 − stems) | issue reporting applies |
| Guide containing `TTTT` | Penalized and flagged (unless `AvoidPolyT = false`) | Poly-T detection is part of the score; `HasPolyT` is reported either way |

### 6.2 Limitations

The composite `Score` is a sequence-composition heuristic with no position weighting; the position-weighted, trained on-target models are the separate Doench 2014 (`CalculateOnTargetDoench2014`) and Doench 2016 Rule Set 2 (`CalculateOnTargetRuleSet2`, also attached to every qualifying designed candidate as `OnTargetScore`) implementations (§2.3), and Rule Set 2 is a trained gradient-boosted-tree model. Chromatin context is not modelled, and `DesignGuideRnas` performs no genome-wide uniqueness check (specificity is the separate off-target API: MIT/CFD aggregates in `Off_Target_Analysis.md`). Neither the heuristic `Score` nor `OnTargetScore` predicts editing outcomes directly.

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
