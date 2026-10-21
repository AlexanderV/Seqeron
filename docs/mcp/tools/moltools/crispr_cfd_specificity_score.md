# crispr_cfd_specificity_score

Published Doench 2016 CFD aggregate guide specificity score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `crispr_cfd_specificity_score` |
| **Method ID** | `CrisprDesigner.CalculateCfdSpecificityScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The **published** aggregate CFD (Cutting Frequency Determination, Doench et al. 2016) guide specificity score, following CRISPOR's `annotateOfftargets` (`calcCfdScore(guideSeq, otSeq)` on each hit, summed and aggregated by `calcMitGuideScore`). By default only the system's own PAM is scanned (see *Alternative PAMs* below). Scans the genome for off-targets with **1..max_mismatches** mismatches ([find_off_targets](find_off_targets.md)), scores each hit with the CFD mismatch/PAM matrices using **that site's own PAM** (read on the protospacer strand, so reverse-strand hits are handled correctly), and aggregates the **raw** per-hit CFD values in `[0, 1]`:

| `aggregation` | Formula | Notes |
|---|---|---|
| `MitStyle` (default) | `100 / (100 + Σ CFD) × 100` | CRISPOR's default, `calcMitGuideScore(sum(cfdScores))`. 100 = no off-targets. |
| `ParkinsonNormalized` | `100 / Σ CFD` | CRISPOR's `doCfdFix` (`fixCfd=1`) variant. Not bounded by 100; **errors when Σ CFD = 0**. |

`ParkinsonNormalized` is written in the reference as `calcCfdGuideScore(sum) = 100 × 100 / sum`, but the same `doCfdFix` flag also rescales every per-hit CFD by ×100 inside `calcCfdScore`; over the raw per-hit values this API takes, the two factors cancel to `100 / Σ raw`.

### Alternative PAMs (`include_alternative_pams`)

CRISPOR runs ONE wider search that also accepts the alternative PAMs of its `offtargetPams` table (`NGG → NAG, NGA`; `NGA → NGG`) and computes `cfdScore = calcCfdScore(guideSeq, otSeq)` for every hit of that search (`annotateOfftargets`, crispor.py L1905-1990). With the default `include_alternative_pams: false` this tool scans the one system's PAM only, so a guide whose only off-target sits on `NAG`/`NGA` scores 100 here but not in the reference. Set it to `true` to route to `CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams` and reproduce the reference's hit set (the same collection as [crispr_mit_specificity_score](crispr_mit_specificity_score.md): main PAM + alternative PAMs, a site matching both counted once). Two properties, read from the source:

- **No 0.2 factor.** `if pam=="NGG" and otSeq[-2:]!="GG": mitScore = mitScore * 0.2` (L1955) multiplies the MIT score only; `cfdScore` is taken afterwards from the full `otSeq` (the CFD PAM table, e.g. `AG` = 0.2593, is what penalises the PAM).
- **The cutoff is on the MIT hit score.** The search pre-filter (L4104-4122, `ALTPAMMINSCORE = 1.0`) drops an alternative-PAM hit whose `calcHitScore(guideSeq, tSeqNoPam)` is `< 1.0` — the MIT/Hsu hit score, even though this aggregate sums CFD values. `min_alternative_pam_hit_score` (default 1.0; `0` = off) is that cutoff.

### On-target interval (`on_target_start` / `on_target_end`)

By default every exact (0-mismatch) match in the genome is dropped as "the on-target". CRISPOR's `annotateOfftargets` (crispor.py L1920-1929) drops only the 0-mismatch hits **inside the on-target interval it is given** (`chrom==inChrom and start >= inStart and end <= inEnd`) and scores any further perfect match as a full off-target (`calcHitScore = 100`; CFD: the PAM score). Pass `on_target_start` (0-based, inclusive) and `on_target_end` (exclusive), forward-strand coordinates of the `genome` string covering the on-target's protospacer **and** PAM, to reproduce that; give both or neither. A hit only partly inside the interval is not skipped; either strand's hit is tested on its forward-strand footprint. CRISPOR's no-interval "auto-ontarget" mode (drop the first 0-mismatch hit it meets, in aligner output order) is not offered: that order is unspecified and, under CFD or with alternative PAMs, the exact hits differ in score.

CFD is defined only for a **20-nt protospacer**, so SaCas9 (21 nt), Cas12a-class and CasX are **refused**. The reference does not score them either: it sets the per-hit sentinel `cfdScore = -1` for a SaCas9 PAM and `cfdScore = 0.0` for a Cpf1 PAM, and the guide-level `guideCfdScore = -1` for a Cpf1 PAM. For SaCas9 its guide-level value is not a sentinel but a vacuous `calcMitGuideScore(sum([])) = 100` — every hit having been dropped by the `if cfdScore != -1` filter — which is exactly the misleading perfect score this tool refuses to return. Prefer this tool over [crispr_specificity_score](crispr_specificity_score.md), which is this library's own unsourced heuristic.

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L2017](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L1941)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `guide_sequence` | string | Yes | Guide RNA (must be 20 nt for the supported systems). |
| `genome` | string | Yes | Genome / reference to scan. |
| `max_mismatches` | integer | No | 0..5 (default 4, the reference's default). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`); must be a 20-nt PAM-after-target system. |
| `aggregation` | enum | No | `MitStyle` (default) or `ParkinsonNormalized`. |
| `include_alternative_pams` | boolean | No | Also scan the reference's alternative PAMs (SpCas9: NAG, NGA; NGA: NGG), each hit CFD-scored by its own PAM, no 0.2 factor (default `false`). |
| `min_alternative_pam_hit_score` | number | No | With `include_alternative_pams`: minimum MIT/Hsu hit score (0..100) for an alternative-PAM hit to count (reference `ALTPAMMINSCORE`; default `1.0`; `0` = no filtering). |
| `on_target_start` | integer | No | Optional on-target interval start (0-based, inclusive; with `on_target_end`). Exact matches whose protospacer+PAM footprint lies inside are skipped; other perfect matches are scored as off-targets. Default: omitted = every exact match dropped. |
| `on_target_end` | integer | No | Optional on-target interval end (0-based, exclusive). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `specificity` | number | Aggregate CFD specificity (0..100 for `MitStyle`; unbounded for `ParkinsonNormalized`). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Guide sequence cannot be null or empty |
| 1002 | Genome cannot be null or empty |
| 1003 | Maximum mismatches must be in the range 0..5 |
| 4001 | Guide length does not match the system's guide length |
| 4002 | The CFD model is defined only for a 20-nt protospacer (SaCas9 / Cas12a-class / CasX are refused) |
| 4003 | The ParkinsonNormalized CFD aggregate is undefined for a CFD sum of 0 |

## Examples

### Example 1: reverse-strand off-target with a NAG PAM

Guide `GACGCATAAAGATGAGACGC` vs genome `CTTGCGTCTCATCTTTATGCGTA` with `system_type: SpCas9_NAG`: one reverse-strand off-target `TACGCATAAAGATGAGACGC` with PAM `AAG`. Per-hit CFD = `rG:dA@1 (1.0) × PAM "AG" (0.259259259)` → `MitStyle = 99.74141115651946`, `ParkinsonNormalized = 385.7142861`.

### Example 2: NAG-only off-target, main PAM vs alternative PAMs

Guide `GACGCATAAAGATGAGACGC` vs genome `GACGCTTAAAGATGAGACGC` + `TAG` + `TTTTTT`: one 1-mismatch off-target (index 5) on the alternative PAM `TAG`, MIT hit score 60.5 (>= 1.0, kept), per-hit CFD (CRISPOR's `calcCfdScore`) = mismatch factor 0.714285714 × PAM `AG` 0.259259259 = 0.18518518492592592. Default (`include_alternative_pams: false`) → `100.0`; `true` → `100/(100 + 0.18518518492592592) × 100 = 99.81515711670932`.

### Example 3: duplicated perfect site with an on-target interval

Guide `GACGCATAAAGATGAGACGC` vs genome `GACGCATAAAGATGAGACGC`+`AGG`+`TTTT`+`GACGCATAAAGATGAGACGC`+`TGG`+`TTTTTT`: default → `100.0`; `on_target_start: 0, on_target_end: 23` → the second perfect site is an off-target with CFD = 1.0 (`NGG` PAM) → `100/(100+1)×100 = 99.00990099009901`.

### Example 4: on-target only → `specificity: 100.0` (`MitStyle`); `ParkinsonNormalized` errors (Σ CFD = 0).

## Performance

- **Time Complexity:** O(genome × guide).
- **Space Complexity:** O(number of off-targets).

## See Also

- [crispr_mit_specificity_score](crispr_mit_specificity_score.md), [find_off_targets](find_off_targets.md), [crispr_specificity_score](crispr_specificity_score.md)
