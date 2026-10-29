# crispr_mit_specificity_score

Published MIT / Hsu 2013 aggregate guide specificity score.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `crispr_mit_specificity_score` |
| **Method ID** | `CrisprDesigner.CalculateMitSpecificityScore` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The **published** MIT / Hsu 2013 aggregate guide specificity score — the "Guide score" of crispr.mit.edu, implemented in CRISPOR as `calcMitGuideScore`. Scans the genome for off-targets with **1..max_mismatches** mismatches ([find_off_targets](find_off_targets.md)), scores each hit with the Hsu position-weight model, and returns `100 / (100 + Σ hitScores) × 100` in `0..100` (100 = no off-targets). The value is **unrounded** — the reference rounds to an integer.

The per-hit model is chosen by system, exactly as the reference's `annotateOfftargets` does: a **SaCas9** PAM (`NNGRRT`/`NNNRRT`) with a 21-nt guide is scored with the SaCas9-specific trained model of **Tycko et al. 2018** (`calcSaHitScore`), every other PAM-after-target system with the Hsu model. The reference's own `calcHitScore` also has a last-20 fudge for a non-SaCas9 guide longer than 20 nt and an `A`-pad for a 19-nt guide; the library keeps both branches for parity, but with the current system table neither can run (every PAM-after-target system has a 20-nt guide, or the SaCas9 21-nt guide that goes to the Tycko model).

### Alternative PAMs (`include_alternative_pams`)

CRISPOR runs ONE wider search that also accepts the alternative PAMs of its `offtargetPams` table (`NGG → NAG, NGA`; `NNGRRT → NNGRRN`; `TTTV → TTTN`) and then downweights such a hit's MIT score by **0.2** (`annotateOfftargets`: `if pam=="NGG" and otSeq[-2:]!="GG": mitScore = mitScore * 0.2` — "alternative PAMs represent only ~10% of all cleaveage events... so that these off-targets are not ranked among the top but still appear in the list somewhat"). With the default `include_alternative_pams: false` this tool scans the one system's PAM only, so it reports a guide as **more specific** than the reference. Set it to `true` to reproduce the reference's behaviour. Two properties of the rule are worth knowing:

- It is gated on the **guide's** PAM being `NGG`, so an `SpCas9_NAG` guide is never downweighted.
- Only the MIT score is downweighted; the CFD aggregate has its own alternative-PAM entry point ([crispr_cfd_specificity_score](crispr_cfd_specificity_score.md), `include_alternative_pams`) which scores each hit with its own PAM and applies no 0.2 factor — its PAM table already penalises non-NGG PAMs.

`SpCas9` scans both `NAG` and `NGA` (systems `SpCas9_NAG` / `SpCas9_NGA`, hits downweighted by 0.2); `SaCas9` scans `NNGRRN` (system `SaCas9_NNGRRN`; hits enter UNweighted because the reference's rule is gated on a main `NGG`; a site matching both `NNGRRT` and `NNGRRN` is counted once); Cas12a is refused. `min_alternative_pam_hit_score` is the reference's `ALTPAMMINSCORE` / `--minAltPamScore` (`crispor.py` L318, L4108-4122): alternative-PAM hits whose unweighted hit score is below it are dropped. Its default is CRISPOR's own `1.0` (B08 round 3 G6; it was `0` = no filtering before — pass `0` to keep every alternative-PAM hit). The C# API also exposes the rule over an arbitrary PAM-annotated hit set (`CrisprDesigner.CalculateMitSpecificityScore(guide, IEnumerable<PamAnnotatedOffTarget>, system)` with `CrisprDesigner.FindOffTargetsWithPam`), where an `NGA` hit IS counted and downweighted.

Cas12a/Cpf1-class and CasX systems (PAM **before** target) are **refused**: the reference reports no MIT off-target score for them at all (`mitScore = 0.0`, `guideScore = -1`). Prefer this tool over [crispr_specificity_score](crispr_specificity_score.md), which is this library's own unsourced heuristic.

### On-target (`on_target_start` / `on_target_end`) and repeat guides (`max_occurrences`)

CRISPOR's `annotateOfftargets` (crispor.py L1917-1928) drops exactly **one kind** of 0-mismatch hit, the guide's own on-target, and scores every other perfect match as a full off-target (`calcHitScore = 100`; CFD: the PAM score). Without an interval this is its "auto-ontarget" mode (`inChrom=='' and foundOneOntarget==False`): the **first** exact match is the on-target. The reference's "first" is the aligner's output order; here it is the scan order (forward-strand PAMs by ascending position, then reverse-strand PAMs) — exact matches of one PAM class score alike, so the score does not depend on it. A guide with 2 / 5 / 50 / 500 perfect sites therefore scores `100/N` = 50 / 20 / 2 / 0.2 (MIT; CRISPOR prints 50 / 20 / 2 / 0 after integer rounding — checked by executing its `annotateOfftargets`). Pass `on_target_start` (0-based, inclusive) and `on_target_end` (exclusive), forward-strand coordinates of the `genome` string covering the on-target's protospacer **and** PAM, to name the on-target locus instead: exact matches inside the interval are dropped, every other one is scored. Give both or neither. A hit only partly inside the interval is not skipped; either strand's hit is tested on its forward-strand footprint. With `include_alternative_pams` the first exact match of the **main-PAM** scan is the on-target; exact alternative-PAM matches are kept.

`max_occurrences` (default 60000, CRISPOR's `MAXOCC`) reproduces its repeat rule (`crispor.py` L1906-1910): a guide with more alignments than this — windows of the guide's length on either strand, **whatever the PAM**, within `max_mismatches` substitutions (BWA's `X0 + X1`; checked against BWA 0.7.17 with the reference's own command line) — is suppressed and scores `0`.


## Core Documentation Reference

- Source: [CrisprDesigner.cs#L1662](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L1662)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `guide_sequence` | string | Yes | Guide RNA (length must match the system). |
| `genome` | string | Yes | Genome / reference to scan. |
| `max_mismatches` | integer | No | 0..5 (default 4, the reference's default). |
| `system_type` | enum | No | CRISPR system (default `SpCas9`); must be PAM-after-target. |
| `include_alternative_pams` | boolean | No | Also scan the reference's alternative PAMs (SpCas9: NAG, NGA, downweighted by 0.2; SaCas9: NNGRRN, unweighted) (default `false`). |
| `min_alternative_pam_hit_score` | number | No | With `include_alternative_pams`: minimum unweighted hit score (0..100) for an alternative-PAM hit to count (reference `ALTPAMMINSCORE`; default `1.0`, the reference's own default; `0` = no filtering). |
| `on_target_start` | integer | No | Optional on-target interval start (0-based, inclusive; with `on_target_end`). Exact matches whose protospacer+PAM footprint lies inside are skipped; other perfect matches are scored as off-targets. Default: omitted = the first exact match is the on-target (CRISPOR auto-ontarget). |
| `on_target_end` | integer | No | Optional on-target interval end (0-based, exclusive). |
| `max_occurrences` | integer | No | Repeat limit (default `60000`, CRISPOR `MAXOCC`): a guide with more alignments than this (guide-length windows, either strand, any PAM, within `max_mismatches`) scores `0`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `specificity` | number | MIT/Hsu aggregate specificity 0..100 (100 = no off-targets). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Guide sequence cannot be null or empty |
| 1002 | Genome cannot be null or empty |
| 1003 | Maximum mismatches must be in the range 0..5 |
| 4001 | Guide length does not match the system's guide length |
| 4002 | The MIT/Hsu 2013 off-target model is not defined for PAM-before-target (Cas12a-class / CasX) systems |

## Examples

### Example 1: one PAM-distal off-target → 50

Guide `GACGCATAAAGATGAGACGC` vs genome `TACGCATAAAGATGAGACGCAGG`: one off-target differing at the PAM-distal position 1 (Hsu weight 0) → single-hit score 100 → `100/(100+100)×100 = 50.0`.

### Example 2: SaCas9 21-nt guide

Guide `ACGTACGTACGTACGTACGTA` vs genome `ACGTACGTACGTACGTACGTCAAGAAT` with `system_type: SaCas9`: one off-target differing at the PAM-proximal position → Tycko 2018 `calcSaHitScore` = 4.2613097364785695 → `95.91285612347565`.

### Example 3: alternative-PAM hit folded in

Guide `GACGCATAAAGATGAGACGC` vs genome `GACGCTTAAAGATGAGACGCAGGTTTTGACGCATAAAGATGAGACTCTAG` holds two off-targets: one on the main `AGG` PAM (`calcHitScore` 60.5) and one on an alternative `TAG` PAM (31.4999…). With `include_alternative_pams: false` → `100/(100+60.5)×100 = 62.30529595015576`; with `true` → `100/(100+60.5+0.2×31.4999…)×100 = 59.95203836930455`. Counting the alternative-PAM hit at full weight would give 52.0833…, so the 0.2 factor is visible in the value.

### Example 4: duplicated perfect site with an on-target interval

Guide `GACGCATAAAGATGAGACGC` vs genome `GACGCATAAAGATGAGACGC`+`AGG`+`TTTT`+`GACGCATAAAGATGAGACGC`+`TGG`+`TTTTTT`: default (auto-ontarget: the first perfect site is the on-target) → the second is an off-target with `calcHitScore = 100` → `100/(100+100)×100 = 50.0`; `on_target_start: 0, on_target_end: 23` names the first site → `50.0`; `on_target_start: 27, on_target_end: 50` names the second → `50.0` as well; an interval holding both (`0, 50`) → `100.0`.

### Example 5: on-target only → `specificity: 100.0`.

## Performance

- **Time Complexity:** O(genome × guide).
- **Space Complexity:** O(number of off-targets).

## See Also

- [crispr_cfd_specificity_score](crispr_cfd_specificity_score.md), [find_off_targets](find_off_targets.md), [crispr_specificity_score](crispr_specificity_score.md)
