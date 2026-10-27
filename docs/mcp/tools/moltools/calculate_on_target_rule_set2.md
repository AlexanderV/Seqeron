# calculate_on_target_rule_set2

Doench et al. 2016 "Rule Set 2" / Azimuth on-target efficacy score for an SpCas9 guide (conventionally 0..1).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `calculate_on_target_rule_set2` |
| **Method ID** | `CrisprDesigner.CalculateOnTargetRuleSet2` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Doench et al. 2016 "Rule Set 2" / Azimuth on-target efficacy score for an SpCas9 guide, conventionally in 0..1 (higher = predicted more active) — the "Doench '16" efficiency score reported by CRISPOR. Rule Set 2 is a trained gradient-boosted-tree model, reproduced here (sklearn-free) from Microsoft Research's Azimuth model (`V3_model_nopos`, or `V3_model_full` with gene context).

Input is the model's 30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt PAM + 3 nt downstream, A/C/G/T only (case-insensitive). The PAM is **not audited**: Azimuth's own default (`pam_audit=True`) rejects a non-GG 30-mer, but CRISPOR calls it with `pam_audit=False` ("do not check for NGG PAM", `crisporEffScores.py` L1174-1178) so that NAG and NGA guides are scored too; the model was trained on NGG, so such a score is an extrapolation.

Pass `amino_acid_cut_position` and `percent_peptide` **together** to use Azimuth's gene-context (full) model instead of the sequence-only one; `percent_peptide = -1` is Azimuth's "no gene position" sentinel and is answered with the sequence-only model, as in the reference (`azimuth/model_comparison.py` selects `V3_model_nopos` when `percent_peptide == -1`). Any other `percent_peptide` outside `[0, 100]` is refused. Call to rank guides by the published on-target activity model.

Source: Doench et al. 2016, Nat Biotechnol 34:184 (PMID 26780180); trained model from https://github.com/MicrosoftResearch/Azimuth (BSD-3-Clause); provenance and reproduction in `scripts/azimuth/extract_azimuth_model.py`. For the published linear Rule Set 1 model use [calculate_on_target_doench2014](calculate_on_target_doench2014.md).

## Core Documentation Reference

- Source: [CrisprDesigner.cs#L1404](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs#L1404)
- Algorithm doc: [Guide_RNA_Design.md](../../../algorithms/MolTools/Guide_RNA_Design.md) §2.3 "Published on-target efficacy models"

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `context_30mer` | string | Yes | 30-nt context: 4 nt upstream + 20 nt protospacer + 3 nt PAM (trained on NGG, not audited) + 3 nt downstream. |
| `amino_acid_cut_position` | integer \| null | No | Optional amino-acid position of the cut site in the target protein, `>= 0` (requires `percent_peptide`); enables the gene-context model. Ignored when `percent_peptide = -1`. |
| `percent_peptide` | number \| null | No | Optional cut position as a percentage 0..100 along the coding sequence, or -1 for Azimuth's "no gene position" sentinel (which falls back to the sequence-only model); requires `amino_acid_cut_position`. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `score` | number | Predicted on-target activity, conventionally 0..1 (higher = more active). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Context 30-mer cannot be null or empty. |
| 1002 | amino_acid_cut_position and percent_peptide must be supplied together (gene-context model) or both omitted. |
| 1003 | Rule Set 2 requires a 30-nt context (4 upstream + 20 protospacer + 3 PAM + 3 downstream); got {length}. |
| 1004 | Rule Set 2 context must contain only A/C/G/T; found '{char}'. |
| 1006 | Rule Set 2 percent-peptide must be in [0, 100], or -1 for Azimuth's "no gene position" sentinel (which selects the sequence-only model). |
| 1007 | Rule Set 2 amino-acid cut position must be >= 0 (a residue count; Azimuth's own fixture spans 0..499). |

## Examples

### Example 1: sequence-only model

`context_30mer: "AAAAAAAAAAAAAAAGCTAACAGCAGGAGT"` → `score: 0.5290374377090171` (verified Azimuth oracle row, reference `0.529037`).

### Example 2: gene-context (full) model

Same context with `amino_acid_cut_position: 10`, `percent_peptide: 25.0` → `score: 0.5720089783500271`.

### Example 3: the `-1` sentinel

Same context with `amino_acid_cut_position: 10`, `percent_peptide: -1` → `score: 0.5290374377090171` (identical to Example 1, as in the reference).

### Example 4: half of the gene context

`amino_acid_cut_position: 10` without `percent_peptide` → error 1002.

## Performance

- **Time Complexity:** O(trees × depth) — fixed 30-nt input; featurization is constant-size.
- **Space Complexity:** O(features) per call; the trained model is loaded once and cached.

## See Also

- [calculate_on_target_doench2014](calculate_on_target_doench2014.md), [design_guide_rnas](design_guide_rnas.md), [evaluate_guide_rna](evaluate_guide_rna.md)
