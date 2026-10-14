# chou_fasman_assign

Discrete Chou-Fasman (1978) secondary-structure assignment: one H/E/T/C state per residue.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `chou_fasman_assign` |
| **Method ID** | `SequenceStatistics.PredictSecondaryStructureChouFasman` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Assigns every residue of a protein to **α-helix (`H`)**, **β-strand (`E`)**, **β-turn (`T`)** or
**coil (`C`)** with the Chou & Fasman (1978) prediction rules and the 1978 conformational
parameters (Pα, Pβ, Pt) and bend frequencies f(i)..f(i+3):

1. **Nucleation** — helix: any 6-residue window with ≥ 4 helix formers; sheet: any 5-residue window
   with ≥ 3 sheet formers.
2. **Extension** — each nucleus grows in both directions while the tetrapeptide formed by the new
   residue and its three segment neighbours has mean propensity ≥ 1.00.
3. **Acceptance** — helix if ⟨Pα⟩ > 1.03 and ⟨Pα⟩ > ⟨Pβ⟩; strand if ⟨Pβ⟩ > 1.05 and ⟨Pβ⟩ > ⟨Pα⟩.
4. **Overlap** — each run covered by both is helix if ⟨Pα⟩ > ⟨Pβ⟩ over the run, otherwise strand.
5. **β-turn** — tetrapeptide i..i+3 with p(t) = f(i)·f(i+1)·f(i+2)·f(i+3) > 7.5×10⁻⁵, ⟨Pt⟩ > 1.00 and
   ⟨Pα⟩ < ⟨Pt⟩ > ⟨Pβ⟩; turns override H/E.

`ruleSet` selects the nucleation rules:

- **`Standard`** (default) — a former is a residue with P > 1.00 (Chen, Gu & Huang 2006 statement of the rules).
- **`Graded1978`** — the graded former/breaker classes of Chou & Fasman (1978) (Prevelige & Fasman 1989
  Table II): Hα/hα count 1, Iα ½, at most one breaker (bα/Bα) per helix nucleus, sheet likewise with
  Hβ/hβ and bβ/Bβ; Pro only in the first three helix residues, Pro/Arg/Trp boundary reassignments.

Residues other than the 20 standard amino acids have no parameters: they are never formers, windows
containing them are ineligible, and they stay `C`. Input is case-insensitive. Arithmetic is exact
(integer propensity units). Typical Q3 accuracy of the method is ~50–60 %. For a windowed propensity
profile instead of an assignment use [predict_chou_fasman](predict_chou_fasman.md).

## Core Documentation Reference

- Source: [SequenceStatistics.cs#L2100](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs#L2100)
- Algorithm: `docs/algorithms/Statistics/Secondary_Structure_Prediction.md`
- Evidence: `docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `proteinSequence` | string | Yes | Protein sequence, one-letter code (min length 1, case-insensitive) |
| `ruleSet` | string | No | `Standard` (default) or `Graded1978` (case-insensitive) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `assignment` | string | One state code per residue (same length as the input): `H`, `E`, `T` or `C` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1002 | Unknown Chou-Fasman rule set '…'. Expected one of: Standard, Graded1978. |

## Examples

### Example 1: Helix followed by a strand

**User Prompt:**
> Predict the Chou-Fasman secondary structure of EEEEEEVVVVVV.

**Expected Tool Call:**
```json
{
  "tool": "chou_fasman_assign",
  "arguments": { "proteinSequence": "EEEEEEVVVVVV" }
}
```

**Response:**
```json
{ "assignment": "HHHHEEEEEEEE" }
```
The E₆ helix extends over the whole chain, the strand covers residues 4..11; on the overlap run
EEVVVVVV ΣPα 9.38 < ΣPβ 10.94, so the strand wins and residues 0..3 stay helix (hand trace).

### Example 2: Graded 1978 classes

**Expected Tool Call:**
```json
{
  "tool": "chou_fasman_assign",
  "arguments": { "proteinSequence": "AAAHHS", "ruleSet": "Graded1978" }
}
```

**Response:**
```json
{ "assignment": "HHHHHH" }
```
A A A are Hα formers and H H are Iα (½ each): 4 formers in 6, no breaker → helix. The `Standard`
rule set counts only P > 1.00 formers (3) and returns `CCCCCC`.

### Example 3: β-turn

**Input:** `{ "proteinSequence": "NPDG" }` → `{ "assignment": "TTTT" }` (p(t) = 1.32×10⁻³, ⟨Pt⟩ = 1.525).

## Performance

- **Time Complexity:** O(n) per nucleus scan; extension and overlap are linear in the segment length.
- **Space Complexity:** O(n).

## See Also

- [predict_chou_fasman](predict_chou_fasman.md) — windowed Chou-Fasman propensity profile
- [hydrophobicity_profile](hydrophobicity_profile.md) — Kyte-Doolittle hydropathy
