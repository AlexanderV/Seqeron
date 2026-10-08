# validate_probe

Validate a probe's specificity against reference sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MolTools |
| **Tool Name** | `validate_probe` |
| **Method ID** | `ProbeDesigner.ValidateProbe` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

1. **Hits.** Scans each reference sequence with ungapped `max_mismatches`-tolerant (Hamming) matching (canonical `ApproximateMatcher.FindWithMismatches`, the given strand) and counts the total hits, intended site included; more than one hit is an issue. `specificityScore` is a library uniqueness score — **0 hits → 0.0**, **N hits → 1/N** — not a published metric.
2. **Self-structure.** For probes of ≤ 60 nt made of A/C/G/T, Primer3's thermodynamic hybridization-probe screen: ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin Tm at the Primer3 probe conditions (50 nM oligo, 50 mM monovalent, no Mg²⁺/dNTP) must not exceed 47 °C (PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH); values equal primer3-py `calc_homodimer` / `calc_end_stability` / `calc_hairpin`. Longer (thal.c THAL_MAX_ALIGN = 60) or non-ACGT probes use Primer3's alignment-mode internal-oligo self-dimer screen (dpal `self_any` / `self_end` > 12.00, PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END, no length limit) and a sequence-only inverted-repeat hairpin stem.
3. **Cross-hybridization (optional).** With `non_target_sequences`, Kane et al. (2000) criteria on both strands of each non-target: overall identity (identical columns of the best local alignment — Smith–Waterman–Gotoh, BLAST+ blastn scoring +2/−3, gap 5/2 — over the probe length) **> 75 %**, or a contiguous identical stretch (longest common substring) **> 15 nt**. Each aligned site also reports the ntthal duplex Tm of the probe with the site's complementary strand (primer3-py `calc_heterodimer`, ≤ 60-nt ACGT probes), flagged when above the optional `max_duplex_tm` (Rouillard et al. 2003, OligoArray 2.0).

`isValid` is true when no issue was recorded. Matching is case-insensitive.

## Core Documentation Reference

- Source: [ProbeDesigner.cs#L1366](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1366)
- Algorithm doc: [Probe_Validation.md](../../../algorithms/MolTools/Probe_Validation.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `probe_sequence` | string | Yes | Probe sequence to validate (non-null; empty → invalid). |
| `reference_sequences` | string[] | Yes | Reference sequences to scan (non-null). |
| `max_mismatches` | integer | No | Maximum allowed mismatches (≥ 0, default 3). |
| `self_complementarity_threshold` | number | No | Legacy fold-back-fraction limit (default 0.3); kept for compatibility, no longer used by the screen. |
| `non_target_sequences` | string[] | No | Known non-target sequences for the Kane criteria (default none). |
| `max_non_target_identity` | number | No | Kane identity threshold in [0,1]; flagged when strictly above (default 0.75). |
| `max_contiguous_match` | integer | No | Kane contiguous-identity threshold in nt; flagged when strictly longer (≥ 0, default 15). |
| `max_duplex_tm` | number | No | Optional OligoArray 2.0-style threshold (°C) on each non-target site's ntthal duplex Tm with the probe; flagged when strictly above (default none). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `isValid` | boolean | True when no issue was recorded. |
| `specificityScore` | number | Library uniqueness score `[0,1]` (0 hits → 0, N → 1/N). |
| `offTargetHits` | integer | Total approximate-match hits across references (intended site included). |
| `selfComplementarity` | number | Fold-back fraction `[0,1]` (informational library metric). |
| `selfAny` / `selfEnd` | number | Primer3 alignment-mode internal-oligo self_any / self_end (dpal; the fallback self-dimer criterion, limit 12.00). |
| `hasSecondaryStructure` | boolean | ntthal hairpin Tm > 47 °C (or the fallback stem screen). |
| `issues` | string[] | Reported issues. |
| `thermodynamicScreen` | boolean | True when the Primer3 ntthal screen was applied. |
| `selfDimerTm` / `selfEndDimerTm` / `hairpinTm` | number \| null | ntthal Tm values in °C (null with the fallback screen). |
| `crossHybridization` | object[] | Per non-target strand: `nonTargetIndex`, `reverseComplementStrand`, `identity`, `identicalColumns`, `alignmentScore`, `longestContiguousMatch`, `exceedsIdentityThreshold`, `exceedsContiguousThreshold`, `exceedsDuplexTmThreshold`, `crossHybridizes`, `siteStart`/`siteEnd` (aligned site in the assessed strand), `duplexTm` (ntthal Tm of the probe on the site's complementary strand, primer3-py `calc_heterodimer`; null for > 60-nt / non-ACGT). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Probe sequence cannot be null |
| 1002 | Reference sequences cannot be null |
| 1003 | Maximum mismatches cannot be negative |
| 1004 | Non-target identity threshold must be in [0, 1] |
| 1005 | Contiguous-match threshold cannot be negative |

## Examples

### Example 1: Unique probe → specificity 1.0

Probe `ATCGATCGATCGATCGATCG` in `NNNNNATCGATCGATCGATCGATCGNNNN` with `max_mismatches = 0`: 1 hit → specificity 1.0.

**Input:** `{ "probe_sequence": "ATCGATCGATCGATCGATCG", "reference_sequences": ["NNNNNATCGATCGATCGATCGATCGNNNN"], "max_mismatches": 0 }`

**Response (abridged):** `{ "offTargetHits": 1, "specificityScore": 1.0 }`

### Example 2: Repetitive probe → 1/N specificity

10-mer poly-A in a 34-mer poly-A: 25 exact-match positions → specificity `1/25 = 0.04`.

### Example 3: Kane cross-hybridization

Probe `TATGCCTCCGGTACATCAACTACAGTTAGCCTTAAGAGAAAAATCCCAAA` against a non-target carrying the probe with a substitution at every 5th position: best local alignment 40/50 identical → identity 0.80 > 0.75 → issue `Cross-hybridization risk with non-target 0: identity 80%, …`, `isValid = false`.

## See Also

- [design_probes](design_probes.md), [design_probes_primer3](design_probes_primer3.md), [design_tiling_probes](design_tiling_probes.md), [crispr_specificity_score](crispr_specificity_score.md)
