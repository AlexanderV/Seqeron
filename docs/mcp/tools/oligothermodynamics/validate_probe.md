# validate_probe

Validate a probe's specificity against reference sequences.

## Overview

| Property | Value |
|----------|-------|
| **Server** | OligoThermodynamics |
| **Tool Name** | `validate_probe` |
| **Method ID** | `ThermodynamicProbeDesigner.ValidateProbe` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

1. **Hits.** Scans each reference sequence with ungapped `max_mismatches`-tolerant (Hamming) matching (canonical `ApproximateMatcher.FindWithMismatches`; the given strand by default, and with `both_strands = true` also the probe's reverse complement — its sites on the other strand of a double-stranded reference, as blastn `-strand both`; a reference position hit in both orientations, i.e. a reverse-palindromic probe or site, is one site and is cross-hybridizing when either orientation meets a Kane criterion) and counts the total hits, intended site included (`offTargetHits`). Each hit is judged by the Kane et al. (2000) criteria on its ungapped diagonal — identity (L − mismatches) / L **> `max_non_target_identity`** (0.75) or a run of identical positions **> `max_contiguous_match`** (15 nt); those sites are counted in `crossHybridizingHits`, and **more than one** is an issue (the probe can bind besides its intended site). With the default radius of 3 mismatches every hit of a ≥ 13-nt probe meets the identity criterion, so the decision equals the former "> 1 hit" rule there. `specificityScore` (**0 hits → 0.0**, **N hits → 1/N**) and the default radius 3 are **library conventions** — no published 1/N specificity score or 3-mismatch hybridization tolerance exists (Kane 2000, OligoArray 2.0 and Primer3 PRIMER_INTERNAL_MAX_LIBRARY_MISHYB judge a site by identity / duplex stability / alignment score); the score is reported only and does not enter `isValid`. To reach every ungapped site above 75 % identity pass `max_mismatches` = ⌈L/4⌉ − 1.
2. **Self-structure.** For A/C/G/T probes of ≤ `thermodynamic_screen_max_length` nt (default 60 = Primer3's THAL_MAX_ALIGN; opt-in up to 10 000 = thal.c compiled with `-DTHAL_MAX_ALIGN`) (with `thermodynamic_screen`, default true), Primer3's thermodynamic hybridization-probe screen: ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin Tm at the reaction conditions `monovalent_mm` / `divalent_mm` / `dntp_mm` / `dna_conc_nm` (default the Primer3 probe conditions 50 mM / 0 / 0 / 50 nM; Primer3 `_pr_data_control` legality: monovalent and oligo > 0, Mg²⁺ and dNTP ≥ 0) must not exceed `max_structure_tm` (default 47 °C, PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH); values equal primer3-py `calc_homodimer` / `calc_end_stability` / `calc_hairpin`. Longer or non-ACGT probes, and `thermodynamic_screen = false`, use Primer3's alignment-mode internal-oligo self-dimer screen (dpal `self_any` / `self_end` > `max_self_any` / `max_self_end`, default 12.00, PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END, no length limit) and the sequence-only stem-loop screen of `hairpin_potential` (`ThermodynamicPrimerDesigner.HasHairpinPotential`: an exactly complementary stem of ≥ 4 bp — library convention — closing a loop of ≥ 3 nt, thal.c `min_hrpn_loop`; almost every random probe of > 60 nt contains one, so its flag is reported in `hasSecondaryStructure` and `warnings` but is **not an issue** and does not decide `isValid` — audit round 5, A5-5; raise `thermodynamic_screen_max_length` for a sourced hairpin-Tm decision on long probes).
3. **Cross-hybridization (optional).** With `non_target_sequences`, Kane et al. (2000) criteria on both strands of each non-target: overall identity (identical columns of the best local alignment — Smith–Waterman–Gotoh, BLAST+ blastn scoring +2/−3, gap 5/2 — over the probe length) **> 75 %**, or a contiguous identical stretch (longest common substring) **> 15 nt**. Each aligned site also reports the ntthal duplex Tm of the probe with the site's complementary strand (primer3-py `calc_heterodimer` at the reaction conditions, ACGT probe and site; computed when the probe or the site is ≤ `thermodynamic_screen_max_length` nt, thal.c `thal_check_errors`), flagged when above the optional `max_duplex_tm` (Rouillard et al. 2003, OligoArray 2.0).

`isValid` is true when no issue was recorded (`warnings` do not count). Matching is case-insensitive.

## Core Documentation Reference

- Source: [Seqeron.Genomics.OligoThermodynamics/ThermodynamicProbeDesigner.Core.cs#L1129](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.OligoThermodynamics/ThermodynamicProbeDesigner.Core.cs#L1129)
- Algorithm doc: [Probe_Validation.md](../../../algorithms/MolTools/Probe_Validation.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `probe_sequence` | string | Yes | Probe sequence to validate. |
| `reference_sequences` | string[] | Yes | Reference sequences to scan for off-target hits. |
| `max_mismatches` | integer | No | Search radius of the ungapped reference scan: maximum mismatches per site (default 3, a library convention; a found site is an issue only when it meets the Kane criteria — pass ceil(L/4) - 1 to reach every site above 75 % identity). |
| `self_complementarity_threshold` | number | No | Legacy fold-back-fraction limit (default 0.3); kept for compatibility, no longer used by the screen (the fallback self-dimer limit is Primer3's PRIMER_INTERNAL_MAX_SELF_ANY/_END = 12). |
| `non_target_sequences` | string[] | No | Optional known non-target sequences for the Kane et al. (2000) cross-hybridization criteria (both strands). |
| `max_non_target_identity` | number | No | Kane identity threshold in [0,1]; a reference site or non-target strand with identity strictly above it is flagged (default 0.75). |
| `max_contiguous_match` | integer | No | Kane contiguous-identity threshold in nt; a reference site or non-target strand with a longer identical stretch is flagged (default 15). |
| `max_duplex_tm` | number | No | Optional OligoArray-style threshold (°C): a non-target site whose ntthal duplex Tm with the probe is above it is flagged (default none). |
| `monovalent_mm` | number | No | PRIMER_INTERNAL_SALT_MONOVALENT: monovalent cation concentration in mM (> 0) for the ntthal screen and the non-target duplex Tm (default 50). |
| `divalent_mm` | number | No | PRIMER_INTERNAL_SALT_DIVALENT: Mg²⁺ concentration in mM (>= 0; default 0). (default: 0) |
| `dntp_mm` | number | No | PRIMER_INTERNAL_DNTP_CONC: dNTP concentration in mM (>= 0; default 0). (default: 0) |
| `dna_conc_nm` | number | No | PRIMER_INTERNAL_DNA_CONC: probe (oligo) concentration in nM (> 0; default 50). (default: 50) |
| `thermodynamic_screen` | boolean | No | true (default) = Primer3 ntthal self-dimer / 3' self-dimer / hairpin Tm screen for A/C/G/T probes <= thermodynamic_screen_max_length nt (default 60) (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1); false = the fallback screens for every probe (Primer3 alignment-mode self_any / self_end + the sequence-only stem-loop screen of hairpin_potential: an exactly complementary stem >= 4 bp closing a loop >= 3 nt). |
| `max_structure_tm` | number | No | PRIMER_INTERNAL_MAX_SELF_ANY_TH = _SELF_END_TH = _HAIRPIN_TH: maximum ntthal self-dimer / 3' self-dimer / hairpin Tm in °C of the thermodynamic screen (default 47). |
| `max_self_any` | number | No | PRIMER_INTERNAL_MAX_SELF_ANY: maximum Primer3 alignment-mode self_any of the fallback screen (default 12). |
| `max_self_end` | number | No | PRIMER_INTERNAL_MAX_SELF_END: maximum Primer3 alignment-mode self_end of the fallback screen (default 12). |
| `thermodynamic_screen_max_length` | integer | No | THAL_MAX_ALIGN of ntthal (ProbeParameters.ThermodynamicScreenMaxLength): longest A/C/G/T probe given the thermodynamic self-structure screen, and the length limit of the non-target duplex Tm (computed when the probe or the site is at most this long). Default 60 = Primer3; opt-in 61..10000 = thal.c compiled with -DTHAL_MAX_ALIGN (cost grows as n^2). (default: 60) |
| `both_strands` | boolean | No | true = also scan the reference sequences for the probe's reverse complement (the probe's sites on the other strand of a double-stranded reference, as blastn -strand both); a position hit in both orientations (reverse-palindromic probe or site) counts once. Default false = the given strand only. (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `isValid` | boolean | True when no issue was recorded. |
| `specificityScore` | number | Library uniqueness score `[0,1]` (0 hits → 0, N → 1/N); library convention, reported only. |
| `offTargetHits` | integer | Total approximate-match hits across references (intended site included). |
| `crossHybridizingHits` | integer | Hits meeting a Kane criterion (identity > `max_non_target_identity` or > `max_contiguous_match` contiguous identical nt); > 1 is an issue. |
| `selfComplementarity` | number | Fold-back fraction `[0,1]` (informational library metric). |
| `selfAny` / `selfEnd` | number | Primer3 alignment-mode internal-oligo self_any / self_end (dpal; the fallback self-dimer criterion, limit 12.00). |
| `hasSecondaryStructure` | boolean | ntthal hairpin Tm > `max_structure_tm` (an issue), or the fallback stem-loop screen: exact ≥ 4-bp stem, loop ≥ 3 nt (a warning only). |
| `issues` | string[] | Reported issues (sourced criteria; any issue → `isValid` false). |
| `warnings` | string[] | Reported findings that do not decide `isValid`: the fallback stem-loop flag "Potential secondary structure formation" (library convention; A5-5). |
| `thermodynamicScreen` | boolean | True when the Primer3 ntthal screen was applied. |
| `selfDimerTm` / `selfEndDimerTm` / `hairpinTm` | number \| null | ntthal Tm values in °C (null with the fallback screen). |
| `crossHybridization` | object[] | Per non-target strand: `nonTargetIndex`, `reverseComplementStrand`, `identity`, `identicalColumns`, `alignmentScore`, `longestContiguousMatch`, `exceedsIdentityThreshold`, `exceedsContiguousThreshold`, `exceedsDuplexTmThreshold`, `crossHybridizes`, `siteStart`/`siteEnd` (aligned site in the assessed strand), `duplexTm` (ntthal Tm of the probe on the site's complementary strand, primer3-py `calc_heterodimer`; null when probe and site are both longer than `thermodynamic_screen_max_length` or non-ACGT). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Probe sequence cannot be null |
| 1002 | Reference sequences cannot be null |
| 1003 | Maximum mismatches cannot be negative |
| 1004 | Non-target identity threshold must be in [0, 1] |
| 1005 | Contiguous-match threshold cannot be negative |
| 1006 | Illegal value for salt or DNA concentration (Primer3: monovalent and oligo concentration > 0) |
| 1007 | Illegal value for divalent salt or dNTP concentration (Primer3: ≥ 0) |

## Examples

### Example 1: Unique probe → specificity 1.0

Probe `ATCGATCGATCGATCGATCG` in `NNNNNATCGATCGATCGATCGATCGNNNN` with `max_mismatches = 0`: 1 hit → specificity 1.0.

**Input:** `{ "probe_sequence": "ATCGATCGATCGATCGATCG", "reference_sequences": ["NNNNNATCGATCGATCGATCGATCGNNNN"], "max_mismatches": 0 }`

**Response (abridged):** `{ "offTargetHits": 1, "specificityScore": 1.0 }`

### Example 2: Repetitive probe → 1/N specificity

10-mer poly-A in a 34-mer poly-A: 25 exact-match positions → specificity `1/25 = 0.04`.

### Example 3: Kane cross-hybridization

Probe `TATGCCTCCGGTACATCAACTACAGTTAGCCTTAAGAGAAAAATCCCAAA` against a non-target carrying the probe with a substitution at every 5th position: best local alignment 40/50 identical → identity 0.80 > 0.75 → issue `Cross-hybridization risk with non-target 0: identity 80%, …`, `isValid = false`.

### Example 4: Stated reaction conditions

`GCGCGCGCGCGCGCGCGCGC` with `monovalent_mm = 100`, `divalent_mm = 2`, `dntp_mm = 0.2`, `dna_conc_nm = 250`: ntthal self-dimer / 3′ self-dimer / hairpin Tm 88.57933031369095 / 88.57933031369095 / 93.4845818217363 °C (primer3-py 2.3.1 at those conditions; 78.86 / 78.86 / 87.30 °C at the defaults) → `isValid = false`.

### Example 5: Reference sites judged by the Kane criteria

12-mer `GATCCGACGCTA` against `TTTTTGATCCGACGCTATTTTTGCTCCTACGGTATTTTT`: 2 hits (exact; 3 mismatches = 9/12 = 0.75, not > 0.75) → `crossHybridizingHits = 1`, no off-target issue, `specificityScore = 0.5`.

### Example 6: Reverse-complement site with `both_strands`

20-mer `GATCCGACGCTATATGCCGT` against a reference holding the exact site and the reverse complement of the probe with 2 substitutions (positions 3, 11) and with 5 substitutions: given strand → 1 hit; `both_strands = true` → 2 hits (the rc site has 18/20 = 0.90 > 0.75) → `crossHybridizingHits = 2`, issue `2 potential off-target sites …`, `isValid = false`, `specificityScore = 0.5` (independent Python brute force over both strands; with `max_mismatches = 5` → 3 hits, 2 Kane sites — the 5-substitution site is 15/20 = 0.75, longest run 3).

## See Also

- [design_probes](design_probes.md), [design_probes_primer3](design_probes_primer3.md), [design_tiling_probes](design_tiling_probes.md), [crispr_specificity_score](crispr_specificity_score.md)
