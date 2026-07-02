# seqeron-mirna — tool map (all 14)

Server: **annotation**. One backing class: `MiRnaAnalyzer.*`.
This skill is **not** in `domain-map.json`, so it has **no** generated `_generated/tools.md` —
**this curated map is the index.** Verify schemas in `docs/mcp/tools/annotation/<tool>.md`.

> Coordinates: site/target positions are **0-based inclusive** into the mRNA. The seed is
> **positions 2–8** (0-based `Substring(1,7)`, so `seedStart=1`, `seedEnd=7`). Free energies are
> **kcal/mol at 37 °C**, Turner 2004 nearest-neighbour. Inputs are case-insensitive; DNA `T` is
> treated as RNA `U`. Always confirm exact I/O in the tool doc.

## miRNA record / seed

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `create_mirna` | `MiRnaAnalyzer.CreateMiRna` | Build a `MiRna` record: upper-case + T→U normalise, extract seed 2–8 (`seedStart=1`, `seedEnd=7`). | `create_mirna.md` |
| `mirna_seed_sequence` | `MiRnaAnalyzer.GetSeedSequence` | Seed only (positions 2–8, upper-cased). Requires ≥ 8 nt; `""` if shorter. | `mirna_seed_sequence.md` |
| `rna_reverse_complement` | `MiRnaAnalyzer.GetReverseComplement` | RNA reverse complement (A↔U, G↔C; T→A; unknown→N). Derives a seed's target pattern (`GAGGUAG`→`CUACCUC`). | `rna_reverse_complement.md` |

## Target-site prediction (⚠ MIRNA-TARGET-001)

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_mirna_target_sites` | `MiRnaAnalyzer.FindTargetSites` | Scan mRNA (T→U internal) for canonical sites keyed on seed RC; classify Bartel/TargetScan hierarchy (8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer), score ≥ `minScore` (default 0.5). ⚠ `score` = site-type + pairing, **not** context++. | `find_mirna_target_sites.md` |
| `analyze_target_context` | `MiRnaAnalyzer.AnalyzeTargetContext` | AU content of a flanking window + near-start/near-end flags → `contextScore∈[0,1]` (`0.5·auContent` + `0.3` if mid-transcript). ⚠ **partial** context proxy, not context++. | `analyze_target_context.md` |
| `site_accessibility` | `MiRnaAnalyzer.CalculateSiteAccessibility` | Local secondary-structure density in a ±50 nt window → `accessibility=max(0,1−density·10)`. Higher = more open/targetable. (Not itself guarded.) | `site_accessibility.md` |
| `align_mirna_to_target` | `MiRnaAnalyzer.AlignMiRnaToTarget` | Antiparallel ungapped duplex: miRNA `i` vs `target[len−1−i]`; per-position `|` WC / `:` G-U / ` ` mismatch; `matches`, `guWobbles`, Turner `freeEnergy`. (Not guarded.) | `align_mirna_to_target.md` |

## Seed families / similarity / variants (not guarded)

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `group_by_seed_family` | `MiRnaAnalyzer.GroupBySeedFamily` | Partition miRNAs into families by **exact** seed equality → `families[{seedFamily,members}]`. | `group_by_seed_family.md` |
| `compare_seed_regions` | `MiRnaAnalyzer.CompareSeedRegions` | Position-wise seed compare → `matches`, `mismatches` (Hamming + length diff), `isSameFamily` (seeds identical). | `compare_seed_regions.md` |
| `find_similar_mirnas` | `MiRnaAnalyzer.FindSimilarMiRnas` | Database miRNAs within `maxMismatches` (default 1) seed Hamming distance; query excluded by name. | `find_similar_mirnas.md` |
| `generate_seed_variants` | `MiRnaAnalyzer.GenerateSeedVariants` | Original seed + every single-position A/C/G/U substitution → `1 + 3L` sequences. `includeWobble` reserved (currently no effect). | `generate_seed_variants.md` |

## Pre-miRNA hairpins (⚠ MIRNA-CLEAVAGE-001, MIRNA-PRECURSOR-001)

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_pre_mirna_hairpins` | `MiRnaAnalyzer.FindPreMiRnaHairpins` | Scan for stem-loop windows (loop 3–25 nt, total in `[minHairpinLength=55, maxHairpinLength=120]`, `matureLength=22`) → span, mature/star arm, dot-bracket, Turner energy. Simplified consecutive-pairing model (bulged real pre-miRNAs may be missed). ⚠ star-arm span approximate; **no** read-support score. | `find_pre_mirna_hairpins.md` |

## Base pairing (miRNA copies — distinct from `RnaSecondaryStructure.*`)

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `can_pair` | `MiRnaAnalyzer.CanPair` | `true` for A-U, G-C, or G-U wobble (case-insensitive; T→U). | `can_pair.md` |
| `is_wobble_pair` | `MiRnaAnalyzer.IsWobblePair` | `true` only for G-U / U-G (Watson-Crick → `false`). | `is_wobble_pair.md` |

## Envelope — three guarded units

Authoritative: [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md).

- **MIRNA-TARGET-001** (MinimumMode `Permissive`, guarded) — full context++ score not out-of-the-box;
  `find_mirna_target_sites` score = site-type + pairing, `analyze_target_context` = partial AU/position
  proxy. Supply 3'UTR set (`ComputeTa3Utr`→`TA_3UTR`) + alignment/tree/sigmoid for `PCT` + `SPS`/ORF
  features, or cite TargetScan. See STOP rule in [`../SKILL.md`](../SKILL.md).
- **MIRNA-CLEAVAGE-001** (MinimumMode `Permissive`, guarded) — 5p cut reproduces miRBase exactly, but
  the miRNA\*(3p)/star-arm span is an approximate 2-nt-3′-overhang offset. Supply miRBase MIMAT
  coordinates or small-RNA-seq read pileups for exact boundaries.
- **MIRNA-PRECURSOR-001** (documented-only, no runtime guard) — miRDeep2 read-stacking / read-support
  score not implemented; `find_pre_mirna_hairpins` yields a structural candidate only. Use miRDeep2
  with the caller's reads for read support.

## Not a tool here (route elsewhere)

- **General RNA folding / MFE / stem-loops / pseudoknots / dot-bracket** → [`seqeron-rna-structure`](../../seqeron-rna-structure/SKILL.md)
  (`RnaSecondaryStructure.*`). The `can_pair`/`is_wobble_pair` here are the **miRNA** copies.
- **mRNA 3'UTR / ORFs / gene structure / splicing** → [`bio-annotation`](../../bio-annotation/SKILL.md).
- **Full context++ / TargetScan tables, small-RNA-seq read support** → out of scope (supply the data
  or use the named external tool — see Envelope).
