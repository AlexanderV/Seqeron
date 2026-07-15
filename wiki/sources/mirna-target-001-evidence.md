---
type: source
title: "Evidence: MIRNA-TARGET-001 (miRNA target-site prediction — site-type classification + context++ scoring)"
tags: [validation, mirna]
doc_path: docs/Evidence/MIRNA-TARGET-001-Evidence.md
sources:
  - docs/Evidence/MIRNA-TARGET-001-Evidence.md
source_commit: aa11631f0f0b525bc218f877ce18b6e69d373542
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: MIRNA-TARGET-001

The validation-evidence artifact for test unit **MIRNA-TARGET-001** — **Target Site Prediction**,
canonical seed-based target-site discovery, site-type classification, and the opt-in TargetScan
context++ score (`MiRnaAnalyzer.FindTargetSites` / `ScoreTargetSiteContextPlusPlus` and helpers).
This is the **fourth and final ingested unit of the MiRNA family**, completing it, and one instance of
the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the
synthesizing concept is [[mirna-target-site-prediction]]. [[test-unit-registry]] tracks the unit.

The unit validates the **antiparallel seed-complement scan** (the reverse complement of the miRNA seed
positions 2-8 is sought on the mRNA), classification into the **Bartel/TargetScan hierarchy**
(8mer > 7mer-m8 > 7mer-A1 > 6mer, plus offset 6mer), and the fully source-traced **context++**
regression scorer.

## What this file records

- **Online sources (primary literature + TargetScan distribution), retrieved verbatim this session:**
  - **Bartel (2009)** *Cell* (PMID 19167326) — seed-site hierarchy 8mer > 7mer-m8 > 7mer-A1 > 6mer;
    seed positions 2-8; antiparallel binding.
  - **Lewis et al. (2005)** *Cell* (PMID 15652477) — "conserved seed pairing, often flanked by
    adenosines"; an adenosine opposite position 1 correlates with efficacy.
  - **Grimson et al. (2007)** *Mol Cell* (PMID 17612493) — site-type efficacy weights; AU-rich context
    favors targeting.
  - **Agarwal et al. (2015)** *eLife* 4:e05005 (PMID 26267216) — the **context++** model: one multiple
    linear regression **per site type**; verbatim fitted coefficients from `Agarwal_2015_parameters.txt`
    and feature/scaling logic from `targetscan_70_context_scores.pl` (both `curl`'d from the nsoranzo
    TargetScan mirror).
  - **TargetScan 8.0 FAQ** — centered sites removed as not reliably functional; 3 canonical + offset 6mer.
  - **Friedman et al. (2009)** (PMID 18955434) — PCT (probability of conserved targeting); branch-length
    score (Bls); most mammalian mRNAs are miRNA targets. PCT sigmoid via `targetscan_70_BL_PCT.pl`.
  - **Garcia et al. (2011)** *Nat Struct Mol Biol* 18:1139 — TA_3UTR = number of non-overlapping 3'UTR
    8mer/7mer-m8/7mer-A1 sites; stored as **log10(count)** in `TA_SPS_by_seed_region.txt`.
  - **McCaskill (1990)** + **Lorenz et al. (2011)** ViennaRNA + RNAplfold man page — partition-function
    unpaired-probability basis for the SA (structural accessibility) feature.
  - **miRBase** — authoritative miRNA sequences for the fixtures (hsa-let-7a-5p, hsa-miR-122-5p, hsa-miR-21-5p).

- **Key definitions captured:** the five site types (8mer = pos 2-8 + A opposite pos 1; 7mer-m8 = 2-8;
  7mer-A1 = 2-7 + A opposite pos 1; 6mer = 2-7; offset 6mer = 3-8); antiparallel orientation (miRNA
  3'→5' vs mRNA 5'→3', pos 1 at the 3' end of the target site); the alignment glyphs `|` WC / `:` G:U
  wobble / space mismatch.

- **Reference oracles:** hsa-let-7a-5p seed `GAGGUAG` → seed RC `CUACCUC`; 8mer site `...CUACCUCA...`,
  7mer-m8 `...CUACCUC...`, 7mer-A1 `...UACCUCA...`, 6mer `...UACCUC...`. Score monotonicity
  `8mer ≥ 7mer-m8 ≥ 7mer-A1 ≥ 6mer ≥ offset 6mer`. Hand-derived context++ contributions (all-G-flank
  8mer/7mer-m8/7mer-A1 partial CS, the miR-21 6mer case, TA = log10(5) = 0.698970004336019, the SA
  partition-function `GAAAC`/`Z = 1 + e^(−5.4/RT)` pin, the PCT branch-length worked tree
  `((A:1.0,B:2.0):0.5,(C:1.5,D:3.0):4.0);` → Bls {A,B}=3.0 / {A,C}=7.0 / {A,B,C,D}=12.0 / {A}=0.0).

- **Documented edge cases:** empty/null → no sites; mRNA shorter than seed → no sites; multiple and
  overlapping sites each reported; DNA `T`→`U` normalization; case-insensitive; `minScore` filter
  (`>1.0` suppresses all, `<0.0` admits all scored); at one position 8mer outranks 7mer-m8; all scores
  in `[0,1]`; well-paired duplex FreeEnergy negative; SeedMatchLength 8/7/6.

## Deviations and assumptions

**Intentionally simplified** (documented in the algorithm doc, not open correctness gaps):

- The default `Score` is a **heuristic** normalization of site class + mismatch burden + extra pairing
  (base 1.0 / 0.52 / 0.32 / 0.15 / 0.10; +0.05 for >10 matches; −0.01 per mismatch; clamp `[0,1]`), not a
  calibrated repression probability. `FreeEnergy` is a heuristic duplex energy, not a nearest-neighbor
  fold; the perl-style running total omits the final alignment position.
- `TargetSiteType.Centered` and `Supplementary` exist in the enum but `FindTargetSites` never emits them
  (accepted deviation) — discovery is limited to canonical seed classes + offset 6mer.
- The opt-in **context++** score is a **partial** CS: `SPS`, `Len_ORF`, `ORF8m`, and the PCT per-family
  sigmoid parameters `b0..b3` are **caller-supplied** (data-/parameter-blocked; live in TargetScan's
  compiled, citation-required tables — not published as numbers, so not bundled/invented). `TA_3UTR`
  **is** now computable from a supplied 3'UTR set (`ComputeTa3Utr` = log10 N). Whatever stays residual is
  listed in `OmittedFeatures`. Only the published equations + the permissive Agarwal coefficients are bundled.

**No source contradictions** — Bartel 2009, Lewis 2005, Grimson 2007, Agarwal 2015, Garcia 2011,
Friedman 2009, TargetScan, and miRBase are mutually consistent on the site ladder and the context++ model.
