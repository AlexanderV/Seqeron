---
type: concept
title: "miRNA target-site prediction (site-type classification + context++ scoring)"
tags: [mirna, algorithm]
sources:
  - docs/Evidence/MIRNA-TARGET-001-Evidence.md
  - docs/algorithms/MiRNA/Target_Site_Prediction.md
source_commit: aa11631f0f0b525bc218f877ce18b6e69d373542
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: mirna-target-001-evidence
      evidence: "Test Unit ID: MIRNA-TARGET-001 ... Algorithm Group: MiRNA ... Target Site Prediction (FindTargetSites)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:seed-sequence-analysis
      source: mirna-target-001-evidence
      evidence: "The seed region is positions 2-8 ... The reverse complement of the seed is sought in the mRNA — the seed extracted by MIRNA-SEED-001 is the INPUT that determines targeting"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:rna-base-pairing
      source: mirna-target-001-evidence
      evidence: "FindTargetSites uses GetReverseComplement (seed → target motif) and AlignMiRnaToTarget (antiparallel duplex, | WC / : G:U wobble); Watson-Crick A:U,G:C + G:U wobble pairing"
      confidence: high
      status: current
---

# miRNA target-site prediction (site-type classification + context++ scoring)

Finding **mRNA target sites** for a mature miRNA and classifying each into the standard animal-miRNA
**site-type hierarchy** — then scoring efficacy. This is the **fourth and final ingested unit of the
MiRNA family** (test unit **MIRNA-TARGET-001**, `MiRnaAnalyzer.FindTargetSites`), **completing** the
family: the record is [[mirna-target-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern. It is the consumer end of the family
— it **depends on** [[seed-sequence-analysis]] (the extracted 2-8 seed **is** the targeting
determinant) and on [[rna-base-pairing]] (`GetReverseComplement` builds the target motif;
`AlignMiRnaToTarget` builds the antiparallel duplex). The sibling precursor unit is
[[pre-mirna-hairpin-detection]].

The core biology (Bartel 2009; Lewis et al. 2005): an animal miRNA recognises targets by antiparallel
pairing between its 5' **seed** and complementary sites (chiefly in 3'UTRs). The **reverse complement
of the seed (positions 2-8) is sought on the mRNA**; miRNA position 1 sits at the **3' end** of the
target site.

## 1. Site-type hierarchy (Bartel 2009 / TargetScan)

Sites are ranked by efficacy (8mer strongest). Let `seedRC = revcomp(m₂…m₈)`; the finder derives the
**6mer core** `seedRC[2..7]` and the **offset-6mer pattern** `seedRC[1..6]`:

| Type | Pairing rule | SeedMatchLength |
|------|--------------|-----------------|
| **8mer** | miRNA pos 2-8 match **+** `A` opposite pos 1 | 8 |
| **7mer-m8** | miRNA pos 2-8 match | 7 |
| **7mer-A1** | miRNA pos 2-7 match **+** `A` opposite pos 1 | 7 |
| **6mer** | miRNA pos 2-7 match | 6 |
| **Offset 6mer** | miRNA pos 3-8 match (marginal efficacy) | 6 |

Worked oracle (**hsa-let-7a-5p**, seed `GAGGUAG` → `seedRC = CUACCUC`): 8mer `…CUACCUCA…`,
7mer-m8 `…CUACCUC…`, 7mer-A1 `…UACCUCA…`, 6mer `…UACCUC…`. Centered sites are removed by TargetScan 8
as not reliably functional; the enum carries `Centered`/`Supplementary` but `FindTargetSites` never
emits them (accepted deviation).

## 2. The finder — `FindTargetSites(mRnaSequence, miRna, minScore = 0.5)`

Two-pass, deterministic, `T`→`U`-normalised, case-insensitive, coordinates **0-based inclusive**:

1. **Pass 1** scans for exact 6mer-core matches and **upgrades** each hit: upstream position-8 base
   present **and** downstream A1 present → 8mer; only position-8 → 7mer-m8; only A1 → 7mer-A1; else 6mer.
2. Each pass-1 hit gets a full antiparallel duplex (`AlignMiRnaToTarget`) and a heuristic score; kept if
   `Score ≥ minScore`, and its coordinates are marked covered.
3. **Pass 2** scans for offset-6mer matches, **skipping** any that overlap a retained pass-1 site or that
   would form a full seed match at the same location — so **higher-priority canonical classes suppress
   overlapping offset-6mer calls** (INV-04).

**Default score (heuristic, `[0,1]`):** base `1.0` (8mer) / `0.52` (7mer-m8) / `0.32` (7mer-A1) /
`0.15` (6mer) / `0.10` (offset 6mer), `+0.05` when the duplex has >10 matches, `−0.01` per mismatch,
then clamped. This **ranks hits inside the repository but is not a calibrated repression probability**;
it is monotone `8mer ≥ 7mer-m8 ≥ 7mer-A1 ≥ 6mer ≥ offset 6mer`. `FreeEnergy` is likewise a heuristic
duplex energy (stacked WC `−2.0`, isolated WC `−1.0`, G:U `−0.5`, mismatch `+0.5`), **not** a
nearest-neighbor fold — negative for a well-paired duplex, sign-reliable only.

## 3. Opt-in TargetScan context++ score (Agarwal et al. 2015)

`ScoreTargetSiteContextPlusPlus` provides the **source-fitted** alternative to the heuristic score,
leaving the default `Score` untouched. context++ is a multiple-linear regression with a **distinct
coefficient set per site type**; the score is the **sum of per-feature contributions**
`contribution = coeff × scaledScore` (`getAgarwalContribution`), with the **verbatim** coefficients
from `Agarwal_2015_parameters.txt` and the scaling/DP/indicator logic ported from
`targetscan_70_context_scores.pl`:

- **Continuous features** (`TA_3UTR, SPS, Local_AU, 3P_score, SA, Len_ORF, Len_3UTR, Min_dist, PCT`) are
  **min-max scaled** `(raw − min)/(max − min)` per site type (not clamped to `[0,1]`).
- **Binary indicators** (`sRNA1{A,C,G}`, `sRNA8{A,C,G}`, `Site8{A,C,G}`, `Off6m`, `ORF8m`) are used **raw**.
- **Intercept** (the site-type term) is added once.

**Realised from the miRNA + the supplied 3'UTR alone:** Intercept, **Local_AU** (position-weighted A/U
fraction over the 30 nt up/downstream), **sRNA1/sRNA8** position-identity indicators (skipped when that
nt is U), **Site8** target-position-8 indicators (defined only for 7mer-A1 / 6mer), **3P_score** (3'
supplementary pairing — max over single-gap offsets of summed runs of ≥2 pairs, min-max scaled with
min=1/max=3.5), **Min_dist** (`log10` distance to nearest 3'UTR end), **Len_3UTR** (`log10` UTR length),
and **Off6m** (offset-6mer count). Two more are **computed** rather than caller-supplied:

- **TA_3UTR** — target-site abundance from a supplied 3'UTR set: `TA = log10(N)`, `N` = total
  non-overlapping **8mer + 7mer-m8 + 7mer-A1** sites of the seed across the UTRs (`ComputeTa3Utr` /
  `CountSeedSites3Utr`, per Garcia 2011; TargetScan stores this as log10 in `TA_SPS_by_seed_region.txt`).
  Hand-checked oracle: five UTRs → N = 5 → TA = 0.698970004336019.
- **SA** — 14-nt-window unpaired probability from the **Turner-2004 McCaskill partition function**
  (`RnaSecondaryStructure.CalculateRegionUnpairedProbability` = `Z_open / Z`, local ±40-nt / 80-nt
  window, then `log10` + min-max scaled) — the accessibility feature, computed exactly per
  `getSA_contribution` (Agarwal Fig 4A: log10 unpaired probability of a 14-nt region centred on the
  match to miRNA nucleotides 7-8).

**Still caller-supplied (data/parameter-blocked):** **SPS** (Garcia 2011 seed-region table lookup),
**Len_ORF** / **ORF8m** (need the transcript ORF), and the **PCT sigmoid parameters** `b0..b3` (per
miRNA family, in TargetScan's compiled citation-required tables — not published as numbers, so not
bundled). **PCT itself** is computed when the caller supplies a phylogenetic tree + conserved-species
set: the library derives the **Friedman 2009 branch-length score (Bls)** = total branch length of the
minimal subtree connecting the conserved species, then maps it via the published logistic
`PCT(Bls) = b0 + b1/(1 + e^(−b2·Bls + b3))` truncated at 0 (`targetscan_70_BL_PCT.pl`). Worked tree
`((A:1.0,B:2.0):0.5,(C:1.5,D:3.0):4.0);` → Bls {A,B}=3.0, {A,C}=7.0, {A,B,C,D}=12.0, {A}=0.0.

Whatever remains residual is reported in `OmittedFeatures`, so `ContextScorePartial` is honestly a
**partial** context++ score, not the published headline CS. Worked 8mer partial CS (let-7a,
`GGGGG+CUACCUCA+GGGGG`, all-G flanks ⇒ Local_AU fraction 0): Intercept −0.589 + Local_AU +0.15461 +
sRNA8G +0.015 + 3P_score +0.016 + Min_dist −0.04976 + Len_3UTR −0.28304 + Off6m −0.020 =
**−0.7561913315126536** (verbatim Agarwal coefficients).

## Invariants and edge cases

- **INV:** every reported site has `0 ≤ Score ≤ 1`; `SeedMatchLength ∈ {8,7,6}`; `Start`/`End` 0-based
  inclusive; duplex is antiparallel; pass-1 canonical sites suppress overlapping offset-6mer.
- Empty/null mRNA, empty `miRna.Sequence`, or seed RC shorter than 7 nt → **no sites** (non-throwing).
- DNA target → `T` treated as `U`. `minScore > 1.0` → no sites; `minScore < 0.0` → every scored
  candidate emitted. Multiple/overlapping seed-RC occurrences are each reported.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for **canonical seed-complement site
discovery + site-type classification**, plus a **fully source-traced but partial** context++ efficacy
score. It models the seed-driven canonical classes only: **bulged, centered, supplementary-only, and
noncanonical** sites are out of scope, as are transcript isoforms and (in the default finder) AU
context / conservation / accessibility. The default `Score` and `FreeEnergy` are ranking aids, not
biophysical or transcriptome-scale measures; the context++ headline CS requires the data-/parameter-
blocked features (use the official TargetScan distribution, or supply them). **No source
contradictions** — Bartel 2009, Lewis 2005, Grimson 2007, Agarwal 2015, Garcia 2011, Friedman 2009,
TargetScan, and miRBase agree on the site ladder and the context++ model; the recorded items are the
intentional heuristic-score / partial-CS simplifications and the unemitted `Centered`/`Supplementary`
enum members.
