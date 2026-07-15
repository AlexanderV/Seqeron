---
type: source
title: "Evidence: META-BETA-001 (beta diversity — Bray-Curtis dissimilarity + Jaccard distance)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-BETA-001-Evidence.md
sources:
  - docs/Evidence/META-BETA-001-Evidence.md
source_commit: 1eb99ff7a5a3f4819f81562f353676b522f5769f
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-BETA-001

The validation-evidence artifact for test unit **META-BETA-001** — **beta diversity**, the
*between-sample* dissimilarity computed by `MetagenomicsAnalyzer.CalculateBetaDiversity`. Second
ingested unit of the Metagenomics family and the direct sibling of the within-sample
[[meta-alpha-001-evidence|META-ALPHA-001]]. One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The method is synthesized in its own
concept, [[beta-diversity]]; [[test-unit-registry]] tracks the unit.

## What this file records

Two classical dissimilarity measures validated together (both deterministic functions of two
taxon→abundance sample profiles), plus a placeholder `UniFracDistance` field:

- **Online sources (mutually consistent, no contradictions):**
  - **Whittaker (1960)** *Vegetation of the Siskiyou Mountains* — original beta-diversity concept.
  - **Bray & Curtis (1957)** *An Ordination of the Upland Forest Communities of Southern Wisconsin* —
    the Bray-Curtis dissimilarity formula and its properties.
  - **Jaccard (1901)** *Étude comparative de la distribution florale…* — the Jaccard index.
  - Wikipedia (Beta diversity / Bray–Curtis dissimilarity / Jaccard index), retrieved 2026-03-09.

- **Extracted formulas & properties:**
  - **Bray-Curtis:** `BC_jk = 1 − 2·C_jk / (S_j + S_k)` = `1 − 2·Σ min(N_ij,N_ik) / Σ(N_ij + N_ik)`;
    range `[0,1]` (0 = identical, 1 = no shared species); **uses abundance**; **NOT a true metric**
    (violates the triangle inequality).
  - **Jaccard distance:** `J = 1 − |A∩B| / |A∪B|`; binary form `d_J = (M01+M10)/(M01+M10+M11)`;
    range `[0,1]`; **presence/absence only**; **is a true distance metric** (satisfies the triangle
    inequality); defined for finite non-empty sample sets.

- **Source-verified invariants:** range `[0,1]` for both; symmetry `D(A,B)=D(B,A)`; identity
  `D(A,A)=0`; disjoint samples → distance `1` for both; **abundance sensitivity** — Bray-Curtis uses
  abundance values, Jaccard uses presence/absence only.

- **Worked oracles (Wikipedia Bray-Curtis fish-tank example):** Tank 1 `{Goldfish:6, Guppy:7,
  Rainbow:4}` (S=17) vs Tank 2 `{Goldfish:10, Guppy:0, Rainbow:6}` (S=16); shared min-sum
  `C = 6+0+4 = 10` → **BC = 1 − 20/33 = 13/33 ≈ 0.3939**. Jaccard on the same data: shared
  `{Goldfish, Rainbow} = 2`, union `= 3` (Guppy absent in Tank 2, count 0) → **J = 1 − 2/3 = 1/3 ≈
  0.3333**. The 2026-03-09 coverage pass also added the exact SymmetryProperty oracle **BC = 3/5,
  Jaccard = 2/3**.

## Implementation notes (from the Evidence file)

`CalculateBetaDiversity` compares the union of taxon keys from both samples in one pass; a taxon is
**present** iff its abundance is strictly `> 0`, and missing keys default to `0`. Bray-Curtis applies
the formula **directly to the supplied abundance values (no internal normalization)** — the caller
decides whether they are raw counts or relative abundances. Jaccard first converts abundances to
presence/absence via the `> 0` rule. The result record also carries `SharedSpecies`, `UniqueToSample1`,
`UniqueToSample2` counts and preserves the two sample-name labels verbatim.

## Design decisions and deviations

- **Empty-sample convention:** both formulas divide by a total; when both samples are empty the
  denominator is `0` and the value is undefined. The implementation returns **`0`** for two identical
  empty samples (zero distance for identical inputs) and — per the algorithm doc's edge-case table —
  **`1`** when one sample is empty and the other has present taxa (maximum distance), consistent with
  the `scipy.spatial.distance` convention. (The Evidence file frames it as "returns 0 for two identical
  empty samples and 1.0 when one sample is empty"; the algorithm doc's guard states all-zero/empty →
  BC 0 / Jaccard 0.)
- **UniFrac placeholder:** `UniFracDistance` is hard-coded to `0` — the phylogenetic UniFrac measure
  requires a tree that is outside this method's contract. The single accepted **deviation**
  (algorithm doc §5.4): the field exists on the `BetaDiversity` record for future extension only.

## Coverage classification (2026-03-09)

11 scenarios covered unchanged; **2 weak → strengthened** (BothSamplesEmpty gained explicit BC=0 /
Jaccard=0 assertions; SymmetryProperty gained exact hand-computed BC=3/5 & Jaccard=2/3 with species
counts and name preservation); **1 duplicate → removed** (DistanceRangeConstraints' four scenarios
subsumed by IdenticalSamples / DifferentSingleSpecies / SymmetryProperty / OneSampleEmpty; range
`[0,1]` implicitly verified by 12 exact-value tests); **0 missing**.

No source contradictions — the Bray & Curtis 1957 and Jaccard 1901 primary definitions and the
Wikipedia formulations are the standard, mutually consistent measures.
