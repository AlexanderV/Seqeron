---
type: concept
title: "Intrinsic-disorder prediction (TOP-IDP sliding window, PredictDisorder)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/DISORDER-PRED-001-Evidence.md
  - docs/Evidence/DISORDER-PROPENSITY-001-Evidence.md
  - docs/Evidence/DISORDER-REGION-001-Evidence.md
  - docs/algorithms/ProteinPred/Disorder_Prediction.md
source_commit: 98b44f1a8112227eb70a11c589272ca8ce62e7af
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: disorder-pred-001-evidence
      evidence: "Test Unit ID: DISORDER-PRED-001 ... Algorithm: Disorder Prediction (Intrinsically Disordered Protein Prediction)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:protein-low-complexity-seg
      source: disorder-pred-001-evidence
      evidence: "Both are protein disorder / features family units operating on compositional signals; SEG measures Shannon-entropy low complexity while PredictDisorder scores TOP-IDP disorder propensity — distinct but sibling algorithms flagged as separate concepts on the SEG page"
      confidence: high
      status: current
---

# Intrinsic-disorder prediction (TOP-IDP sliding window, PredictDisorder)

Estimating **intrinsic disorder** — regions of a protein that lack a single stable 3-D structure
under physiological conditions and instead populate dynamic conformational ensembles (Dunker et al.
2001) — from **amino-acid composition alone**. Seqeron's `DisorderPredictor.PredictDisorder` scores
each residue as the **sliding-window average of the normalized TOP-IDP propensity scale** (Campen et
al. 2008) and labels residues at or above a cutoff as disordered. Validated under test unit
**DISORDER-PRED-001**; the validation record is [[disorder-pred-001-evidence]] and
[[test-unit-registry]] tracks the unit. The **raw per-residue propensity primitives** underneath the
predictor (`GetDisorderPropensity` + the Dunker classification) are separately validated as
**DISORDER-PROPENSITY-001** ([[disorder-propensity-001-evidence]]) — see the primitives section
below. See [[algorithm-validation-evidence]] for the artifact pattern.

This is the **third ingested unit of the protein disorder / features family** (DISORDER-LC / MORF /
PRED / PROPENSITY / REGION) and the **shared `PredictDisorder` anchor** the family was expected to
grow: [[morf-prediction-dip-in-disorder|MoRF prediction]] reads *this* per-residue disorder profile
to find its ordered "dip within disorder", and **disordered-region detection** aggregates *this*
profile into contiguous runs — validated as **DISORDER-REGION-001**
([[disorder-region-001-evidence]]) and written up in the region-detection section below. It is a
**distinct algorithm** from [[protein-low-complexity-seg|SEG
low-complexity]] (which measures *compositional Shannon entropy*, not disorder propensity) — low-
complexity regions overlap with, but are not identical to, intrinsically disordered regions.

## The TOP-IDP scale and per-residue score

The **TOP-IDP scale** (Campen et al. 2008, PMC2676888) is an amino-acid scale *optimized* to
discriminate order from disorder — derived by surveying 517 scales and applying simulated annealing
(ARV 0.761, an 11% improvement over the best prior scale). For residue `i`, the disorder score is the
mean of the min-max-normalized propensity over a local window `Wᵢ`:

```
Sᵢ = (1/|Wᵢ|) · Σ_{c ∈ Wᵢ} (p(c) − p_min) / (p_max − p_min)
```

with `p_min = −0.884` (W, most order-promoting) and `p_max = 0.987` (P, most disorder-promoting), so
the normalization divisor is `p_max − p_min = 1.871` and every score lands in **`[0, 1]`** (INV-01).
A residue is **disordered** when `Sᵢ ≥` the cutoff.

TOP-IDP Table 2 values (order → disorder) and their normalized equivalents:

| Residue | TOP-IDP raw | Normalized | Dunker class |
|---------|-------------|------------|--------------|
| W | −0.884 | 0.000 | order-promoting |
| F | −0.697 | 0.100 | order-promoting |
| Y | −0.510 | 0.200 | order-promoting |
| I | −0.486 | 0.213 | order-promoting |
| L | −0.326 | 0.298 | order-promoting |
| V | −0.121 | 0.408 | order-promoting |
| N | 0.007 | 0.476 | order-promoting |
| C | 0.020 | 0.483 | order-promoting |
| A | 0.060 | 0.505 | disorder-promoting |
| G | 0.166 | 0.561 | disorder-promoting |
| R | 0.180 | 0.569 | disorder-promoting |
| Q | 0.318 | 0.642 | disorder-promoting |
| S | 0.341 | 0.655 | disorder-promoting |
| K | 0.586 | 0.786 | disorder-promoting |
| E | 0.736 | 0.866 | disorder-promoting |
| P | 0.987 | 1.000 | disorder-promoting |

(M and T are the ambiguous residues; see the classification below.)

## Parameters and the 0.542 cutoff

| Parameter | Default | Role / source |
|-----------|---------|---------------|
| `windowSize` | **21** | width of the TOP-IDP averaging window (Campen's evaluation window); edge windows are **clipped** to the sequence bounds, not padded, so terminal residues score from < 21 positions |
| `disorderThreshold` | **0.542** | disorder cutoff — the TOP-IDP **maximum-likelihood** decision threshold (Campen et al. 2008); `IsDisordered = Sᵢ ≥ 0.542` |
| `minRegionLength` | **5** | minimum contiguous disordered run emitted into `DisorderedRegions` (region detail: `Disordered_Region_Detection.md`) |

At single-residue resolution the 0.542 cutoff falls between **A** (normalized 0.505 → ordered) and
**G** (0.561 → disordered). Note the sibling [[morf-prediction-dip-in-disorder|MoRF]] unit uses the
**0.5** order/disorder threshold (PMC2570644) against the same normalized score — different published
thresholds for different purposes, **not** a contradiction.

## Dunker (2001) residue classification

Three disjoint sets covering all 20 residues, exposed as public properties and confirmed by Campen:

- **Disorder-promoting** — {A, R, G, Q, S, P, E, K} (polar/charged + Ala).
- **Order-promoting** — {W, C, F, I, Y, V, L, N} (bulky hydrophobics + Cys, Asn).
- **Ambiguous/borderline** — {D, H, M, T} (explicitly **NOT** disorder-promoting).

Pro promotes disorder via its rigid cyclic side chain (α-helix breaker); Gly via its side-chain-less
flexibility. `CalculateHydropathy` is a bundled Kyte & Doolittle (1982) utility returning a window's
mean hydropathy, supporting the Uversky et al. (2000) charge–hydropathy view of disorder.

## Per-residue propensity primitives (DISORDER-PROPENSITY-001)

Beneath the windowed predictor sit four **per-residue lookup primitives**, validated separately as
test unit **DISORDER-PROPENSITY-001** ([[disorder-propensity-001-evidence]]):

- `GetDisorderPropensity(residue)` — returns the **raw, un-normalized TOP-IDP Table 2 value**
  (`W = −0.884` … `P = +0.987`), **not** the `[0, 1]` normalized `p(c)` the windowed `Sᵢ` uses. Case
  is folded (input upper-cased first), and any residue outside the 20-residue scale (B, J, O, U, X, Z,
  gaps) returns **`0.0`** — a `GetValueOrDefault(..., 0)` implementation contract, *not* a
  source-defined value.
- `IsDisorderPromoting(residue)` — `true` iff the residue is in the Dunker disorder-promoting set;
  `false` for both order-promoting and the ambiguous `{D, H, M, T}`.
- `DisorderPromotingAminoAcids` = `{A, E, G, K, P, Q, R, S}` and `OrderPromotingAminoAcids` =
  `{C, F, I, L, N, V, W, Y}` — the two public sets (8 members each); the two plus ambiguous
  `{D, H, M, T}` are pairwise disjoint and cover all 20 residues (8 + 8 + 4).

**Do not conflate the two value spaces:** `GetDisorderPropensity` exposes the *raw* scale value,
whereas `PredictDisorder`'s `Sᵢ` averages the *min-max-normalized* value `(p − (−0.884))/1.871`. Anchor
residues coincide (W = raw min / normalized 0, P = raw max / normalized 1), but interior values differ.
The DISORDER-PROPENSITY-001 evidence flags a **ranking-vs-value discrepancy**: the rendered ranking
string places `…Q, K, S, E, P`, yet the Table 2 values give `S = 0.341 < K = 0.586` (so by value
`Q, S, K, E, P`) — the **numeric values are authoritative**, the ranking string is a presentation-order
artifact with no correctness impact.

## Result object

`PredictDisorder` returns a `DisorderPredictionResult`: the **uppercased** `Sequence`;
`ResiduePredictions` (one per residue — 0-based `Position`, `Residue`, `DisorderScore`,
`IsDisordered`); `DisorderedRegions` (contiguous runs ≥ `minRegionLength`); `OverallDisorderContent`
(fraction of residues flagged disordered); and `MeanDisorderScore` (mean of all residue scores).

## Disordered-region detection (DISORDER-REGION-001)

The **aggregation layer** over the per-residue profile: it collapses a scored `ResiduePredictions`
list into `DisorderedRegions` and labels each region. Validated as test unit **DISORDER-REGION-001**
([[disorder-region-001-evidence]]) — the **fifth unit of the protein-disorder family**, and exactly
the growth this anchor anticipated. It introduces no new per-residue math; it re-uses *this* TOP-IDP
profile (cutoff 0.542, window 21).

**Contiguous-run aggregation.** A region is a maximal run of residues with `IsDisordered == true` of
length ≥ `minRegionLength` (default **5**). Each region carries `Start`/`End`, and a `MeanScore` =
mean of its constituent residue `DisorderScore`s. Boundary contract (the oracle set): empty
predictions → **no regions**; all-ordered → **no regions**; all-disordered with length ≥ minLength →
**exactly one region spanning the sequence**; an isolated disordered residue shorter than minLength →
**no region**; a region reaching the sequence end must be captured (**no off-by-one trailing bug**).
Because the underlying window blurs order↔disorder transitions by ±½-window, region boundaries in
mixed sequences may shift by up to `window/2`.

**Region classification — two schemes.** Every region is labelled by a **default first-principles
`RegionType`** and may optionally be re-labelled by the **sourced MobiDB-lite flavor scheme**:

| Scheme | Trigger | Labels | Threshold |
|--------|---------|--------|-----------|
| **Default `RegionType`** (always on) | dominant residue fraction over the region | Proline-rich (P) · Acidic (D/E) · Basic (K/R) · Ser/Thr-rich (S/T) · else **Long IDR** (len > 30, van der Lee 2014) · else **Standard IDR** | fraction **> 0.25** (internal ~5×-random heuristic — **NOT** Das & Pappu NCPR; see contradiction) |
| **`ClassifyRegionFlavorMobiDbLite`** (opt-in) | composition-only, Necci et al. 2020 v3 source | charge: Polyampholyte / Positive- / Negative-Polyelectrolyte / WeaklyCharged; then composition: CysteineRich → ProlineRich → GlycineRich → SEG-low-complexity → Polar `{S,T,N,Q}` | charge FCR/NCPR at **0.35**; enrichment at **≥ 0.32** (inclusive); 9-residue window, sub-region ≥ 9 |

Each region also carries a rescaled **`Confidence` ∈ [0,1]** (INV) — a declared first-principles
heuristic; **MobiDB-lite defines no per-residue confidence**, so `Confidence` is unchanged when the
opt-in flavor scheme is used, and so are the region boundaries. van der Lee et al. 2014 supplies the
compositional IDR subtypes (proline-rich/acidic/basic/Ser-Thr-rich) and the short-vs-long (≥30)
functional split the default labels reproduce; Dunker 2001 supplies the >30-residue long-IDR
significance.

**Classification oracles.** Homopolymers pin the labels: 30×P → 1 proline-rich region; 30×E → acidic;
K/R 40-mer → basic; 31×S → Ser/Thr-rich (note **T alone** normalizes to ≈0.504, *below* the 0.542
cutoff, so it is not disordered on its own); 30×W → 0 regions. MobiDB-lite hand-traced anchors:
`RKDERKDE`→Polyampholyte, `RKRKPPPPPP`→PositivePolyelectrolyte (charge beats P), `RKRKRKR`+13×A →
FCR = 0.35 **not** > 0.35 → WeaklyCharged, `CCCCPPPPAA`→CysteineRich (C first in priority), 8×C+17×A →
C = 0.32 (inclusive) → CysteineRich, 7×C+18×A → 0.28 → WeaklyCharged.

**Contradiction flagged:** the artifact warns that **Das & Pappu 2013's 0.25 boundary is NCPR** (a
globule/coil conformational-state threshold), **not** a compositional-enrichment threshold — so the
default `RegionType` 0.25 cutoff is an internal heuristic, *not* citable to Das & Pappu. MobiDB-lite's
own 0.35 / 0.32 thresholds *are* source-exact. Full trace, the 12-row flavor table and coverage list
are on [[disorder-region-001-evidence]].

## Canonical oracle and corner cases

Homopolymer anchors pin the normalization endpoints and cutoff (interior residues, away from clipped
edges):

- `WWWW…` (30×W) → interior score **0.0** (normalized minimum) → ordered.
- `PPPP…` (30×P) → interior score **1.0** (normalized maximum), `OverallDisorderContent = 1.0`.
- `EEEE…` (30×E) → interior score **≈ 0.866**, above the 0.542 cutoff → disordered.

Other cases: **null / empty → empty result** (zero summary stats, no throw); **case-insensitive**
(input uppercased first); **non-canonical residues** (X/B/Z…) are preserved in `Residue` but skipped
in the average — a window with *no* recognized residues scores **0.0**; `minRegionLength` larger than
every disordered run → empty `DisorderedRegions` even with residues above threshold. Complexity is
**O(n·w)** time / O(n) space (each residue recomputes its window average).

## Deviations, assumptions, and scope

The **evidence artifact records no assumptions** — every parameter (TOP-IDP Table 2, cutoff 0.542,
Dunker sets, Kyte–Doolittle hydropathy) is traced to peer-reviewed sources. The **implementation** is
an explicitly **simplified, single-feature TOP-IDP heuristic**: it adds no evolutionary profiles,
predicted secondary structure, or trained-model features and is **not competitive with modern
predictors** (the source remarks name IUPred2A, MobiDB-lite) — useful for coarse in-repository
screening only. The one algorithm-doc assumption (ASM-01 / deviation #1) is that local TOP-IDP
composition suffices for first-pass disorder ranking, with non-canonical residues skipped in the
average. A [[research-grade-limitations|research-grade]] implementation.

## References

Campen A. et al. (2008) *Protein Pept Lett* 15(9):956–963 (PMC2676888, PMID 18991772, TOP-IDP scale
+ 0.542 cutoff); Dunker A.K. et al. (2001) *J Mol Graph Model* 19(1):26–59 (order/disorder/ambiguous
sets); Kyte J. & Doolittle R.F. (1982) *J Mol Biol* 157(1):105–132 (hydropathy); Uversky V.N. et al.
(2000) *Proteins* 41(3):415–427 (charge–hydropathy). Full citations in [[disorder-pred-001-evidence]]
(do not duplicate here).
