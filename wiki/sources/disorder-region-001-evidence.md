---
type: source
title: "Evidence: DISORDER-REGION-001 (disordered-region detection + IDR classification)"
tags: [validation, analysis]
doc_path: docs/Evidence/DISORDER-REGION-001-Evidence.md
sources:
  - docs/Evidence/DISORDER-REGION-001-Evidence.md
source_commit: 98b44f1a8112227eb70a11c589272ca8ce62e7af
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: DISORDER-REGION-001

The validation-evidence artifact for test unit **DISORDER-REGION-001** — **disordered-region
detection**, the **contiguous-run aggregation + region classification layer** that sits on top of the
per-residue [[intrinsic-disorder-prediction-top-idp|`PredictDisorder`]] TOP-IDP profile. This is the
**fifth ingested unit of the protein disorder / features family** (after
[[disorder-lc-001-evidence|DISORDER-LC-001]], [[disorder-morf-001-evidence|DISORDER-MORF-001]],
[[disorder-pred-001-evidence|DISORDER-PRED-001]] and
[[disorder-propensity-001-evidence|DISORDER-PROPENSITY-001]]) and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The method itself — how the
per-residue disorder profile is collapsed into regions and how each region is labelled — is written up
in the **"Disordered-region detection"** section of the anchor concept
[[intrinsic-disorder-prediction-top-idp]]; this file records the source trace, the classification
thresholds and the worked oracles. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (five, with authority ranks):**
  - **Campen et al. 2008** "TOP-IDP-Scale" (*Protein Pept Lett* 15(9):956–963, PMC2676888, PMID
    18991772, **rank 1 / primary**) — reused for the scale, cutoff 0.542 and the region idea itself:
    *"contiguous windows predicted as disordered constitute a disordered region."* Window size 21.
  - **Dunker et al. 2001** "Intrinsically disordered protein" (*J Mol Graph Model* 19(1):26–59, PMID
    11381529, **rank 1**) — the disorder/order/ambiguous residue sets, and the biological fact that
    **long (>30-residue) disordered regions** are functionally significant (~33% of eukaryotic
    proteins), with region-function types (molecular-recognition, flexible linkers, entropic
    springs/bristles).
  - **van der Lee et al. 2014** "Classification of intrinsically disordered regions and proteins"
    (*Chem Rev* 114(13):6589–6631, PMC4095912, PMID 24773235, **rank 1 / review**) — the **IDR
    functional/compositional subtypes** the default `RegionType` label reproduces: **proline-rich,
    acidic, basic, Ser/Thr-rich**; compositional bias / homo-repeats (polyQ, polyE) within IDRs; and
    **length-based classification** — short (<30) vs long (≥30) IDRs differ in function (long IDRs
    more likely functionally autonomous). Region boundaries = the disorder→order transition.
  - **Necci et al. 2020/2021** "MobiDB-lite 3.0: fast consensus annotation of intrinsic disorder
    flavors" (*Bioinformatics* 36(22-23):5533–5534, DOI 10.1093/bioinformatics/btaa1045, PMID
    33325498, **rank 1 paper + rank 3 version-pinned reference impl** `BioComputingUP/MobiDB-lite`
    branch `v3`) — a **citable, deterministic disorder-flavor scheme** ingested as the **opt-in
    alternative** `ClassifyRegionFlavorMobiDbLite` (see below).
  - **Wikipedia** "Intrinsically disordered proteins" (**rank 4**, citing Dunker 2001 / Campen 2008 /
    van der Lee 2014) — prevalence of long disorder by kingdom (2.0% archaea / 4.2% eubacteria /
    33.0% eukaryotes, DISOPRED2 / Ward 2004); low-complexity ↔ disorder overlap.
- **Corner cases / failure modes:** short sequences (< window 21) have reduced accuracy (boundary
  effects); order→disorder transitions blur boundaries by ±half-window; **empty predictions → no
  regions**; **all-ordered → no regions**; **all-disordered & length ≥ minLength → exactly one region
  spanning the sequence**; isolated disordered residues (< minLength) → no region; a **trailing
  region** reaching the sequence end must be captured (no off-by-one); classification ties / no
  fraction > 0.25 → fallback label.
- **Datasets / oracles:** homopolymer classification anchors — 30×P → 1 proline-rich region
  (normalized TOP-IDP 1.0); 30×E → acidic (0.866); K/R 40-mer → basic (K fraction 0.5); 31×S →
  Ser/Thr-rich (0.655) with the note that **T alone** normalizes to ≈0.504 (**below** the 0.542
  cutoff); 30×W → 0 regions (ordered, 0.0); mixed W₁₀·P₃₀·W₁₀ → ≥1 central region (boundary may shift
  by window/2). Plus a **hand-traced MobiDB-lite flavor table** (12 rows) computed directly from the
  v3 source.
- **Recommended coverage (18 items):** MUST — empty→no regions; all-ordered→no regions; contiguous
  disorder→one region; region Start/End correct; region MeanScore = mean of constituent residue
  scores; regions < minLength excluded; trailing region captured; the four composition labels
  (proline-rich / acidic / basic / Ser-Thr-rich at fraction > 0.25); Long IDR (length > 30, no bias);
  Standard IDR fallback; confidence ∈ [0,1]. SHOULD — multiple regions split by an ordered segment;
  region at sequence start; region length exactly == minLength. COULD — classification priority when
  multiple biases co-occur (internal order, no published source).

## Region classification — the two schemes

**1. Default `RegionType` heuristic (first-principles, always on).** After a contiguous run ≥
`minRegionLength` is emitted, it is labelled by the dominant residue fraction over the region:
**Proline-rich** (P > 0.25), **Acidic** (D/E > 0.25), **Basic** (K/R > 0.25), **Ser/Thr-rich** (S/T >
0.25), else **Long IDR** (length > 30, van der Lee) or **Standard IDR** (fallback). The **0.25
threshold is an internal ~5×-random heuristic**, NOT the Das & Pappu 2013 value — see the
contradiction note. Each region also carries a first-principles rescaled **`Confidence` ∈ [0,1]**.

**2. `ClassifyRegionFlavorMobiDbLite` (opt-in, source-exact).** Necci et al. 2020's deterministic
scheme, hand-traced from the `v3` `states.py`/`consensus.py`. It is **composition-only over the
region's residues** and maps cleanly onto the profile; boundaries and the default `RegionType` /
`Confidence` are **unchanged** when it is used.
- **Charge classes** (`get_disorder_class`, Das & Pappu 2013 "diagram of states"): translate `RK→+`,
  `DE→−`; `f₊=count(+)/L`, `f₋=count(−)/L`, `FCR=f₊+f₋`, `NCPR=|f₊−f₋|`. **If FCR > 0.35**:
  **Polyampholyte** if `NCPR ≤ 0.35` (or both f₊,f₋ > 0.35), else **PositivePolyelectrolyte** if
  f₊ > 0.35, else **NegativePolyelectrolyte** if f₋ > 0.35; **else WeaklyCharged**.
- **Composition classes** (only if charge = WeaklyCharged), in **priority order**: **CysteineRich** →
  **ProlineRich** → **GlycineRich** → SEG low-complexity → **Polar** (`{S,T,N,Q}`). Enrichment
  fires at fraction **≥ 0.32 (inclusive)** over the region.
- **Window/length:** 9-residue sliding window; sub-regions reported if ≥ 9 residues long.
- Hand-traced oracles (from the artifact, no code run): `RKDERKDE`→Polyampholyte;
  `RKRKRKRKRR`→PositivePolyelectrolyte; `DEDEDEDEDD`→NegativePolyelectrolyte;
  `RKRKPPPPPP`→PositivePolyelectrolyte (charge beats P); `RKRKRKR`+13×A → FCR = 0.35 (NOT > 0.35) →
  WeaklyCharged; `CCCCAAAAAA`→CysteineRich; `PPPPAAAAAA`→ProlineRich; `GGGGAAAAAA`→GlycineRich;
  `SSTTNNQQAA`→Polar; `CCCCPPPPAA`→CysteineRich (C first in priority); 8×C+17×A → C = 0.32 (≥,
  inclusive) → CysteineRich; 7×C+18×A → C = 0.28 < 0.32 → WeaklyCharged. **MobiDB-lite defines no
  per-residue confidence**, so the repository `Confidence` remains a declared heuristic.

## Deviations, assumptions, contradictions

**Contradiction flagged in-source (reference #6):** the artifact explicitly warns that **Das & Pappu
2013's 0.25 boundary is NCPR** (net charge per residue, a globule/coil *conformational-state*
threshold), **NOT a compositional enrichment threshold** for IDR classification. The default
`RegionType` 0.25 composition cutoff is therefore an **internal first-principles heuristic (~5×
random)**, not a value from Das & Pappu — do not cite the paper for it. The distinct MobiDB-lite
charge thresholds (FCR/NCPR at 0.35, enrichment at 0.32) *are* source-exact.

**Assumptions:** region boundaries and the per-residue profile inherit DISORDER-PRED-001 (validated
TOP-IDP, no assumptions there). The default `RegionType` composition thresholds and the rescaled
`Confidence` are declared **first-principles heuristics**; the MobiDB-lite flavors are the sourced,
opt-in alternative added 2026-06-24. No source contradictions beyond the 0.25-NCPR clarification.
