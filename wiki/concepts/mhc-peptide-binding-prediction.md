---
type: concept
title: "MHC-peptide binding prediction + binder classification (thresholds / BIMAS+SMM matrix / MHCflurry)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-MHC-001-Evidence.md
source_commit: ad78834fc39b17b44fb5fb87ba380d83ef7a3eee
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-mhc-001-evidence
      evidence: "Test Unit ID: ONCO-MHC-001, Algorithm: MHC-Peptide Binding (length filtering + affinity/%rank thresholds; matrix-based prediction)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:hla-nomenclature-and-allele-specific-loh
      source: onco-mhc-001-evidence
      evidence: "The MHCflurry predictor keys each allele to a 37-residue MHC-pocket pseudosequence by HLA allele name (e.g. HLA-A*02:01 -> YFAMYGEKVAHTHVDTLYGVRYDHYYTWAVLAYTWYA); HLA LOH removes the presentation platform for the neoantigens restricted to a lost allele."
      confidence: medium
      status: current
    - predicate: relates_to
      object: concept:immune-infiltration-deconvolution
      source: onco-mhc-001-evidence
      evidence: "The evidence file explicitly analogizes its caller-supplied trained matrix / non-embedded weights packaging boundary to 'CIBERSORT LM22 in ONCO-IMMUNE-001' (twice)."
      confidence: medium
      status: current
    - predicate: relates_to
      object: concept:neoantigen-peptide-generation
      source: onco-mhc-001-evidence
      evidence: "This is the affinity gate downstream of neoantigen-peptide generation (ONCO-NEO-001, GenerateNeoantigenPeptides); the resolved class I 8-14 window propagates to the neoantigen generator default."
      confidence: high
      status: current
---

# MHC-peptide binding prediction + binder classification

The **twenty-second ingested Oncology unit** (**ONCO-MHC-001**) and the wiki's first
**MHC / peptide-presentation prediction** method — the immuno-oncology step that decides which
peptides an HLA (MHC) molecule presents, the affinity gate at the heart of **neoantigen candidate
scoring**. The literature-traced record is [[onco-mhc-001-evidence]]; [[test-unit-registry]] tracks
the unit and [[algorithm-validation-evidence]] describes the evidence-artifact pattern.

The unit spans three layers that share one output space — a predicted **binding affinity (IC50, nM)**
or percentile (**%Rank**) reduced to a **binder tier** (`Strong` / `Weak` / `NonBinder`):

1. **Classification** — given an affinity or percentile plus the peptide's length, apply MHC-class
   thresholds.
2. **Matrix-based prediction** — turn a peptide + a caller-supplied position-weight matrix into an
   IC50 (BIMAS product rule; SMM additive rule).
3. **Neural-network prediction** — a ported **MHCflurry 2.0** pan-allele class-I affinity network.

## 1. Binder classification (thresholds)

### Affinity (IC50) tiers — `ClassifyBindingAffinity`

From Sette et al. (1994)'s ≈ 500 nM CTL-response threshold and the IEDB affinity tiers:

| IC50 (nM) | Tier |
|-----------|------|
| **< 50** | Strong binder |
| **< 500** (and ≥ 50) | Weak / intermediate binder |
| **≥ 500** | NonBinder |

### Percentile (%Rank) tiers — `ClassifyBindingRank`

From NetMHCpan-4.1 / NetMHCIIpan-4.0 (Reynisson et al. 2020) — the tiers differ by MHC **class**:

| %Rank | Class I | Class II |
|-------|---------|----------|
| Strong binder | **< 0.5** | **< 2** |
| Weak binder | **< 2** | **< 10** |

**All cutoffs are strict `<`.** A peptide *exactly* at a cutoff falls to the looser tier: IC50 = 50 nM is
**not** strong (→ Weak), IC50 = 500 nM is **not** a binder (→ NonBinder); %Rank = 0.5 is **not** strong
(→ Weak, since 0.5 < 2), %Rank = 2.0 is **not** weak (→ NonBinder). %Rank is a percentile in [0, 100];
IC50 is a positive concentration — negative / NaN / ∞ IC50 and negative / NaN / > 100 %Rank are rejected.

### Peptide-length validity by class — `IsValidPeptideLength`

| MHC class | Accepted peptide length | Source |
|-----------|-------------------------|--------|
| **Class I** | **8–14 aa** (default 8–11) | NetMHCpan-4.1 class I window (8/9/…/14-mer options) |
| **Class II** | **13–25 aa** (9-mer binding core) | IEDB class II tool description |

The class I ceiling was widened **11 → 14** (`MhcClassIMaxPeptideLength`) to match the full NetMHCpan-4.1
window; the change propagates to the neoantigen-peptide generator default (**ONCO-NEO-001**,
`GenerateNeoantigenPeptides`). The MHCflurry encoder supports lengths **5–15** (`max_length = 15`), so it
scores 12/13/14-mers unchanged; callers may pass `maxLength` up to 15 to reach the encoding ceiling.

## 2. Matrix-based prediction

Both variants take a peptide and a **position-specific scoring matrix** (one row per peptide position;
`Predict…` throws `ArgumentException` if `peptide.Length ≠ matrix rows`). An **unlisted / ambiguous
residue** is neutral — coefficient **1.0** for BIMAS (no effect on a product), contribution **0.0** for
SMM (no effect on a sum).

- **BIMAS product rule** (Parker, Bednarek & Coligan 1994 / BIMAS docs): the score is a **half-time of
  dissociation** — running score starts 1.0 and is multiplied by each position's coefficient, then by a
  final constant. Per allele: **20 × 9 = 180 coefficients + 1 final constant**. Accuracy ≈ factor of 5.
  Worked: `LMV` = 2.0 · 3.0 · 1.5 · const 10 = **90**; `AAA` (all unlisted) = 1·1·1 · 10 = **10**.
- **SMM additive rule** (Peters & Sette 2005) with the IEDB **log50k linearisation**:

  ```
  log50k = 1 − log(IC50)/log(50000)   ⇒   IC50 = 50000^(1 − score)
  ```

  Exact anchors: score 0 → **50000 nM** (NonBinder), 0.5 → **√50000 = 223.6068 nM** (Weak), 0.9 →
  **2.9505 nM** (Strong), 1.0 → **1 nM** (Strong); round-trip IC50 = 500 nM → score **0.4256252**. SMM
  IC50 is always finite and **> 0** for any finite score, so it always satisfies the classifier's
  `IC50 > 0` precondition. `PredictAndClassifySmm` chains predict → classify.

## 3. MHCflurry 2.0 pan-allele class-I neural network

A source-verified port of **MHCflurry 2.0** (O'Donnell et al. 2020; Apache-2.0 source + `models_class1_pan`
20200610 weights):

- **Encoding** — 21-symbol alphabet `ACDEFGHIKLMNPQRSTVWYX`, **BLOSUM62** vectors. The peptide uses
  `left_pad_centered_right_pad` (`max_length = 15`): placed left-aligned, centred, and right-aligned in a
  3 × 15 = 45-position layout → 45 × 21 = **945** values. Each allele → a fixed **37-residue MHC-pocket
  pseudosequence** (e.g. `HLA-A*02:01 → YFAMYGEKVAHTHVDTLYGVRYDHYYTWAVLAYTWYA`) → 37 × 21 = **777** values.
- **Network** — peptide (945) and allele (777) Flattened + concatenated → **1722**, then `tanh` `Dense`
  layers and a final `Dense(1, sigmoid)`. Two ensemble topologies: `feedforward` and
  `with-skip-connections` (each hidden layer i ≥ 1 receives `concat(prev_prev_input, prev_activation)`);
  the port reads `topology` per network and wires skips accordingly (ignoring it gives a dimension
  mismatch). The allele **embedding is empty in the `.npz`** — the representation is supplied at predict
  time from the bundled pseudosequence table.
- **Output** — `to_ic50(x) = 50000^(1 − x)` (the same transform as SMM). The 10-network ensemble combines
  by **geometric mean** of the per-network IC50s (`exp(mean(log ic50))`).

**Oracle** (`mhcflurry` 2.1.5; C# port reproduces both single-net and full-ensemble IC50s to **< 0.03 %**):
SIINFEKL/HLA-A\*02:01 (self, non-binder) ≈ 11.5 µM; known HLA-A\*02:01 binders GILGFVFTL (flu M1) ≈ 19 nM,
NLVPMVATV (CMV pp65) ≈ 17 nM, ELAGIGILTV (MART-1) ≈ 84–119 nM. Length-12/13/14 variants (≈ 25–33 µM)
confirm the widened 8–14 window scores unchanged.

## Packaging boundary — trained data is caller-supplied

Like several oncology units, ONCO-MHC-001 embeds **algorithms, not proprietary trained data**. No
redistributable, cross-verifiable trained HLA coefficient matrix was obtainable (BIMAS files served by a
defunct CGI; Parker 1994's table paywalled; IEDB SMM matrices non-commercial), so the library ships the
**scoring rules** + a `LoadScoringMatrix` loader and the caller supplies matrix values under their own
licence. The full 80 MB MHCflurry ensemble is likewise not embedded (near-incompressible float32, repo
health) — one member (~4.6 MB) is embedded for CI parity, the rest load via `LoadWeightPack`; only the
~0.7 MB Apache-2.0 pseudosequence table is bundled. The evidence file **explicitly analogizes** this to
the caller-supplied **CIBERSORT LM22** matrix of [[immune-infiltration-deconvolution]] (ONCO-IMMUNE-001).

## Relation to the immuno-oncology / neoantigen area

This is the **affinity gate** downstream of [[neoantigen-peptide-generation|neoantigen-peptide generation]]
(**ONCO-NEO-001**, `GenerateNeoantigenPeptides`) and upstream of neoantigen ranking: a mutant peptide is a
candidate only if it is a valid-length, sufficiently strong binder to one of the tumour's HLA alleles.
[[hla-nomenclature-and-allele-specific-loh]] (ONCO-HLA-001) is the presentation-platform sibling: HLA LOH
**removes** an allele, so the neoantigens this predictor scores against that allele can no longer be
presented — immune escape. [[immune-infiltration-deconvolution]] (ONCO-IMMUNE-001) quantifies the immune
*cell content* of the microenvironment rather than the antigen-presentation machinery.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for MHC-binding classification and prediction,
**not for clinical or diagnostic use**. The threshold, matrix, and MHCflurry algorithms are all fully
source-traced; the **trained** matrices/weights are caller-supplied (a licensing/size packaging boundary,
not a correctness assumption). The class I 8–14 window is resolved against NetMHCpan-4.1. No source
contradictions — the NetMHCpan %Rank tiers, the Sette/IEDB IC50 tiers, the BIMAS/Parker/SMM matrix rules,
and the MHCflurry network cover complementary methods and agree on the shared affinity/percentile framing.
