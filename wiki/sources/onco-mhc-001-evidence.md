---
type: source
title: "Evidence: ONCO-MHC-001 (MHC-peptide binding — length/affinity/%rank classification + matrix + MHCflurry neural net)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-MHC-001-Evidence.md
sources:
  - docs/Evidence/ONCO-MHC-001-Evidence.md
source_commit: ad78834fc39b17b44fb5fb87ba380d83ef7a3eee
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-MHC-001

The validation-evidence artifact for test unit **ONCO-MHC-001** — **MHC-Peptide Binding**
(length filtering + affinity/%rank classification thresholds; matrix-based prediction; and a
ported MHCflurry 2.0 neural-network affinity predictor). The **twenty-second ingested unit of
the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The method is synthesized in its own
concept, [[mhc-peptide-binding-prediction]]; [[test-unit-registry]] tracks the unit.

## What this file records

The unit has three layers: (1) **classification** of a peptide given an affinity (IC50) or
percentile (%Rank) plus peptide-length validity by MHC class; (2) **matrix-based prediction**
(BIMAS product rule + SMM additive rule → IC50); (3) a **ported MHCflurry 2.0** pan-allele
class-I neural network.

- **Online sources (mutually consistent; complementary methods, no contradictions):**
  - **Reynisson et al. (2020)** *Nucleic Acids Res* 48(W1):W449–W454 — **NetMHCpan-4.1 /
    NetMHCIIpan-4.0**: verbatim **%Rank** binder tiers — class I **SB < 0.5 %, WB < 2 %**; class II
    **SB < 2 %, WB < 10 %**; class I peptide **length 8–14** (default 8–11) (rank 1).
  - **Sette et al. (1994)** *J Immunol* 153(12):5586–5592 — the **≈ 500 nM affinity threshold**
    (preferably ≤ 50 nM) for a peptide to elicit a CTL response; the biological basis of the 50 nM
    (strong) / 500 nM (binder) IC50 cutoffs (rank 1).
  - **IEDB** threshold help — **IC50 tiers**: < 50 nM high affinity, < 500 nM intermediate, < 5000 nM
    low; 500 nM = strong-binder demarcation (rank 5; page 403 → verbatim from search snippet).
  - **Roomp, Antes & Lengauer (2010)** *BMC Bioinformatics* 11:90 — independently corroborates the
    **500 nM binder/non-binder demarcation** (rank 1).
  - **IEDB** class II tool description — class II peptides **13–25 aa**, **9-mer binding core** (rank 5).
  - **BIMAS** HLA peptide-motif scoring docs (Parker; archived NIH/CIT) — the **product rule**: running
    score starts 1.0, is multiplied by a per-position coefficient (20×9 = 180 coefficients + 1 final
    constant per allele), unlisted residue coefficient = **1.0** (neutral); score = estimated **half-time
    of dissociation** (rank 5, restating Parker 1994).
  - **Parker, Bednarek & Coligan (1994)** *J Immunol* 152(1):163–175 — primary HLA-A2 coefficient table
    (**180 coefficients = 20 aa × 9 positions**, product-of-coefficients rule, accuracy ≈ factor of 5);
    table itself paywalled (rank 1).
  - **Peters & Sette (2005)** *BMC Bioinformatics* 6:132 + **IEDB log50k** — the **SMM** additive matrix
    and the score↔IC50 linearisation `log50k = 1 − log(IC50)/log(50000)` ⇒ **`IC50 = 50000^(1 − score)`**
    (score 0 → 50000 nM, 1 → 1 nM, 0.5 → √50000 = 223.6068 nM) (rank 1/3).
  - **O'Donnell et al. (2020)** *Cell Systems* 11(1):42–48 + **MHCflurry** source (Apache-2.0) —
    **MHCflurry 2.0** class-I pan-allele affinity NN: BLOSUM62 `left_pad_centered_right_pad` peptide
    encoding (45×21 = 945), **37-residue allele pseudosequence** (37×21 = 777, bundled), Flatten +
    concatenate → 1722, `tanh` Dense stack (`feedforward` + `with-skip-connections`), final
    `Dense(1, sigmoid)`, output **`to_ic50(x) = 50000^(1−x)`**, ensemble = **geometric mean** of
    per-network IC50s (rank 1 paper / rank 3 reference impl `mhcflurry` 2.1.5, `models_class1_pan`
    20200610).

- **Implemented surface:** `ClassifyBindingAffinity` (IC50 tiers), `ClassifyBindingRank` (%Rank tiers,
  class I/II), `IsValidPeptideLength` (class I 8–14, class II 13–25), a combined length-gate + affinity
  helper; matrix prediction — BIMAS `Predict…` (product rule) and SMM (additive + `IC50 = 50000^(1−score)`),
  `PredictAndClassifySmm`, a `LoadScoringMatrix` loader (`CONST=` / `RESIDUE=VALUE`); and the MHCflurry
  port (peptide/allele BLOSUM62 encoders, forward pass with per-network `topology` skip-wiring,
  `to_ic50`, geometric-mean combiner, `LoadWeightPack`).

- **Documented corner cases / failure modes:** peptide length outside encoder [5,15] →
  `ArgumentOutOfRangeException`; non-canonical residue → the `X` index (neutral); unknown allele →
  `KeyNotFoundException`; BIMAS/SMM peptide length must equal the matrix position count →
  `ArgumentException`; SMM IC50 always finite and > 0 (50000^x > 0); unlisted residue → BIMAS coefficient
  1.0 / SMM contribution 0.0 (both no-effect); **strict `<`** boundary semantics (IC50 exactly 50 →
  not strong; %Rank exactly 0.5 → not strong; exactly 500 nM / 2.0 % → not a binder); the empty
  embedding array in `weights_*.npz` (allele representation supplied at predict time).

- **Datasets (deterministic worked oracles):**
  - **Class I %Rank:** 0.4 → Strong, 0.5 → Weak, 1.0 → Weak, 2.0 → NonBinder, 5.0 → NonBinder.
  - **IC50 (nM):** 10 → Strong, 50 → Weak, 200 → Weak, 500 → NonBinder, 1000 → NonBinder.
  - **Peptide length:** class I 7 invalid / 9 valid / 14 valid / 15 invalid; class II 12 invalid /
    15 valid.
  - **SMM transform anchors:** score 0 → 50000 (NonBinder), 0.5 → 223.6068 (Weak), 0.9 → 2.9505
    (Strong), 1.0 → 1 (Strong); round-trip IC50 = 500 → score 0.4256252.
  - **BIMAS product:** `LMV` coeffs 2·3·1.5 × const 10 → T½ 90; `AAA` (all unlisted) 1·1·1 × 10 → 10.
  - **MHCflurry oracle** (`mhcflurry` 2.1.5, single-net + full 10-net geometric mean; C# port matches to
    **< 0.03 %** rel error): SIINFEKL/HLA-A\*02:01 (self, non-binder) 11483.2 / 11927.2 nM; GILGFVFTL
    (flu M1) 19.12 / 19.96; NLVPMVATV (CMV pp65) 17.54 / 16.57; ELAGIGILTV (MART-1) 119.05 / 83.55;
    length-12/13/14 (`GILGFVFTLAAA…`) 25274.9 / 32389.1 / 32972.2 nM confirm the widened 8–14 window.

- **Coverage recommendations:** MUST — IC50 & %Rank classification at/around both cutoffs; class II %Rank
  (2/10); length validity at class boundaries; invalid inputs rejected (negative/NaN/∞ IC50,
  negative/NaN/>100 %Rank); SMM transform anchor IC50s exact; `PredictAndClassifySmm` chains
  predict→classify; BIMAS product with unlisted = 1.0; loader parse + malformed → `FormatException`.

## Deviations and assumptions

- **Caller-supplied trained matrix (Framework):** no redistributable, cross-verifiable trained HLA
  coefficient matrix was obtainable — BIMAS files are served by a now-defunct dynamic CGI (unarchived),
  Parker 1994's 180-value table is paywalled, and IEDB SMM matrices carry a non-commercial /
  no-redistribution licence. The library embeds only the published **scoring rules** (BIMAS product; SMM
  `IC50 = 50000^(1−score)`) plus `LoadScoringMatrix`; the caller supplies the values under their own
  licence. **Explicitly mirrors the caller-supplied CIBERSORT LM22 matrix of ONCO-IMMUNE-001**
  ([[immune-infiltration-deconvolution]]).
- **Full MHCflurry ensemble not embedded (size):** the 10-network `models_class1_pan` weights total
  ≈ 80 MB of near-incompressible float32, declined for repo health. The encoders, the 37-residue
  pseudosequence table (Apache-2.0, ~0.7 MB, bundled), the forward-pass engine, `to_ic50`, and the
  geometric-mean combiner are all ported; **one** ensemble member (~4.6 MB) is embedded for CI parity;
  the full ensemble loads from a caller-supplied weight pack via `LoadWeightPack`. Not a correctness
  assumption — a packaging boundary, again analogous to ONCO-IMMUNE-001's LM22.
- **RESOLVED — class I length window = 8–14:** `MhcClassIMaxPeptideLength` widened 11 → 14 (the full
  NetMHCpan-4.1 class I window; MHCflurry encoder supports 5–15 via `max_length = 15`), propagating to
  `IsValidPeptideLength` and the `GenerateNeoantigenPeptides` default (ONCO-NEO-001). The narrower
  pVACtools 8–11 default is no longer assumed; lengths 12–14 reproduce the live oracle within RelTol 1e-3,
  8–11 byte-identical.
- **Strict-inequality boundary reading:** the cutoffs are stated as strict `<` in the sources, so a
  peptide exactly at a cutoff falls to the looser tier — the tests encode this literal reading.

No source contradictions — the NetMHCpan %Rank tiers, the Sette/IEDB IC50 tiers, the BIMAS/Parker/SMM
matrix rules, and the MHCflurry neural net cover complementary prediction/classification methods and
agree on the shared affinity/percentile framing.
