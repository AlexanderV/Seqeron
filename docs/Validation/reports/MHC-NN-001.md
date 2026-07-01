# Validation Report: MHC-NN-001 — MHCflurry Pan-Allele NN Binding Affinity

- **Validated:** 2026-06-25   **Area:** Oncology
- **Canonical method(s):** `MhcflurryAffinityPredictor.EncodePeptide` / `EncodePseudosequence` /
  `GetPseudosequence` / `ToIc50` / `LoadWeightPack` / `Network.ForwardRaw` / `Network.PredictIc50` /
  `PredictIc50(ensemble)` / `PredictIc50WithPseudosequence` / `PredictAndClassify`, chained into
  `OncologyAnalyzer.ClassifyBindingAffinity`.
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Authoritative sources opened this session

1. O'Donnell, Rubinsteyn, Hodes, Carrasco, Hammerbacher, "MHCflurry: Open-Source Class I MHC Binding
   Affinity Prediction", Cell Systems 7(1):129-132.e4 (2018).
2. O'Donnell, Rubinsteyn, Laserson, "MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented
   Peptides by Incorporating Antigen Processing", Cell Systems 11(1):42-48.e7 (2020),
   doi:10.1016/j.cels.2020.06.010.
3. openvax/mhcflurry source (Apache-2.0), modules retrieved this session:
   `regression_target.py`, `encodable_sequences.py`, `amino_acid.py`, `ensemble_centrality.py`,
   `class1_neural_network.py`.
4. The installed `mhcflurry==2.1.5` package + the fully fetched `models_class1_pan` weight pack
   (release 20200610, `models.combined`, 10 PAN-CLASS1 networks) — used as the live differential oracle.

## Stage A — Description

- **IC50 transform.** `regression_target.to_ic50`: `ic50 = max_ic50 ** (1.0 - x)` with `max_ic50 = 50000.0`;
  inverse `from_ic50`: `x = 1 - ln(ic50)/ln(max_ic50)` clamped to [0,1]. Matches the code's
  `ToIc50 = 50000^(1-x)` and the round-trip test exactly.
- **Peptide encoding.** `encodable_sequences.left_pad_centered_right_pad`: shape `max_length*3` (=45)
  positions, three copies (left-aligned, centred, right-aligned); `center_left_padding =
  floor((max_length-length)/2)`, `center_left_offset = max_length + center_left_padding`;
  `unknown_character = "X"`; arbitrary minimum length 5. Confirmed `max_length=15`, BLOSUM62 from the
  bundled member's hyperparameters.
- **Amino-acid order / BLOSUM62.** `amino_acid.py`: `AMINO_ACIDS = ACDEFGHIKLMNPQRSTVWYX` (X appended last,
  index 20); `vector_encoding_length("BLOSUM62") = 21`; matrix reindexed to that order. The C# `Blosum62`
  table and `AminoAcidOrder` reproduce this; the bundled X-row is all-zero except the X channel (=1).
- **Ensemble combiner.** `ensemble_centrality.CENTRALITY_MEASURES["mean"]` operates in log space →
  geometric mean of per-network IC50s. Matches `exp(mean(log(ic50)))`.
- **Topology.** `class1_neural_network.make_network`: feedforward; tanh hidden, sigmoid output; batch-norm
  default False; for `with-skip-connections`, layer 1 input = concat(merged, out[0]), layer i≥2 input =
  concat(out[i-2], out[i-1]); output layer has no skip. Matches the C# `Network.ForwardRaw` wiring.
- **Edge semantics.** Unknown allele → lookup failure; unsupported residue handling is governed by
  `allow_unsupported_amino_acids` (default False → KeyError; True → map to X). See divergence note below.

Stage A verdict: **PASS** — every formula, constant, ordering and the ensemble/topology semantics match the
primary source and the openvax implementation verbatim.

## Stage B — Implementation

Code path reviewed: `MhcflurryAffinityPredictor.cs` lines 45–681; chain into
`OncologyAnalyzer.ClassifyBindingAffinity` (line 8197, strong<50 nM / weak<500 nM).

### Independent oracle reproduction (this session)

The full `models_class1_pan` weight pack was available locally. The bundled
`tests/.../TestData/Mhcflurry/mhcflurry_single_net.bin` was decoded and compared byte-for-byte against the
official `weights_PAN-CLASS1-1-3ed9fb2d2dcc9803.npz` (the smallest member, feedforward [512,512]):

- Header: magic `MHCF`, version 1, 1 network, topology 0, 3 dense layers (1722→512, 512→512, 512→1).
- `max|kernel diff| = 0.000e+00`, `max|bias diff| = 0.000e+00` for all three layers, 0 trailing bytes.
  The empty `(0,777)` embedding array from the npz is correctly dropped (pseudosequence supplied externally).

An **independent NumPy forward pass** over those exact weights (BLOSUM62 + left_pad_centered_right_pad +
tanh/tanh/sigmoid + `50000^(1-x)`) reproduced the test's golden IC50s to 6 decimals:

| Peptide | Allele | NumPy oracle IC50 (nM) | C# test golden | match |
|---|---|---|---|---|
| SIINFEKL | HLA-A*02:01 | 11483.195201 | 11483.195201 | ✓ |
| GILGFVFTL | HLA-A*02:01 | 19.123150 | 19.123150 | ✓ |
| NLVPMVATV | HLA-A*02:01 | 17.542640 | 17.542640 | ✓ |
| ELAGIGILTV | HLA-A*02:01 | 119.054961 | 119.054961 | ✓ |
| AAAWYLWEV | HLA-A*02:01 | 16.559303 | 16.559303 | ✓ |
| SIINFEKL | HLA-B*07:02 | 28830.796646 | 28830.796646 | ✓ |
| SLYNTVATL | HLA-A*02:01 | 28.972028 | 28.972028 | ✓ |
| CINGVCWTV | HLA-A*02:01 | 92.105940 | 92.105940 | ✓ |

The C# `PredictIc50` reproduces these within RelTol 1e-3 (observed gap well under the claimed <0.03%; the
golden values themselves are bit-equal to 6 dp to the NumPy oracle). The MHCflurry `Class1AffinityPredictor`
high-level API (its model-selected single pan network) was also run live and produced the same-order values
(e.g. SIINFEKL/A*02:01 ≈ 11927 nM, GILGFVFTL/A*02:01 ≈ 20 nM), confirming the pipeline shape.

### with-skip-connections coverage

The bundled member is feedforward, so the densenet skip wiring was not exercised by the original suite. A
tiny synthetic `with-skip-connections` pack (3 hidden + output, deterministic NumPy weights) was generated;
the independent NumPy forward pass following the verified wiring gives raw output **0.501112284692050**
(IC50 220.9319 nM). The embedded `mhcflurry_skip_net.bin` is byte-identical (md5) to that generated pack,
and `Network.ForwardRaw` reproduces the raw output within 1e-6 — proving the skip-connection branch.

### Test-quality audit & gaps closed

The original 19-test fixture's golden values trace to the real model (not code echoes) and cover the
left_pad_centered_right_pad layout, the pseudosequence table, the IC50 transform, single-net parity,
ranking, the geometric-mean ensemble, the predict→classify chain, and loader edge cases. Gaps found and
**fixed this session** (+6 tests, 25 total):

- `EncodePseudosequence` had no direct test → added position-major BLOSUM62 + diagonal-dominance check.
- Non-canonical residue (X-fallback) and case-insensitivity were untested → added.
- `ToIc50` inverse/monotonicity (vs `from_ic50`) untested → added round-trip + strictly-decreasing checks.
- `Network.PredictIc50` (instance) and `PredictIc50WithPseudosequence` overloads untested → added agreement
  test (and confirms allele lookup == direct pseudosequence).
- `with-skip-connections` topology unexercised → added the synthetic-pack parity test above.
- `ForwardRaw` wrong-width guard untested → added.

### Divergence (documented, not a defect)

MHCflurry's default `allow_unsupported_amino_acids=False` **raises** on residues outside the 20 standard AAs;
the C# `IndexOfResidue` instead maps any unsupported residue to the X (unknown) vector — i.e. it implements
the `allow_unsupported_amino_acids=True` path. This is explicitly documented in the source comment and is a
deliberately more lenient input policy; it does not affect parity on any valid peptide/allele input and the
X-fallback itself is the correct MHCflurry "unknown" vector. Recorded as a note, not a defect.

### Boundary (acceptable, documented)

The full 10-network `models_class1_pan` ensemble is **caller-supplied** via `LoadWeightPack` (the pack is
~80 MB and the members have heterogeneous topologies/widths). The bundled single member is verified
bit-exact and its forward pass + IC50 transform reproduce the published oracle exactly; the geometric-mean
combiner is verified via the duplicated-network identity. This caller-supplied weight pack is an accepted
boundary, not a limitation.

Stage B verdict: **PASS** — code faithfully realises the validated description; all expected values trace to
the real MHCflurry weights / published transform; the full suite is green (0 failed, 0 warnings on changed
files).

## Verdict & follow-ups

✅ **CLEAN.** No defect. Description and implementation independently confirmed against O'Donnell et al.
(2018/2020) and the openvax/mhcflurry source + live `models_class1_pan` oracle. Test coverage extended from
19 to 25 tests to cover every public method/overload and the skip-connection topology. One documented
divergence (lenient unsupported-residue handling) and one documented boundary (caller-supplied full ensemble)
— neither is LIMITED.
