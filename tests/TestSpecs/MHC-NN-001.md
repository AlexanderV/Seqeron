# Test Specification: MHC-NN-001

**Test Unit ID:** MHC-NN-001
**Area:** Oncology
**Algorithm:** MHCflurry Pan-Allele NN Binding Affinity
**Status:** ☑ Validated (Stage A ✅ / Stage B ✅ / CLEAN) — 2026-06-25
**Last Updated:** 2026-06-25

> Independently re-validated under the two-stage protocol against O'Donnell et al. (2018/2020) MHCflurry,
> the openvax/mhcflurry source, and a live `mhcflurry==2.1.5` + `models_class1_pan` differential oracle.
> See `docs/Validation/reports/MHC-NN-001.md`.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | O'Donnell et al. (2018) "MHCflurry", Cell Systems 7(1):129-132.e4 |
| 2 | O'Donnell, Rubinsteyn, Laserson (2020) "MHCflurry 2.0", Cell Systems 11(1):42-48.e7, doi:10.1016/j.cels.2020.06.010 |
| 3 | openvax/mhcflurry (Apache-2.0): regression_target.py, encodable_sequences.py, amino_acid.py, ensemble_centrality.py, class1_neural_network.py |
| 4 | mhcflurry 2.1.5 + models_class1_pan (release 20200610) live oracle |

## 2. Canonical Method(s)

`EncodePeptide`, `EncodePseudosequence`, `GetPseudosequence`, `ToIc50`, `LoadWeightPack`,
`Network.ForwardRaw`, `Network.PredictIc50`, `PredictIc50(ensemble)`, `PredictIc50WithPseudosequence`,
`PredictAndClassify`.

- **Source file:** `MhcflurryAffinityPredictor.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MhcflurryAffinityPredictor_PredictIc50_Tests.cs` (25 tests)

## 3. Contract / Invariants

- R: IC50 = 50000^(1-x) > 0; strictly decreasing in raw output; inverse of from_ic50.
- R: peptide encoding length 945 (3×15×21), allele 777 (37×21), input 1722; padding/unknown → X vector.
- R: ensemble combiner = geometric mean of per-network IC50s (exp(mean(log·))); duplicated net = single net.
- D: deterministic given weights; lowercase folded to canonical AA; unsupported residue → X (lenient).

## 4. Cross-check / Differential Oracle

- **Reference:** mhcflurry 2.1.5 (`models_class1_pan`), member `PAN-CLASS1-1-3ed9fb2d2dcc9803` (feedforward
  [512,512]); bundled `mhcflurry_single_net.bin` verified byte-identical to the official npz.
- **Numbers (NumPy oracle == C# golden, 6 dp):** SIINFEKL/A*02:01 11483.195201; GILGFVFTL/A*02:01 19.123150;
  NLVPMVATV 17.542640; ELAGIGILTV 119.054961; AAAWYLWEV 16.559303; SIINFEKL/B*07:02 28830.796646;
  SLYNTVATL 28.972028; CINGVCWTV 92.105940. C# parity within RelTol 1e-3 (observed <0.03%).

## 5. Validation Checklist

- [x] Stage A: sources retrieved this session; formula/constants/ordering/ensemble/topology confirmed.
- [x] Stage B: implementation reviewed; bundled weights bit-exact to official npz; NumPy oracle reproduced.
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
- [x] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.

## 6. Notes

- Documented divergence: C# maps unsupported residues to X (MHCflurry `allow_unsupported_amino_acids=True`
  path) rather than raising as MHCflurry's default; benign, no parity impact on valid input.
- Documented boundary: full 10-network ensemble is caller-supplied via `LoadWeightPack` (large, heterogeneous
  topologies). Single-member parity + geometric-mean combiner are fully verified — accepted, not LIMITED.
