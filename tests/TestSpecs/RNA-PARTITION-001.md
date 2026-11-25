# Test Specification: RNA-PARTITION-001

**Test Unit ID:** RNA-PARTITION-001
**Area:** RnaStructure
**Algorithm:** RNA Partition Function (McCaskill) on the Turner 2004 model and Boltzmann Structure Probability
**Status:** ☑ Complete (re-validated 2026-11-25, B12)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | McCaskill JS (1990) *Biopolymers* 29:1105-1119 | 1 | https://doi.org/10.1002/bip.360290621 (PMID 1695107) | 2026-06-14 |
| 2 | Will S, MIT 18.417 — McCaskill inside / base-pair-probability slides | 1 | https://math.mit.edu/classes/18.417/Slides/mccaskill.pdf, `mccaskill2.pdf` | 2026-06-14 |
| 3 | Mathews DH et al. (2004) *PNAS* 101:7287; NNDB Turner 2004 | 1 | https://rna.urmc.rochester.edu/NNDB/turner04/ | 2026-11-25 |
| 4 | ViennaRNA 2.7.2 sdist: `partfunc/`, `probabilities/probs_structures.c` (`pr_structure`, `mean_bp_distance`, `centroid`), `params/params.c` (`RESCALE_dG`, `pf_scale`), `params/default.c` + `intl*dH.h` (enthalpies) | 1 (reference implementation) | PyPI `viennarna-2.7.2.tar.gz`; Python `RNA` 2.7.2 | 2026-11-25 |
| 5 | Lorenz R et al. (2011) ViennaRNA Package 2.0, *Algorithms Mol Biol* 6:26 | 1 | https://doi.org/10.1186/1748-7188-6-26 | 2026-11-25 |

### 1.2 Key Evidence Points

1. `Z = Σ_S exp(−E(S)/RT)` over all pseudoknot-free structures, E from the Turner 2004 loop model — Sources 1, 3, 5.
2. `P(i,j) = Σ_{S∋(i,j)} exp(−E(S)/RT) / Z` via inside + outside recursion — Sources 1, 2.
3. `G = −RT ln Z` (`pf()`), `P(S) = exp((G − E(S))/RT)` (`pr_structure`), `⟨d⟩ = 2Σ p(1−p)` (`mean_bp_distance`), centroid = pairs with `p > 0.5` — Source 4.
4. `RT = 1.98717e-3·T` kcal/mol (`GASCONST`); default T = 37 °C = 310.15 K — Source 4.
5. Temperature: `ΔG(T) = ΔH − (ΔH − ΔG₃₇)·T/310.15` for every table entry with the Turner 2004 enthalpies; dangles / exterior+multiloop mismatches clamped ≤ 0 — Source 4 (`params.c`).
6. `pf_scale` per nucleotide (`sfact = 1.07`, from the MFE) keeps long sequences in the double range — Source 4.

### 1.3 Documented Corner Cases

- Hairpin loops < 3 nt forbidden; a sequence that cannot pair has `Z = 1` (Sources 2, 3).
- N / IUPAC / gap characters never pair (ViennaRNA `vrna_nucleotide_encode`); T ≡ U.
- `Z` may exceed the double range (`> e^709`): `ln Z` / `G` stay finite.

### 1.4 Known Failure Modes / Pitfalls

1. Ambiguous decomposition over-counts structures (Source 2) — tested by exhaustive enumeration.
2. ViennaRNA's `pf()` weights exterior dangles with `SMOOTH()` (≤ 0.017 kcal/mol) and its MFE truncates rescaled parameters to 0.01 kcal/mol at T ≠ 37 °C — locked tolerances (Evidence §9).
3. `params_load` then `RNA.md()`: pass the `md` explicitly to `fold_compound` (oracle hygiene).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculatePartitionFunction(string, int, double)` | RnaSecondaryStructure | **Canonical** | McCaskill O(n³) on the Turner model; Z, P(i,j), G, unpaired, MFE + P(MFE), diversity, centroid |
| `CalculateStructureProbability(double, double, double)` | RnaSecondaryStructure | **Canonical** | `exp((G − E)/RT)` |
| `CalculateUnpairedProbabilities(string, int, double)` | RnaSecondaryStructure | **Delegate** | projection of the canonical result (RNA-ACCESS-001 tests) |
| `GenerateRandomRna(int, double)` | RnaSecondaryStructure | **Delegate** | forwards to seeded overload with `new Random()` |
| `GenerateRandomRna(int, Random, double)` | RnaSecondaryStructure | **Canonical** | deterministic random generation (seeded) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `Z ≥ 1`; `G = −RT ln Z ≤ MFE` | Yes | Source 1 (open chain); MFE is one term |
| INV-2 | Every `P[i,j] ∈ [0,1]`; matrix rows sum ≤ 1 | Yes | Source 1, 2 |
| INV-3 | `P_unpaired(i) + Σ_j P(i,j) = 1` | Yes | Source 5 |
| INV-4 | `P(MFE structure) = exp((G − E_mfe)/RT)`, `E_mfe = CalculateFreeEnergy(MfeStructure)` | Yes | Source 4 |
| INV-5 | `Z`, all `P`, `MFE`, `P(MFE)` = exhaustive enumeration with the loop-by-loop structure evaluator | Yes | independent oracle |
| INV-6 | `CalculateStructureProbability(E, E) = 1`, `= exp(−ΔE/RT)`, no overflow for |E| ≫ RT | Yes | Source 1, 4 |
| INV-7 | `GenerateRandomRna` deterministic for a fixed seed; alphabet A/C/G/U | Yes | convention |

---

## 4. Test Cases

### 4.1 MUST Tests

| ID | Test Case | Expected Outcome | Test (`RnaSecondaryStructure_PartitionFunction_Tests`) |
|----|-----------|------------------|--------------------------------------------------------|
| M1 | no admissible pair / span ≤ 3 (`AAAA`, `GC`, `GAAC`, `G`) | Z = 1, all unpaired | `CalculatePartitionFunction_NoAdmissiblePair_ReturnsOne` |
| M3 | exhaustive enumeration, 28 sequences × 3 temperatures | Z, G, MFE, every P within 1e-9 | `CalculatePartitionFunction_EqualsExhaustiveEnumeration` |
| M5 | single-pair `GAAAAC` | Z = 1 + w, P = w/(1+w) | `..._SingleAdmissiblePair_ProbabilityIsBoltzmann` |
| M6 | ViennaRNA 2.7.2 locked values (17 cases: 37/20/60 °C; hairpin, bulge/interior, multiloop, tRNA-like, AU-rich) | G ±0.02, MFE, P(MFE), diversity, unpaired, top pairs, centroid | `CalculatePartitionFunction_MatchesViennaRna` |
| M6b | sanity `GGGAAAACCC` | `(((....)))`, MFE −2.02, G −2.1529 | `..._SanityHairpin_GGGAAAACCC` |
| M6d | normalisation for 40 random sequences 10–120 nt | INV-2, INV-3, INV-1 | `..._Normalisation_HoldsForRandomSequences` |
| M7 | P(MFE) formula and MFE structure evaluation | INV-4 | `..._MfeProbability_IsBoltzmannFactorOfMfeStructure` |
| M8/M9 | Boltzmann probability | 1.0; `exp(−1/RT)` = 0.19735 | `CalculateStructureProbability_*` |
| M10 | `null` | `ArgumentNullException` | `..._NullSequence_Throws` |
| M11 | empty | Z = 1, G = 0, no pairs | `..._EmptySequence_ReturnsOneNoPairs` |
| M12 | seeded determinism | identical strings | `GenerateRandomRna_SameSeed_IsDeterministicAndValidAlphabet` |
| M13 | diversity / centroid definitions | by definition | `..._DiversityAndCentroid_FollowTheirDefinitions` |
| M14 | T / lower case equal to RNA | equal Z | `..._TAndLowerCase_AreReadAsUpperCaseRna` |
| M15 | N / IUPAC never pair | no pair at N, Z = 1 for all-N | `..._UnknownBases_NeverPair` |
| M16 | `minLoopSize` clamp / honour | clamp to 3; no loop < 5 for 5 | `..._MinLoopSize_ClampedAndHonoured` |
| M17 | invalid temperature | `ArgumentOutOfRangeException` (0, −5, NaN, ∞) | `..._InvalidTemperature_Throws` |
| M18 | rescaling formula | `ΔG(T) = ΔH − (ΔH − ΔG₃₇)·T/310.15` | `TurnerModel_Stack_FollowsEnthalpyRescaling` |
| M19 | free energy at 20/60 °C | ViennaRNA eval −3.39 / −0.14 | `CalculateFreeEnergy_AtOtherTemperatures_MatchesViennaRna` |
| M20 | long sequence n = 500 | no NaN/overflow, normalised | `..._Length500_IsStableAndNormalised` |

### 4.2 SHOULD Tests

| ID | Test Case | Expected | Test |
|----|-----------|----------|------|
| S1 | `CalculateUnpairedProbabilities` = same engine | equal fields | `CalculateUnpairedProbabilities_EqualsPartitionFunctionFields` |
| S2 | G increases with T for a stable hairpin; default T = 37 °C | monotone | `..._FreeEnergy_IncreasesWithTemperature_AndDefaultIs37C` |
| S3 | extreme temperatures 1 mK … 10⁶ K | finite G, no NaN | `..._ExtremeTemperatures_AreFinite` |
| S4 | large energies in `CalculateStructureProbability` | no overflow | `CalculateStructureProbability_LargeEnergies_NoOverflow` |
| S5 | `GenerateRandomRna` boundaries / null `Random` / default overload | length, alphabet, exception | `GenerateRandomRna_*` |

### 4.3 Heavy tier (Fuzzing / Properties / Metamorphic / Combinatorial / Algebraic)

n = 1000 and n = 2000 stability and speed (`RnaPartitionFuzzTests`), temperature edges, validation, longer-stem monotonicity, worked example,
random-fold invariants; combinatorial sequence × temperature grid (`RnaStructureCombinatorialTests`).

---

## 5. Audit of Existing Tests

The 2026-06 tests locked the **fixed-energy-per-pair** model (Z counts 16 / 20, `E_bp` monotonicity, hand-enumerated P). That model was removed
(Evidence §8, finding F-PF1); its combinatorial content is retained by M3 (exhaustive enumeration on the Turner model), INV-1, INV-2.

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M20, S1–S5 | ✅ Covered | `Unit/Analysis/RnaSecondaryStructure_PartitionFunction_Tests.cs` |
| `CalculateUnpairedProbabilities` analytic pins | ✅ Re-locked to GASCONST 1.98717 | `RnaSecondaryStructure_UnpairedProbabilities_Tests.cs` |

---

## 6. Test Data

ViennaRNA 2.7.2 values: Evidence §9 (script `pfcmp3.py`, scratchpad). Enumeration oracle: in-test (`BruteForce`).

---

## 7. Open Questions / Decisions

- ViennaRNA default parameter file vs NNDB-precision: NNDB values used (Evidence §9.3); default-file residual ≈ 0.1 kcal/mol mean.
- `T ≠ 37 °C` MFE is unrounded (ViennaRNA truncates): tolerance 0.15 in the locked cases.
