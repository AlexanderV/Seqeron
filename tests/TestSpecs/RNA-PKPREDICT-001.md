# Test Specification: RNA-PKPREDICT-001

**Test Unit ID:** RNA-PKPREDICT-001
**Area:** RnaStructure
**Algorithm:** Pseudoknot structure prediction — at most one canonical H-type knot (pknotsRG class csr-PK, exact DP over the shared Turner 2004 model)
**Status:** ☑ Complete (re-reviewed 2026-11-25, Stage A / Stage B)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| 1 | Reeder & Giegerich (2004), BMC Bioinformatics 5:104 — canonical simple recursive pseudoknots, rules 1–3, O(n⁴) | 1 | https://doi.org/10.1186/1471-2105-5-104 | 2026-11-25 |
| 2 | pknotsRG source: `pknotsRG-mfe.lhs` (`pknot`, `getbeta`), `Energy.lhs` (pkinit 900, npp 30, pbp 0), `Algebras.lhs` (`pk'`, `pkml` 600) | 3 | https://github.com/jensreeder/pknotsRG | 2026-11-25 |
| 3 | ViennaRNA 2.7.2 `RNA.E_Hairpin/E_IntLoop/E_MLstem/E_ExtLoop/fold` + NNDB-precision parameter file (oracle primitives) | 1 | PyPI sdist `ViennaRNA` | 2026-11-25 |
| 4 | Antczak et al. (2018) — crossing i<k<j<l | 1 | https://doi.org/10.1093/bioinformatics/btx783 | 2026-06-23 |
| 5 | Su et al. (1999) Nat Struct Biol 6:285 / PDB 437D (BWYV sequence) | 1 | https://www.rcsb.org/structure/437D | 2026-06-23 |

### 1.2 Key Evidence Points

1. Knot = a·u·b·v·a'·w·b'; equal-length bulge-free maximal helices (≥ 2 pairs); rule 3 shortens b where helices compete for v; one nt after a and two nt before b' unpaired; v may be empty — Sources 1, 2.
2. Energy: 9.0 + 3·0.3 + stacks + 0.3/unpaired loop nt + helix ends (shared d2 model); 6.0 (`pkml`) for a knot as a multiloop / knot-loop component — Source 2.
3. The result is the exact optimum of that space; flanks and loops fold with the shared Turner 2004 engine; never above the plain MFE — Sources 1, 2.
4. Oracle: independent grammar implementation on ViennaRNA primitives agrees to 0.00 on 1 971 sequences; exhaustive table-free enumeration agrees on 70 small sequences (Evidence §Stage B).

### 1.3 Documented Corner Cases

- Knot only when strictly better than the pseudoknot-free optimum (tie → plain MFE structure).
- Sequences with NNDB-only terms ViennaRNA lacks (special C bulge, 1-nt bulge degeneracy, GGUC stack, special GU closure, all-C loops) are not oracle-comparable (excluded from the lock rows).
- BWYV: model optimum contains a knot (−10.51 vs −7.90); not claimed equal to the crystal pairing.

### 1.4 Known Failure Modes / Pitfalls (fixed 2026-11-25)

1. Flanks 5' of a / 3' of b' left unfolded and unscored → knots with a flanking hairpin lost (`GGGCGAAAGCCC…` −15.98 → −21.70).
2. Knot energy without bridging nt / helix ends / loop context (designed H-type −8.76 → −13.56).
3. v ≥ 1 and dropped overlap candidates (`UCAUGUUAGCUCAUGAAUGGCU` −4.22 → −5.93).

---

## 2. Canonical Methods Under Test

| Method | Type | Notes |
|--------|------|-------|
| `PredictStructurePseudoknot(string, int, double)` | Canonical | ≤ 1 knot in the exterior loop; engine shared with RNA-PKRECURSIVE-001 |
| `rna_pseudoknot_structure` (MCP, `allowNestedKnots=false`) | Wrapper | `RnaPseudoknotStructureTests` |

---

## 3. Invariants

| ID | Invariant |
|----|-----------|
| INV-PK-01 | ΔG ≤ plain MFE |
| INV-PK-02 | HasPseudoknot ⇔ crossing pair in the pair set |
| INV-PK-03 | each position paired ≤ once, hairpins ≥ 3, dot-bracket round trip |
| INV-PK-04 | no knot ⇒ structure and ΔG = plain MFE |
| INV-PK-05 | recursive ≤ single ≤ MFE |
| INV-PK-06 | DP = exhaustive enumeration of the search space |

---

## 4. Test Cases (`Unit/Analysis/RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs`)

| ID | Test | Expected |
|----|------|----------|
| M1 | designed H-type `GGGGAACCCCAACCCCAAGGGG` | pairs (0,15)…(9,18)…, `((((..[[[[..))))..]]]]`, −13.56 (MFE −7.84) |
| M2 | knot strictly beats plain MFE | ΔG < MFE − 1 |
| M3 | hairpin + knot `GGGCGAAAGCCCAAGGGGAACCCCAACCCCAAGGGG` | `((((....))))..((((..[[[[..))))..]]]]`, −21.70 (MFE −15.98) |
| M4 | empty middle loop `UCAUGUUAGCUCAUGAAUGGCU` | `(((((..[[[[)))))..]]]]`, −5.93 (MFE −4.22) |
| M5 | pknotsRG README sequence, BWYV | −12.31, −10.51 (oracle) |
| T1 | oracle rows (29) | energy, dot-bracket, no-knot ⇒ MFE structure, validity |
| T2 | exhaustive enumeration, random + H-type-biased 11–16 nt (70) and 4 designed | DP = enumeration (0.006) |
| I1 | random sweep (80, 11–40 nt) | ΔG ≤ MFE, no knot ⇒ MFE structure, valid, recursive ≤ single |
| I2 | two knots possible | single places exactly one (4 `[`), recursive two (8) |
| E1 | DNA / lower case; minLoopSize 0 → 3; temperature (310.15 = default, 298.15 lower); parallel calls; 150 nt | parity / clamp / monotone / deterministic / completes |
| E2 | null, empty, < 11 nt, hairpin, N / IUPAC | empty / plain MFE / valid |

MCP: `Seqeron.Mcp.Analysis.Tests/RnaPseudoknotStructureTests.cs` (schema + binding), contract list updated.

---

## 5. Audit of Existing Tests

Replaced `RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs` (old expectations: knot-only energy −8.76, BWYV "no knot", "single method cannot place a second knot"). Heavy-tier files
(Fuzzing `RnaPseudoknotPredictFuzzTests`, Properties, Metamorphic, Combinatorial) kept; their nested-knot datasets updated to oracle-verified grammar cases. Shared helpers:
`PknotsRgBruteForce.cs`, `PknotsRgOracleData.cs`.
