# Test Specification: RNA-PKRECURSIVE-001

**Test Unit ID:** RNA-PKRECURSIVE-001
**Area:** RnaStructure
**Algorithm:** Recursive pknotsRG grammar — several / nested canonical H-type knots (exact DP, shared Turner 2004 model)
**Status:** ☑ Complete (re-reviewed 2026-11-25, Stage A / Stage B)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | URL | Accessed |
|---|--------|---------------|-----|----------|
| 1 | Reeder & Giegerich (2004) BMC Bioinformatics 5:104 — loops u, v, w fold "including simple recursive pseudoknots" | 1 | https://doi.org/10.1186/1471-2105-5-104 | 2026-11-25 |
| 2 | Reeder, Steffen & Giegerich (2007) NAR 35:W320 — knot value competes per interval | 1 | https://doi.org/10.1093/nar/gkm258 | 2026-06-23 |
| 3 | pknotsRG source `pknotsRG-mfe.lhs` (`block`, `comps`, `knot`), `Algebras.lhs` (`pkml` = 600) | 3 | https://github.com/jensreeder/pknotsRG | 2026-11-25 |
| 4 | ViennaRNA 2.7.2 primitives (oracle) | 1 | PyPI `ViennaRNA` | 2026-11-25 |
| 5 | Antczak et al. (2018) crossing condition | 1 | https://doi.org/10.1093/bioinformatics/btx783 | 2026-06-23 |

### 1.2 Key Evidence Points

1. Knots may appear side by side (exterior loop), as multiloop branches (+6.0), and as components of knot loops (+6.0) — Sources 1–3.
2. Same energy model, canonization and exactness as RNA-PKPREDICT-001; one engine serves both predictors.
3. Oracle (naive grammar on ViennaRNA primitives incl. its own multiloop-knot V recursion): 0 energy differences over 1 971 comparable sequences; recursive < single for 98 of 809 designed sequences.

### 1.3 Corner Cases / Pitfalls

- A single helix cannot directly enclose a knot (not in the grammar); the former "over-arching helix" production is gone (its datasets tested a non-pknotsRG class).
- Recursive = single-knot optimum whenever no second / nested knot fits (< ≈ 22 nt) — asserted against the exhaustive enumeration.
- Two strong knots are not always the optimum (cross-region nested folds can win) — engineered cases only.

---

## 2. Canonical Methods Under Test

| Method | Type |
|--------|------|
| `PredictStructurePseudoknotRecursive(string, int, double)` | Canonical |
| `rna_pseudoknot_structure` (MCP, `allowNestedKnots=true`) | Wrapper |

---

## 3. Invariants

INV-1 recursive ≤ single ≤ MFE; INV-2 `HasPseudoknot` ⇔ crossing pair; INV-3 valid structure (each position ≤ 1 pair, hairpins ≥ 3, dot-bracket round trip); INV-4 no knot ⇒ plain MFE structure; INV-5 recursion terminates and is deterministic / thread safe.

---

## 4. Test Cases (`Unit/Analysis/RnaSecondaryStructure_PredictStructurePseudoknotRecursive_Tests.cs`)

| ID | Test | Expected |
|----|------|----------|
| T1 | oracle rows (29) | recursive energy and dot-bracket |
| M1 | twin designed knots (47 nt) | `((((..[[[[..))))..]]]]...((((..[[[[..))))..]]]]`, −27.52 (single −24.46) |
| M2 | two knots of different size (44 nt) | −22.38 (single −18.85) |
| M3 | knot inside knot loop v (38 nt) | `(((.[[[((((..[[[[..))))...]]]])))..]]]`, −12.15 (single −11.56) |
| M4 | knot as multiloop branch (64 nt) | −46.66 (single = MFE −43.46) |
| M5 | hairpin | no knot |
| E1 | recursive = exhaustive single-knot optimum, 11–16 nt (40) | equal |
| I1 | random sweep (60, 11–55 nt) | recursive ≤ single ≤ MFE, valid, no knot ⇒ MFE structure |
| E2 | null/empty/short, DNA + lower case, minLoopSize 0 → 3, temperature, concurrency, 150 nt | as RNA-PKPREDICT-001 |

---

## 5. Audit of Existing Tests

Replaced the previous file (over-arching A·U helix datasets, `RecursivePkFolder` internals). Property / Metamorphic / Combinatorial / Fuzzing files updated where they used the old nested dataset.
