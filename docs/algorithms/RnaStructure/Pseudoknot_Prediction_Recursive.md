# RNA Pseudoknot Prediction — Recursive pknotsRG Grammar (several / nested knots)

| Field | Value |
|-------|-------|
| Algorithm Group | RnaStructure |
| Test Unit ID | RNA-PKRECURSIVE-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Implemented — exact optimum of the full pknotsRG csr-PK grammar over the shared Turner 2004 model |
| Last Reviewed | 2026-11-25 |

## 1. Overview

`PredictStructurePseudoknotRecursive` folds the whole sequence with the **recursive** pknotsRG grammar [1][2]: several knots side
by side in the exterior loop, knots as branches of multiloops, and knots inside the loops u, v, w of other knots
("the unpaired strands u, v, w … fold internally in an arbitrary way, including simple recursive pseudoknots" [1]; the
pseudoknot value "competes with values of unknotted foldings for the interval" [2]). It is the same engine, energy model,
canonization and verification as [Pseudoknot_Prediction](Pseudoknot_Prediction.md) — only the single-knot restriction is
removed — so everything there (formula §2.3, parameters, deviations §5.3, oracle §5.2) applies. Never worse than
`PredictStructurePseudoknot` and `CalculateMfeStructure`.

## 2. Formal Basis

Interval tables over [i,j] (increasing span): V (closed pair, shared Zuker engine), WM (multiloop part), W (exterior prefix), and the knot layer

```
K(i,j)   best canonical knot spanning exactly [i,j]                      (O(n²) candidates (k,l) per interval)
C(i,j)   best knot-loop segment of [i,j]:  0.3 + C(i+1,j)                 (i unpaired)
                                       |  V(i,q) + Ext(i,q) + C(q+1,j)    (stem component)
                                       |  K(i,q) + 6.0 + C(q+1,j)         (knot component — the recursion)
WM(i,j) ∋ K(i,j) + 6.0                                                   (knot as a multiloop branch, pknotsRG pkml)
W(j)    ∋ W(p−1) + K(p,j)                                                (knot in the exterior loop)
K(i,j)  = min_{k,l} 9.0 + 0.9 + stacks(a) + stacks(b) + ends + C(u') + C(v) + C(w')
```

Termination: every table entry depends on strictly smaller intervals (the loops of a knot are proper sub-intervals of [i,j]), so
the recursion is a finite DP, O(n⁴) time and O(n²) space. Grammar limits (class csr-PK): helices of equal length without bulges,
maximal, ≥ 2 pairs; kissing hairpins, chained helices and a knot directly enclosed by a single helix (interior loop) are outside the
grammar and are not generated.

## 3. Contract

Signature `PredictStructurePseudoknotRecursive(string rnaSequence, int minLoopSize = 3, double temperature = 310.15)` →
`PseudoknotStructure`; inputs, outputs, edge cases as in [Pseudoknot_Prediction §3, §6](Pseudoknot_Prediction.md).
`HasPseudoknot` ⇔ the pair set has a crossing pair. MCP: `rna_pseudoknot_structure` (`allowNestedKnots = true`).

## 4. Examples (oracle-verified, ViennaRNA-primitive implementation of the grammar)

| Sequence | single (≤ 1 knot) | recursive | what the recursion adds |
|----------|-------------------|-----------|-------------------------|
| `GGGGAACCCCAACCCCAAGGGGAAAGGGGAACCCCAACCCCAAGGGG` | −24.46 | **−27.52** `((((..[[[[..))))..]]]]...((((..[[[[..))))..]]]]` | two knots side by side |
| `CCGCUUGGAGUGCGGAAACUCCUUUCGGAAGCCGUCCGUACGGC` | −18.85 | **−22.38** | two knots of different size |
| `GCCAGCGGCCUUACGGGAUAGGCAUACCCGGGCAACGC` | −11.56 | **−12.15** `(((.[[[((((..[[[[..))))...]]]])))..]]]` | a knot in loop v of another knot |
| `CGCCCGGGCAAGGCGGGAAACCGCCAAGCGGGAGCCCCCCCCCGCCUGGGGGCAAGCCCGGGCG` | −43.46 (= MFE) | **−46.66** `(((((((((..(((((....)))))..(((((.[[[[[[.)))))..]]]]]]..)))))))))` | a knot as a multiloop branch (6.0) |

Among 809 designed sequences the recursive optimum is strictly below the single-knot optimum for 98; all 95 oracle-comparable ones agree to 0.00.

## 5. Implementation Notes

[RnaSecondaryStructure.Pseudoknot.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.Pseudoknot.cs) (`KnotLayer` with `NestedKnots = true`).
**Replaced:** the former `RecursivePkFolder` (interval memo calling `CalculateMfeStructure` per interval, anchored left knots, an
"enclosing helix over a knotted interior" production that is not part of the pknotsRG grammar, no multiloop knots, 6.9 s at 80 nt) and its
duplicate knot scorer/loop scorer; one engine now serves both predictors on the shared `FillDp` / `TracebackMfe`.

## 6. References

1. Reeder J, Giegerich R. (2004). BMC Bioinformatics 5:104. https://doi.org/10.1186/1471-2105-5-104
2. Reeder J, Steffen P, Giegerich R. (2007). pknotsRG: RNA pseudoknot folding including near-optimal structures and sliding windows. Nucleic Acids Res 35:W320–W324. https://doi.org/10.1093/nar/gkm258
3. pknotsRG source (github.com/jensreeder/pknotsRG): `pknotsRG-mfe.lhs` (`block`, `comps`, `knot`), `Algebras.lhs` (`pkml`).
4. Antczak M, et al. (2018). Bioinformatics 34:1304–1312 (crossing condition i < k < j < l).
