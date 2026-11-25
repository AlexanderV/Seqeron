# RNA Pseudoknot Structure Prediction (canonical simple recursive pseudoknots, pknotsRG class)

| Field | Value |
|-------|-------|
| Algorithm Group | RnaStructure |
| Test Unit ID | RNA-PKPREDICT-001 (single knot) — shares the engine of [RNA-PKRECURSIVE-001](Pseudoknot_Prediction_Recursive.md) |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Implemented — exact optimum of the pknotsRG csr-PK search space over the shared Turner 2004 model (energy parameters differ from pknotsRG's, see §5.3) |
| Last Reviewed | 2026-11-25 |

## 1. Overview

`PredictStructurePseudoknot` returns the minimum-free-energy structure that may contain **at most one canonical H-type
pseudoknot**, the knot standing in the exterior loop, its three loops folding pseudoknot-free, the rest of the sequence
(before, after, inside the loops) folding with the shared Turner 2004 engine. It is the exact optimum of that search space
(no heuristic, no enumeration cut-offs): the pknotsRG dynamic programme of Reeder & Giegerich [1][2] with the knot count
limited to one. `PredictStructurePseudoknotRecursive` is the same engine without the limit (all knots, knots in
multiloops and in knot loops). Both are never worse than `CalculateMfeStructure`.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A pseudoknot is a crossing arrangement of base pairs: (i,j) and (k,l) cross iff i < k < j < l [5]. The H-type knot is two
crossing helices with three loops. pseudoknot-free folders (Zuker–Stiegler) cannot produce it; general pseudoknot DPs are
O(n⁶) (Rivas & Eddy 1999) and pknotsRG restricts the class so that O(n⁴) time / O(n²) space suffice [1].

### 2.2 Search space (Reeder & Giegerich 2004 [1], pknotsRG-mfe.lhs `pknot` [2])

```
knot = a · u · b · v · a' · w · b'        (5' → 3', a pairs a', b pairs b')
```

Boundaries i (start of a), k (end of a'), l (start of b), j (end of b'); α = |a|, β = |b|.

| Rule | Statement | Implementation |
|------|-----------|----------------|
| Rule 1 | both strands of a helix have equal length, no bulges (class csr-PK; the relaxed class with a 1-nt bulge is not included, §5.3) | helices are chains of consecutive pairs (i+t, k−t), (l+t, j−t) |
| Rule 2 | helices have maximal extent from their outer pair | α = length of the maximal pair chain from (i,k); β from (l,j) |
| Rule 3 | where two maximal helices compete for the bases of v the boundary is fixed arbitrarily | b is shortened (pknotsRG `getbeta`): β ≤ k − α − l + 1 (v ≥ 0) and β ≤ j − k − 2 (w ≥ 2) |
| size | each helix ≥ 2 pairs | α, β ≥ 2 |
| bridging | one nt after a and two nt before b' stay unpaired (stereochemistry; `front` / `back`) | u' = [i+α+1, l−1] after the explicit unpaired nt i+α; w' = [k+1, j−β−2] before the explicit nts j−β−1, j−β; v = [l+β, k−α] may be empty |

Loops u', v, w' fold as **segments**: a sequence of closed stems and unpaired nucleotides (pknotsRG `comps`); in the
recursive variant a segment component may also be a knot. The four free boundaries i, l, k, j give O(n⁴) candidates; the
loop segments are table look-ups.

### 2.3 Energy model

Every nested term is the shared `TurnerModel` (Turner 2004 NNDB-precision parameters, ViennaRNA `-d2`, temperature
rescaling). The pknotsRG pseudoknot parameters are used unchanged [2]:

```
ΔG(knot) = 9.0                                         pkinit  (Energy.lhs: pkinit = 900)
         + 3 · 0.3                                     the three bridging unpaired nt (pk' : 3*npp)
         + Σ stack(a) + Σ stack(b)                     Turner 2004 stacks incl. the GGUC/CUGG three-stack; 0.0 per pair (pbp = 0)
         + C(u') + C(v) + C(w')                        knot-loop segments
         + Ext(i,k) + Ext(l,j)                         d2 exterior-stem term of the outer pairs of a and b (terminal AU/GU + mismatch)
         + Inner(a5last, a3first) + Inner(b5last, b3first)   d2 mismatch + terminal AU/GU of the inner pairs read from the enclosed loop
C(x..y)  = best sequence of components, 0.3 per unpaired nt (npp = 30),
           stem (p,q): V(p,q) + Ext(p,q);   knot (recursive variant only): K(p,q) + 6.0   (pkml = 600)
```

`V(p,q)` is the Zuker closed-pair energy of the shared engine (hairpin / interior / multiloop). The knot is a component of
the exterior loop at cost K, of a multiloop at cost K + 6.0 (instead of the branch term; recursive variant only).
`PredictStructurePseudoknot`: at most one knot; knot loops and multiloops are knot-free.

### 2.4 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-PK-01 | 1 M NaCl, 37 °C default; Turner terms rescale with `temperature`, the pseudoknot constants (9.0 / 0.3 / 6.0) do not | knot/no-knot decision shifts with conditions not modelled |
| ASM-PK-02 | secondary-structure thermodynamics only; tertiary interactions (triplexes, Mg²⁺) are not modelled | tertiary-stabilised knots (BWYV) need not be the model optimum |
| ASM-PK-03 | canonical class csr-PK only (§2.2) | kissing hairpins, bulged/unequal helices, chained helices are not generated |

### 2.5 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-PK-01 | ΔG ≤ plain MFE | the pseudoknot-free options are part of every minimisation |
| INV-PK-02 | `HasPseudoknot` ⇔ the returned pair set has a crossing pair | a knot's helices cross by construction (a5 < b5 < a3 < b3); a knot-free result is the pseudoknot-free traceback |
| INV-PK-03 | every position pairs at most once; hairpins ≥ 3 nt | segments and helices occupy disjoint index runs |
| INV-PK-04 | no knot unless it lowers ΔG | the traceback prefers the pseudoknot-free option on ties; result = plain MFE structure then |
| INV-PK-05 | recursive ≤ single ≤ MFE | each is the optimum of a larger space |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| rnaSequence | string | required | RNA or DNA (T read as U), case-insensitive; other characters never pair |
| minLoopSize | int | 3 | minimum hairpin loop (values < 3 are clamped to 3) |
| temperature | double | 310.15 | kelvin; Turner terms rescaled with their enthalpies |

### 3.2 Output

`PseudoknotStructure(Sequence, DotBracket, BasePairs, FreeEnergy, HasPseudoknot)`; `DotBracket` is rendered by `ToDotBracket`
(ViennaRNA `vrna_db_from_ptable` families: crossing helices on `[]`, `{}`, …); `BasePairs` sorted by 5' position; `FreeEnergy`
rounded to 0.01 at 37 °C. null / empty input → empty structure (ΔG 0); n < 11 → plain MFE (no knot fits).

## 4. Algorithm

The knot layer (`KnotLayer` in `RnaSecondaryStructure.Pseudoknot.cs`) is hooked into the shared `FillDp` / `TracebackMfe`
(no second energy engine): for each interval [i,j], in increasing span, after V(i,j):

1. `K(i,j)` = min over (k, l) of the formula of §2.3 (O(n²) per interval, helix chains and helix energies precomputed once, O(n²) memory each).
2. `C(i,j)` = knot-loop segment value (first component, then the rest).
3. WM (multiloop) gets the option K + 6.0 and the exterior W the option W(p−1) + K(p,j) (recursive variant); the single-knot
   variant keeps a second exterior prefix array "at most one knot".
4. Traceback re-derives the minimising option (K, C items), adds the helix pairs, recurses into the segments.

Time O(n⁴) (four free boundaries; measured ≈ 0.25 s per fold at n = 200 in a Debug build, 0.14 s for n = 80 versus 6.9 s of the
previous interval recursion), space O(n²).

## 5. Implementation Notes

### 5.1 Location

[RnaSecondaryStructure.Pseudoknot.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.Pseudoknot.cs)
(`PredictStructurePseudoknot`, `PredictStructurePseudoknotRecursive`, `KnotLayer`); hooks in `FillDp` / `TracebackMfe`
(`RnaSecondaryStructure.cs`); `TurnerModel.InnerStemEnd` (helix end seen from the loop it encloses). MCP: `rna_pseudoknot_structure`.

### 5.2 Verification

Independent oracle (Python, ViennaRNA 2.7.2 `E_Hairpin` / `E_IntLoop` / `E_MLstem` / `E_ExtLoop` with the NNDB-precision parameter file, naive
memoised grammar, explicit enumeration of every knot): 1 971 sequences (357 + 459 + 452 random/designed ≤ 34 nt, 809 designed nesting cases
≤ 45 nt, 2 multiloop-knot cases of 62–64 nt) compared in both modes: 0 energy differences (to 0.00); 2 differences, both explained by NNDB-only terms ViennaRNA
lacks (1-nt bulge degeneracy, special-C bulge inside a loop segment); 106 sequences excluded because their plain MFE already contains such a term.
In-tree: 70-sequence exhaustive table-free enumeration (all flank structures × all canonical knots × all loop structures) = DP to 0.006.

### 5.3 Conformance to pknotsRG and deviations

**Implemented as in pknotsRG [1][2]:** grammar `knot`, canonization rules 1–3 (csr-PK), pkinit 9.0, npp 0.3 per unpaired knot-loop nt,
pbp 0.0, pkml 6.0, bridging nucleotides (1 + 2), v may be empty, knots in exterior / multiloop / knot loops, exactness (the DP minimises over the whole space).

**Deviations (all documented, none silent):**

| # | pknotsRG | Here | Reason |
|---|----------|------|--------|
| 1 | Turner 1999 loop/stack parameters (Energy.lhs, Intloop*.lhs) | shared Turner 2004 `TurnerModel` | campaign rule: one energy model; numbers are therefore not comparable with pknotsRG output |
| 2 | stem component in a knot loop: +0.4 (`ul`) + dangles | d2 exterior-stem term, no branch term | the 2004 model has no 0.4 branch term for an unclosed loop; ViennaRNA d2 treats every helix end alike |
| 3 | dangles on the bridging nucleotides (`dangles`) | d2 mismatch + terminal AU on the four helix ends | same, uniform d2 |
| 4 | coaxial stacking term `emptymid` when v = ∅ | none | neither the shared model nor the ViennaRNA default has coaxial stacking |
| 5 | pknotsRG-mfe program also admits a single 1-nt bulge in a helix ("relaxed csr-PK", 11 of 212 PseudoBase knots) | strict rule 1(b) | not implemented: needs bulge-position enumeration and the seven `getbeta` overlap cases; recorded as follow-up |
| 6 | helix = minimum-energy stack among the non-bulged alternatives | maximal helix (paper rule 2) | identical unless an additional pair raises the energy (positive GU stacks) |
| 7 | knots only as exterior / multiloop / knot-loop components | same; a knot directly enclosed by one helix (interior loop) is not generated | not in the grammar (the previous recursive code did generate it, see Recursive doc) |

### 5.4 Limitations

Tertiary interactions are not modelled (BWYV, the 28-nt frameshift knot, is predicted with S1 = 5 bp, S2 = 4 bp, ΔG −10.51 vs −7.90 for the
pseudoknot-free MFE; the published crystal pairing was not compared). Knot classes outside csr-PK (§ASM-PK-03) cannot be returned.
Sequences ≳ 500 nt are slow (O(n⁴)).

## 6. Edge Cases

| Case | Behaviour |
|------|-----------|
| null / empty | empty structure, ΔG 0 |
| < 11 nt | plain MFE (shortest knot: 2·2 + 2·2 + 1 + 0 + 2 = 11) |
| strong hairpin / no profitable knot | HasPseudoknot = false, structure and ΔG = plain MFE |
| DNA T, lower case | identical to RNA spelling |
| N / IUPAC | never pair |

## 7. Examples

```csharp
var pk = RnaSecondaryStructure.PredictStructurePseudoknot("GGGGAACCCCAACCCCAAGGGG");
// pk.DotBracket == "((((..[[[[..))))..]]]]", pk.FreeEnergy == -13.56 (plain MFE -7.84)
var flanked = RnaSecondaryStructure.PredictStructurePseudoknot("GGGCGAAAGCCCAAGGGGAACCCCAACCCCAAGGGG");
// "((((....))))..((((..[[[[..))))..]]]]", -21.70  (hairpin + knot; before the fix: no knot, -15.98)
```

Related: [Minimum_Free_Energy](Minimum_Free_Energy.md), [Pseudoknot_Detection](Pseudoknot_Detection.md), [Dot_Bracket_Notation](Dot_Bracket_Notation.md),
tests [RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs),
evidence [RNA-PKPREDICT-001-Evidence.md](../../Evidence/RNA-PKPREDICT-001-Evidence.md).

## 8. References

1. Reeder J, Giegerich R. (2004). Design, implementation and evaluation of a practical pseudoknot folding algorithm based on thermodynamics. BMC Bioinformatics 5:104. https://doi.org/10.1186/1471-2105-5-104
2. pknotsRG source: README, Foldingspace.lhs, Energy.lhs (pkinit 900, npp 30, pbp 0), Algebras.lhs (`energyalgext`: pkml 600, `pk'`), pknotsRG-mfe.lhs (`pknot`, `getbeta`). https://github.com/jensreeder/pknotsRG
3. Rivas E, Eddy SR. (1999). A dynamic programming algorithm for RNA structure prediction including pseudoknots. J Mol Biol 285:2053–2068.
4. Su L, Chen L, Egli M, Berger JM, Rich A. (1999). Minor groove RNA triplex in the crystal structure of a ribosomal frameshifting viral pseudoknot. Nat Struct Biol 6:285–292 (PDB 437D).
5. Antczak M, et al. (2018). New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation. Bioinformatics 34:1304–1312.
