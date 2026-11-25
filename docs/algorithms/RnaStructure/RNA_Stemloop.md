# RNA Stem-Loop Detection

| Field | Value |
|-------|-------|
| Algorithm Group | RNA Structure |
| Test Unit ID | RNA-STEMLOOP-001 |
| Related Projects | N/A |
| Implementation Status | Complete (perfect-stem enumeration + opt-in bulge/internal-loop stems) |
| Last Reviewed | 2026-11-25 (B12 RNA-STEMLOOP-001) |

## 1. Overview

Stem-loop detection scans a single RNA sequence for hairpin motifs formed by a paired stem closing an unpaired terminal loop. It enumerates every hairpin loop of the allowed size closed by a canonical pair and reports the stem that closes it: by default a perfect helix grown outward while the next pair can form; optionally (`maxInternalLoopSize > 0`) the minimum-free-energy *unbranched* stem, which may contain bulges and internal loops. Every candidate is scored as an isolated structure with the same Turner 2004 loop functions as the MFE folder (ViennaRNA `eval_structure` counterpart). It is a motif enumerator, not a global structure optimizer.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Stem-loops consist of a paired stem and a terminal loop. In RNA, standard admissible pairs include Watson-Crick `A-U` and `G-C` pairs plus wobble `G-U` pairs. Hairpin loops shorter than three nucleotides are sterically disfavored and are treated as impossible in standard RNA secondary-structure models.

### 2.2 Core Model

A stem-loop is a hairpin loop L = [p+1 .. q−1] (|L| in [minLoopSize, maxLoopSize], (p,q) a canonical pair) plus the chain of nested pairs closing it.

* **Perfect stem** (`maxInternalLoopSize = 0`, default): the NNDB *helix* — consecutive stacked pairs (p,q), (p−1,q+1), … grown outward while the next outer pair can form (Wikipedia *Stem-loop*: complementary regions "base-pair to form a double helix that ends in an unpaired loop"). Dropped if shorter than `minStemLength`.
* **Unbranched stem** (`maxInternalLoopSize = M > 0`): the minimum-free-energy chain of pairs from any outer pair (i,j) down to (p,q) whose consecutive pairs enclose only degree-2 loops (stack, bulge, internal loop) with n1 + n2 ≤ M — Zuker & Stiegler's V(i,j) recursion without the multiloop branch, ViennaRNA MAXLOOP semantics. Computed once per sequence by an "outside" DP

$$
G(k,l) = \min\Big\{\,\mathrm{AU/GU}(k,l),\ \min_{(i,j)} G(i,j) + I(i,j;k,l),\ G(k-3,l+3) + \mathrm{GGUC}\Big\}
$$

  (I = the shared `InteriorTerm`: Turner 2004 stacks, bulges, int11/int21/int22 and generic internal loops); the stem of hairpin (p,q) is the argmin chain of G(p,q) (the hairpin term is fixed by (p,q)). Dropped if it has fewer than `minStemLength` pairs.

`TotalFreeEnergy` is the free energy of the motif as an isolated structure — the sub-sequence [Start..End] with its own dot-bracket — evaluated by the shared structure evaluator (`EvaluateStructureEnergy`, the ViennaRNA `eval_structure` counterpart): stacks, bulge/internal loops, hairpin loop (`HairpinTerm`) and the terminal AU/GU penalty of both helix ends (no dangles: the motif has no neighbours). For a perfect stem this is the NNDB worked-example sum ΔG = ΔG(helix) + ΔG(hairpin loop), with `Stem.FreeEnergy` = `CalculateStemEnergy` (stacks + per-end AU/GU penalty), except for the four A-U-closed special hexaloops whose measured value already contains the closure penalty (ViennaRNA `vrna_E_hairpin` returns `Hexaloop_E`). For an unbranched stem `Stem.FreeEnergy` = total − hairpin-loop term.

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Stem-loop motifs can be identified by local complementary extension around a loop candidate | More global structural context is ignored during motif detection |
| ASM-02 | Wobble pairing is acceptable when `allowWobble = true` | Candidate sets change when wobble pairs are biologically inappropriate for the use case |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Returned stem-loops have `Loop.Size >= 3` | `minLoopSize` is clamped to `3` before scanning |
| INV-02 | Returned stem-loops satisfy `Stem.Length >= minStemLength` | Candidates shorter than the threshold are discarded |
| INV-03 | Returned loop type is always `Hairpin` | The constructor uses `LoopType.Hairpin` for this API |
| INV-04 | `DotBracketNotation.Length = End - Start + 1` | The dot-bracket helper allocates exactly the detected span length |
| INV-05 | `TotalFreeEnergy = EvaluateStructureEnergy(seq[Start..End], DotBracketNotation)` | The total is computed by the shared evaluator |
| INV-06 | With `maxInternalLoopSize = M > 0`, `TotalFreeEnergy` ≤ the energy of any other unbranched stem closing the same hairpin loop | Outside DP minimum (exact; checked vs ViennaRNA loop functions) |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `rnaSequence` | `string` | required | RNA sequence to analyze | Empty or too-short input yields no results |
| `minStemLength` | `int` | `3` | Minimum number of paired bases in a returned stem | Candidates below this threshold are rejected |
| `minLoopSize` | `int` | `3` | Minimum hairpin loop size | Values below `3` are clamped to `3` |
| `maxLoopSize` | `int` | `10` | Maximum scanned loop size | Restricts candidate enumeration |
| `allowWobble` | `bool` | `true` | Whether `G-U` / `U-G` wobble pairs may be stem pairs | When `false`, wobble pairs stop extension / are excluded from chains |
| `maxInternalLoopSize` | `int` | `0` | 0 = perfect stems; M > 0 = MFE unbranched stems with bulges/internal loops of ≤ M unpaired nt | Negative → `ArgumentOutOfRangeException` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Start` | `int` | Inclusive 5' start of the detected stem-loop span |
| `End` | `int` | Inclusive 3' end of the detected stem-loop span |
| `Stem` | `Stem` | Stem coordinates, base pairs, and stem free energy |
| `Loop` | `Loop` | Hairpin-loop coordinates, size, and sequence |
| `TotalFreeEnergy` | `double` | Sum of stem and loop free energies |
| `DotBracketNotation` | `string` | Local dot-bracket string covering only the detected stem-loop span |

### 3.3 Preconditions and Validation

The input sequence is uppercased internally and DNA `T` is read as `U`; any other symbol (N, IUPAC codes, gaps) never pairs (ViennaRNA default pair set) and, as a hairpin first-mismatch base, takes ViennaRNA's `N` value (maximum over A/C/G/U of the mismatch table). If `minLoopSize < 3`, it is clamped to `3` (NNDB: hairpin loops < 3 nt are prohibited). When the sequence is empty or shorter than `minStemLength * 2 + minLoopSize`, the method yields no candidates. Stem extension stops at the first invalid base pair, or at the first wobble pair when `allowWobble` is `false`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Clamp `minLoopSize` to `3` and reject too-short sequences.
2. Enumerate loop-start positions across the sequence.
3. Enumerate loop sizes between `minLoopSize` and `maxLoopSize` for each start.
4. Perfect stems: extend the helix outward while the next pair can form. Unbranched stems: trace the argmin chain of the outside DP G(p,q) (filled once, O(n²·M²)).
5. Discard candidates whose stem has fewer than `minStemLength` pairs.
6. Build `Stem`, `Loop`, local dot-bracket and `TotalFreeEnergy` = `EvaluateStructureEnergy(motif, dot-bracket)`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Pairs come from the shared `PairLookup` table (`CanPair`/`GetBasePairType`). Base pairs are listed outermost (5') first. All energies come from the shared Turner 2004 functions of the MFE engine (`HairpinTerm`, `InteriorTerm`, GGUC/CUGG three-stack, AU/GU end penalty); nothing is re-implemented.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindStemLoops` (perfect) | `O(n · L · s)` | `O(k)` | `n` = length, `L` = loop-size range, `s` = motif span (extension + evaluation), `k` = results |
| `FindStemLoops` (`M > 0`) | `O(n² · M²)` | `O(n²)` | one outside DP per call, then `O(s)` per candidate |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.FindStemLoops(string, int, int, int, bool)`: Enumerates candidate stem-loop motifs.
- `RnaSecondaryStructure.DetectPseudoknots(IReadOnlyList<BasePair>)`: Post hoc crossing-pair detector for already assembled base-pair sets.

### 5.2 Current Behavior

The detector returns all candidates that satisfy the local constraints; it does not suppress overlaps or attempt to choose a single best motif. Each result includes both structural coordinates and a locally scored energy. Pseudoknots are not predicted directly from sequence in this method, but the same class exposes `DetectPseudoknots` for already assembled base-pair sets.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Minimum hairpin-loop size of three nucleotides (NNDB).
- ViennaRNA default pair set (A·U, G·C, G·U; T ≡ U), optional wobble.
- Turner 2004 energies identical to ViennaRNA 2.7.2 `eval_structure` on the isolated motif (NNDB-precision parameters), plus the sourced NNDB-only terms (special GU closure, all-C loops, GGUC/CUGG three-stack, 1-nt bulge degeneracy / special C) — 156,286 perfect-stem candidates and 35,155 unbranched-stem optima checked, 0 mismatches (Evidence §11).
- Stems with bulges and internal loops (opt-in) as the MFE unbranched chain (Zuker–Stiegler V without multiloops).

**Design scope (not a simplification):**

- Every hairpin loop yields one candidate and overlapping candidates are all returned (enumeration, by definition); a non-overlapping / globally optimal set is the job of `CalculateMfeStructure` / `PredictStructureMfe`.
- A stem never contains a multiloop (it would no longer be a single stem-loop).

**Not implemented:**

- Sequence-level pseudoknot prediction in this method; **users should rely on:** `RnaSecondaryStructure.DetectPseudoknots(...)` for post hoc detection or `RnaSecondaryStructure.PredictStructure(...)` for higher-level assembly.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or too-short sequence | No results | No valid loop-plus-stem span can be formed |
| `allowWobble = false` | Wobble pairs terminate stem extension | Explicit branch in the extension loop |
| Overlapping candidates | All qualifying candidates are yielded | The method performs no overlap filtering |

### 6.2 Limitations

This API is restricted to single hairpin motifs (no multiloops, no pseudoknots). Energies of isolated motifs exclude exterior dangles. Hairpin loops > 30 nt use the NNDB Jacobson–Stockmayer formula 6.4 + 1.75·RT·ln(n/9) unrounded; ViennaRNA truncates `lxc·ln(n/30)` to integer dcal/mol, so such loops can differ from `RNA.eval_structure` by ≤ 0.01 kcal/mol. The unbranched mode needs O(n²) memory.

## 7. Examples and Related Material

- [RNA-STEMLOOP-001](../../../tests/TestSpecs/RNA-STEMLOOP-001.md) documents the repository's stem-loop test specification.
- [RNA_Secondary_Structure.md](./RNA_Secondary_Structure.md) documents the higher-level structure-prediction API that consumes these candidates.

## 8. References

1. Wikipedia contributors. Stem-loop. Wikipedia. https://en.wikipedia.org/wiki/Stem-loop
2. Wikipedia contributors. Tetraloop. Wikipedia. https://en.wikipedia.org/wiki/Tetraloop
3. Wikipedia contributors. Pseudoknot. Wikipedia. https://en.wikipedia.org/wiki/Pseudoknot
4. Woese, C. R., S. Winker, and R. R. Gutell. 1990. Architecture of ribosomal RNA: constraints on the sequence of tetra-loops. Proceedings of the National Academy of Sciences 87(21):8467-8471.
5. Heus, H. A., and A. Pardi. 1991. Structural features that give rise to the unusual stability of RNA hairpins containing GNRA loops. Science 253(5016):191-194.
6. Svoboda, P., and A. Cara. 2006. Hairpin RNA: A secondary structure of primary importance. Cellular and Molecular Life Sciences 63(7):901-908.
7. Turner, D. H., and D. H. Mathews. NNDB: Turner 2004 hairpin loops (incl. worked examples §5.3) and special hairpin loops. https://rna.urmc.rochester.edu/NNDB/turner04/hairpin.html
8. Zuker, M., and P. Stiegler. 1981. Optimal computer folding of large RNA sequences using thermodynamics and auxiliary information. Nucleic Acids Research 9(1):133-148.
9. Lorenz, R., et al. 2011. ViennaRNA Package 2.0. Algorithms for Molecular Biology 6:26 (version 2.7.2 source: `eval/eval_hairpin.c`, `params/default.c`, `sequences/alphabet.c`).
