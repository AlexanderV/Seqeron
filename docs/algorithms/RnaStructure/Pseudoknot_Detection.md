# Pseudoknot Detection

| Field | Value |
|-------|-------|
| Algorithm Group | RnaStructure |
| Test Unit ID | RNA-PSEUDOKNOT-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-11-25 |

## 1. Overview

Pseudoknot detection identifies pairs of RNA base pairs that *cross* — a non-nested arrangement that classical (nested) secondary-structure models cannot represent. Given a set of base pairs, the algorithm reports every crossing pair-of-pairs. It is an exact, deterministic combinatorial test (no thermodynamics, no scoring): the output depends only on nucleotide positions. Use it to flag pseudoknotted regions in a predicted or experimentally-derived structure. The companion `AssignPseudoknotOrders` decomposes the pairs into nested layers (pseudoknot orders: order 0 = a maximum nested subset, i.e. the pairs kept by optimal pseudoknot removal) and `ToDotBracketByPseudoknotOrder` / `ToDotBracket` render the pseudoknotted structure in dot-bracket-letter notation.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

RNA secondary structure is the set of base pairs of a folded RNA. Most structures are *nested*: pairs are either disjoint or fully contained, so the structure can be drawn as non-crossing arcs (or as balanced parentheses). A **pseudoknot** is "a nucleic acid secondary structure containing at least two stem-loop structures in which half of one stem is intercalated between the two halves of another stem" [4]; equivalently, "base pairs occur that 'overlap' one another in sequence position" [4]. Such structures are *not well nested* and cannot be expressed with a single bracket pair [1][4].

### 2.2 Core Model

Let two base pairs be written with open < close positions as (i, j) and (k, l). The pairs **cross** (form a pseudoknot) iff

> i < k < j < l

[1][3]. Equivalent verbatim form from Antczak et al. (2018): for a pair (i, i′) there exists another (j, j′) with i < j < i′ < j′ [1]. The complementary arrangements are not pseudoknots:

- **Nested:** i < k < l < j (one pair fully inside the other) [1][4].
- **Disjoint:** j < k (non-overlapping ranges).

### 2.2a Pseudoknot order and bracket families

*Pseudoknot order* is "a minimum number of base pair set decompositions resulting in a nested structure"; order 0 → `()`, 1 → `[]`, 2 → `{}`, 3 → `<>`, 4.. → `Aa`, `Bb`, … [1]. `AssignPseudoknotOrders` follows the reference implementation biotite `structure.pseudoknots` [3] (source code read: "removes base pairs until the remaining base pairs are completely nested … the pseudoknot order of the removed base pairs is incremented and the procedure is repeated … Base pairs are removed in a way that maximizes the number of remaining base pairs", the Smit et al. dynamic programming [2]): each layer is an exact maximum-cardinality nested subset of the remaining pairs. biotite returns *all* optimal solutions; this method returns one, keeping the 5'-most opening pair at the first difference (lexicographically greatest subset by opening position).

The iterated maximum is an upper bound on the minimum of [1]: the minimum number of nested layers is the chromatic number of the crossing (circle) graph, NP-complete to compute in general [5], and maximising layer 0 can cost a layer — pairs (0,4),(1,12),(2,6),(3,10),(8,13),(9,14) need 4 layers here although 3 colours suffice (exhaustive search; 139/3,000 random sets of 2–9 pairs). RNApdbee [1] likewise offers several heuristics (EG, FCFS, DP, hybrid) because the two criteria conflict.

`ToDotBracket` assigns families exactly as ViennaRNA 2.7.2 `vrna_db_from_ptable` [6]: repeated 5'→3' passes, a pair takes the current family when it nests inside the pairs opened earlier in the pass, otherwise it is deferred; families `([{<` then `A..Z`. This first-fit rendering can differ from the order rendering (a larger 3' helix crossing a smaller 5' one: order string `[.(((.].)))`, ViennaRNA `(.[[[.).]]]`). Both parse back to the same pairs with `ParseDotBracket`.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Each base pair has two distinct positions interpreted as (open, close) after min/max normalization. | A degenerate pair (i = j) cannot cross anything; it is silently never part of a pseudoknot. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported pseudoknot's two pairs satisfy i < k < j < l. | Direct guard in the crossing test [1][3]. |
| INV-02 | Nested pairs (i < k < l < j) produce no pseudoknot. | They fail j < l. [1][4] |
| INV-03 | Disjoint pairs (j < k) produce no pseudoknot. | They fail k < j. |
| INV-04 | Fewer than two base pairs → empty result. | A pseudoknot needs two crossing pairs [2]. |
| INV-05 | Detection is deterministic and order-independent over the same pair set. | Pure double loop with symmetric (open-first) ordering of each pair-of-pairs. |

### 2.5 Comparison with Related Methods

| Aspect | Pseudoknot Detection (this) | Nested folding (e.g. Nussinov / MFE) |
|--------|------------------------------|--------------------------------------|
| Input | a fixed base-pair set | a sequence |
| Output | crossing pair-of-pairs | a nested base-pair set |
| Pseudoknots | detects them | excludes them by construction |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| basePairs | `IReadOnlyList<BasePair>` | required | Base pairs to scan; positions are 0-based. | `null` allowed (→ empty result); positions need not be pre-ordered. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (each) | `Pseudoknot` | `Start1<End1` and `Start2<End2` with Start1<Start2; `Start1<Start2<End1<End2` (crossing). `CrossingPairs` holds the two original `BasePair`s. |
| `AssignPseudoknotOrders` | `IReadOnlyList<int>` | order of each input pair (index-aligned); pairs of equal order never cross; order 0 is a maximum nested subset. |
| `ToDotBracket(length, pairs)` | `string` | ViennaRNA `vrna_db_from_ptable` rendering; `ParseDotBracket` returns the same pairs. |
| `ToDotBracketByPseudoknotOrder(length, pairs)` | `string` | order *k* written with family *k* of `() [] {} <> Aa …` [1]. |

### 3.3 Preconditions and Validation

`null` or a set of fewer than two pairs returns an empty sequence (no exception). Positions are 0-based; each pair is normalized to (open < close) before comparison, so a pair stored as (close, open) is treated identically. No alphabet/case handling is involved (the test is purely positional). The method is lazy (`yield`) and never throws; pairs that share an endpoint (not a secondary structure) or have i = j never satisfy the strict inequalities.

`AssignPseudoknotOrders`, `ToDotBracket` and `ToDotBracketByPseudoknotOrder` require a secondary structure and throw `ArgumentException` for a self pair, a negative position (or ≥ `length`), or a position used by two pairs (duplicates included) — ViennaRNA's pair table cannot hold such input and `vrna_db_from_ptable` returns NULL for an inconsistent table. A structure needing more than 30 families throws (ViennaRNA logs a warning and silently drops the remaining pairs).

## 4. Algorithm

### 4.1 High-Level Steps

1. If the set is `null` or has fewer than two pairs, return empty.
2. For each unordered pair-of-pairs (a, b):
   1. Normalize each base pair to (open < close).
   2. Order the two pairs so the one opening first is (i, j).
   3. If i < k < j < l, emit a `Pseudoknot` carrying both original pairs.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Single decision rule, the crossing condition i < k < j < l [1][3]. No scoring tables or thermodynamic parameters are involved.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| DetectPseudoknots | O(n²) | O(1) extra (plus O(p) output) | n = number of base pairs; p = number of crossing pairs reported (worst case Θ(n²), so the scan is output-optimal in the worst case). |
| AssignPseudoknotOrders | O(K·m²) | O(m²) | m = 2 × pairs, K = number of orders; interval DP per layer. |
| ToDotBracket | O(F·L) | O(L) | L = length, F ≤ 30 families (passes). |

No substring/pattern search is performed (the input is a set of integer-position pairs, not a sequence), so the repository suffix tree is not applicable.

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.DetectPseudoknots(IReadOnlyList<BasePair>)`: returns one `Pseudoknot` per crossing pair-of-pairs (MCP `detect_pseudoknots`).
- `RnaSecondaryStructure.AssignPseudoknotOrders(IReadOnlyList<BasePair>)`, `ToDotBracketByPseudoknotOrder(int, IReadOnlyList<BasePair>)`, `ToDotBracket(int, IEnumerable<(int,int)>)` (MCP `pseudoknot_orders`). `ToDotBracket` is also the single renderer used by the MFE traceback, `PredictStructure`, `FindStemLoops` and both pseudoknot predictors (their former explicit "stem 2 → `[]`" path removed; outputs unchanged).

### 5.2 Current Behavior

The method is a lazy `O(n²)` all-pairs scan. Each pair is normalized to (open < close) and the two pairs are reordered by opening position, making the result independent of how endpoints were stored and of input ordering. Reused infrastructure: none — the suffix tree was evaluated and rejected because the input is a positional base-pair set, not a sequence to search.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Crossing condition i < k < j < l identifying a pseudoknot between two base pairs [1][3].
- Nested (i<k<l<j) and disjoint (j<k) arrangements correctly excluded [1][4].

**Intentionally simplified:**

- (none)

- Pseudoknot order assignment per biotite / Smit et al. [2][3] and DBL rendering [1]; ViennaRNA `vrna_db_from_ptable` family assignment [6].

**Not implemented:**

- Exact minimum layer count (chromatic number of the circle graph, NP-complete [5]); the sourced iterated-maximum procedure is used, see §2.2a.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Each crossing pair-of-pairs reported separately | Contract | A region with many mutually-crossing pairs yields several `Pseudoknot` results; grouping into layers is `AssignPseudoknotOrders` | accepted | Conflict relation of [1][2]. |
| 2 | > 30 families rejected | Deviation | ViennaRNA drops the excess pairs with a warning; we throw so the rendering is never lossy | accepted | `ParseDotBracket` reads exactly these 30 families. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` input | empty result | INV-04 / contract |
| empty or single pair | empty result | < 2 pairs cannot cross [2] |
| nested pairs (0,5)+(1,4) | empty result | INV-02 [1][4] |
| disjoint pairs (0,2)+(3,5) | empty result | INV-03 |
| pair stored as (close,open) | normalized; same result as (open,close) | open<close normalization |

### 6.2 Limitations

`DetectPseudoknots` reports the pairwise crossing relation; layers/orders come from `AssignPseudoknotOrders` (order-0 layer = a maximum nested subset, i.e. optimal pseudoknot removal). The number of orders may exceed the theoretical minimum (§2.2a). No H-type / kissing-hairpin classification is made. The input base pairs must already be determined (no folding).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

// H-type pseudoknot "([)]": pairs (0,2) and (1,3)
var pairs = new List<BasePair>
{
    new(0, 2, 'A', 'U', BasePairType.WatsonCrick),
    new(1, 3, 'G', 'C', BasePairType.WatsonCrick),
};
var knots = DetectPseudoknots(pairs).ToList();
// knots.Count == 1; knots[0].Start1==0, End1==2, Start2==1, End2==3
```

**Numerical / biological walk-through:**

For `([)]` the two pairs are (i,j)=(0,2) and (k,l)=(1,3). Check 0 < 1 < 2 < 3 → all four inequalities hold → crossing → one pseudoknot [1]. For nested (0,5)+(1,4): 0<1 but 5<4 is false → no pseudoknot.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_DetectPseudoknots_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_DetectPseudoknots_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [RNA-PSEUDOKNOT-001-Evidence.md](../../../docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md)
- Related algorithms: [RNA Base Pairing](../RnaStructure/RNA_Base_Pairing.md)

## 8. References

1. Antczak M, Popenda M, Zok T, Zurkowski M, Adamiak RW, Szachniuk M. 2018. New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation. *Bioinformatics* 34(8):1304–1312. https://academic.oup.com/bioinformatics/article/34/8/1304/4721780
2. Smit S, Rother K, Heringa J, Knight R. 2008. From knotted to nested RNA structures: a variety of computational methods for pseudoknot removal. *RNA* 14(3):410–416. https://rnajournal.cshlp.org/content/14/3/410
3. biotite.structure.pseudoknots. Biotite documentation. https://www.biotite-python.org/latest/apidoc/biotite.structure.pseudoknots.html
4. Pseudoknot. Wikipedia (citing Rivas E, Eddy SR. 1999). https://en.wikipedia.org/wiki/Pseudoknot
5. Garey MR, Johnson DS, Miller GL, Papadimitriou CH. 1980. The complexity of coloring circular arcs and chords. *SIAM J. Algebraic Discrete Methods* 1(2):216–227.
6. ViennaRNA 2.7.2 source, `src/ViennaRNA/structures/structure_dotbracket.c` (`vrna_db_from_ptable`) and `structure_pairtable.c` (`vrna_ptable_from_string`, `vrna_pt_pk_remove`), PyPI sdist.
