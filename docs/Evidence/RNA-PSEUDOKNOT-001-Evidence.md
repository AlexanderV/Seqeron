# Evidence Artifact: RNA-PSEUDOKNOT-001

**Test Unit ID:** RNA-PSEUDOKNOT-001
**Algorithm:** Pseudoknot Detection (crossing base pairs)
**Date Collected:** 2026-06-14

---

## Online Sources

### Antczak et al. (2018) — "New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation"

**URL:** https://academic.oup.com/bioinformatics/article/34/8/1304/4721780
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *Bioinformatics*)
**Retrieved how:** WebSearch "Antczak 2018 RNA pseudoknot order dot-bracket crossing base pairs definition" → WebFetch of the Oxford Academic article URL above.

**Key Extracted Points:**

1. **Crossing/conflict definition (verbatim from retrieved text):** Two base pairs are crossed/conflicted when "for any pair (*i, i'*) there exists another one, (*j, j'*), such that *i < j < i' < j'*." (Renaming to this unit's symbols: pairs (i,j) and (k,l) cross when **i < k < j < l**.)
2. **Pseudoknot order:** defined as "a minimum number of base pair set decompositions resulting in a nested structure." If removing one double-stranded region leaves a structure without conflicts, that pseudoknot has order 1.
3. **Dot-bracket-letter (DBL) notation by region order:** order 0 → `( )`, order 1 → `[ ]`, order 2 → `{ }`, order 3 → `< >`, orders 4–8 → letters `A/a`,`B/b`,`C/c`,`D/d`,`E/e`. H-type pseudoknot encoded `( [ ) ]`.

### Smit & Knight et al. (2008) — "From knotted to nested RNA structures: a variety of computational methods for pseudoknot removal"

**URL:** https://rnajournal.cshlp.org/content/14/3/410
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *RNA*)
**Retrieved how:** WebSearch "Smit 2008 RNA pseudoknot detection dynamic programming from knotted to nested" → result block from rnajournal.cshlp.org.

**Key Extracted Points:**

1. **Pseudoknot presence requires crossing pairs:** pseudoknotted structures contain base pairs that cannot be drawn without crossing; obtaining a nested (pseudoknot-free) structure requires disregarding some base pairs of the knotted structure.
2. This algorithm family (pseudoknot removal / order assignment) is the origin of the reference-implementation approach used by biotite (see below).

### biotite.structure.pseudoknots (reference implementation)

**URL:** https://www.biotite-python.org/latest/apidoc/biotite.structure.pseudoknots.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation in an established library)
**Retrieved how:** WebSearch result → WebFetch of the biotite API page above.

**Key Extracted Points:**

1. **Crossing condition:** base pairs are "crossing/knotted when they cannot be arranged in a purely nested (non-crossing) configuration"; pairs conflict when their positions interleave (one pair's open lies inside the other pair while its close lies outside).
2. **Order assignment / layering:** pseudoknot order = "the minimum number of base pair set decompositions resulting in a nested structure"; nested pairs receive order 0, knotted pairs receive order 1+ (algorithm by Smit et al.).

### Wikipedia — "Pseudoknot"

**URL:** https://en.wikipedia.org/wiki/Pseudoknot
**Accessed:** 2026-06-14
**Authority rank:** 4 (cites primaries Rivas & Eddy 1999; Antczak et al. 2018)
**Retrieved how:** WebFetch of the Wikipedia article URL above.

**Key Extracted Points:**

1. **Qualitative definition (verbatim):** "A pseudoknot is a nucleic acid secondary structure containing at least two stem-loop structures in which half of one stem is intercalated between the two halves of another stem."
2. **Non-nested wording (verbatim):** "The base pairing in pseudoknots is not well nested; that is, base pairs occur that 'overlap' one another in sequence position."
3. **Cited primaries:** Rivas E, Eddy SR (1999) dynamic-programming algorithm for pseudoknots; Antczak et al. (2018) dot-bracket representation.

---

## Documented Corner Cases and Failure Modes

### From Antczak et al. (2018) / biotite

1. **Nested pairs are NOT pseudoknots:** a pair fully contained in another (i < k < l < j) is nested → no pseudoknot.
2. **Disjoint pairs are NOT pseudoknots:** non-overlapping ranges (j < k) are side-by-side → no pseudoknot.
3. **A pseudoknot needs ≥ 2 base pairs:** an empty pair set or a single pair can never cross → no pseudoknot.

### From Smit & Knight (2008)

1. **Order direction is symmetric:** crossing is a relation between two pairs; which pair is "the pseudoknot" is a labeling choice for removal/order, but the *presence* of a crossing is symmetric.

---

## Test Datasets

### Dataset: H-type pseudoknot (Antczak 2018 §DBL example "( [ ) ]")

**Source:** Antczak et al. (2018), DBL example — the minimal H-type pseudoknot `([)]`.

Sequence positions (0-based) and the two base pairs of `([)]`:

| Symbol index | 0 | 1 | 2 | 3 |
|--------------|---|---|---|---|
| DBL          | ( | [ | ) | ] |

| Pair | open | close |
|------|------|-------|
| P1 `(...)` | 0 | 2 |
| P2 `[...]` | 1 | 3 |

Check i<k<j<l with (i,j)=(0,2), (k,l)=(1,3): 0 < 1 < 2 < 3 → **crossing → exactly one pseudoknot detected.**

### Dataset: Nested pairs (negative control)

**Source:** Wikipedia "not well nested" wording; Antczak nested condition i<k<l<j.

Pairs (0,5) and (1,4): 0 < 1 < 4 < 5 → fully nested → **no pseudoknot.**

### Dataset: Disjoint pairs (negative control)

Pairs (0,2) and (3,5): ranges [0,2] and [3,5] do not overlap → **no pseudoknot.**

---

## Assumptions

1. **ASSUMPTION: Empty / single-pair input returns no pseudoknots.** Not stated as an explicit corner case in the sources, but follows directly from the crossing definition: a pseudoknot requires two pairs that cross, so fewer than two pairs cannot produce one. Treated as a derived (not invented) consequence, not a free parameter.
2. **ASSUMPTION: Each crossing pair-of-pairs is reported as one Pseudoknot result.** The sources define *crossing as a binary relation between two base pairs*; the method reports each crossing relation. Grouping multiple mutually-crossing pairs into a single higher-order pseudoknot (DBL order assignment, Antczak 2018) is a richer feature that is Not Implemented here and documented as such.

---

## Recommendations for Test Coverage

1. **MUST Test:** H-type `([)]` pairs (0,2)+(1,3) → exactly one pseudoknot whose crossing pairs are those two. — Evidence: Antczak (2018) crossing condition i<k<j<l and `([)]` example.
2. **MUST Test:** Nested pairs (0,5)+(1,4) → zero pseudoknots. — Evidence: Antczak nested condition i<k<l<j; Wikipedia "not well nested".
3. **MUST Test:** Disjoint pairs (0,2)+(3,5) → zero pseudoknots. — Evidence: crossing requires range overlap.
4. **MUST Test:** Pair endpoints stored in either order (Position1>Position2) are normalized via min/max before the crossing test. — Evidence: crossing is defined on open<close positions.
5. **SHOULD Test:** Single pair and empty input → zero pseudoknots. — Rationale: derived corner case (≥2 pairs needed).
6. **SHOULD Test:** Determinism / order independence — same input set always yields the same pseudoknots. — Rationale: invariant for a pure combinatorial scan.
7. **COULD Test:** Property test — for random pair sets, every reported pseudoknot's two pairs satisfy i<k<j<l (O(n²) invariant). — Rationale: O(n²) unit requires a property-based invariant test.

---

## References

1. Antczak M, Popenda M, Zok T, Zurkowski M, Adamiak RW, Szachniuk M (2018). New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation. *Bioinformatics* 34(8):1304–1312. https://academic.oup.com/bioinformatics/article/34/8/1304/4721780
2. Smit S, Rother K, Heringa J, Knight R (2008). From knotted to nested RNA structures: a variety of computational methods for pseudoknot removal. *RNA* 14(3):410–416. https://rnajournal.cshlp.org/content/14/3/410
3. biotite.structure.pseudoknots — Biotite documentation. https://www.biotite-python.org/latest/apidoc/biotite.structure.pseudoknots.html (accessed 2026-06-14)
4. Pseudoknot — Wikipedia (citing Rivas & Eddy 1999; Antczak et al. 2018). https://en.wikipedia.org/wiki/Pseudoknot (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
