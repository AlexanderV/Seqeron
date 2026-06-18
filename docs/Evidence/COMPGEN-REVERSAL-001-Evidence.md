# Evidence Artifact: COMPGEN-REVERSAL-001

**Test Unit ID:** COMPGEN-REVERSAL-001
**Algorithm:** Reversal Distance (breakpoint-based lower bound)
**Date Collected:** 2026-06-14

---

## Online Sources

### Bafna V, Pevzner PA (1998). Sorting by Transpositions. SIAM J. Discrete Math. 11(2):224–240.

**URL:** https://www.ic.unicamp.br/~meidanis/courses/mo640/2008s2/textos/Bafna-Pevzner-1998.pdf
**Accessed:** 2026-06-14
**Retrieved by:** WebSearch "Bafna Pevzner sorting by reversals breakpoint distance lower bound" → fetched the PDF (binary), then extracted text with `pdftotext`.
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Breakpoint definition (verbatim, §2):** "For all 0 ≤ i ≤ n, the pair (π_i, π_{i+1}) is a breakpoint if π_{i+1} ≠ π_i + 1." The permutation is considered extended.
2. **Identity has zero breakpoints (verbatim):** "Observe that the identity permutation is the only permutation with 0 breakpoints, and therefore, sorting a permutation corresponds to decreasing the number of breakpoints."
3. **Operation/breakpoint lower bound principle (verbatim, for transpositions):** "a transposition can decrease the number of breakpoints by at most 3, implying a trivial lower bound of d(π) ≥ #breakpoints / 3." The same construction yields the reversal bound d(π) ≥ #breakpoints / 2 (a reversal cuts two adjacencies — see Hunter Lecture 16 below).
4. **b(π) notation:** the number of breakpoints of π; sorting reduces b(π) to 0.

---

### Hunter College CS, Computational Biology — Lecture 16: Genome rearrangements, sorting by reversals.

**URL:** https://www.cs.hunter.cuny.edu/~saad/courses/compbio/lectures/lecture16.pdf
**Accessed:** 2026-06-14
**Retrieved by:** WebSearch "sorting by reversals breakpoint at most 2 lower bound lecture notes" → fetched the PDF (binary), then extracted text with `pdftotext`.
**Authority rank:** 2 (course material from established curriculum, restating primary results of Bafna–Pevzner / Hannenhalli–Pevzner)

**Key Extracted Points:**

1. **Extended permutation (verbatim):** "We define the extended version of α as the following signed permutation: (α(0), α(1), …, α(n), α(n+1)) = (0, α(1), …, α(n), n+1)."
2. **Breakpoint definition (verbatim, signed):** "If (x, y) appear in (extended) α but neither (x, y) nor (−y, −x) appear in (extended) β, then (x, y) is a breakpoint of α with respect to β."
3. **Worked example (verbatim, exact numbers):** "α = (−2, −3, +1, +6, −5, −4) and β = (+1, +2, +3, +4, +5, +6). Then the extended versions … are (0, −2, −3, +1, +6, −5, −4, +7) … The breakpoints of α with respect to β are: (0, −2), (−2, −3), (−3, +1), (+1, +6), (+6, −5), (−4, +7). Note that (−5, −4) is not a breakpoint since (4, 5) appears in β." → **b(α) = 6** in the signed model.
4. **Reversal removes at most two breakpoints (verbatim):** "a reversal ρ = [i, j] can reduce the number of breakpoints by at most two since it separates between α(i−1) and α(i) and between α(j) and α(j+1) … b(α) − b(αρ) ≤ 2."
5. **Lower bound (verbatim derivation):** "b(α) − b(αρ_1 … ρ_t) ≤ 2t. But b(αρ_1 … ρ_t) = b(β) = 0. Therefore b(α) ≤ 2t. But d(α) ≥ t; therefore, d(α) ≥ b(α)/2. This lower bound is not very tight."
6. **Reversal distance is symmetric (verbatim):** "d_β(α) = d_α(β)"; β can be taken as the identity without loss of generality.

---

### Hübotter J (2020). On Sorting by Reversals (survey). + search-result confirmations.

**URL:** https://jonhue.github.io/min-sbr/paper.pdf
**Accessed:** 2026-06-14
**Retrieved by:** WebSearch "reversal distance breakpoint lower bound signed permutation worked example" and "unsigned permutation breakpoint definition reversal distance lower bound" → URL surfaced; the PDF body was compressed and not text-extractable, but the search-engine summaries (corroborating multiple independent results) confirmed the unsigned statement below.
**Authority rank:** 4 (survey citing primaries; used only to corroborate the unsigned variant, which is also derivable from Bafna–Pevzner §2)

**Key Extracted Points:**

1. **Unsigned breakpoint (verbatim from corroborating search results):** "An unsigned reversal breakpoint exists between a pair of consecutive elements (π_i, π_{i+1}) if |π_{i+1} − π_i| ≠ 1, for 0 ≤ i ≤ n." This is the unsigned specialization of Bafna–Pevzner §2 (drop the +1 direction, use |Δ|).
2. **Unsigned lower bound (verbatim from corroborating search results):** "A reversal can eliminate at most two breakpoints, which implies b(π)/2 ≤ d_r(π)."

---

### Bergeron A, Mixtacki J, Stoye J (2009). The Inversion Distance Problem. In *Mathematics of Evolution and Phylogeny* (course chapter).

**URL:** https://gi.cebitec.uni-bielefeld.de/_media/teaching/2018winter/cg/inversionbergeron.pdf
**Accessed:** 2026-06-14
**Retrieved by:** WebSearch "unsigned permutation breakpoint definition reversal distance lower bound Pevzner" → fetched PDF (binary), extracted text with `pdftotext`.
**Authority rank:** 1 (textbook chapter)

**Key Extracted Points:**

1. **Adjacency vs breakpoint (verbatim):** "P1 has two adjacencies, −2 · −1 and 6 · 7. All other points of P1 are breakpoints." (Confirms: consecutive elements that are NOT a consecutive pair form a breakpoint; the rest are adjacencies.)
2. **Breakpoint graph attribution (verbatim):** "Bafna and Pevzner [2] extended the breakpoint graph to signed [permutations]. The most common version of the breakpoint graph is based on an unsigned [permutation]."

---

## Documented Corner Cases and Failure Modes

### From Bafna & Pevzner (1998)

1. **Identity permutation:** the only permutation with 0 breakpoints ⇒ distance 0.
2. **Lower bound is not exact:** "a permutation with few breakpoints may be more distant from the identity permutation than one with more breakpoints" — the breakpoint bound is a *lower bound*, not the true reversal distance.

### From Hunter Lecture 16

1. **A reversal that removes two breakpoints need not make progress:** "a reversal ρ that decreases the number of breakpoints by two is not necessarily one that makes progress" — confirms the bound is loose; exact distance needs the cycle/hurdle refinement (Hannenhalli–Pevzner).

---

## Test Datasets

### Dataset: Hunter Lecture 16 worked permutation (unsigned specialization)

**Source:** Hunter College CompBio Lecture 16 (worked example), applied to the unsigned breakpoint model implemented here.

The lecture's signed example uses α = (−2, −3, +1, +6, −5, −4) with target identity. This unit's `CalculateReversalDistance` operates on **unsigned** gene-order indices, so the unsigned specialization uses the magnitudes [2,3,1,6,5,4] vs identity [1,2,3,4,5,6].

| Parameter | Value |
|-----------|-------|
| Input perm1 | [2, 3, 1, 6, 5, 4] |
| Input perm2 (target) | [1, 2, 3, 4, 5, 6] |
| Extended unsigned | [0, 2, 3, 1, 6, 5, 4, 7] |
| Unsigned breakpoints | 0→2 (\|Δ\|=2 BP), 2→3 (adj), 3→1 (BP), 1→6 (BP), 6→5 (adj), 5→4 (adj), 4→7 (BP) ⇒ **b = 4** |
| Lower bound d ≥ b/2 | 2 |
| Returned value ⌈b/2⌉ | **2** |

### Dataset: Fully reversed sequence

**Source:** Direct application of breakpoint definition (Bafna–Pevzner §2).

| Parameter | Value |
|-----------|-------|
| Input perm1 | [4, 3, 2, 1] |
| Input perm2 (target) | [1, 2, 3, 4] |
| Extended unsigned | [0, 4, 3, 2, 1, 5] |
| Unsigned breakpoints | 0→4 (BP), 4→3 adj, 3→2 adj, 2→1 adj, 1→5 (BP) ⇒ **b = 2** |
| Lower bound ⌈b/2⌉ | **1** (a single full reversal sorts it) |

### Dataset: Identity

| Parameter | Value |
|-----------|-------|
| perm1 = perm2 | [1, 2, 3, 4, 5] |
| Breakpoints | 0 |
| Distance | 0 |

---

## Assumptions

1. **ASSUMPTION: Unsigned model boundary rounding** — The lower bound theorem states d ≥ b/2 (real-valued). For an integer return, the smallest integer satisfying the bound is ⌈b/2⌉. The implementation returns ⌈b/2⌉ = `(b + 1) / 2` (integer arithmetic). This is the canonical integer breakpoint lower bound; it is the tightest integer guaranteed by the theorem, not an invented value.
2. **ASSUMPTION: Unequal-length inputs throw** — The sources define reversal distance only between two permutations of the same marker set (same n). The implementation throws `ArgumentException` when lengths differ. Not separately specified by the sources; it is the only well-defined behavior.

---

## Recommendations for Test Coverage

1. **MUST Test:** Identity permutations → distance 0 — Evidence: Bafna–Pevzner §2 ("identity is the only permutation with 0 breakpoints").
2. **MUST Test:** Hunter worked example (unsigned specialization) [2,3,1,6,5,4] → 2 — Evidence: Hunter Lecture 16 worked example + unsigned breakpoint definition.
3. **MUST Test:** Fully reversed [4,3,2,1] → 1 — Evidence: breakpoint definition (b=2, ⌈b/2⌉=1).
4. **MUST Test:** Lower-bound property d ≤ true distance and d ≥ 0 — Evidence: Hunter Lecture 16 ("d(α) ≥ b(α)/2").
5. **SHOULD Test:** Single-element and empty inputs → 0 — Rationale: n ≤ 1 has no internal adjacency / no breakpoint.
6. **SHOULD Test:** Unequal lengths throw `ArgumentException` — Rationale: distance undefined across different marker sets.
7. **COULD Test:** Symmetry d(α,β) = d(β,α) — Rationale: Hunter Lecture 16 ("d_β(α) = d_α(β)").

---

## References

1. Bafna V, Pevzner PA (1998). Sorting by Transpositions. *SIAM Journal on Discrete Mathematics* 11(2):224–240. https://www.ic.unicamp.br/~meidanis/courses/mo640/2008s2/textos/Bafna-Pevzner-1998.pdf
2. Hunter College Computational Biology — Lecture 16: Genome rearrangements, sorting by reversals. https://www.cs.hunter.cuny.edu/~saad/courses/compbio/lectures/lecture16.pdf
3. Hübotter J (2020). On Sorting by Reversals. https://jonhue.github.io/min-sbr/paper.pdf
4. Bergeron A, Mixtacki J, Stoye J (2009). The Inversion Distance Problem. https://gi.cebitec.uni-bielefeld.de/_media/teaching/2018winter/cg/inversionbergeron.pdf

---

## Change History

- **2026-06-14**: Initial documentation.
