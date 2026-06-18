# Evidence Artifact: COMPGEN-REARR-001

**Test Unit ID:** COMPGEN-REARR-001
**Algorithm:** Genome Rearrangement Detection by Breakpoints (signed gene-order comparison)
**Date Collected:** 2026-06-13

---

## Online Sources

### Hunter College Computational Biology — Lecture 16 "Genome rearrangements, sorting by reversals"

**URL:** https://www.cs.hunter.cuny.edu/~saad/courses/compbio/lectures/lecture16.pdf
**Accessed:** 2026-06-13
**Authority rank:** 1 (course exposition of Hannenhalli–Pevzner / Bafna–Pevzner theory; primary-derived)
**Retrieved how:** WebSearch query `breakpoint distance genome rearrangement definition adjacency signed permutation formula` → fetched the PDF; WebFetch could not decode the compressed PDF, so the binary was extracted to text with `pdftotext` (poppler) and read with the Read/Bash tools. Facts below are copied from the extracted text (lines 254–312).

**Key Extracted Points:**

1. **Signed permutation (verbatim):** "A signed permutation α over the set L = 1,2,...,n is a permutation of L such that α(i) = +a or −a, where a ∈ L." Example given: `α = (+3, −2, −1, +4, −5)`.
2. **Reversal definition (verbatim):** `α' = α[i,j] = (α(1),...,α(i−1), −α(j),...,−α(i), α(j+1),...,α(n))` — a reversal of a contiguous block reverses the segment **and negates the sign** of each element in it.
3. **Extended permutation (verbatim):** "We define the extended version of α as ... `(α(0), α(1),...,α(n), α(n+1)) = (0, α(1),...,α(n), n+1)`" — i.e. prepend `0` and append `n+1`.
4. **Breakpoint definition (verbatim):** "If (x, y) appear in (extended) α but neither (x, y) nor (−y, −x) appear in (extended) β, then (x, y) is a breakpoint of α with respect to β."
5. **Worked example (verbatim numbers):** `α = (−2, −3, +1, +6, −5, −4)`, `β = (+1, +2, +3, +4, +5, +6)`. Extended: `α = (0, −2, −3, +1, +6, −5, −4, +7)`. "The breakpoints of α with respect to β are: (0,−2), (−2,−3), (−3,+1), (+1,+6), (+6,−5), (−4,+7)." **b(α) = 6.** "Note that (−5,−4) is not a breakpoint since (4,5) appears in β." (Here `(−5,−4)` as a consecutive pair in α corresponds to the descending adjacency that maps to identity adjacency `(4,5)` after accounting for sign, hence not a breakpoint.)
6. **b(β) = 0 (verbatim):** "Note that b_β(β) = 0." The identity has no breakpoints.
7. **Breakpoint lower bound (verbatim derivation):** a reversal changes the breakpoint count by at most 2: "b(α) − b(αρ) ≤ 2"; therefore over t reversals `b(α) ≤ 2t` and since `d(α) ≥ t`, **`d(α) ≥ b(α)/2`** — the breakpoint distance is a lower bound for reversal distance.
8. **Why breakpoints matter (verbatim):** "if (x, y) is a breakpoint, then in order to transform α into β, some reversal must separate between x and y."

### Tannier, Zheng & Sankoff — "On the Complexity of Rearrangement Problems under the Breakpoint Distance" (via PMC full text)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3887456/
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper; cites Tannier et al. 2009 for the breakpoint distance)
**Retrieved how:** WebSearch query `"breakpoint distance" genome rearrangement number of conserved adjacencies n - 1 - adjacencies formula telomere` → opened the PMC HTML page with WebFetch and extracted the definitions.

**Key Extracted Points:**

1. **Breakpoint (verbatim):** "When two genes (or conserved segments or markers) are adjacent in one genome, but not in the other, we call this position a breakpoint."
2. **Adjacency (verbatim):** "An edge between extremities x and y, called adjacency, indicates that x and y are adjacent in the genome." Each gene `g` has two extremities, a head `g^h` and a tail `g^t`.
3. **Breakpoint distance formula (verbatim):** `d(π₁, π₂) = n − sim(π₁, π₂)`, where `n` = number of genes and `sim(π₁, π₂)` = "the number of common adjacencies" (plus half the common telomeric adjacencies in models with telomeres).
4. **Telomeric adjacency (verbatim):** "A telomeric adjacency xT_x indicates that x is an end of a linear chromosome." Telomere extensions correspond to the extended-permutation endpoints `0` and `n+1`.

### Bafna & Pevzner — "Sorting by Transpositions" (SIAM J. Discrete Math. 11(2):224–240, 1998; via UNICAMP course copy)

**URL:** https://www.ic.unicamp.br/~meidanis/courses/mo640/2008s2/textos/Bafna-Pevzner-1998.pdf
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed paper, the canonical transposition reference)
**Retrieved how:** WebSearch query `Bafna Pevzner genome rearrangements sorting reversals transpositions breakpoint definition adjacency 1995` → fetched PDF; decoded with `pdftotext` (poppler) and read pages 1–4. Facts copied from the extracted text.

**Key Extracted Points:**

1. **Rearrangement operations (verbatim):** "Genomes evolve by inversions and transpositions as well as by more simple operations of deletion, insertion, and duplication of fragments." This enumerates the operation classes modelled by `RearrangementType`.
2. **Permutation model & extension (verbatim):** "We assume that the order of genes in a genome is represented by a permutation π = π₁π₂,...,πₙ. Extend the permutation to include π₀ = 0 and π_{n+1} = n + 1."
3. **Transposition (verbatim):** "ρ(i, j, k) ... 'inserts' an interval [i, j−1] of π between π_{k−1} and π_k ... has the effect of moving genes π_i, π_{i+1},...,π_{j−1} to a new location in a genome." A transposition **moves a block, preserving its internal orientation** (no sign change), distinguishing it from a reversal.
4. **Inversion is a reversal (verbatim):** the paper treats "inversions and transpositions of long fragments"; inversions correspond to the reversal operation (sign-flipping) of the sorting-by-reversals problem.

---

## Documented Corner Cases and Failure Modes

### From Hunter lecture / Bafna–Pevzner

1. **Identity / collinear order:** `b(β) = 0` — two genomes with identical signed gene order have zero breakpoints, hence no rearrangements detected.
2. **Sign-only adjacency:** a descending consecutive pair whose unsigned values are consecutive in the target (e.g. `(−5, −4)` ↔ identity `(4,5)`) is **not** a breakpoint; the criterion must test both `(x,y)` and `(−y,−x)`.
3. **Boundary breakpoints:** the extended endpoints `(0, π₁)` and `(πₙ, n+1)` can themselves be breakpoints (see worked example: `(0,−2)` and `(−4,+7)` are counted).

### From Tannier et al. (PMC3887456)

1. **Telomeres:** ends of linear chromosomes are telomeric adjacencies; in the unichromosomal permutation model they are the extended endpoints `0` and `n+1`.
2. **No common adjacencies:** when `sim = 0`, `d_BP = n` (every adjacency is a breakpoint).

---

## Test Datasets

### Dataset: Hunter lecture worked example (b(α) = 6)

**Source:** Hunter Lecture 16, lines 287–289 (extracted text).

| Parameter | Value |
|-----------|-------|
| α (signed) | (−2, −3, +1, +6, −5, −4) |
| β (target, identity) | (+1, +2, +3, +4, +5, +6) |
| Extended α | (0, −2, −3, +1, +6, −5, −4, +7) |
| Breakpoints | (0,−2), (−2,−3), (−3,+1), (+1,+6), (+6,−5), (−4,+7) |
| **b(α)** | **6** |
| Non-breakpoint | (−5,−4) — because (4,5) ∈ β |
| Reversal-distance lower bound d ≥ b/2 | d ≥ 3 |

### Dataset: Identity permutation (b = 0)

**Source:** Hunter Lecture 16, "b_β(β) = 0".

| Parameter | Value |
|-----------|-------|
| α | (+1, +2, +3, +4, +5) |
| β | (+1, +2, +3, +4, +5) |
| **b(α)** | **0** |

### Dataset: Single inversion (one reversed block)

**Source:** reversal definition (Hunter, line 256) — `α[i,j]` reverses a segment and negates signs.

Target `β = (+1,+2,+3,+4,+5)`. Reverse the block [2,4]: `α = (+1, −4, −3, −2, +5)`.
Extended α = `(0, +1, −4, −3, −2, +5, +6)`. Adjacencies: `(0,+1)` ok (→1 adjacency 0,1); `(+1,−4)` breakpoint; `(−4,−3)` → maps to identity adjacency (3,4) via `(−y,−x)=(+3,+4)`? pair is `(−4,−3)`; `(−y,−x) = (+3,+4)` which **is** in β → not a breakpoint; `(−3,−2)` → `(+2,+3)` ∈ β → not a breakpoint; `(−2,+5)` breakpoint; `(+5,+6)` ok. **b(α) = 2.**

| Parameter | Value |
|-----------|-------|
| α | (+1, −4, −3, −2, +5) |
| **b(α)** | **2** |
| d ≥ b/2 | d ≥ 1 (one reversal suffices) |

---

## Assumptions

1. **ASSUMPTION: Anchors supplied as an ordered ortholog mapping.** The breakpoint definition is over a permutation of common markers. This unit derives the permutation from two gene lists plus an `orthologMap` (genome-1 gene id → genome-2 gene id), reusing the same input shape as the sibling synteny/ortholog units (COMPGEN-SYNTENY-001, COMPGEN-ORTHO-001). Anchor *generation* is out of scope and delegated to those units; the breakpoint computation itself is unchanged by this separation, so this is not an open correctness gap.
2. **ASSUMPTION: Strand char `'+' / '-'` encodes the sign.** Each `Gene.Strand` is mapped to the permutation sign (`+` → positive, `-` → negative). The reversed-target gene also carries the opposite strand. This is the standard signed-permutation encoding (Hunter, signed permutation definition) and is not an open gap.
3. **ASSUMPTION: Translocation / Deletion / Insertion / Duplication are NOT classified by this method.** These require chromosome identifiers (translocation) or gene-set differences (indel/duplication) that a single in-order signed permutation does not express, and no authoritative single-permutation rule assigns them. `DetectRearrangements`/`ClassifyRearrangement` therefore classify only the two operations definable from a signed in-order permutation: **Inversion** (sign-flip + local order reversal, per the reversal definition) and **Transposition** (block relocation preserving orientation, per Bafna–Pevzner). Other `RearrangementType` values remain in the enum for callers but are out of scope for these two methods. Documented as a "Not implemented" limitation, not an invented behavior.

---

## Recommendations for Test Coverage

1. **MUST Test:** Hunter worked example `α=(−2,−3,+1,+6,−5,−4)` yields exactly 6 breakpoints. — Evidence: Hunter Lecture 16 lines 287–289.
2. **MUST Test:** identity / collinear gene order yields 0 breakpoints (no events). — Evidence: "b_β(β)=0" (Hunter).
3. **MUST Test:** single reversed block `α=(+1,−4,−3,−2,+5)` yields exactly 2 breakpoints. — Evidence: reversal definition (Hunter line 256) + breakpoint criterion.
4. **MUST Test:** descending sign-consecutive pair `(−5,−4)` is NOT a breakpoint (criterion tests `(−y,−x)`). — Evidence: Hunter "(−5,−4) is not a breakpoint since (4,5) appears in β".
5. **MUST Test:** `ClassifyRearrangement` returns Inversion for a sign-flipped local reversal. — Evidence: reversal negates signs (Hunter line 256).
6. **MUST Test:** `ClassifyRearrangement` returns Transposition for a same-strand block relocation. — Evidence: Bafna–Pevzner transposition moves a block preserving orientation.
7. **MUST Test:** fewer than 2 mappable orthologs → no events (a permutation of < 2 markers has no internal adjacency). — Evidence: breakpoints are over consecutive pairs.
8. **MUST Test:** null inputs → ArgumentNullException for every public method. — Evidence: contract / sibling convention.
9. **SHOULD Test:** ortholog whose target is absent in genome2 is skipped, not a crash. — Rationale: robustness, mirrors sibling synteny unit.
10. **SHOULD Test:** the number of detected breakpoint events equals `d_BP = n − (common adjacencies)` consistency on a small case. — Rationale: cross-checks against Tannier formula.
11. **COULD Test:** property — breakpoint count is always in `[0, n+1]` for any input (n+1 internal pairs in the extended permutation). — Rationale: value-bound invariant.
12. **COULD Test:** property — identical inputs give 0 breakpoints regardless of size (idempotence of identity). — Rationale: invariant `b(β)=0`.

---

## References

1. Bafna V, Pevzner PA. (1998). Sorting by Transpositions. *SIAM Journal on Discrete Mathematics* 11(2):224–240. https://doi.org/10.1137/S089548019528280X — course copy retrieved: https://www.ic.unicamp.br/~meidanis/courses/mo640/2008s2/textos/Bafna-Pevzner-1998.pdf
2. Tannier E, Zheng C, Sankoff D. (2009). Multichromosomal median and halving problems under different genomic distances. *BMC Bioinformatics* 10:120 — cited within: "On the Complexity of Rearrangement Problems under the Breakpoint Distance", PMC full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3887456/
3. Hunter College CSCI Computational Biology, Lecture 16: Genome rearrangements, sorting by reversals (exposition of Hannenhalli–Pevzner / Bafna–Pevzner). https://www.cs.hunter.cuny.edu/~saad/courses/compbio/lectures/lecture16.pdf (accessed 2026-06-13).

---

## Change History

- **2026-06-13**: Initial documentation.
