---
type: source
title: "Evidence: COMPGEN-REARR-001 (Genome rearrangement detection by breakpoints, signed gene-order comparison)"
tags: [validation, comparative-genomics]
doc_path: docs/Evidence/COMPGEN-REARR-001-Evidence.md
sources:
  - docs/Evidence/COMPGEN-REARR-001-Evidence.md
source_commit: d4fdff436171d32c7e3ac8000f49cfda9e1e28fb
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: COMPGEN-REARR-001

The validation-evidence artifact for test unit **COMPGEN-REARR-001** — **genome rearrangement
detection by breakpoints** on a *signed gene-order permutation* (Hannenhalli–Pevzner /
Bafna–Pevzner theory). This is a **Comparative-genomics** family Evidence file and one instance of
the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the
algorithm, its definitions, invariants, worked oracles, and corner cases are summarized in the
dedicated concept [[genome-rearrangement-breakpoint-distance]]. It is the **signed-permutation /
breakpoint** counterpart to the block-signal [[synteny-and-rearrangement-detection]] (a different
formulation of the same "what rearrangements separate two genomes?" problem — see that concept). Its
sibling COMPGEN units are [[average-nucleotide-identity]],
[[ortholog-detection-reciprocal-best-hits]], [[conserved-gene-clusters-common-intervals]],
[[dot-plot-word-match]], and [[genome-comparison-core-dispensable]]. See [[test-unit-registry]] for
how units are tracked.

## What this file records

- **Online sources (all authority rank 1):**
  - **Hunter College Computational Biology, Lecture 16** ("Genome rearrangements, sorting by
    reversals"; exposition of Hannenhalli–Pevzner / Bafna–Pevzner) — the verbatim definitions: a
    **signed permutation** `α(i)=±a`; a **reversal** `α[i,j]` that reverses a contiguous block *and
    negates each sign*; the **extended permutation** `(0, α(1),…,α(n), n+1)`; the **breakpoint**
    criterion "if `(x,y)` appears in extended α but neither `(x,y)` nor `(−y,−x)` appears in extended
    β"; and the **reversal-distance lower bound** — a reversal changes the breakpoint count by ≤ 2,
    so `d(α) ≥ b(α)/2`. Worked example `α=(−2,−3,+1,+6,−5,−4)` vs identity → **b(α)=6** (with
    `(−5,−4)` *not* a breakpoint because `(4,5)∈β`); identity → **b(β)=0**. Retrieved via WebSearch →
    PDF extracted with `pdftotext` (poppler), read from the extracted text.
  - **Tannier, Zheng & Sankoff** ("On the Complexity of Rearrangement Problems under the Breakpoint
    Distance", PMC3887456; cites Tannier et al. 2009) — the **breakpoint / adjacency** vocabulary
    (adjacency = an edge between gene extremities head `gʰ`/tail `gᵗ`), the **breakpoint-distance
    formula** `d(π₁,π₂) = n − sim(π₁,π₂)` where `sim` = number of common adjacencies (plus half the
    common telomeric adjacencies), and **telomeric adjacencies** (chromosome ends = the extended
    endpoints `0` and `n+1`).
  - **Bafna & Pevzner** ("Sorting by Transpositions", *SIAM J. Discrete Math.* 11(2):224–240, 1998;
    UNICAMP course copy, `pdftotext`) — the operation classes ("inversions and transpositions as
    well as deletion, insertion, duplication"), the same extended-permutation model, and the
    **transposition** definition `ρ(i,j,k)` that *moves a block preserving its internal orientation*
    (no sign change) — the discriminator between Transposition and the sign-flipping Inversion
    (reversal).

- **Documented corner cases / failure modes:** identity / collinear order → `b=0` (no events);
  a descending sign-consecutive pair `(−5,−4)` mapping to identity adjacency `(4,5)` is **not** a
  breakpoint (the criterion must test both `(x,y)` and `(−y,−x)`); the extended endpoints `(0,π₁)`
  and `(πₙ,n+1)` can themselves be breakpoints; telomeres are the extended endpoints; `sim=0` →
  `d_BP=n` (every adjacency a breakpoint).

- **Datasets (documented oracles):**
  - *Hunter worked example* — `α=(−2,−3,+1,+6,−5,−4)`, β=identity, extended
    `(0,−2,−3,+1,+6,−5,−4,+7)` → 6 breakpoints `(0,−2)(−2,−3)(−3,+1)(+1,+6)(+6,−5)(−4,+7)`,
    `b(α)=6`, `(−5,−4)` excluded, reversal lower bound `d≥3`.
  - *Identity* — `α=β=(+1,+2,+3,+4,+5)` → `b(α)=0`.
  - *Single inversion* — reverse block [2,4] of identity → `α=(+1,−4,−3,−2,+5)`, extended
    `(0,+1,−4,−3,−2,+5,+6)` → `b(α)=2` (`(+1,−4)` and `(−2,+5)` are breakpoints; the internal
    `(−4,−3)`,`(−3,−2)` map to identity adjacencies via `(−y,−x)`), `d≥1`.

- **Test-coverage recommendations:** MUST — Hunter example → exactly 6 breakpoints; identity → 0;
  single reversed block → exactly 2; `(−5,−4)` NOT a breakpoint (criterion tests `(−y,−x)`);
  `ClassifyRearrangement` → **Inversion** for a sign-flipped local reversal and **Transposition**
  for a same-strand block relocation; `< 2` mappable orthologs → no events; null → `ArgumentNullException`.
  SHOULD — ortholog absent in genome2 skipped not crashed; detected-breakpoint count consistent with
  `d_BP = n − common adjacencies`. COULD — breakpoint count always in `[0, n+1]`; identical inputs → 0
  regardless of size.

## Deviations and assumptions

**Three ASSUMPTIONS**, all source-backed scoping decisions, none an open correctness gap:

- **Anchors supplied as an ordered ortholog mapping.** The permutation of common markers is derived
  from two gene lists + an `orthologMap` (the same input shape as the sibling
  [[ortholog-detection-reciprocal-best-hits|ORTHO]]/synteny units); anchor *generation* is delegated
  to those units, leaving the breakpoint computation itself unchanged.
- **Strand char `'+'/'-'` encodes the sign** of each element in the signed permutation (standard
  signed-permutation encoding, Hunter); the reversed-target gene carries the opposite strand.
- **Only Inversion and Transposition are classified.** These are the two operations definable from a
  single in-order signed permutation (Inversion = sign-flip + local order reversal per the reversal
  definition; Transposition = block relocation preserving orientation per Bafna–Pevzner).
  Translocation / Deletion / Insertion / Duplication require chromosome ids or gene-set differences a
  single permutation cannot express and no authoritative single-permutation rule assigns them, so
  they are a documented **"Not implemented"** limitation for these two methods (the enum retains the
  values for callers). This contrasts with the block-signal
  [[synteny-and-rearrangement-detection|CHROM-SYNT]] `DetectRearrangements`, which *does* classify
  Translocation/Deletion/Duplication from adjacent-block coordinates.

No contradictions among sources — Hunter (signed-permutation/breakpoint/reversal-bound), Tannier et
al. (adjacency vocabulary + `d=n−sim`), and Bafna–Pevzner (transposition vs inversion) are mutually
consistent, each governing a distinct part of the theory.
