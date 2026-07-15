---
type: source
title: "Evidence: ANNOT-REPEAT-001 (Repetitive Element Detection & Classification)"
tags: [validation, annotation]
doc_path: docs/Evidence/ANNOT-REPEAT-001-Evidence.md
sources:
  - docs/Evidence/ANNOT-REPEAT-001-Evidence.md
source_commit: 3ecff734f0954cc6cd463f353d3994e99068f93f
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: ANNOT-REPEAT-001

The validation-evidence artifact for test unit **ANNOT-REPEAT-001** (Repetitive Element
Detection and Classification — tandem repeats, inverted repeats, repeat-class assignment). One
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the algorithm and its three sub-problems are summarized in
[[repetitive-element-detection]], the shared anchor for the repeats/tandem family. See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources** (accessed 2026-06-13):
  - **Wikipedia "Tandem repeat"** (rank 4, citing Duitama et al. 2014, MeSH, Jorda et al. 2010)
    — head-to-tail adjacency definition, ≥2 copies minimum, STR 1–6 bp / minisatellite 10–60 bp
    unit sizes, and the non-primitive-unit corner case.
  - **Wikipedia "Inverted repeat"** (rank 4, citing Ussery et al. 2008, Ye et al. 2014) — arm +
    downstream reverse-complement definition, spacer of any length including zero, zero-spacer ⇒
    palindrome.
  - **IUPACpal, Hampson et al. 2021** (BMC Bioinformatics 22:51, PMC7866733, rank 1) — the formal
    IR grammar `W W̄ᴿ` (perfect), `W G W̄ᴿ` (gapped, `|G| ≥ 0`), imperfect within `k` via
    `δ_H(W, W̄ᴿ) ≤ k`, and the four detection parameters (min/max arm, max gap, max mismatches).
  - **RepeatMasker documentation** (Smit, Hubley & Green; Repbase, rank 3) — the repeat-class
    vocabulary (SINE/LINE/LTR/DNA/Satellite/Simple_repeat/Low_complexity/Small RNA/Unknown),
    homology-by-best-Smith-Waterman-Gotoh-match classification against Repbase, and the
    no-match ⇒ Unclassified/Unknown behaviour.
- **Datasets** — tandem `ATTCGATTCGATTCG` → unit `ATTCG`, 3 copies, span `[0,15)`; inverted
  `GAATTC` → arms `GAA`/`TTC` revcomp, gap 0 (palindrome); gapped inverted
  `TTACGAAAAAACGTAA` → left `TTACG`, right `CGTAA` = revcomp, gap 6; motif-size → class table
  (`A`/`CA`/`CAG` = Simple_repeat mono/di/tri, 7+ bp = Satellite/minisatellite).
- **Corner cases** — single copy (copies = 1) is not a tandem repeat; prefer the primitive
  (shortest-period) unit; zero-gap IR = even-length palindrome; odd-length palindrome leaves the
  centre base unpaired; no library match above threshold ⇒ Unknown (never force a class).
- **Recommended coverage** — six MUST tests (the two tandem cases, the two IR cases, exact
  library-class match, no-match fallback), plus SHOULD (primitive-unit `AAAAAA` → `A`;
  null/empty/too-short) and COULD (the `sequence == input[start..end]` integer-copies invariant).

## Assumption (from the artifact)

**`ClassifyRepeat` library matching is exact-substring containment, not Smith-Waterman
homology.** RepeatMasker scores against a curated Repbase library; a full local-alignment +
Repbase library is out of scope for one unit. Seqeron screens the query for library elements
**exactly contained** within it (element ⊆ query, one-directional) and assigns the class of the
**longest** match, falling back to motif-size Simple_repeat classification (else Unknown) when
nothing matches. Only the matching relaxation (exact substring vs. scored alignment) is assumed;
the SINE/LINE/LTR/DNA/Satellite/Simple_repeat/Unknown vocabulary is source-backed. Documented as
a Framework/Simplified limitation, not an invented constant.

No contradictions among the four sources — the Wikipedia and IUPACpal IR definitions are the same
grammar at different formality, and the RepeatMasker class list is the shared vocabulary. The only
deviation is the classification-matching relaxation above, which does not affect the tandem-repeat
or inverted-repeat detection correctness.
