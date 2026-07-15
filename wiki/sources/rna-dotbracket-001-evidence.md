---
type: source
title: "Evidence: RNA-DOTBRACKET-001 (RNA dot-bracket / extended WUSS notation — parse & validate)"
tags: [validation, rna]
doc_path: docs/Evidence/RNA-DOTBRACKET-001-Evidence.md
sources:
  - docs/Evidence/RNA-DOTBRACKET-001-Evidence.md
source_commit: 4c6115f113e0ef75414b050883a31e2594d48d35
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: RNA-DOTBRACKET-001

The validation-evidence artifact for test unit **RNA-DOTBRACKET-001** — **Dot-Bracket (extended
WUSS) Notation parsing and validation** (`ParseDotBracket` / `ValidateDotBracket` on
`RnaSecondaryStructure`). This is the **notation/representation layer of the RNA secondary-structure
family** and one instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern; the synthesizing concept is [[rna-dot-bracket-notation]]. [[test-unit-registry]]
tracks the unit.

The parsed pairs are RNA base pairs, so the unit is adjacent to the base-pairing chemistry primitive
[[rna-base-pairing]], and the balanced-dot-bracket structure it validates is exactly what the
[[pre-mirna-hairpin-detection|pre-miRNA hairpin]] unit emits and the RNA-STRUCT-001 MFE folder reads.

## What this file records

- **Online sources (all rank 3–5, reference-implementation / curated-database docs):**
  - **ViennaRNA — RNA Structure Notations** (readthedocs, rank 3): base pairs = matching `()`,
    unpaired = `.`; extended notation adds `<>`, `{}`, `[]` and uppercase/lowercase letter pairs;
    different bracket families are **not required to nest** (enables pseudoknots); uppercase = 5'
    opener, lowercase = 3' closer.
  - **ViennaRNA — Dot-Bracket Notation** (tbi.univie.ac.at, rank 3): basic `((((....))))`; the three
    **verbatim equivalent** crossing encodings `<<<<[[[[....>>>>]]]]`, `((((AAAA....))))aaaa`,
    `AAAA{{{{....aaaa}}}}`.
  - **ViennaRNA — WUSS / `vrna_db_from_WUSS()`** (rank 3) + **Infernal** (Nawrocki & Eddy 2013,
    rank 3): nested base pairs annotated by `<>`/`()`/`{}`/`[]`, pseudoknots by upper/lowercase
    letter pairs; **any matched pair indicates a base pair, the exact symbol has no meaning** as long
    as partners match; `vrna_db_from_WUSS()` flattens all brackets and treats letter-pair pseudoknots
    as unpaired — confirming each family is an **independent pairing system**.
  - **Rfam glossary** (EBI/Rfam, rank 5): `<>` = simple stem-loop pairs, `()`/`[]`/`{}` = pairs
    enclosing multifurcations; the non-bracket symbols `-` (loops/bulges), `,` (single strand between
    helices), `:` (external single-stranded), `.` (insertions) are **all single-stranded (unpaired)**.

- **Documented corner cases / failure modes:** crossing bracket families are valid (one stack per
  family; a shared stack mis-pairs `([)]`); a closer must match an opener of the **same** family
  (`(]` is malformed); uppercase = opening (5'), matching lowercase = closing (3'); non-bracket WUSS
  symbols (`-`,`,`,`:`,`.`) must be ignored, not errored.

- **Datasets / oracles:**
  - **Parsing** (0-based, 5'→3'): `((((....))))` → (0,11),(1,10),(2,9),(3,8); `([)]` → `(`:(0,2) and
    `[`:(1,3); `<<<<[[[[....>>>>]]]]` and `((((AAAA....))))aaaa` → `<`/`(`:(0–3,15–12) and
    `[`/`A`:(4–7,19–16).
  - **Validation:** `(((...)))` / `(([[]]))` / `([)]` → **true**; `(((...)` (unclosed) / `...)`
    (unopened) / `)(` (reversed) / `(]` (mismatched families) → **false**.

- **Test-coverage recommendations:** MUST — simple nested hairpin parses to exact positions; crossing
  families parse independently; validation accepts balanced/nested + crossing and rejects
  unclosed/unopened/reversed/mismatched. SHOULD — uppercase/lowercase letter pairs (uppercase = 5'
  opener); non-bracket WUSS symbols + dots treated as unpaired. COULD — empty/null handling.

## Deviations and assumptions

Two recorded **API-contract assumptions** (behavior for input the sources leave undefined; neither
invents a numeric value):

1. **Best-effort parse of malformed input** — `ParseDotBracket` on a malformed string yields only the
   pairs it can match and silently drops unmatched closers rather than throwing; the documented
   contract is to test well-formedness with `ValidateDotBracket` first (the sources define behavior
   only for well-formed notation).
2. **Empty/null string** — treated as a valid, pair-free structure (`ValidateDotBracket("") == true`,
   `ParseDotBracket("")` empty); the empty string is unambiguously balanced under the balanced-bracket
   definition, which ViennaRNA does not itself define.

**No source contradictions** — ViennaRNA, WUSS/Infernal, and the Rfam glossary agree on the four
bracket families, letter-pair pseudoknots, independent per-family matching, and the single-stranded
non-bracket symbols.
