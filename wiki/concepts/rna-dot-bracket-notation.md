---
type: concept
title: "RNA dot-bracket (extended WUSS) notation — parsing and validation"
tags: [rna, algorithm]
sources:
  - docs/Evidence/RNA-DOTBRACKET-001-Evidence.md
source_commit: 4c6115f113e0ef75414b050883a31e2594d48d35
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: rna-dotbracket-001-evidence
      evidence: "Test Unit ID: RNA-DOTBRACKET-001 ... Algorithm: Dot-Bracket (extended WUSS) Notation — parsing and validation"
      confidence: high
      status: current
---

# RNA dot-bracket (extended WUSS) notation — parsing and validation

The **notation layer** of the RNA secondary-structure family (test unit **RNA-DOTBRACKET-001**):
parse a structure string into the set of base pairs `(i, j)` it encodes, and validate that the
string is **well-formed** (every bracket has a correctly ordered, same-family partner). This is the
representation on which the folding units operate — the [[pre-mirna-hairpin-detection|pre-miRNA
hairpin detector]] emits a balanced `(`/`.`/`)` dot-bracket and the MFE folder (RNA-STRUCT-001)
reads its hairpin from the real dot-bracket structure. The record is [[rna-dotbracket-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the artifact
pattern. It is distinct from the chemistry of **which bases pair** ([[rna-base-pairing]]) — this unit
is purely about the **bracket string**, not the nucleotides.

## 1. The notation (ViennaRNA / WUSS)

- **Base pairs** are written as **matching pairs of brackets**; **unpaired** residues as dots `.`
  (ViennaRNA). Basic notation uses only `()` — e.g. `((((....))))` is one hairpin of 4 bp.
- **Extended (WUSS) notation** adds three more bracket families — `<>`, `{}`, `[]` — and
  **uppercase/lowercase letter pairs** `Aa`, `Bb`, …. Each family is an **independent pairing
  system**: `<` closes only with `>`, `[` only with `]`, `A` (5', opener) only with `a` (3', closer).
  The **exact choice of symbol carries no meaning** — Infernal/`vrna_db_from_WUSS()` flattens all
  brackets and treats any matched pair of `()`/`[]`/`{}` as a base pair.
- **Pseudoknots (crossing helices):** because families are matched independently, **different
  families are not required to nest** with each other. Crossing is what encodes a pseudoknot —
  `([)]` is two crossing 1-bp helices, and the three ViennaRNA equivalents
  `<<<<[[[[....>>>>]]]]`, `((((AAAA....))))aaaa`, `AAAA{{{{....aaaa}}}}` all encode the **same**
  pair of size-4 crossing helices.
- **Non-bracket WUSS symbols are all single-stranded (unpaired):** Rfam's `-` (internal loops /
  bulges), `,` (single strand between helices), `:` (external single-stranded), and `.` (insertions).
  A parser/validator must **ignore** these, not error on them.

## 2. Parsing (`ParseDotBracket`)

**One stack per bracket family.** Scanning left→right, an opener is pushed onto its family's stack;
a closer pops that family's stack to obtain its partner index, emitting the pair `(i, j)` with
`i < j`. Keeping a **separate stack per family** is what makes pseudoknots parse correctly — a single
shared stack would mis-pair `([)]`. Uppercase letters open, lowercase letters close (`A…a`).

Worked oracles (0-based, 5'→3'):

| Input | Base pairs |
|-------|------------|
| `((((....))))` | (0,11),(1,10),(2,9),(3,8) — outermost nests with outermost |
| `([)]` | `(` family (0,2); `[` family (1,3) — crossing |
| `<<<<[[[[....>>>>]]]]` | `<`: (0,15),(1,14),(2,13),(3,12); `[`: (4,19),(5,18),(6,17),(7,16) |
| `((((AAAA....))))aaaa` | `(`: (0,15),(1,14),(2,13),(3,12); `A`: (4,19),(5,18),(6,17),(7,16) |

## 3. Validation (`ValidateDotBracket`)

Well-formed ⟺ **every family's stack is empty at end of scan and never underflows** (no closer
without a matching earlier opener of the **same** family). A closer must match an opener of the
same family — `(]` is malformed even though one open and one close are present.

| Input | Valid? | Reason |
|-------|--------|--------|
| `(((...)))` | true | balanced, nested |
| `(([[]]))` | true | nested mixed families, each balanced |
| `([)]` | true | crossing families, each family balanced |
| `(((...)` | false | unclosed `(` |
| `...)` | false | `)` with no open partner |
| `)(` | false | `)` before any `(` |
| `(]` | false | `(` unclosed **and** `]` unopened (mismatched families) |

## Invariants and edge cases

- **INV:** each parsed pair has `i < j`; a well-formed string parses to exactly `(#openers)` pairs
  with balanced per-family stacks.
- **INV:** validation is per-family — validity requires each family independently balanced; crossing
  families are valid (pseudoknots), mismatched families are not.
- **Best-effort parse of malformed input** (Assumption 1): `ParseDotBracket` on a malformed string
  returns only the pairs it *can* match and **silently drops** unmatched closers rather than throwing;
  the documented contract is to gate with `ValidateDotBracket` first.
- **Empty/null** (Assumption 2): treated as a valid, pair-free structure —
  `ValidateDotBracket("") == true`, `ParseDotBracket("")` empty (the empty string is unambiguously
  balanced; ViennaRNA does not define it).

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for the **notation** only — it does not
fold, score, or assign energies (the [[rna-free-energy-turner-model|Turner 2004 free-energy terms]]
are a separate layer, summed by the RNA-STRUCT-001 MFE folder). Symbol identity is
semantically ignored (any matched `()`/`[]`/`{}`/letter pair is just "a base pair"), matching the
Infernal `vrna_db_from_WUSS()` flattening convention. **No source contradictions** — ViennaRNA
(readthedocs + tbi.univie.ac.at), the WUSS/`vrna_db_from_WUSS()` docs, Infernal (Nawrocki & Eddy
2013), and the Rfam glossary agree on the four bracket families, letter-pair pseudoknots, and the
single-stranded non-bracket symbols; the only recorded items are the two API-contract assumptions
(malformed best-effort parse; empty/null valid).
