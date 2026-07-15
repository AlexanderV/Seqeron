---
type: source
title: "Evidence: PROTMOTIF-PROSITE-001 (PROSITE pattern matching — ConvertPrositeToRegex + FindMotifByProsite, [G>] C-term-in-brackets)"
tags: [validation, protein, motif]
doc_path: docs/Evidence/PROTMOTIF-PROSITE-001-Evidence.md
sources:
  - docs/Evidence/PROTMOTIF-PROSITE-001-Evidence.md
source_commit: 0908c5f04255fcb3d51c7706d74a689eee481faa
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROTMOTIF-PROSITE-001

The validation-evidence artifact for test unit **PROTMOTIF-PROSITE-001** — **PROSITE Pattern
Matching** (`ConvertPrositeToRegex`, `FindMotifByProsite`). It is a **third Evidence artifact
over the same PROSITE→regex engine** as [[protmotif-find-001-evidence]] (PROTMOTIF-FIND-001)
and [[protmotif-pattern-001-evidence]] (PROTMOTIF-PATTERN-001): same primitives, same PA-line
grammar, no contradictions. The engine model, PROSITE→regex table, overlap mechanism and
scoring live in [[protein-motif-pattern-search]]. One instance of the templated
[[algorithm-validation-evidence|evidence artifact]] pattern; see [[test-unit-registry]].

Its distinct contributions are (1) the **`[G>]` C-terminus-inside-brackets** corner case
(PS00267, PS00539) → regex alternation `(?:G|$)`; (2) a real published oracle — **Human
Transferrin (P02787) scanned against PS00001** via ScanProsite, 2 N-glycosylation sites; and
(3) the ScanProsite extended-syntax / match-mode notes.

## What this file records

- **Online sources (all rank 1–2):**
  - **PROSITE User Manual — PA line spec** (rank 2): the full PA-line grammar — IUPAC letters,
    `x`=any, `[ALT]` class, `{AM}` exclusion, `-` separator, `(n)`/`(n,m)` repetition (**range
    only valid with `x`**, fixed `(n)` valid on any element), `<` N-terminus / `>` C-terminus
    anchors, `>` **inside** brackets (`[G>]` = G-or-end), and a trailing `.` that ends the PA
    line.
  - **ScanProsite documentation** (rank 2): extended syntax (separators `-` may be omitted when
    the pattern has no ambiguous residues, e.g. `MASKE` = `M-A-S-K-E`); three match modes
    (**greedy** default, **overlap**, **include**); `A(2,4)` is not a valid element.
  - **PROSITE entries (rank 2):** PS00001 ASN_GLYCOSYLATION `N-{P}-[ST]-{P}.`; PS00028
    ZINC_FINGER_C2H2_1 `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H.` (max repeat 36 fingers/
    protein).
  - **ScanProsite verified scan (rank 2):** P02787 (TRFE_HUMAN) × PS00001 JSON output — 2 hits.
  - **Hulo et al. 2007** (rank 1, PROSITE database) and **De Castro et al. 2006** (rank 1,
    ScanProsite — overlapping matches supported via match-mode parameters).

- **`[G>]` C-terminus inside brackets (NEW corner case):** in rare patterns `>` appears inside
  a character class meaning "that residue **or** end-of-sequence". `ConvertPrositeToRegex`
  emits a regex alternation `(?:G|$)` rather than a plain class. Verified against
  **PS00267** `F-[IVFY]-G-[LM]-M-[G>].` → `F[IVFY]G[LM]M(?:G|$)` and **PS00539**
  `F-[GSTV]-P-R-L-[G>].` → `F[GSTV]PRL(?:G|$)`. (Added 2026-02-13 per the file's change
  history; "removed all assumptions".)

- **Period terminates parsing (sharpened):** a `.` ends the pattern (User Manual §IV.E) — not
  just a trailing decoration. `R-G-D.` → `RGD`, and `R-G-D.A-B-C` → `RGD` (content after the
  period is dropped).

## Test datasets (exact oracles)

**Syntax conversion** (selection; full table in the source):

| PROSITE | Regex | Source |
|---------|-------|--------|
| `A-x(2,4)-G` | `A.{2,4}G` | User Manual (range on `x`) |
| `<M-x-K` | `^M.K` | N-terminus anchor |
| `A-x-G>` | `A.G$` | C-terminus anchor |
| `[RK](2)-x-[ST]` | `[RK]{2}.[ST]` | PS00004 (fixed count on a class) |
| `R-G-D.A-B-C` | `RGD` | period terminates parsing |
| `F-[IVFY]-G-[LM]-M-[G>].` | `F[IVFY]G[LM]M(?:G\|$)` | PS00267 `[G>]` |
| `F-[GSTV]-P-R-L-[G>].` | `F[GSTV]PRL(?:G\|$)` | PS00539 `[G>]` |

**Matching** (0-based starts): `AAARGDAAA` × `R-G-D` → [3]; `AANASAAANGTAAAA` × `N-{P}-[ST]-{P}`
→ [2, 8]; `AANPSAAA` × `N-{P}-[ST]-{P}` → [] (Pro at the `{P}` position); `AASARKAAA` ×
`[ST]-x-[RK]` → [2]; empty sequence or empty pattern → [].

**Real published data:** Human Transferrin **P02787** (TRFE_HUMAN) × PS00001 `N-{P}-[ST]-{P}`
→ 2 matches, ScanProsite 1-based 432–435 and 630–633 (= 0-based 431–434 and 629–632).

## Recommended coverage

MUST: `ConvertPrositeToRegex` for each atom (`x`, `x(n)`, `x(n,m)`, `[..]`, `{..}`, `A(n)`,
`<`, `>`, trailing `.`, mid-pattern `.` termination) and each real pattern (PS00016/00001/
00028/00008/00018); the **`[G>]` C-term-in-brackets** cases (PS00267, PS00539); empty input →
empty. `FindMotifByProsite` matches PS00001/PS00016 at the exact positions, is
case-insensitive, finds multiple matches, and handles the `[G>]` pattern via **both** the `G`
branch and the C-terminus branch (and rejects it mid-sequence without `G`). SHOULD: Human
Transferrin PS00001 real-protein oracle; N-/C-terminal anchors match only at start/end.

## Relationship to the other PROSITE units

No contradictions with PROTMOTIF-FIND-001 or PROTMOTIF-PATTERN-001. Those two established the
three primitives, information-content scoring, overlapping-lookahead enumeration, and corrected
two impl patterns (PS00007, PS00018). PROSITE-001 revalidates `ConvertPrositeToRegex` /
`FindMotifByProsite` and adds the `[G>]` alternation rule, the mid-pattern period-termination
oracle, and the P02787 real-protein positive control. All three keep 0-based
`MotifMatch.Start`/`End` vs ScanProsite's 1-based coordinates (no correctness effect).
