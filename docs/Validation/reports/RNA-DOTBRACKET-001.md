# Validation Report: RNA-DOTBRACKET-001 — Dot-Bracket (extended WUSS) Notation

- **Validated:** 2026-06-16   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.ParseDotBracket(string)`, `RnaSecondaryStructure.ValidateDotBracket(string)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (not trusting repo citations)

1. **ViennaRNA — RNA Structure Notations** (https://viennarna.readthedocs.io/en/latest/io/rna_structures.html, fetched 2026-06-16). Verbatim: dot-bracket "denotes base pairs by matching pairs of parenthesis `()` and unpaired nucleotides by dots `.`". Extended version "may use additional pairs of brackets, such as `<>`, `{}`, and `[]`, and matching pairs of uppercase/lowercase letters." Different bracket pairs are **not required to be nested** (enables pseudoknots); within a family standard LIFO nesting applies. Three equivalent crossing-helix encodings quoted verbatim: `<<<<[[[[....>>>>]]]]`, `((((AAAA....))))aaaa`, `AAAA{{{{....aaaa}}}}`.
2. **ViennaRNA (web search of the same doc set)** confirmed verbatim: "crossing pairs may be annotated by matching uppercase/lowercase letters from the alphabet A-Z, where **the uppercase letter must be the 5' and the lowercase letter the 3'** nucleotide of the base pair." This independently grounds letter direction (S1) — not taken from the repo Evidence.
3. **ViennaRNA — WUSS notation** (https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/utils/struct/wuss.html, fetched 2026-06-16): `vrna_db_from_WUSS()` "flatten[s] all brackets, and treats pseudo-knots annotated by matching pairs of upper/lowercase letters as unpaired nucleotides" — confirms each family is an independent pairing system.
4. **Rfam Glossary — WUSS** (https://docs.rfam.org/en/latest/glossary.html, fetched 2026-06-16): base-pair symbols `<>` (simple stem loops) and `()`,`[]`,`{}` (enclosing multifurcations); single-stranded symbols `-` (internal loops/bulges), `,` (single strand between helices), `:` (external), `.` (insertions).
5. **Infernal/WUSS spec** (web search, 2026-06-16): full unpaired set incl. `_` (hairpin loops), `-` (bulge/interior), `,` (multifurcation), `:` (external), `.` (insertion), `~` (unaligned). Base-pair depths `<>`,`()`,`[]`,`{}`.

### Formula / model check
- Base pair `(i,j)`, `i<j` = opening symbol at `i`, matching closing symbol of the **same family** at `j`. Confirmed by source 1.
- Four bracket families `()[]{}<>` + uppercase/lowercase letter pairs, each an **independent** LIFO stack → crossing helices (pseudoknots). Confirmed by sources 1, 3.
- Uppercase = 5' opener, lowercase = 3' closer. Confirmed verbatim by source 2.
- All non-bracket/non-letter symbols (`.`,`-`,`,`,`:`,`_`,`~`) are single-stranded. Confirmed by sources 4, 5.
- "Partners must match up" (closer matches an opener of the same family) → `(]` is ill-formed. Confirmed by sources 1, 3.

### Edge-case semantics
- Empty/null = balanced empty structure (valid, zero pairs). Not defined by ViennaRNA but unambiguous under the balanced-bracket definition; documented as an accepted contract assumption. OK.
- Best-effort parse of malformed input (stray closer dropped, no throw) — a contract decision, not a sourced numeric value; flagged as an assumption. OK.

### Independent cross-check (hand-computed from source 1's verbatim examples)
The scientific claim is that the three encodings represent the **same** base-pair set. Hand-tracing per-family LIFO:

| Encoding | Base-pair set (0-based) |
|----------|-------------------------|
| `<<<<[[[[....>>>>]]]]` | `<`: (0,15),(1,14),(2,13),(3,12); `[`: (4,19),(5,18),(6,17),(7,16) |
| `((((AAAA....))))aaaa` | `(`: (0,15),(1,14),(2,13),(3,12); `A`: (4,19),(5,18),(6,17),(7,16) |

Both yield `{(0,15),(1,14),(2,13),(3,12),(4,19),(5,18),(6,17),(7,16)}` — identical, confirming equivalence (M3). `((((....))))` → {(0,11),(1,10),(2,9),(3,8)} (M1). `([)]` → independent stacks give {(0,2),(1,3)} (M2).

### Findings / divergences
None. The description (`docs/algorithms/RnaStructure/Dot_Bracket_Notation.md`), TestSpec, and Evidence all match the external sources exactly. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs:1649–1762`.

- `OpeningToClosing`/`ClosingToOpening` maps cover the four families (lines 1649–1663).
- `ParseDotBracket` (1678): per-symbol `Dictionary<char, Stack<int>>`; opening bracket or uppercase letter → push index; closing bracket → pop its family's stack and yield `(opener,i)`; lowercase letter → map to uppercase, pop, yield; everything else skipped. Best-effort: a closer with an empty/absent stack is silently ignored.
- `ValidateDotBracket` (1724): same scan; a closer with no matching open stack ⇒ `false`; any non-empty stack after the scan ⇒ `false`; else `true`. null/empty ⇒ `true`.

### Formula realised correctly?
Yes. Per-family independent LIFO stacks are exactly the model Stage A validated. Verified by tracing each external example: M1/M2/M3/M6/S1/S2/S3 parse to the hand-computed sourced sets; M4 (`(((...)))`,`(([[]]))`,`([)]`,`....`,`""`,null → true) and M5 (`(((...)`,`...)`,`)(`,`(]` → false) trace correctly. `(]`: `(` pushed onto `(`-stack; `]` maps to opener `[`, no `[`-stack present ⇒ `false` — correctly rejects the mismatched family that a single-counter validator would accept.

### Cross-verification table recomputed vs code (full suite run)
All 10 unit tests pass; values match the externally sourced expectations (not code echoes). Full suite: **6592 passed, 0 failed, 0 skipped**.

### Variant/delegate consistency
N/A — two canonical static methods, no `*Fast`/instance variants. The pre-existing duplicate parse/validate tests in `RnaSecondaryStructureTests.cs` were removed during consolidation; remaining `DotBracket`/`DotBracketNotation` references there belong to other units (stem-loop, MFE) and are out of scope.

### Test quality audit (HARD gate)
- **Sourced, not echoed:** every expected pair set / boolean traces to ViennaRNA verbatim examples or the WUSS "partners must match up" rule, hand-verified above. A deliberately-wrong single-stack parser would fail M2/M3 and a single-counter validator would fail M5 (`(]`) — the tests discriminate against the classic blind spots.
- **No green-washing:** exact `Is.EquivalentTo` full pair sets (no `Contains`/`Greater`/ranges), exact boolean assertions, no skips/ignores, no widened tolerances.
- **Coverage:** every public method and every code branch exercised — open bracket, close bracket, uppercase open, lowercase close, unpaired skip; crossing families; letter/bracket equivalence; WUSS unpaired symbols; best-effort stray closer; empty/all-dots/null; unclosed/unopened/reversed/mismatched-family rejection.
- **Honest green:** FULL unfiltered suite `Failed: 0`; `dotnet build` 0 errors. (4 pre-existing NUnit-analyzer warnings live in unrelated `ApproximateMatcher_EditDistance_Tests.cs`; no file was changed this session.)

Gate result: **PASS.**

### Findings / defects
None. No code or test changes were required.

## Verdict & follow-ups
- **Stage A: PASS.** Description matches ViennaRNA/WUSS/Rfam/Infernal sources exactly.
- **Stage B: PASS.** Implementation faithfully realises the validated model; tests are exact, sourced, and cover all branches/edge cases.
- **End-state: CLEAN** — no defect found; algorithm fully functional. No defect logged.
