# Validation Report: PROTMOTIF-PROSITE-001 — PROSITE Pattern Matching

- **Validated:** 2026-06-24   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.ConvertPrositeToRegex(prositePattern)`,
  `ProteinMotifFinder.FindMotifByProsite(sequence, pattern, name)`
  (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs:245-425)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

This is an independent re-validation (fresh context). The unit was previously validated at
cb113ce (PASS/PASS) with a hardening fix (reject unsupported `*` rather than silently drop).
This session re-derived every PROSITE syntax rule from the authoritative ExPASy manual and
re-confirmed the implementation by hand; no defect found, no code changed.

---

## Stage A — Description

### Sources opened & what they confirm

- **ScanProsite Documentation / PROSITE User Manual §IV.E (PA-line grammar)** —
  https://prosite.expasy.org/scanprosite/scanprosite_doc.html (fetched 2026-06-24). Verbatim
  rules extracted and used to validate the spec/code:
  - "Each element in a pattern is separated from its neighbor by a '-'." → `-` is a separator (removed).
  - "The symbol 'x' is used for a position where any amino acid is accepted." → `x` ≡ `.`.
  - "[ALT] stands for Ala or Leu or Thr." → `[ABC]` ≡ regex `[ABC]` (any-of).
  - "between a pair of curly brackets '{ }' the amino acids that are **not** accepted … {AM}
    stands for all any amino acid **except** Ala and Met." → `{ABC}` ≡ `[^ABC]` (EXCLUSION, not
    inversion of anything else).
  - "x(3) corresponds to x-x-x", "A(3) corresponds to A-A-A" → fixed `(n)` repetition valid on
    any element → regex `{n}`.
  - "x(2,4) corresponds to x-x or x-x-x or x-x-x-x" and "**Ranges can only be used with 'x'**,
    for instance 'A(2,4)' is not a valid pattern element." → range `(n,m)` ≡ `.{n,m}`, PROSITE-legal
    only on `x`.
  - "a pattern restricted to the N- or C-terminal … starts with '<' or ends with '>'." →
    `<`≡`^`, `>`≡`$`.
  - "'>' can also occur inside square brackets for the C-terminal element. 'F-[GSTV]-P-R-L-[G>]'
    is equivalent to 'F-[GSTV]-P-R-L-G' or 'F-[GSTV]-P-R-L>'." → `[G>]` ≡ `(?:G|$)`.
  - A period terminates the PA-line pattern. 1-based residue numbering in ScanProsite output.

### Operator-by-operator verification

| PROSITE element | Meaning (manual, verbatim) | Correct regex | Verified |
|---|---|---|---|
| `-` | element separator | (removed) | ✅ |
| `x` | any amino acid | `.` | ✅ |
| `[ABC]` | one of A/B/C | `[ABC]` | ✅ |
| `{ABC}` | any **except** A/B/C | `[^ABC]` | ✅ |
| `x(n)` / `A(n)` | n consecutive copies | `.{n}` / `A{n}` | ✅ |
| `x(m,n)` | m..n copies (x only) | `.{m,n}` | ✅ |
| `<` | N-terminus | `^` | ✅ |
| `>` | C-terminus | `$` | ✅ |
| `[G>]` | G **or** C-term | `(?:G\|$)` | ✅ |
| `.` | end of pattern | terminate | ✅ |
| `*` (extended query) | not in PA-line grammar | rejected (FormatException) | ✅ |

### Worked examples (hand-computed against the manual)

- **N-glycosylation PS00001 `N-{P}-[ST]-{P}` → `N[^P][ST][^P]`.** On `AANASAAANGTAAA`
  (0-based): index 2 = `N,A,S,A` (A≠P, S∈[ST], A≠P) → match @2; index 8 = `N,G,T,A` → match @8.
  On `AANPSAAA` the P at the `{P}` position blocks the match → 0 matches. Matches manual canonical example.
- **`{AM}` (manual's own example) → `[^AM]`.** Matches any residue except A and M.
- **Zinc finger PS00028 `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` →
  `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H`** — exercises both fixed `x(n)` and ranged `x(m,n)`.
- **`[G>]` PS00267/PS00539 → `(?:G|$)`**: a literal G at that position OR end-of-sequence.

### Findings / divergences

- The implementation is *more permissive* than strict PROSITE: it will also convert an element
  range such as `[RK](2,4)` → `[RK]{2,4}` even though the manual says ranges are valid only on
  `x`. This never mis-handles a valid PROSITE pattern (those never carry an element range) and
  always yields valid .NET regex (INV-2), so it is benign permissiveness, not a defect.

**Stage A: PASS** — the description (operators, exclusion semantics, anchors, `[G>]` special
case, terminator, 1-based convention) matches the authoritative ExPASy source exactly.

---

## Stage B — Implementation

### Code path reviewed

`ConvertPrositeToRegex` (ProteinMotifFinder.cs:257-425); `FindMotifByProsite` (:245-252) →
`FindMotifByPattern` (:178-240).

### Formula realised correctly? (evidence)

- `-` skipped (:269-273); `x` / `x(...)` → `.` / `.{...}` (:274-298); `[...]` preserved, with the
  `[G>]`/`[<…]` special case → `(?:letters|$)` / `(?:^|letters)` (:299-351); `{...}` → `[^...]`
  exclusion correctly negated (:352-367); bare `<`→`^`, `>`→`$` (:368-379); element `(n)`/`(m,n)`
  → `{n}`/`{m,n}` (:380-394); letters uppercased literals (:395-400); `.` breaks the loop
  terminating the pattern (:401-405); any other char (notably `*`, `?`, `+`) → `FormatException`
  naming the offending char/position/pattern and listing the supported grammar (:406-421).
- **Matching & position base:** `FindMotifByPattern` wraps the regex in a lookahead `(?=(regex))`
  (:198) for overlapping discovery (ScanProsite-consistent), folds case via
  `RegexOptions.IgnoreCase` + uppercasing, skips zero-width captures (:226-227), and reports
  `Start`/`End` as **0-based inclusive** (:232-233). A 2 s match timeout guards catastrophic
  backtracking → no matches rather than a hang.

### Cross-verification (recomputed against the running code)

| Input | Expected | Observed | ✅ |
|---|---|---|---|
| `N-{P}-[ST]-{P}` | `N[^P][ST][^P]` | `N[^P][ST][^P]` | ✅ |
| `[ST]-x-[RK]` | `[ST].[RK]` | `[ST].[RK]` | ✅ |
| `{AM}` | `[^AM]` | `[^AM]` | ✅ |
| PS00028 full | `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H` | identical | ✅ |
| `<M-x-K` | `^M.K` | `^M.K` | ✅ |
| `A-x-G>` | `A.G$` | `A.G$` | ✅ |
| `F-[GSTV]-P-R-L-[G>].` | `F[GSTV]PRL(?:G\|$)` | identical | ✅ |
| `R-G-D.A-B-C` (period terminates) | `RGD` | `RGD` | ✅ |
| `<{C}*>`, `A-x*-B` (`*`) | FormatException | FormatException | ✅ |

Match-position cross-checks (via direct runs, not just the test asserts):
PS00001 on `AANASAAANGTAAA` → matches @2 (`NASA`) and @8 (`NGTA`); on `AANPSAAA` → 0;
`{AM}` on `ALM` → matches `L` only; PS00267 `[G>]` → matches via the G branch (`FVGLMG`),
via the C-term branch (`FVGLM` at end), and **rejects mid-sequence** (`FVGLMAAAA` → 0).

### Variant / delegate consistency

`FindMotifByProsite` delegates to `FindMotifByPattern` after `ConvertPrositeToRegex`;
`FindCommonMotifs` (:150-164) iterates the `CommonMotifs` table whose stored `RegexPattern`
strings were spot-checked to equal `ConvertPrositeToRegex(Pattern)` (e.g. PS00001, PS00005,
PS00028, PS00008, PS00018). Consistent.

### Test quality audit

`ProteinMotifFinder_PrositePattern_Tests.cs` — 90 tests run, all green. Assertions check exact
regex strings and exact 0-based match positions (not "no-throw"); the `*`-rejection tests assert
`FormatException` with the offending char in the message; the `[G>]` tests cover all three
branches (G / end / mid-seq reject); S1 verifies Human Transferrin P02787 → exactly 2 PS00001
sites at 0-based 431/629 (= ScanProsite 1-based 432/630).

### Findings / defects

None. Implementation faithfully realises the validated PROSITE PA-line grammar.

---

## Verdict & follow-ups

- **Stage A: PASS. Stage B: PASS. End state: CLEAN.**
- No code changed (validation only). Full suite: **18213 passed / 0 failed** (Prosite filter: 90/90).
- Note (informational, not a defect): the converter accepts element ranges `[X](m,n)` that
  strict PROSITE forbids (ranges are x-only); benign because no valid PROSITE pattern uses them
  and output stays valid regex.
