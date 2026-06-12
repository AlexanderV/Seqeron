# Validation Report: PROTMOTIF-PROSITE-001 — PROSITE Pattern Matching

- **Validated:** 2026-06-12   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.ConvertPrositeToRegex(prositePattern)`,
  `ProteinMotifFinder.FindMotifByProsite(sequence, pattern, name)`
  (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm

- **ExPASy ScanProsite / PROSITE User Manual §IV.E** —
  https://prosite.expasy.org/scanprosite/scanprosite_doc.html (fetched 2026-06-12).
  Confirmed verbatim the PROSITE pattern grammar used by the spec and implementation:
  - Elements separated by `-` (hyphen); standard IUPAC one-letter codes.
  - `x` = any amino acid.
  - `[ABC]` = one of the listed residues (ALLOWED set).
  - `{ABC}` = any residue EXCEPT those listed (EXCLUDED set).
  - `A(n)` = n consecutive copies (`A-A-A` for `A(3)`); `x(n)` likewise.
  - `x(m,n)` = a variable run of m..n residues; **ranges `(m,n)` are valid only with `x`**,
    not with fixed residues (`A(2,4)` is illegal PROSITE).
  - `<` at the start = restriction to the N-terminus; `>` at the end = restriction to the
    C-terminus.
  - **`]>`/`[...>]` special case:** `>` may appear inside brackets for the C-terminal
    element — `F-[GSTV]-P-R-L-[G>]` means "either `F-[GSTV]-P-R-L-G` or `F-[GSTV]-P-R-L>`"
    (i.e. `G` **or** end-of-sequence). Verified directly from the manual text.
  - A period `.` terminates the pattern (PA-line convention).
  - 1-based residue numbering in ScanProsite output.

### Operator-by-operator verification (each mapped to matching)

| PROSITE element | Meaning (manual) | Correct regex | Verified |
|---|---|---|---|
| `-` | separator | (removed) | ✅ |
| `x` | any aa | `.` | ✅ |
| `[ABC]` | one of A/B/C | `[ABC]` | ✅ |
| `{ABC}` | any **except** A/B/C | `[^ABC]` | ✅ (EXCLUDE, not inverted) |
| `x(n)` | exactly n | `.{n}` | ✅ |
| `x(m,n)` | m..n | `.{m,n}` | ✅ |
| `A(n)` / `[AB](n)` | element repeated n | `A{n}` / `[AB]{n}` | ✅ |
| `<` | N-terminus | `^` | ✅ |
| `>` | C-terminus | `$` | ✅ |
| `[G>]` | G or C-term | `(?:G\|$)` | ✅ |
| `.` | end of pattern | terminate | ✅ |

### Worked examples (hand-computed)

- **N-glycosylation `N-{P}-[ST]-{P}` → `N[^P][ST][^P]`.**
  Sequence `AANASAAANGTAAA` (0-based): index 2 = `N,A,S,A` (A≠P, S∈[ST], A≠P) → **match @2**;
  index 8 = `N,G,T,A` (G≠P, T∈[ST], A≠P) → **match @8**. Expected matches {2, 8} — confirms
  spec M15 and the canonical manual example. `AANPSAAA` → the P at the {P} position blocks the
  match → **0 matches** (spec M17). ✅
- **ATP/GTP P-loop PS00017 `[AG]-x(4)-G-K-[ST]` → `[AG].{4}GK[ST]`.**
  6+ residues: first ∈{A,G}, four arbitrary, then G, K, and S/T — matches the published
  Walker-A motif. ✅ (present in `CommonMotifs["PS00017"]`).
- **PKC site PS00005 `[ST]-x-[RK]` → `[ST].[RK]`.** ✅
- **Zinc finger PS00028 `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` →
  `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H`** (exercises both fixed and ranged `x`). ✅

### Findings / divergences

- The Evidence/TestSpec correctly note that ranges `(m,n)` are PROSITE-legal only with `x`.
  The implementation is *more permissive*: it will also convert `[RK](2,4)` → `[RK]{2,4}`.
  This never mis-handles a valid PROSITE pattern and always yields valid .NET regex (INV-2),
  so it is a benign permissiveness, not a defect.

Stage A description (operators, semantics, anchors, terminator, 1-based convention) matches the
authoritative ExPASy source exactly. **Stage A: PASS.**

---

## Stage B — Implementation

### Code path reviewed

`ConvertPrositeToRegex` (ProteinMotifFinder.cs:215-371) and `FindMotifByProsite`
(:203-210) → `FindMotifByPattern` (:161-198).

### Formula realised correctly? (evidence)

- `-` skipped (:227); `x`/`x(...)` → `.`/`.{...}` (:232-256); `[...]` preserved, with the
  `[G>]`/`[<...]` special case producing `(?:letters|$)` / `(?:^|letters)` (:257-309);
  `{...}` → `[^...]` (:310-325) — **exclusion, correctly negated**; bare `<`→`^`, `>`→`$`
  (:326-337); `(n)`/`(m,n)` after an element → `{n}`/`{m,n}` (:338-352); letters uppercased
  literals (:353-358); `.` breaks the loop terminating the pattern (:359-363).
- **Matching & position base:** `FindMotifByPattern` wraps the regex in a lookahead
  `(?=(regex))` (:176) so overlapping occurrences are found (ScanProsite-consistent), case
  is folded via `RegexOptions.IgnoreCase` + uppercasing, and `Start`/`End` are **0-based
  inclusive** (:189-196). Anchors `^`/`$` inside the lookahead assert correctly at string
  boundaries.

### Cross-verification table recomputed vs code (all match)

| Input | Expected regex | Test |
|---|---|---|
| `R-G-D` | `RGD` | M1 ✅ |
| `A-x-G` | `A.G` | M2 ✅ |
| `x(3)-A` | `.{3}A` | M3 ✅ |
| `A-x(2,4)-G` | `A.{2,4}G` | M4 ✅ |
| `[ST]-x-[RK]` | `[ST].[RK]` | M5 ✅ |
| `N-{P}-[ST]-{P}` | `N[^P][ST][^P]` | M6 ✅ |
| `<M-x-K` | `^M.K` | M7 ✅ |
| `A-x-G>` | `A.G$` | M8 ✅ |
| `[RK](2)-x-[ST]` | `[RK]{2}.[ST]` | M9 ✅ |
| PS00028 full | `C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H` | M10 ✅ |
| PS00008 | `G[^EDRKHPFYW].{2}[STAGCN][^P]` | M11 ✅ |
| PS00018 | full EF-hand regex | M12 ✅ |
| `R-G-D.` / `R-G-D.A-B-C` | `RGD` | M14/M14b ✅ |
| `F-[IVFY]-G-[LM]-M-[G>].` | `F[IVFY]G[LM]M(?:G\|$)` | M22 ✅ |
| `F-[GSTV]-P-R-L-[G>].` | `F[GSTV]PRL(?:G\|$)` | M22b ✅ |

Match-position cross-checks: RGD @3 (M16); N-glyc {2,8} (M15); exclusion blocks (M17);
N/C-term anchors only at ends (S2/S3); `[G>]` matches via G branch, via C-term branch, and
rejects mid-sequence (M23/M23b/M23c); **Human Transferrin P02787 → exactly 2 PS00001 sites at
0-based 431 (`NKSD`) and 629 (`NVTD`)** = ScanProsite 1-based 432/630 (S1). All recomputed
against the running code via the test suite.

### Variant/delegate consistency

`FindMotifByProsite` delegates to `FindMotifByPattern`; `FindCommonMotifs`/`FindDomains` reuse
the same regex engine. `CommonMotifs` dictionary entries' precomputed regexes were spot-checked
against their PROSITE patterns (PS00001/4/5/6/7/8/9/16/17/18/28/29) — all consistent.

### Test quality audit

46 tests in `ProteinMotifFinder_PrositePattern_Tests.cs` (35 cited + sibling fixtures).
Assertions check exact regex strings, exact 0-based positions, exact matched substrings, and a
real UniProt sequence against ScanProsite-verified coordinates — not tautologies. All
Stage-A edge cases (empty pattern/sequence, trailing/mid period, exclusion block, both anchors,
`[G>]` three branches) are covered.

### Findings / defects

- **No defects.** The `{ABC}` exclusion is correctly negated (not inverted), ranges map
  correctly, anchors are enforced, positions are 0-based inclusive with no off-by-one.
- **Out-of-scope note (not a defect):** the extended ScanProsite query metacharacter `*`
  (Kleene star, e.g. `<{C}*>`) is not handled — it falls through and is silently dropped.
  This is a *query-form* extension, not part of standard PROSITE PA-line patterns, is outside
  this unit's spec and test scope, and no documented requirement exercises it.

**Stage B: PASS.**

---

## Verdict & follow-ups

- **Stage A: PASS** — every operator matches the ExPASy PROSITE User Manual §IV.E exactly;
  worked examples (N-glyc, P-loop, PKC, zinc finger) reproduce the manual/ScanProsite results.
- **Stage B: PASS** — implementation faithfully realises the grammar; all 46 PROSITE tests pass;
  full suite 4486 passed / 0 failed (baseline preserved).
- **State: CLEAN** — no defect found; no code changes required.

---

## Fix applied (2026-06-12)

**Hardening (follow-up to the out-of-scope note above).** The previously-noted silent drop of the
extended ScanProsite query metacharacter `*` (Kleene star, e.g. `<{C}*>`) is now **rejected, not
silently dropped**, mirroring the "reject, don't silently drop" policy already used by the Newick
parser (`PhylogeneticAnalyzer`, `FormatException` on unsupported/trailing input).

- **Change:** the final fall-through `else` branch in `ConvertPrositeToRegex`
  (ProteinMotifFinder.cs) previously did `i++` — silently skipping any unrecognised character.
  It now throws a **`FormatException`** naming the offending character, its position, and the
  pattern, and listing the supported grammar. `FormatException` matches the type already used for
  malformed-input rejection in the codebase (Newick parser).
- **Scope of the throw:** fires only on genuinely-unsupported tokens (`*`, and other stray
  metacharacters such as `?`/`+`) that reach the fall-through. All **supported** operators are
  handled by earlier branches and are **unchanged** — `[...]`, `{...}`, `x`, `x(n)`, `x(m,n)`,
  element repetition `(n)`/`(m,n)`, `<`, `>`, `-` separators, the `[G>]`/`[<…]` bracket special
  case, and the `.` terminator all produce **byte-for-byte identical** regex output and identical
  match results (verified by regression tests asserting exact regex strings and exact match
  positions).
- **Tests added** (`ProteinMotifFinder_PrositePattern_Tests.cs`): 4 rejection tests
  (`*` in `<{C}*>` and `A-x*-B` via both `ConvertPrositeToRegex` and `FindMotifByProsite`, plus
  stray `?`/`+`) asserting `FormatException` with the offending char in the message; 2 regression
  tests asserting the standard operators (N-glyc `N-{P}-[ST]-{P}`, P-loop `[AG]-x(4)-G-K-[ST]`,
  anchored `<M-x-K>`, `x(m,n)`, element repeat, trailing `.`, `[G>]`) still convert and match
  exactly. The four `*`/stray-metachar tests were **confirmed to FAIL first** against the
  pre-fix code (silently succeeded/mis-parsed); both regression tests passed before and after.
- **Suite:** full project **4495 passed / 0 failed** (baseline 4489 + 6 new); `~Prosite`,
  `~MotifSearch` (34) and `~Domain` (29) filters all green — no regression in the other consumers
  of the PROSITE engine.
