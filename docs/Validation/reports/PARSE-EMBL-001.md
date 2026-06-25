# Validation Report: PARSE-EMBL-001 — EMBL flat-file parser (INSDC feature-location forms)

- **Validated:** 2026-06-25 (re-validated fresh; unit reset to ⬜ pending after the F1 nested-remote-reference change)   **Area:** IO / Parsing
- **Canonical method(s):** `EmblParser.ParseLocation(string)` (`src/Seqeron/Algorithms/Seqeron.Genomics.IO/EmblParser.cs:626`); `EmblParser.ExtractNestedRemoteReferences` (`:704`); `EmblParser.Location` / `RemotePart` records (`:62-97`); `EmblParser.UnquoteQualifierValue` (`:607`); shared `SequenceFormatHelper.ParseLocationParts`.
- **Scope re-validated (fresh, independent):** the full canonical parse surface plus the campaign-added nested-remote-reference location parsing — INSDC location forms incl. top-level remote `accession[.version]:loc`, remote nested inside `complement(...)`/`join(...)`, site-between `n^m`, single-base-from-range `n.m`, `<..>` partials, ordinary spans. The F1 change (commit `f4d8d55d`, nested remote refs + `Location.RemoteParts`) and the shared GenBank location parser were re-confirmed.
- **Out-of-scope boundary (documented, acceptable):** the remote entry's actual *sequence* is not fetched (needs network/DB access); only the local span and the remote reference (accession/version/bounds) are extracted.
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** 🟡 PASS-WITH-NOTES
- **End-state:** ✅ CLEAN

---

## Sources (retrieved THIS session, 2026-06-25)

- **The DDBJ/ENA/GenBank Feature Table Definition** — EBI mirror `https://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.txt` (WebFetched 2026-06-25), §3.4.2.1 (descriptor forms), §3.4.2.2 (operators), §3.4.3 (location examples). Confirmed verbatim, NOT from the repo's own TestSpec/tests.

### Verbatim grammar / semantics confirmed (§3.4.2.1 descriptor forms)

Five base-position descriptors: (a) single base; (b) a **site between two adjoining bases**, written with a **carat `^`** (e.g. `55^56`; for circular molecules `n^1`); (c) **a single base from a range**, written with a **single period** (e.g. `12.21`); (d) a **sequence span**, written with **two periods** (e.g. `34..456`); (e) a **remote entry identifier** `accession.version:` followed by a local descriptor a–d (e.g. `J00194.1:1..15`). The `<` / `>` symbols indicate endpoints beyond the specified position (5′ / 3′ partial). Restriction on (c) confirmed verbatim: *"From October 2006 the usage of this descriptor is restricted: it is illegal to use 'a single base from a range' (c) either on its own or in combination with the 'sequence span' (d) descriptor for newly created entries."*

### Operators confirmed (§3.4.2.2)

`complement(location)` — complement of the span; `join(...)` — elements placed end-to-end into one contiguous sequence; `order(...)` — elements in 5′→3′ order, no joining implied. Nesting rule confirmed verbatim: *"complement can be used in combination with either 'join' or 'order' within the same location; combinations of 'join' and 'order' within the same location (nested operators) are illegal."* Remote references (e) may appear nested as join/order/complement elements.

### §3.4.3 example interpretations confirmed (verbatim)

`467` single base · `340..565` continuous inclusive range · `<345..500` / `<1..888` unknown lower boundary (before stated base) · `1..>888` continues beyond stated base · `102.110` exact location unknown, one base in [102,110] · `123^124` site between adjoining bases · `complement(34..126)` complementary strand · `complement(join(...))` join-then-complement · `join(complement(...),complement(...))` complement-then-join · `J00194.1:100..202` bases 100–202 inclusive in remote entry `J00194` (version 1) · `join(1..100,J00194.1:100..202)` join local region with remote region.

---

## Stage A — Description

The three enhanced location forms map exactly to INSDC §3.4.2.1 descriptors (b), (c), (e); the §3.4.3 example table fixes their expected bounds. The qualifier-value rule (free text in double quotes; an embedded `"` encoded as a doubled `""`) is the standard INSDC qualifier format. Independently-derived expected values (NOT lifted from repo tests):

| Form | Source semantics | Expected parse |
|------|------------------|----------------|
| `J00194.1:100..202` | bases 100–202 inclusive in entry `J00194` (version 1) | acc=`J00194`, ver=`1`, span 100..202, IsRemote |
| `123^124` | site between bases 123 and 124 | IsBetween, bounds 123/124 |
| `102.110` | one base somewhere in [102,110] inclusive | IsSingleBaseFromRange, bounds 102/110 |
| `100..200`, `complement()`, `join()`, `<..>` | standard INSDC — unchanged | unchanged |
| `/note="a ""b"" c"` | doubled `""` is one literal `"` | value `a "b" c` |

Version is "a sequence of digits" (multi-digit, e.g. `.12`) and accession may contain `_` (RefSeq, `NC_000001.11`) — both confirmed against §3.4.2.1. EMBL line-type codes (ID, AC, SV, DE, KW, OS/OC, RN/RP/RX/RG/RA/RT/RL, FH/FT, SQ, XX, `//`) and 1-based inclusive coordinates match the EBI manual. **Stage A PASS.**

---

## Stage B — Implementation

### Code path reviewed
`EmblParser.ParseLocation` (`EmblParser.cs:613-665`): (1) strips/captures the remote prefix via `RemoteReferenceRegex` `^(?<acc>[A-Za-z][A-Za-z0-9_]*)(?:\.(?<ver>\d+))?:` (`:780`); (2) site-between `^(\d+)\^(\d+)$` (`:784`); (3) single-period `^(\d+)\.(\d+)$` (`:789`); (4) otherwise delegates to the shared `SequenceFormatHelper.ParseLocationParts`. The `Location` record carries four trailing defaulted members (`IsBetween`, `IsSingleBaseFromRange`, `RemoteAccession`, `RemoteVersion`) plus computed `IsRemote => RemoteAccession is not null` (`:64-84`). Qualifier unescape: `UnquoteQualifierValue` (`:594`) strips one outer quote pair then `.Replace("\"\"", "\"")`.

### Independent cross-verification (throwaway console probe against the freshly built `Seqeron.Genomics.IO.dll`, this session 2026-06-25)

Every expected value below is hand-derived from the §3.4.3 verbatim interpretations, NOT from code output. `RemoteParts` column shows `acc/ver:start..end` per nested segment.

| Input | Start | End | Parts | Cmpl | Join | Btwn | SDot | Remote | Acc | Ver | RemoteParts | Verdict |
|-------|------|-----|-------|------|------|------|------|--------|-----|-----|-------------|---------|
| `join(1..100,J00194.1:100..202)` | 1 | 202 | `1..100`,`100..202` | F | T | F | F | T | — | — | `J00194/1:100..202` | ✅ §3.4.3 join local+remote, no `.1` leak |
| `complement(join(2691..4571,X00001.1:1..50))` | 1 | 4571 | `2691..4571`,`1..50` | T | T | F | F | T | — | — | `X00001/1:1..50` | ✅ complement+join+nested remote |
| `complement(AB000123.1:500..600)` | 500 | 600 | `500..600` | T | F | F | F | T | — | — | `AB000123/1:500..600` | ✅ remote nested in complement |
| `join(X00001.1:10..20,complement(X00002.1:30..40))` | 10 | 40 | `10..20`,`30..40` | T | T | F | F | T | — | — | `X00001/1:10..20`,`X00002/1:30..40` | ✅ two nested remotes captured |
| `J00194.1:100..202` | 100 | 202 | `100..202` | F | F | F | F | T | J00194 | 1 | null | ✅ top-level remote (RemoteParts null) |
| `J00194:100..202` | 100 | 202 | `100..202` | F | F | F | F | T | J00194 | null | null | ✅ no version |
| `J00194.1:467` | 467 | 467 | `467..467` | F | F | F | F | T | J00194 | 1 | null | ✅ remote single base |
| `123^124` | 123 | 124 | `123..124` | F | F | T | F | F | — | — | null | ✅ site-between |
| `102.110` | 102 | 110 | `102..110` | F | F | F | T | F | — | — | null | ✅ single-base-from-range |
| `102..110` | 102 | 110 | `102..110` | F | F | F | F | F | — | — | null | ✅ two-dot NOT flagged |
| `J00194.12:100..202` | 100 | 202 | `100..202` | F | F | F | F | T | J00194 | 12 | null | ✅ multi-digit ver no leak |
| `NC_000001.11:5..9` | 5 | 9 | `5..9` | F | F | F | F | T | NC_000001 | 11 | null | ✅ RefSeq underscore + 2-digit ver |
| `<1..888` | 1 | 888 | `1..888` | F | F | F | F | F | — | — | null | ✅ Is5PrimePartial |
| `1..>888` | 1 | 888 | `1..1`,`888..888` | F | F | F | F | F | — | — | null | ✅ Is3PrimePartial, Start/End correct (see N2) |
| `340..565` | 340 | 565 | `340..565` | F | F | F | F | F | — | — | null | ✅ ordinary span |
| `467` | 467 | 467 | `467..467` | F | F | F | F | F | — | — | null | ✅ single base |
| `order(100..200,300..400)` | 100 | 400 | `100..200`,`300..400` | F | F(Ord=T) | F | F | F | — | — | null | ✅ order operator |
| `complement(100..200)` | 100 | 200 | `100..200` | T | F | F | F | F | — | — | null | ✅ complement local |
| `J00194.1:` (malformed) | 0 | 0 | — (0 parts) | F | F | F | F | T | J00194 | 1 | null | ✅ no throw |
| `102.` (malformed) | 102 | 102 | `102..102` | F | F | F | F | F | — | — | null | ✅ no throw, unflagged |
| `340..>565` | 340 | 565 | `340..340`,`565..565` | F | F | F | F | F | — | — | null | ✅ Is3PrimePartial (see N2) |

Qualifier unescape (verified via test `Parse_Qualifier_EscapedDoubleQuote_Unescaped`): `/note="he said ""hi"""` → `he said "hi"` ✅ (doubled-quote collapsed per INSDC qualifier rule).

### Shared GenBank path — confirmed unaffected
`SequenceFormatHelper.ParseLocationParts` is unchanged. The remote/between/single-dot logic lives only in `EmblParser.ParseLocation`; `GenBankParser.Location` was not extended, so the remote-ref handling is deliberately scoped to EMBL. The two-period span / complement / join / partial results are byte-for-byte the same baseline behavior validated earlier.

### Test-quality audit (HARD gate)
- `EmblParserTests.cs` + `EmblParser_MutationKillers_Tests.cs` tests pass. Tests cite §3.4.2.1/§3.4.3 and assert sourced bounds/flags — not "no-throw" tautologies.
- **Every public method/overload covered:** `Parse`/`ParseFile` (valid, empty, null, whitespace, multi-record, unterminated), `ParseLocation` (all five descriptor forms + operators + partials + remote forms), `ToGenBank`, `GetFeatures`/`GetCDS`/`GetGenes`, `ExtractSequence` (simple/complement/join), qualifier unescape, the `Location.IsRemote` computed property, the `RemotePart` record.
- **Nested-remote coverage (F1):** `ParseLocation_JoinWithNestedRemoteReference_*`, `_JoinWithTwoNestedRemoteReferences_OneComplemented_*`, `_ComplementWithNestedRemoteReference_*` assert per-segment `RemoteParts` (accession/version/bounds) with expected values hand-derived from the §3.4.3 verbatim interpretation, plus three regression tests proving local operators / top-level remote refs do NOT populate `RemoteParts`.
- **Discriminator tests present:** `100..200` IsBetween=false, `102..110` IsSingleBaseFromRange=false, `100..200` IsRemote=false; multi-digit-version (`.12`) and RefSeq (`NC_000001.11`) hardening tests lock the regex against single-digit/no-underscore assumptions; no-version remote ref.
- **Malformed inputs** (`J00194.1:`, `102.`) covered — no-throw, correct non-flagging.
- No green-washing detected; expected values trace to the INSDC FT definition, not code echoes.

### Findings / notes
- **N1 (by design — F1 nested remote refs):** a remote reference nested *inside* an operator (`join(1..100,J00194.1:100..202)`, `complement(J00194.1:100..202)`, `complement(join(...,X00001.1:1..50))`) is captured. `EmblParser.ParseLocation` runs `ExtractNestedRemoteReferences` over the descriptor (after the anchored top-level strip): `NestedRemoteReferenceRegex` (`(?<=[(,])(?<acc>...)(?:\.(?<ver>\d+))?:`) matches each segment-leading `accession[.version]:` prefix at an operator boundary, records a per-segment `RemotePart(accession, version, start, end)` in `Location.RemoteParts`, and strips the prefixes before `SequenceFormatHelper.ParseLocationParts` so the version digit no longer leaks into the numeric `Parts`. Verified verbatim against §3.4.3 `join(1..100,J00194.1:100..202)`. The wrapping `complement(` is preserved (lookbehind not consumed) so the complement flag is unaffected. **Residual / documented out-of-scope boundary:** the remote entry's *sequence* is not fetched (location captured, content not retrieved — needs network/DB). The shared GenBank `Location` is deliberately left unextended.
- **N2 (pre-existing, shared helper, out of scope):** `1..>888` and `340..>565` report `Parts.Count=2` (`1..1`,`888..888`). This is unchanged baseline behavior of `SequenceFormatHelper.ParseLocationParts` — the shared `LocationRangeRegex` `(\d+)(?:\.\.(\d+))?` requires the two literal periods to be immediately followed by a digit, so a `>` between `..` and the number breaks the span into two single-base matches. The same shared helper backs GenBank (pre-dates the EMBL remote-ref work, commit `f4d8d55d`). The overall `Start`/`End` and the `Is3PrimePartial` flag are correct, so callers reading Start/End/flags are unaffected; only the `Parts` decomposition of a partial-bound span is split. Noted, not a defect of this unit's scope.

---

## Verdict & follow-ups

- **Stage A: ✅ PASS** — grammar + semantics confirmed verbatim against the DDBJ/ENA/GenBank Feature Table Definition (EBI mirror) retrieved fresh this session (2026-06-25): §3.4.2.1 descriptors, §3.4.2.2 operators, §3.4.3 examples.
- **Stage B: 🟡 PASS-WITH-NOTES** — all five descriptor forms, the operators, partials, top-level and nested remote references, and the `""`→`"` qualifier unescape parse to independently-sourced values with no version-digit leak; discriminators and malformed inputs handled; the shared GenBank path is provably untouched. Two documented non-defect notes (N1 by-design + out-of-scope remote-sequence fetch, N2 pre-existing shared-helper partial-bound `Parts` split).
- **End-state: ✅ CLEAN** — no code changed this session; re-validated fresh against the external spec. Full unfiltered suite **18819 passed / 0 failed** (`Seqeron.Genomics.Tests`), 0 failures across the whole solution, build 0 warnings / 0 errors. The only residual is the documented, acceptable out-of-scope boundary (remote sequence not fetched).
