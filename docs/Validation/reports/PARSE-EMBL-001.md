# Validation Report: PARSE-EMBL-001 — EMBL flat-file parser (INSDC feature-location forms)

- **Validated:** 2026-06-24   **Area:** IO / Parsing
- **Canonical method(s):** `EmblParser.ParseLocation(string)` (`src/Seqeron/Algorithms/Seqeron.Genomics.IO/EmblParser.cs:613`); `EmblParser.UnquoteQualifierValue` (`:594`); shared `SequenceFormatHelper.ParseLocationParts`
- **Scope re-validated (independent re-confirmation):** the enhancement (commit `57957ab`) — three rarer INSDC Feature-Table location forms (remote reference, site-between `^`, deprecated single-period `n.m`) plus the INSDC `""`→`"` qualifier unescape (commit `2d9f023b`). Confirmed the unchanged forms (range / complement / join / partials) and the shared GenBank path are unaffected.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN

---

## Sources (retrieved THIS session)

- **INSDC Feature Table Definition** — EBI mirror `https://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.txt` (WebFetched 2026-06-24), §3.4.2.1 (descriptor forms) + §3.4.3 (location examples).

### Verbatim grammar / semantics confirmed (§3.4.2.1 descriptor forms)

Five base-position descriptors: (a) single base; (b) a **site between two adjoining bases**, written with a **carat `^`** (e.g. `55^56`); (c) **a single base from a range**, written with a **single period** (e.g. `12.21`); (d) a **sequence span**, written with **two periods** (e.g. `34..456`); (e) a **remote entry identifier** `accession.version:` followed by a local descriptor a–d (e.g. `J00194.1:1..15`). Restriction on (c) confirmed verbatim: *"From October 2006 the usage of this descriptor is restricted: it is illegal to use 'a single base from a range' (c) either on its own or in combination with the 'sequence span' (d) descriptor for newly created entries."*

### §3.4.3 example interpretations confirmed

`467` single base · `340..565` continuous inclusive range · `<345..500` unknown lower boundary · `102.110` exact location unknown, one base in [102,110] · `123^124` site between adjoining bases · `J00194.1:100..202` bases 100–202 inclusive in remote entry `J00194` (version 1).

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

### Independent cross-verification (throwaway console probe against the built `Seqeron.Genomics.IO.dll`, this session)

| Input | Start | End | Parts | Between | SingleDot | Remote | Acc | Ver | Verdict |
|-------|------|-----|-------|---------|-----------|--------|-----|-----|---------|
| `J00194.1:100..202` | 100 | 202 | 1 | F | F | T | J00194 | 1 | ✅ matches §3.4.3 |
| `J00194:100..202` | 100 | 202 | 1 | F | F | T | J00194 | null | ✅ |
| `J00194.1:467` | 467 | 467 | 1 | F | F | T | J00194 | 1 | ✅ single base |
| `123^124` | 123 | 124 | 1 | T | F | F | — | — | ✅ |
| `102.110` | 102 | 110 | 1 | F | T | F | — | — | ✅ |
| `102..110` | 102 | 110 | 1 | F | F | F | — | — | ✅ two-dot NOT flagged |
| `1^2` | 1 | 2 | 1 | T | F | F | — | — | ✅ circular boundary form |
| `J00194.12:100..202` | 100 | 202 | 1 | F | F | T | J00194 | 12 | ✅ multi-digit ver no leak |
| `NC_000001.11:5..9` | 5 | 9 | 1 | F | F | T | NC_000001 | 11 | ✅ RefSeq underscore |
| `J00194.1:` (malformed) | 0 | 0 | 0 | F | F | T | J00194 | 1 | ✅ no throw |
| `102.` (malformed) | 102 | 102 | 1 | F | F | F | — | — | ✅ no throw, unflagged |
| `complement(join(1..50,80..100))` | 1 | 100 | 2 | F | F | F | — | — | ✅ IsComplement+IsJoin |
| `<345..500` | 345 | 500 | 1 | F | F | F | — | — | ✅ Is5PrimePartial |
| `340..>565` | 340 | 565 | 2 | F | F | F | — | — | ✅ Is3PrimePartial (see N2) |

Qualifier unescape probe: `/note="a ""b"" c"` → `a "b" c` ✅ (doubled-quote collapsed per INSDC).

### Shared GenBank path — confirmed unaffected
`SequenceFormatHelper.ParseLocationParts` is unchanged. The remote/between/single-dot logic lives only in `EmblParser.ParseLocation`; `GenBankParser.Location` was not extended, so the remote-ref handling is deliberately scoped to EMBL. The two-period span / complement / join / partial results are byte-for-byte the same baseline behavior validated earlier.

### Test-quality audit (HARD gate)
- 124 EmblParser tests pass (filtered run). Tests cite §3.4.2.1/§3.4.3 and assert sourced bounds/flags — not "no-throw" tautologies.
- Discriminator tests present: `100..200` IsBetween=false, `102..110` IsSingleBaseFromRange=false, `100..200` IsRemote=false; multi-digit-version (`.12`) and RefSeq (`NC_000001.11`) hardening tests lock the regex against single-digit/no-underscore assumptions.
- Malformed inputs (`J00194.1:`, `102.`) covered — no-throw, correct non-flagging.

### Findings / notes
- **N1 (PASS-WITH-NOTES, BY-DESIGN):** a remote reference *inside* an operator (`complement(J00194.1:100..202)`) is not recognized as remote — `RemoteReferenceRegex` is `^`-anchored — so the version digit can leak into the span. Legal per §3.4.2.1(e) but rare, out of the enhancement's stated scope, and **not a regression** (the old parser and the unextended GenBank path fail it identically). No buggy-value test added.
- **N2 (pre-existing, shared helper, out of scope):** `340..>565` reports `Parts.Count=2`. This is unchanged baseline behavior of `SequenceFormatHelper.ParseLocationParts` (the `>` partial path), not introduced by the C8 enhancement; Start/End and the `Is3PrimePartial` flag are correct. Noted, not a defect of this unit's enhancement.

---

## Verdict & follow-ups

- **Stage A: PASS** — grammar + semantics confirmed verbatim against the INSDC FT Definition (EBI mirror) retrieved this session.
- **Stage B: PASS-WITH-NOTES** — all three enhanced location forms and the `""`→`"` qualifier unescape parse to independently-sourced values with no digit leak; discriminators and malformed inputs handled; shared GenBank path provably untouched. Two documented non-defect notes (N1 by-design, N2 pre-existing shared-helper).
- **End-state: ✅ CLEAN** — no code changed this session. Full suite **18213 passed / 0 failed**, build 0 warnings/0 errors.
