# Validation Report: PARSE-EMBL-001 — EMBL flat-file parser (INSDC feature-location forms)

- **Validated:** 2026-06-17   **Area:** IO / Parsing
- **Canonical method(s):** `EmblParser.ParseLocation(string)` (`src/Seqeron/Algorithms/Seqeron.Genomics.IO/EmblParser.cs:606`); shared `SequenceFormatHelper.ParseLocationParts` (`src/Seqeron/Algorithms/Seqeron.Genomics.IO/SequenceFormatHelper.cs:25`)
- **Scope re-validated:** the C8 enhancement (commit `57957ab`) — three rarer INSDC Feature-Table location forms: remote references, the site-between operator `^`, the deprecated single-period range. Plus confirmation that the unchanged forms and the shared GenBank path are unaffected.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES
- **End-state:** ✅ CLEAN

---

## Sources (retrieved THIS session)

- **INSDC Feature Table Definition** — EBI mirror `https://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.txt` (WebFetched 2026-06-17).
- **INSDC** — `https://www.insdc.org/submitting-standards/feature-table/` (WebFetched 2026-06-17, corroborating section 3.4.2.1 / 3.4.3).

### Verbatim grammar / semantics (INSDC FT Definition §3.4.2.1 descriptor forms)

> (a) "a single base number"
> (b) "a site between two indicated adjoining bases"
> (c) "a single base chosen from within a specified range of bases"
> (d) "the base numbers delimiting a sequence span"
> (e) "a remote entry identifier followed by a local location descriptor (i.e., a-d)"

Notation rules (verbatim): sites between bases use a **"carat (^)"**; a single-base-from-range uses **"a single period"** (e.g. `12.21`); sequence spans use **"two periods"** (e.g. `34..456`). The single-period form is restricted: *"From October 2006 the usage of this descriptor is restricted: it is illegal to use 'a single base from a range' (c) either on its own or in combination."* The site-between permitted formats: *"n^n+1 (for example 55^56), or, for circular molecules, n^1."*

### Verbatim §3.4.3 location examples

> `467` — "Points to a single base in the presented sequence"
> `340..565` — "Points to a continuous range of bases bounded by and including the starting and ending bases"
> `102.110` — "Indicates that the exact location is unknown but that it is one of the bases between bases 102 and 110, inclusive"
> `123^124` — "Points to a site between bases 123 and 124"
> `J00194.1:100..202` — "Points to bases 100 to 202, inclusive, in the entry (in this database) with primary accession number 'J00194'"

---

## Stage A — Description

The three enhanced forms map exactly to INSDC §3.4.2.1 descriptors (b), (c), (e), and the example
table §3.4.3 fixes their expected bounds. Independently confirmed expected values (NOT lifted from
the repo tests):

| Form | Source semantics | Expected parse |
|------|------------------|----------------|
| `J00194.1:100..202` | bases 100–202 inclusive in entry `J00194` (version 1) | acc=`J00194`, ver=`1`, span 100..202, IsRemote |
| `123^124` | site between bases 123 and 124 | IsBetween, bounds 123/124 |
| `102.110` | one base somewhere in [102,110] inclusive | IsSingleBaseFromRange, bounds 102/110 |
| `100..200`, `complement()`, `join()`, `<..>`, `""`→`"` | standard INSDC — unchanged | unchanged |

The version is "a sequence of digits" (not a single digit) and the accession may contain an
underscore (RefSeq, e.g. `NC_000001.11`) — both verified against §3.4.2.1. **Stage A PASS.**

---

## Stage B — Implementation

### Code path reviewed
`EmblParser.ParseLocation` (`EmblParser.cs:606-658`): (1) strips/captures the remote prefix via
`RemoteReferenceRegex` `^(?<acc>[A-Za-z][A-Za-z0-9_]*)(?:\.(?<ver>\d+))?:` (line 779); (2) detects
site-between `^(\d+)\^(\d+)$` (line 783); (3) detects single-dot `^(\d+)\.(\d+)$` (line 788);
(4) otherwise delegates to the shared `SequenceFormatHelper.ParseLocationParts`. The `Location` record
gains four trailing members (`IsBetween`, `IsSingleBaseFromRange`, `RemoteAccession`, `RemoteVersion`,
all defaulted) plus a computed `IsRemote => RemoteAccession is not null` (lines 64-84).

### Independent cross-verification (throwaway console probe against the built `Seqeron.Genomics.IO.dll`)

| Input | Start | End | Parts | Between | SingleDot | Remote | Acc | Ver | Verdict |
|-------|------|-----|-------|---------|-----------|--------|-----|-----|---------|
| `J00194.1:100..202` | 100 | 202 | 1 | F | F | T | J00194 | 1 | ✅ matches §3.4.3 |
| `J00194:100..202` | 100 | 202 | 1 | F | F | T | J00194 | null | ✅ |
| `J00194.1:467` | 467 | 467 | — | F | F | T | J00194 | 1 | ✅ single base |
| `123^124` | 123 | 124 | 1 | T | F | F | — | — | ✅ |
| `102.110` | 102 | 110 | 1 | F | T | F | — | — | ✅ |
| `102..110` | 102 | 110 | 1 | F | F | F | — | — | ✅ two-dot NOT flagged |
| `1^2` | 1 | 2 | 1 | T | F | F | — | — | ✅ |
| `J00194.12:100..202` | 100 | 202 | 1 | F | F | T | J00194 | 12 | ✅ multi-digit ver no leak |
| `NC_000001.11:5..9` | 5 | 9 | 1 | F | F | T | NC_000001 | 11 | ✅ RefSeq underscore |
| `J00194.1:` (malformed) | 0 | 0 | 0 | F | F | T | J00194 | 1 | ✅ no throw |
| `102.` (malformed) | 102 | 102 | 1 | F | F | F | — | — | ✅ no throw, unflagged |

### OLD-parser mis-handling demonstrated (commit `57957ab^`)

The pre-enhancement `ParseLocation` had no remote/between/single-dot logic; it fed the whole string
to the shared regex `(\d+)(?:\.\.(\d+))?`:
- `J00194.1:100..202` → matched `00194`(=194), `1`, `100..202` → 3 parts, **Start leaked to 1** (since `overallStart = min`). Now: 1 part, Start=100, version captured.
- `123^124` → two single-base parts (123),(124), **no flag** to distinguish a span. Now: `IsBetween=true`.
- `102.110` → `..` required, so single `.` split into (102),(110) two parts, **no flag**. Now: `IsSingleBaseFromRange=true`, 1 part (102,110).

### Shared GenBank path — confirmed unaffected
`SequenceFormatHelper.ParseLocationParts` is byte-for-byte unchanged; `GenBankParser.ParseLocation`
still produces `100..200`→1 part, `complement(1..50)`→span, `join(...)`→2 parts. GenBank's `Location`
record was *not* extended, so `GenBankParser.ParseLocation("J00194.1:100..202")` still yields the old
leaky Start=1 / 3-parts result — consistent with the enhancement deliberately scoping remote-ref
handling to `EmblParser` only.

### Test-quality audit (HARD gate)
- **Sourced exact assertions** — every new test cites §3.4.2.1/§3.4.3 and asserts the sourced bounds/flags; no "no-throw"-only or tautological assertions.
- **Discriminator (false-flag) tests** — `100..200` IsBetween=false, `102..110` IsSingleBaseFromRange=false, `100..200` IsRemote=false: a wrong flag fails a test.
- **Malformed variants covered** — `J00194.1:` and `102.` (no-throw, correct non-flagging).
- **Mutation-verified this session** — forcing `IsBetween:false` killed 2 tests; forcing `remoteVersion=null` killed 4 tests. The tests genuinely lock behavior.
- **Validator hardening added (test-only, 0 code change)** — `ParseLocation_RemoteReference_MultiDigitVersion_DoesNotLeakIntoSpan` (`.12`) and `ParseLocation_RemoteReference_RefSeqUnderscoreAccession_CapturesAccessionVersion` (`NC_000001.11`) lock the regex against a single-digit-version / no-underscore assumption.

### Findings / notes
- **N1 (PASS-WITH-NOTES, BY-DESIGN limitation):** A remote reference *inside* an operator —
  `complement(J00194.1:100..202)` / `join(J00194.1:1..9, ...)` — is **not** recognized as remote
  (`RemoteReferenceRegex` is anchored `^`), so its accession/version are not captured and the version
  digit leaks into the span (Start→1). This is **legal** per §3.4.2.1(e) but **rare**, was **out of
  the enhancement's stated scope**, and is **not a regression** (the old parser failed it identically,
  and the unextended GenBank path still does). Logged as BY-DESIGN; no buggy-value test added.

---

## Verdict & follow-ups

- **Stage A: PASS.** Grammar + semantics confirmed verbatim against the INSDC FT Definition (EBI mirror) retrieved this session.
- **Stage B: PASS-WITH-NOTES.** All three enhanced forms parse to the independently-sourced values with no digit leak; discriminators and malformed inputs handled; shared GenBank path provably untouched. Tests are strict, sourced, and mutation-killed. One documented BY-DESIGN limitation (operator-wrapped remote refs, N1).
- **End-state: ✅ CLEAN.** 2 hardening tests added (test-only); full suite **6771 passed / 0 failed**, build 0 errors, clean tree.
