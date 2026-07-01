# Validation Report: PARSE-EMBL-001 — EMBL flat-file parser + remote-aware location-sequence assembly (INSDC feature-location forms)

- **Validated:** 2026-06-26 (re-validated fresh; unit reset to ⬜ pending after the limitation-fix `d40318fb` ADDED remote-aware location-sequence assembly)   **Area:** IO / Parsing
- **Canonical method(s):**
  - **New/changed (primary focus):** `FeatureLocationHelper.ResolveLocationSequence(string rawLocation, string localSequence, RemoteSequenceResolver?)` (`src/Seqeron/Algorithms/Seqeron.Genomics.IO/FeatureLocationHelper.cs:157`); the `EmblParser.Location` overload (`:194`); `EmblParser.ResolveLocationSequence(EmblRecord, Location, RemoteSequenceResolver?)` (`EmblParser.cs:866`); the `delegate string? RemoteSequenceResolver(string accession, int? version)` (`FeatureLocationHelper.cs:92`); helpers `EnumerateSegments`/`ResolveSegment`/`SliceSpan` (`:207/236/281`).
  - **Pre-existing (confirmed still holds):** `EmblParser.ParseLocation` (`EmblParser.cs:626`), `ExtractNestedRemoteReferences`, `Location.RemoteParts`, `UnquoteQualifierValue`, shared `SequenceFormatHelper.ParseLocationParts`.
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **End-state:** ✅ CLEAN

---

## Sources (retrieved THIS session, 2026-06-26)

- **The DDBJ/ENA/GenBank Feature Table Definition (INSDC, v11.x)** — `https://www.insdc.org/submitting-standards/feature-table/` (WebFetched 2026-06-26), §3.4 Location / §3.5 Operators.
- **DDBJ mirror** — `https://www.ddbj.nig.ac.jp/ddbj/feature-table-e.html` (WebFetched 2026-06-26), cross-checked for the §3.5 complement/join equivalence worked example.

Both retrieved fresh this session; NOT taken from the repo's own TestSpec/Evidence/tests.

### Verbatim semantics confirmed (both mirrors agree)

- **Single base** `467` — points to a single base. **Span `n..m`** — "Points to a continuous range of bases bounded by **and including** the starting and ending bases" → 1-based **inclusive**.
- **Site between bases** `n^m` (carat); **single base from a range** `n.m` (single period, restricted since Oct 2006).
- **Remote entry** — "the accession-number and sequence version of the remote entry, followed by a colon ':', followed by a location descriptor which applies to that entry's sequence" (e.g. `J00194.1:100..202`).
- **`<` / `>` partials** — "indicate that an end point is beyond the specified base number"; the stated number is the only available coordinate.
- **`complement(location)`** — "read the complement of the presented strand in its 5'-to-3' direction" → reverse-complement.
- **`join(...)`** — "elements should be joined (placed end-to-end) to form one contiguous sequence" (listed order). **`order(...)`** — same 5'→3' order, joining not implied (identical concatenation for sequence extraction).
- **Equivalence (verbatim worked example, §3.5):** `complement(join(2691..4571,4918..5163))` **==** `join(complement(4918..5163),complement(2691..4571))` — both place the feature on the complementary strand; the inner segment order is **reversed** when the complement is distributed.

---

## Stage A — Description

The new assembly contract maps exactly to INSDC §3.4/§3.5: 1-based inclusive slicing of every span (local or remote), listed-order concatenation under `join`/`order`, reverse-complement under `complement(...)` with the documented order-reversal equivalence, remote spans fetched by accession+version and sliced identically, `<`/`>` markers using the stated number verbatim. The offline-first design (no network I/O; a caller-supplied `RemoteSequenceResolver` performs any fetch, the library does the assembly) is an architectural choice, not a spec divergence — the INSDC grammar says nothing about *how* the remote sequence is obtained, only how the location is interpreted once obtained. **Stage A PASS.**

---

## Stage B — Implementation

### Code path reviewed
`FeatureLocationHelper.ResolveLocationSequence` (`FeatureLocationHelper.cs:157-188`): (1) trims; empty/null → `""`. (2) `OuterComplementRegex` `^complement\((?<inner>.*)\)$` detects a whole-location complement and strips it. (3) `EnumerateSegments` removes a `join|order(...)` wrapper and depth-splits the inner on top-level commas (parenthesis-aware, so per-segment `complement(...)` stays intact). (4) each segment → `ResolveSegment`: a per-element `complement(...)` is RC'd individually; a `SegmentRemoteRegex` `^[A-Za-z][A-Za-z0-9_]*(\.\d+)?:` match routes to `remoteResolver?.Invoke(acc, ver)`, else to `localSequence`; `SliceSpan` slices 1-based inclusive via `SequenceFormatHelper.LocationRangeRegex` (`(\d+)(?:\.\.(\d+))?`, single number → single base), clamping to bounds. (5) whole assembly RC'd if outer-complement. Reverse-complement is `DnaSequence.GetReverseComplementString` (reverse + per-base complement, IUPAC-aware). The whole-result RC realises the §3.5 equivalence (reversing the concatenation reverses segment order AND complements each base in one operation) without explicitly reversing the segment list. **No `System.Net`/`HttpClient`/`Socket` anywhere in `Seqeron.Genomics.IO` (grep: NONE).**

### Independent cross-verification (throwaway console probe vs the freshly built `Seqeron.Genomics.IO.dll`, this session 2026-06-26)

Probe used **fixtures distinct from the repo tests** to avoid echoing them — Local `ATGCATGCATGC`; resolver REM=`AAAAACCCCC`, Z=`GATCG`. Every expected string is hand-derived from the §3.4/§3.5 verbatim semantics, NOT from code output. **All 17 reproduced exactly.**

| Case | Input | Hand-derivation | Expected | Code |
|------|-------|-----------------|----------|------|
| local-only join | `join(1..4,9..12)` | 1..4=ATGC, 9..12=ATGC; resolver NOT invoked | `ATGCATGC`, invoked=false | ✅ |
| remote 1-based slice | `REM.1:3..8` | REM bases 3-8 = A A A C C C | `AAACCC` | ✅ |
| mixed join order | `join(1..4,REM.1:6..10)` | ATGC + REM 6-10 (CCCCC) | `ATGCCCCCC` | ✅ |
| complement(remote) | `complement(Z.1:1..5)` | RC(GATCG)=cmpl CTAGC, rev CGATC | `CGATC` | ✅ |
| complement(join local+remote) | `complement(join(1..4,REM.1:1..3))` | inner ATGC+AAA=ATGCAAA; RC = cmpl TACGTTT rev TTTGCAT | `TTTGCAT` | ✅ |
| equivalence (other form) | `join(complement(REM.1:1..3),complement(1..4))` | RC(AAA)=TTT + RC(ATGC)=GCAT | `TTTGCAT` (== prior) | ✅ |
| null resolver | `join(1..4,REM.1:1..5)` (resolver=null) | remote→"", local 1..4=ATGC | `ATGC` | ✅ |
| resolver returns null | `join(1..4,UNKNOWN.1:1..5)` | unknown acc→"", local kept | `ATGC` | ✅ |
| `<` partial | `<1..4` | marker ignored, bases 1..4 | `ATGC` | ✅ |
| `>` partial | `join(1..4,>9..12)` | ATGC + bases 9..12 (ATGC) | `ATGCATGC` | ✅ |
| single base | `3` | base 3 = G | `G` | ✅ |
| order operator | `order(1..4,9..12)` | same concat as join | `ATGCATGC` | ✅ |
| version pass-through | `REM.7:1..3` | resolver sees acc=REM, ver=7 | acc=REM, ver=7 | ✅ |
| no-version → null | `REM:1..3` | resolver sees ver=null | ver=null | ✅ |
| empty / null location | `""`, `null` | contract parity with ExtractSequence | `""`, `""` | ✅ |
| Location overload | `EmblParser.ParseLocation("join(1..4,REM.1:6..10)")` | RawLocation reused | `ATGCCCCCC` | ✅ |
| spec equivalence (large) | `complement(join(5..10,13..15))` vs `join(complement(13..15),complement(5..10))` on `A⁴C⁶A²G³` | inner CCCCCC+GGG=CCCCCCGGG; RC=CCCGGGGGG | both `CCCGGGGGG` | ✅ |

### Pre-existing surface — confirmed still holds
The `d40318fb` diff to `EmblParser.cs` is **purely additive** (one new `ResolveLocationSequence` overload delegating to `FeatureLocationHelper`). `ParseLocation`, `ExtractNestedRemoteReferences`, `Location.RemoteParts`, `UnquoteQualifierValue`, and `SequenceFormatHelper.ParseLocationParts` are byte-for-byte unchanged. `RawLocation` is preserved in every `Location` constructor path (used by the overload). The full `EmblParser*` test suite (bare span, `n^m`, `n.m`, `<..>`, complement/join/order, remote `accession.version:loc`, multi-digit version, RefSeq underscore, no-version, `""`→`"`) is part of the green suite.

### Test-quality audit (HARD gate)
- `FeatureLocationHelper_ResolveLocationSequence_Tests.cs` — 17 RESV cases assert **exact hand-derived strings** (`ACGGTA`, `GCCCCCAAAA`, `CCAATACGT`, `GTTT`, …) with base-by-base derivation in the comments, not code echoes or "no-throw" tautologies. Deterministic, no shared RNG.
- Coverage: every operator (`join`, `order`, `complement`, `complement(join(...))`), remote slice, complement-of-remote, the INSDC equivalence (asserted equal AND equal to a hand value, twice — small + large), null & null-returning resolver, `<`/`>` partials, single base, empty/null location, both raw-string and `EmblParser.Location` overloads, version present/absent pass-through.
- My independent probe (different fixtures) reproduced 17 hand-derived strings → tests are not masking a fixture-specific quirk. No green-washing.

### Findings / notes
- **N1 (by design — offline-first):** the library performs NO network I/O. A remote span's sequence is obtained via the caller-supplied `RemoteSequenceResolver`; when absent/returning null the remote element contributes the empty string and local spans are still spliced in place (mirrors the clamping, non-throwing local-extraction contract). This is the documented design (LIMITATIONS §3 reworded to by-design), not a defect.
- **N2 (out of scope for the assembly path):** `SliceSpan` parses via `LocationRangeRegex`, so a degenerate `n^m`/`n.m` descriptor passed to `ResolveLocationSequence` would be read as a single base / two-period span rather than a between-site/uncertain base. These descriptors denote a *site* or an *uncertain single base*, not an extractable contiguous span, and are out of scope for the resolve contract (which targets spans / `join` / `order` / `complement` / remote refs). The structured `ParseLocation` surface flags them correctly (`IsBetween`/`IsSingleBaseFromRange`). Not a defect.

---

## Verdict & follow-ups

- **Stage A: ✅ PASS** — INSDC §3.4/§3.5 location semantics + the complement/join order-reversal equivalence confirmed verbatim against insdc.org and the DDBJ mirror, both retrieved fresh this session (2026-06-26).
- **Stage B: ✅ PASS** — `ResolveLocationSequence` assembles to 17 independently hand-derived strings (distinct fixtures), honouring listed order, 1-based inclusive slicing, whole-location and per-segment reverse-complement, the §3.5 equivalence, remote slicing with version pass-through, partials, and null-resolver degradation. No network I/O in the core. Pre-existing `ParseLocation` surface provably unchanged and still green.
- **End-state: ✅ CLEAN** — no code changed this session; re-validated fresh against the external spec. Full unfiltered `dotnet test Seqeron.sln -c Debug` → `Seqeron.Genomics.Tests` **18861 passed / 0 failed**, 0 failures across the solution, exit 0. The only residual is the documented, acceptable offline-first boundary (caller supplies the remote fetch).
