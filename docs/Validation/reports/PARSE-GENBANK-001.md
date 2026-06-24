# Validation Report: PARSE-GENBANK-001 — GenBank Flat File Parser

- **Validated:** 2026-06-24   **Area:** FileIO
- **Canonical method(s):** `GenBankParser.Parse(string)`, `GenBankParser.ParseFile(string)`, `GenBankParser.ParseLocation(string)`; helpers `SequenceFormatHelper.ParseLocationParts`, `GenBankParser.FinalizeQualifierValue`/`UnquoteQualifierValue`, `FeatureLocationHelper.ExtractSequence`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (one real defect found and fully fixed this session)
- **End state:** CLEAN (defect fixed, tests added, full suite green)

## Stage A — Description

### Sources opened & what they confirm
- **INSDC Feature Table Definition** — https://www.insdc.org/submitting-standards/feature-table/ (fetched this session):
  - Single base `467`; range `340..565` = "continuous range of bases **bounded by and including** the starting and ending bases" → **1-based, inclusive**.
  - `complement(location)` = "the complement of the presented sequence in the span" → reverse complement.
  - `join(loc,loc,…)` = "placed end-to-end to form one contiguous sequence"; `order(loc,…)` = "found in the specified order (5'→3'), but nothing is implied about the reasonableness of joining them" → distinct from `join`.
  - `<`/`>` = "an end point is beyond the specified base number" (`<345..500` begins before 345 = 5' partial); indicate boundary uncertainty, **not** strand.
  - Qualifier value: free text "enclosed in double quotation marks"; embedded `"` "must be escaped by placing a second double quotation mark immediately before it (`""`)"; worked example `/note="This is an example of ""escaped"" quotation marks"`.
  - Multi-line value: "a single set of quotation marks is required at the beginning and at the end of the text"; continuation lines begin at column 22.
- **NCBI GenBank Sample Record (U49845)** — LOCUS fields (name, length+`bp`, molecule type, topology, division `PLN`, date `DD-MMM-YYYY`); section keyword order; ORIGIN = left position number + 6 groups of 10 lowercase bases (numbers/spaces stripped). (Per Evidence doc + prior report; unchanged.)
- **Biopython `Bio.GenBank.Scanner` / `__init__.py`** (reference implementation, fetched this session) — the authoritative resolution of the multi-line ambiguity the INSDC text leaves open:
  - `parse_feature()` joins continuation lines with `"\n"`; `feature_qualifier()` then does `q_value.replace("\n", " ")` → **one space inserted at each wrap point**.
  - `_BaseGenBankConsumer.remove_space_keys = ["translation"]` → for `/translation`, **all** spaces are removed so the amino-acid string reassembles contiguously.

### Formula / convention check
1-based inclusive coordinates; `complement` = reverse complement; `join`/`order` distinct; `<`/`>` partial markers; `""`→`"` unescaping; ORIGIN number/space stripping; `//` record terminator — all match INSDC/NCBI. Multi-line qualifier reconstruction rule (wrap→single space; `/translation`→no space) matches Biopython.

### Independent cross-check (by hand)
- `complement(join(150..170,180..200))`: parts (150,170),(180,200) 1-based inclusive → 0-based slices [149,170)+[179,200); concatenate then reverse-complement. Matches `FeatureLocationHelper`.
- Qualifier `/note="…""escaped""…"`: outer pair removed, `""`→`"`. Matches `UnquoteQualifierValue`.
- Multi-line `/translation` `…GHIKLMNPQRS` ⏎ `TVWYACDEFG…`: correct value joins with NO space → `…GHIKLMNPQRSTVWYACDEFG…`.

### Findings / divergences
Description: none. The TestSpec/Evidence are correct; they were silent on the multi-line wrap rule, which Biopython resolves (recorded above).

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.IO/GenBankParser.cs` (`Parse`, `ExtractSections`, `ParseLocusLine`, `ParseFeatures`, `FinalizeQualifierValue`/`UnquoteQualifierValue`, `ParseSequence`, `ParseLocation`)
- `src/Seqeron/Algorithms/Seqeron.Genomics.IO/SequenceFormatHelper.cs` (`ParseLocationParts`, `LocationRangeRegex`)
- `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FeatureLocationHelper.cs` (`ExtractSequence`)

### Realised correctly? (evidence)
- Record split on `"\n//"`, keep `LOCUS`-prefixed fragments; section extraction by column-12 header with ORIGIN short-header handling; REFERENCE blocks concatenated. ✓
- LOCUS: length `parts[2]`; molecule type / topology / division (known-code set) / date (`DD-MMM-yyyy`/`yy`) positionally. ✓
- Location: 1-based inclusive; `complement` via StartsWith (GenBank), `join`/`order` via Contains, partial via `<`/`>`, ranges `(\d+)(?:\.\.(\d+))?`, Start=min/End=max. ✓
- ORIGIN: keep `char.IsLetter`, uppercase. ✓
- Feature extraction: `realStart=partStart-1`, `realEnd=partEnd` (1-based→0-based half-open), join concatenation, complement = `GetReverseComplementString` (IUPAC-aware). ✓

### Defect found and fixed — multi-line qualifier reconstruction
**Defect:** the old continuation branch did `qualifiers[lastQual] += " " + line.Trim().Trim('"')`. Two faults, both unexercised by existing tests (every fixture qualifier fit on one line):
1. For `/translation` it **inserts a spurious space** at each wrap point, corrupting the amino-acid string (e.g. `…GHIKLMNPQRS TVWYACDEFG…`). This also corrupts `TranslateCDS`, which returns the `/translation` qualifier verbatim.
2. For any quoted multi-line value the **opening `"`** of the first fragment was never removed (single-line `UnquoteQualifierValue` was bypassed once a continuation appeared) and per-line `Trim('"')` over-strips — yielding values like `"This is a long…` with a leading quote.

**Fix:** `ParseFeatures` now buffers the raw (possibly multi-line) value per qualifier and finalizes it once via new `FinalizeQualifierValue`: wrap points (`\n`) → single space (Biopython `feature_qualifier`), then `UnquoteQualifierValue` (outer pair + `""`→`"`), then for keys in `NoSpaceQualifierKeys = {"translation"}` strip all whitespace (Biopython `remove_space_keys`). Single-line and escaped-quote paths are unchanged.

### Cross-verification recomputed vs code (new tests, all green)
- `Parse_MultiLineTranslation_ReassemblesWithoutSpuriousSpace` — wrapped `/translation` → `MKLLVVPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQR` (contiguous, joins exactly at wrap). RED before fix.
- `Parse_MultiLineFreeText_ConcatenatesWithoutInsertedSpace` — wrapped `/note` → one space at wrap, no leading quote. RED before fix (leading `"` + length off-by-one).

### Variant/delegate consistency
`ParseFile`→`Parse`; EMBL shares `SequenceFormatHelper.ParseLocationParts` and an analogous `FeatureLocationHelper` overload (out of scope, unaffected). `TranslateCDS` benefits automatically (returns the now-correct `/translation`).

### Test quality audit
`GenBankParserTests.cs` assertions use exact values for LOCUS fields, sequence equality, keyword lists, location Start/End, partial/order flags, references, extracted subsequences, and now multi-line qualifier reconstruction. The pre-existing suite had a real coverage hole (no qualifier wrapped across lines); now closed.

## Verdict & follow-ups
Stage A PASS, Stage B PASS-WITH-NOTES. One genuine defect (multi-line qualifier reconstruction: spurious space in `/translation`, unstripped opening quote in any wrapped quoted value) found, fixed to the Biopython/INSDC reference behavior, and locked with two new RED→GREEN tests. **State: CLEAN.**

- Build: `Seqeron.Genomics.Tests` — succeeded, 0 warnings.
- GenBank-filtered tests: 110 passed, 0 failed.
- Full unfiltered project: **18213 passed, 0 failed**.

### Code changed
- `src/Seqeron/Algorithms/Seqeron.Genomics.IO/GenBankParser.cs` — buffer multi-line qualifier values; new `FinalizeQualifierValue` (wrap→space, unquote, `/translation` whitespace-strip via `NoSpaceQualifierKeys`).
- `tests/Seqeron/Seqeron.Genomics.Tests/GenBankParserTests.cs` — added 2 tests + `RecordWithWrappedQualifiers` fixture.
