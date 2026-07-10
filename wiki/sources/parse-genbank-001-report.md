---
type: source
title: "Validation report: PARSE-GENBANK-001 (GenBank flat-file parser — LOCUS/feature location/qualifiers/ORIGIN)"
tags: [validation, file-io, governance]
doc_path: docs/Validation/reports/PARSE-GENBANK-001.md
sources:
  - docs/Validation/reports/PARSE-GENBANK-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: PARSE-GENBANK-001

The two-stage **validation write-up** for test unit **PARSE-GENBANK-001** (GenBank flat-file
parser — LOCUS parsing, INSDC feature-location descriptors, multi-line qualifier reconstruction,
and ORIGIN sequence extraction), validated 2026-06-24. This is the *report* artifact that feeds
one row of the [[validation-ledger]]; it records the validator's **verdict** on both the algorithm
description and the shipped code. The two-stage methodology is the [[validation-protocol]]; the
location grammar it validates is summarized in [[insdc-feature-location]]. Distinct from the
pre-implementation [[parse-genbank-001-evidence]] artifact.

## Verdict

**Stage A: PASS · Stage B: PASS-WITH-NOTES · State: CLEAN.** The description is faithful to the
INSDC/NCBI standard (no defect). One genuine **code** defect — **multi-line qualifier
reconstruction** — was found and fully fixed to the Biopython/INSDC reference behaviour and locked
with two new RED→GREEN tests. GenBank-filtered tests **110 passed / 0 failed**; full unfiltered
project **18213 passed / 0 failed**; build 0 warnings.

## Stage A — description (algorithm faithfulness)

- Canonical methods: `GenBankParser.Parse(string)`, `ParseFile(string)`, `ParseLocation(string)`;
  helpers `SequenceFormatHelper.ParseLocationParts`, `GenBankParser.FinalizeQualifierValue` /
  `UnquoteQualifierValue`, `FeatureLocationHelper.ExtractSequence`.
- Sources opened this session: **INSDC Feature Table Definition** (fetched) — single base `467`;
  range `340..565` = **1-based, inclusive**; `complement(location)` = reverse complement;
  `join(...)` end-to-end contiguous vs `order(...)` (no adjacency implied); `<`/`>` = end point
  beyond the stated base (boundary uncertainty, not strand); qualifier value in double quotes with
  embedded `""` → `"`; multi-line value with one quote pair, continuation at column 22. **NCBI
  Sample Record U49845** — LOCUS fields, section order, ORIGIN layout. **Biopython
  `Bio.GenBank.Scanner`** (fetched) — the authoritative resolution of the multi-line ambiguity:
  `parse_feature()` joins continuation lines with `\n`, `feature_qualifier()` does
  `replace("\n"," ")` (one space per wrap point), and `remove_space_keys = ["translation"]` strips
  **all** spaces for `/translation` so the amino-acid string reassembles contiguously.
- Formula / convention check: 1-based inclusive coordinates; `complement` = reverse complement;
  `join`/`order` distinct; `<`/`>` partial markers; `""`→`"` unescaping; ORIGIN number/space
  stripping; `//` record terminator — all match INSDC/NCBI. Multi-line reconstruction (wrap→single
  space; `/translation`→no space) matches Biopython.
- Independent cross-check (by hand): `complement(join(150..170,180..200))` → 0-based slices
  [149,170)+[179,200) then reverse-complement (matches `FeatureLocationHelper`); `/note` with
  `""escaped""` → outer pair removed, `""`→`"` (matches `UnquoteQualifierValue`); multi-line
  `/translation` joins with **no** space.
- **Stage A finding: PASS, none.** The TestSpec/Evidence are correct; they were merely silent on the
  multi-line wrap rule, which Biopython resolves (recorded above).

## Stage B — implementation (code review + cross-check)

- Code path: `GenBankParser.cs` (`Parse`, `ExtractSections`, `ParseLocusLine`, `ParseFeatures`,
  `FinalizeQualifierValue` / `UnquoteQualifierValue`, `ParseSequence`, `ParseLocation`);
  `SequenceFormatHelper.cs` (`ParseLocationParts`, `LocationRangeRegex`); `FeatureLocationHelper.cs`
  (`ExtractSequence`).
- Realised correctly: record split on `"\n//"` keeping `LOCUS`-prefixed fragments; section
  extraction by column-12 header with ORIGIN short-header handling; LOCUS length / molecule type /
  topology / division / date parsed positionally; location 1-based inclusive with `complement` via
  StartsWith, `join`/`order` via Contains, partial via `<`/`>`; ORIGIN keeps `char.IsLetter`
  uppercased; feature extraction `realStart=partStart-1`, `realEnd=partEnd` (1-based→0-based
  half-open), join concatenation, complement = IUPAC-aware `GetReverseComplementString`. ✓
- **DEFECT found & fixed — multi-line qualifier reconstruction.** The old continuation branch did
  `qualifiers[lastQual] += " " + line.Trim().Trim('"')`, with two faults (both unexercised because
  every existing fixture qualifier fit on one line): (1) for `/translation` it inserted a spurious
  space at each wrap point, corrupting the amino-acid string (and `TranslateCDS`, which returns
  `/translation` verbatim); (2) for any quoted multi-line value the opening `"` of the first fragment
  was never removed and per-line `Trim('"')` over-strips, yielding a leading-quote value. **Fix:**
  `ParseFeatures` now buffers the raw multi-line value per qualifier and finalizes once via new
  `FinalizeQualifierValue` — wrap points (`\n`) → single space (Biopython `feature_qualifier`), then
  `UnquoteQualifierValue` (outer pair + `""`→`"`), then for keys in `NoSpaceQualifierKeys =
  {"translation"}` strip all whitespace (Biopython `remove_space_keys`). Single-line and escaped-quote
  paths unchanged.
- Cross-verification (new tests, all green): `Parse_MultiLineTranslation_ReassemblesWithoutSpuriousSpace`
  → contiguous amino-acid string joining exactly at the wrap (RED before fix);
  `Parse_MultiLineFreeText_ConcatenatesWithoutInsertedSpace` → one space at wrap, no leading quote
  (RED before fix).
- Variant/delegate consistency: `ParseFile`→`Parse`; EMBL shares `ParseLocationParts` and an
  analogous `FeatureLocationHelper` overload (out of scope, unaffected); `TranslateCDS` benefits
  automatically.
- Test-quality audit: exact-value asserts for LOCUS fields, sequence equality, keyword lists,
  location Start/End, partial/order flags, references, extracted subsequences, and now multi-line
  qualifier reconstruction — the pre-existing real coverage hole (no qualifier wrapped across lines)
  is now closed.

## Findings

- **One genuine code defect (multi-line qualifier reconstruction: spurious space in `/translation`,
  unstripped opening quote in any wrapped quoted value), found and fixed** to the Biopython/INSDC
  reference behaviour and locked with two new RED→GREEN tests. State CLEAN.
- **Code changed:** `GenBankParser.cs` — buffer multi-line qualifier values + new
  `FinalizeQualifierValue` (wrap→space, unquote, `/translation` whitespace-strip via
  `NoSpaceQualifierKeys`); `GenBankParserTests.cs` — added 2 tests + `RecordWithWrappedQualifiers`
  fixture.

See the full report at `docs/Validation/reports/PARSE-GENBANK-001.md`.
