# Validation Report: PARSE-GFF-001 — GFF3/GTF Format Parser (file → features)

- **Validated:** 2026-06-24   **Area:** FileIO
- **Canonical method(s):** `GffParser.Parse(content/reader)`, `GffParser.ParseFile(path)`; supporting `ParseLine`, `ParseAttributes`, `UnescapeGff`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.IO/GffParser.cs`
- **Test files:** `tests/Seqeron/Seqeron.Genomics.Tests/GffParserTests.cs`, `GffParser_MutationKillers_Tests.cs` (55 tests total)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

> Scope: this unit is GFF3/GTF **reading/parsing** (file text → `GffRecord` features).
> Distinct from ANNOT-GFF-001 (record building / writing). `ToGff3`/`WriteToStream`/`EscapeGff`
> are exercised here only for round-trip context.

## Stage A — Description

### Sources opened & what they confirm
- **Sequence Ontology GFF3 Specification v1.26** (Lincoln Stein, 18 Aug 2020),
  https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md — fetched and
  read this session. Every claim below quoted/confirmed verbatim.
- Wikipedia "General feature format" and UCSC Genome Browser FAQ (FAQformat) — corroborate the
  9-column layout, 1-based inclusive coordinates, strand/frame conventions, and GTF
  `key "value";` attribute syntax (cited in the Evidence doc; consistent with the SO spec).

### Spec facts verified EXACTLY against GFF3 v1.26 (this session)
1. **Nine TAB-delimited columns**, in order: `seqid, source, type, start, end, score, strand,
   phase, attributes`. ✓ (spec: "comprise the nine tab-delimited columns in sequence")
2. **Coordinates are 1-based and INCLUSIVE.** Spec verbatim: *"The start and end coordinates of
   the feature are given in positive 1-based integer coordinates"* and *"Start is always less
   than or equal to end."* Therefore **feature length = end − start + 1** (contrast BED, which is
   0-based half-open, length = end − start). ✓
3. **strand** ∈ {`+`, `-`, `.`, `?`}; `.` = non-stranded, `?` = relevant-but-unknown. ✓
4. **phase** ∈ {`0`, `1`, `2`}, *"REQUIRED for all CDS features"*; `.` for non-CDS. ✓
5. **score**: undefined fields = `.`. ✓
6. **Column 9**: `tag=value` pairs separated by `;`; reserved chars percent-encoded per RFC 3986
   (`%3B` `;`, `%3D` `=`, `%26` `&`, `%2C` `,`, plus tab/newline/CR/`%25`/controls) and must be
   percent-**decoded** on read. ✓
7. **Attribute names are case-sensitive** — spec verbatim: *"attribute names are case sensitive.
   'Parent' is not the same as 'parent'."* ✓
8. **Parent multi-valued** — comma-separated (`Parent=AF2312,AB2812,abc-3`). ✓
9. **Directives** begin `##`; **comments** begin `#`. ✓

### Hand-constructed worked example (1-based inclusive)
Canonical SO CDS line (tab-separated): `ctg123 . CDS 1201 1500 . + 0 ID=cds00001;Parent=mRNA00001`
- start=**1201**, end=**1500**, score `.`→null, strand `+`, phase `0`,
  attributes `{ID=cds00001, Parent=mRNA00001}`.
- Length = 1500 − 1201 + 1 = **300** bp = **100 codons** (divisible by 3, consistent with phase 0). ✓
- Gene line `... gene 1000 9000 ...` → length 9000 − 1000 + 1 = **8001** bp (NOT 8000), confirming
  1-based inclusive. ✓

### Findings / divergences
None. Description (TestSpec + Evidence) matches GFF3 v1.26 exactly.

## Stage B — Implementation

### Code path reviewed
`GffParser.cs`: `Parse` (84–117) → `ParseLine` (119–152) → `ParseAttributes` (154–197) →
`UnescapeGff` (199–202). Source unchanged since last validation
(`git log` head for this file: `82869eff` coverage-classification, `5464a7fc` folder refactor).

### Conformance checks
- **9-column TAB split**: `line.Split('\t')` (121); columns 0–7 → seqid/source/type/start/end/
  score/strand/phase; column 8 → attributes. ✓
- **1-based INCLUSIVE preserved verbatim**: start/end parsed as raw `int` (129–132); no off-by-one
  at parse time, so length = `End − Start + 1`. The 0-based half-open conversion lives only in
  `ExtractSequence` (`Start − 1` .. `End`, 495–496) — the correct place to slice. ✓
- **score**: `.` → null; else invariant-culture float parse (134–139). ✓
- **strand**: first char of column 7; `+`/`-`/`.`/`?` retained (141). ✓
- **phase**: `.` → null; else int parse (143–145). ✓
- **GFF3 attributes**: split on `;`, then split each on the **first** `=` via `IndexOf('=')`
  (183–193). Key and value are percent-DECODED **after** the structural split, so an encoded
  `%3B`/`%3D` inside a value is never mistaken for a separator. ✓
- **Percent-DECODE**: `Uri.UnescapeDataString` — single-pass RFC 3986 (`%3B`→`;`, `%253B`→`%3B`,
  no double-decode). ✓
- **Case sensitivity**: GFF3 dict uses `StringComparer.Ordinal` (case-sensitive per spec);
  GTF/GFF2 uses `OrdinalIgnoreCase` (156–158). ✓
- **Multi-valued Parent**: kept as raw comma-joined string on parse; `BuildGeneModels` splits on
  `,` to assign children to all parents (285–294). ✓
- **Directives / comments / blank lines**: `##gff-version` sets format (digit after pragma is
  `'3'` → GFF3 else GFF2, 98–105); other `##` skipped; `#` comments skipped; blank lines skipped. ✓
- **Version mis-detection guard**: `##gff-version 2.3` → GFF2 (checks first char is `'3'`). ✓

### Cross-verification recomputed vs code (tests assert exact sourced values)
| Case | Input | Expected (spec) | Test |
|------|-------|-----------------|------|
| 1-based inclusive length | `gene 1 100`, `exon 50 75` | 100, 26 (`End−Start+1`) | `Parse_1BasedCoordinates_Validated` |
| Phase | CDS 0/1/2, gene `.` | 0,1,2,null | `Parse_Phase_ParsedCorrectly` |
| Strand | `+ - . ?` | all retained | `Parse_Strand_AllValidValues` |
| Score `.` / value | undefined / 99.5 | null / 99.5 | `Parse_NoScore_ScoreIsNull`, `Parse_WithScore_ScoreIsParsed` |
| Percent decode | `Name=Test%3BGene` | `Test;Gene` | `Parse_SpecialCharacters_Unescaped` |
| No double-decode | `Name=Test%253BGene` | `Test%3BGene` | `Parse_DoubleEncodedPercent_DecodedCorrectly` |
| Multi-Parent | `Parent=mRNA1,mRNA2` | raw, split in models | `Parse_MultipleParentValues`, `BuildGeneModels_MultiParentExons…` |
| Case-sensitive attrs | `ID` ≠ `id` | only `ID` present | `Parse_AttributeCaseSensitive` |
| Comments/blank/directive | `#`/blank/`##` | skipped | `Parse_SkipsComments`, `Parse_SkipsEmptyLines`, `Parse_DetectsGFF3Version` |
| GFF2 not GFF3 | `##gff-version 2.3` | GFF2 (case-insensitive) | `Parse_VersionDetection_DoesNotMisdetectGFF2AsGFF3` |
| Malformed | `< 8` fields | line skipped | `Parse_MalformedLine_Skips` |
| GTF attrs | `key "value";` | unquoted value, case-insensitive | `Parse_GTFAttributes_*` |

### Test quality audit
Assertions check exact sourced values (start/end ints, `End−Start+1` lengths, decoded strings,
null score/phase), are deterministic, and cover every Stage-A edge case. Earlier weak
`GreaterThanOrEqualTo`/palindrome assertions were already replaced with exact `EqualTo`/non-
palindromic checks (per TestSpec). 55 tests pass.

### Findings / minor notes (not defects)
- **Lenient column count**: `ParseLine` rejects only `< 8` fields (122); a line with exactly 8
  (no attributes column) parses with an empty attribute set. Tolerant of real-world files that
  omit the trailing column; never corrupts a field. Not a correctness defect.
- **Strand stored as first char**: a multi-character strand field keeps only `[0]`. In practice
  strand is a single character; the valid set is fully covered.

## Verdict & follow-ups
- Stage A: **PASS** — description matches GFF3 v1.26 exactly (9 TAB columns; 1-based inclusive,
  length = end − start + 1; strand {+,−,.,?}; phase {0,1,2,.} required for CDS; column-9
  `key=value;` with RFC-3986 percent-decoding; case-sensitive tags; multi-valued Parent).
- Stage B: **PASS** — code realises the validated description; coordinates preserved verbatim
  with the 0-based conversion isolated in `ExtractSequence`; percent-decode applied after the
  structural split; case-sensitivity and version detection correct.
- **End state: CLEAN.** No code changed. Filtered suite: 55 passed / 0 failed.
