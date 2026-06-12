# Validation Report: PARSE-EMBL-001 — EMBL flat-file format parser

- **Validated:** 2026-06-12   **Area:** FileIO/Parsing
- **Canonical method(s):** `EmblParser.Parse`, `ParseFile`, `ParseLocation`, `ToGenBank`, `GetFeatures`/`GetCDS`/`GetGenes`, `ExtractSequence` (`src/Seqeron/Algorithms/Seqeron.Genomics.IO/EmblParser.cs`); location parsing delegated to `SequenceFormatHelper.ParseLocationParts`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN (no defect found)

## Stage A — Description

### Sources opened
- **EBI/ENA EMBL User Manual** (flat-file format spec), `https://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt` — confirms:
  - **Line-type prefixes** are two letters followed by three blanks, content begins at column 6. Key types confirmed: ID, AC, DE, KW, OS, FT, SQ, plus `//` entry terminator. (Also DT, OC, OG, RN/RP/RX/RG/RA/RT/RL, DR, CC, FH, XX.)
  - **ID line token order**: `ID  <accession>; SV <version>; <topology>; <molecule type>; <data class>; <taxonomic division>; <length> BP.` Example: `ID   CD789012; SV 4; linear; genomic DNA; HTG; MAM; 500 BP.`
  - **SQ header** carries length + base composition counts: `SQ   Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;`
  - **Sequence data lines**: 60 bases per line, in groups of 10 separated by a blank, beginning at position 6; columns 73–80 contain a **right-justified running base count** (number of the last base on each line) that **must be removed** during parsing.
- **INSDC DDBJ/ENA/GenBank Feature Table location syntax** (insdc.org) — confirms locations are **1-based inclusive**; `34..456` includes both endpoints; single base `467`; site-between `55^56` via caret; partial `<345..500` / `1..>888`; operators `complement()`, `join()` (end-to-end contiguous), `order()` (order without joining). `complement` may wrap `join`/`order`; `join`/`order` do not nest together. EMBL FT uses the **same** INSDC location syntax as GenBank.

### Formula / convention check
EMBL is a line-oriented flat file, not a numeric formula. All structural claims in the TestSpec and Evidence doc (line prefixes, ID field order, SQ composition, sequence-line stripping of numbers/spaces, INSDC 1-based inclusive locations) match the authoritative sources verbatim. The Evidence doc's controlled vocabularies (10 data classes incl. STS; 15 divisions; 11 INSDC mol_type values) match the cited sources.

### Edge-case semantics
Null/empty/whitespace → empty collection (defined). `KW   .` → empty keyword list. Multiple records split on the `//` terminator. Multi-line continuation for DE/KW/OS/OC/RA/RT/RL and FT qualifiers. All defined and sourced.

### Independent cross-check
Parsed the EBI X56734-style worked record (length set to 60 for a fully spelled-out sequence line): ID fields, SQ stripping, FT `complement(14..49)`, and RP/RX reference capture all reproduced (numbers below).

**Stage A findings:** none. Description is faithful to EBI EMBL manual + INSDC v11.x.

## Stage B — Implementation

### Code path reviewed
- `Parse` (EmblParser.cs:144) splits on `"\n//"`, requires each chunk to `StartsWith("ID")`.
- `ParseIdLine` (286–336): splits on `;`; accession = first token's first word; remaining tokens classified by **controlled-vocabulary membership** (topology `linear`/`circular`, mol_type set, data-class set, division set) and `…BP` for length / `SV …` for version. Robust to token position; rejects out-of-vocabulary values (bare `DNA`, `UNK`).
- `ParseSequence` (657–683): after the `SQ` line, appends only `char.IsLetter(c)` uppercased until `//`. This **strips spaces, leading indentation, and the right-justified running base count** (digits are not letters). The SQ header line itself is skipped via `continue`.
- Locations via `SequenceFormatHelper.ParseLocationParts` (1-based values kept verbatim; `..` range regex; complement/join/order flags; `<`/`>` partial flags; overall Start/End = min/max of parts).
- `ParseFeaturesFromLines` (443): column-aware FT parsing; qualifier values parsed with `IndexOf('=',1)` so `/`-containing values (e.g. `UniProtKB/Swiss-Prot:...`) are not truncated.

### Formula realised correctly?
Yes. 1-based inclusive coordinates preserved; complement/join/order/partial/single-base/site-between handled; SQ running-count and whitespace stripped per spec.

### Cross-verification table (recomputed vs code, worked record)
| Field | Expected (spec) | Code output |
|---|---|---|
| Accession / SV / Topology | X56734 / 1 / linear | X56734 / 1 / linear ✓ |
| Mol / DataClass / Division / Len | mRNA / STD / PLN / 60 | mRNA / STD / PLN / 60 ✓ |
| Sequence | 60 letters, no digits/spaces, uppercase | len 60, HasDigit=False, HasSpace=False, uppercase ✓ |
| FT source | 1..60, IsComplement=False | start 1 end 60 compl False ✓ |
| FT CDS | complement(14..49), IsComplement=True | start 14 end 49 compl True ✓ |
| Ref Positions (RP) | 1-1859 | 1-1859 ✓ |
| Ref CrossRef (RX) | DOI + PUBMED 1907511 | "DOI; 10.1007/BF00039495.; PUBMED; 1907511." ✓ |

Location unit checks reproduced: `100..200`→(100,200); `complement(100..200)`→compl,(100,200); `join(1..50,60..100)`→2 parts,(1,100); `order(...)`→IsOrder; `<1..200`→5′partial; `100..>500`→3′partial; `467`→(467,467); `123^124`→(123,124); `complement(join(...))`→compl+join. ExtractSequence complement/join verified (reverse-complement and concatenation).

### Variant/delegate consistency
`GetCDS`/`GetGenes` are `GetFeatures` filters; `ToGenBank` maps every field 1:1 (Location struct copied, references mapped). EMBL and GenBank share `ParseLocationParts` (EMBL uses `Contains("complement(")`, GenBank `StartsWith`); both give identical results on the tested cases.

### Test quality audit
106 tests in `EmblParserTests.cs` (102 spec'd + 4 parameterized expansions visible) assert **exact** values: exact sequence strings, exact field values, exact location Start/End and flags, controlled-vocabulary accept/reject (bare `DNA`→empty, `UNK`→empty), `/`-in-value preservation, RP/RG capture, multi-record sequences. Not tautological.

### Findings / defects
None. Documented limitations LIM-1…LIM-5 (site-between has no dedicated flag but `^` preserved in RawLocation; remote refs rare/unsupported; deprecated single-dot range; `""`→`"` not unescaped; dead `QualifierRegex` path) are acceptable, accurately scoped, and do not corrupt sequence or mis-parse standard locations. The earlier commit (`6ad5837`) that "fixed EmblParser tests and perf threshold" left the suite green; no fragility observed.

## Verdict & follow-ups
- **Stage A: PASS** — line-type prefixes, ID field order, SQ composition + running-count stripping, and INSDC 1-based-inclusive complement/join locations all confirmed against EBI EMBL manual + INSDC feature-table syntax.
- **Stage B: PASS** — code realises the validated spec; worked EBI record and all location/sequence edge cases recompute correctly; running base count and whitespace correctly stripped; FT complement/join 1-based coordinates correct.
- **State: CLEAN.** No code changed. Build green; 106 Embl-filtered tests pass.

## Fix applied (2026-06-12)

Resolves **LIM-4** (escaped doubled-quote not unescaped) and **LIM-5** (dead `QualifierRegex` path) from the limitations table above.

### LIM-4 — INSDC doubled-quote `""` → `"` unescape
- **Source confirmation:** INSDC Feature Table Definition, qualifier value format — free text is enclosed in double quotation marks and a literal `"` embedded in the value is encoded as **two consecutive** quotation marks (`""`). `""` is therefore the correct escape for an embedded quote.
- **Defect:** `FinishQualifier` previously did `value.ToString().Trim().Trim('"')`. `Trim('"')` strips *all* leading/trailing quote characters, so a value such as `"he said ""hi"""` was both over-stripped (trailing `"""`) and left with the inner `""` un-collapsed — producing `he said ""hi` instead of `he said "hi"`.
- **Fix:** new `UnquoteQualifierValue` helper (EmblParser.cs). When a value is enclosed in double quotes, exactly one outer pair is removed and embedded `""` pairs are collapsed to `"` via `Replace("\"\"", "\"")`. Applied **only** to quoted qualifier values inside `FinishQualifier`; unquoted values, locations, and other line types are untouched.

### LIM-5 — dead `QualifierRegex` removed
- Verified via `grep` that `QualifierRegex()` was reached only through `ParseQualifierString`, itself reached only through the unused `ParseFeatures(IEnumerable<string>)` method (the live path is `ParseFeaturesFromLines`). Removed all three (`QualifierRegex` `[GeneratedRegex]`, `ParseQualifierString`, `ParseFeatures`). `LengthRegex` and `ReferenceNumberRegex` remain in use.

### GenBank
`GenBankParser.ParseFeatures` shares the identical INSDC qualifier syntax and had the same `.Trim('"')` defect. Applied the analogous `UnquoteQualifierValue` helper there (single-line quoted values). See `PARSE-GENBANK-001.md`.

### Tests
- EmblParserTests.cs: added `Parse_Qualifier_EscapedDoubleQuote_Unescaped`, `Parse_Qualifier_NoEscapedQuotes_ValueUnchanged`, `Parse_Qualifier_MultipleEscapedQuotePairs_AllUnescaped`. The unescape tests were confirmed RED against the pre-fix code first.
- All Embl-filtered tests pass; full `Seqeron.Genomics.Tests` suite green at **4503** passed (4498 baseline + 5 new EMBL/GenBank tests), 0 failed. No existing EMBL test regressed.
