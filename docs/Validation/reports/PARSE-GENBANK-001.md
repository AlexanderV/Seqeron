# Validation Report: PARSE-GENBANK-001 — GenBank Flat File Parser

- **Validated:** 2026-06-12   **Area:** FileIO
- **Canonical method(s):** `GenBankParser.Parse(string)`, `GenBankParser.ParseFile(string)`, `GenBankParser.ParseLocation(string)`; helpers `SequenceFormatHelper.ParseLocationParts`, `FeatureLocationHelper.ExtractSequence`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN (no defect found)

## Stage A — Description

### Sources opened & what they confirm
- **NCBI GenBank Sample Record (U49845)** — https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
  - LOCUS line `LOCUS  SCU49845  5028 bp  DNA  PLN  21-JUN-1999`: fields = locus name, length+`bp`, molecule type, (topology), division (`PLN` = plant/fungal/algal), modification date `DD-MMM-YYYY`.
  - Section keyword order: DEFINITION → ACCESSION → VERSION → KEYWORDS → SOURCE → ORGANISM → REFERENCE → FEATURES → ORIGIN → `//` terminator.
  - FEATURES table: feature key + location + qualifiers in `/key="value"` form (e.g. `CDS 687..3158 /gene="AXL2" /product="Axl2p"`).
  - Location syntax: span `687..3158`, partial 5' `<1..206`, complement `complement(3300..4037)`, 3' partial uses `>`.
  - ORIGIN block: line number at left, then 6 groups of 10 bases (60 bases/line); sequence shown lowercase. Parser must strip numbers and spaces.
  - LOCUS declares 5028 bp for U49845 (the ORIGIN reconstructs to 5028 bases; the WebFetch summary's "4981" was a fetch artifact, not the record length).
- **INSDC Feature Table Definition** — https://www.insdc.org/submitting-standards/feature-table/
  - §3.4.2: "first base (5' end) … is base 1" → **1-based**. §3.4.2.1: spans `n..m` are "bounded by and including the starting and ending bases" → **inclusive**.
  - `complement(location)`: read the complement of the span in 5'→3' direction (i.e. reverse complement).
  - `join(...)`: elements placed end-to-end to form one contiguous sequence.
  - `order(...)`: elements in 5'→3' order, *nothing implied about reasonableness of joining* — distinct from `join`.
  - `<` / `>`: endpoint lies beyond the specified base (`<345` starts before 345 = 5' partial; `1..>888` continues past 888 = 3' partial).

### Findings / divergences
None. Every claim in the TestSpec and Evidence doc (section keywords, 1-based-inclusive coords, complement = reverse complement, join vs order, partial `<`/`>`, ORIGIN number/space stripping, `//` terminator) matches the authoritative NCBI/INSDC sources.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.IO/GenBankParser.cs`
- `src/Seqeron/Algorithms/Seqeron.Genomics.IO/SequenceFormatHelper.cs` (`ParseLocationParts`, `LocationRangeRegex`)
- `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FeatureLocationHelper.cs` (`ExtractSequence`)
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs:149` (`GetReverseComplementString`)

### Realised correctly? (evidence)
- **Record splitting / `//` terminator:** `Parse` splits on `"\n//"` and keeps fragments starting with `LOCUS` (`GenBankParser.cs:97-108`). Multi-record verified (REC1/REC2/REC3).
- **Section extraction:** column-12 header detection with ORIGIN short-header handling; REFERENCE blocks concatenated (`ExtractSections`, `SaveSection`). FEATURES preserves leading whitespace for column-based key/qualifier parsing.
- **LOCUS fields:** length parsed from `parts[2]`; molecule type / topology / division / date detected positionally with a known-division-code set and `DD-MMM-yyyy`/`DD-MMM-yy` date formats — matches NCBI fields M1.1–M1.6.
- **ORIGIN reconstruction (`ParseSequence`, line 536):** keeps only `char.IsLetter` and uppercases → numbers, spaces, newlines stripped. Confirms M4.2–M4.4 and INSDC 7.4.1 IUPAC retention.
- **Location parsing:** 1-based inclusive; `complement` via StartsWith (GenBank convention), `join`/`order` via Contains, partial via `<`/`>` Contains; ranges via `(\d+)(?:\.\.(\d+))?`; overall Start = min part, End = max part. Matches INSDC §3.4.2.
- **Feature extraction (`FeatureLocationHelper`):** `realStart = partStart-1`, `realEnd = partEnd` (1-based→0-based half-open), join concatenates parts, complement applies `GetReverseComplementString` (true reverse + complement, IUPAC-aware). Matches INSDC complement/join semantics.

### Worked examples recomputed against the actual code
Driver compiled against the real `Seqeron.Genomics.IO` assembly:
- ORIGIN `acgtacgtac gt` (numbered) → `Sequence="ACGTACGTACGT"` (len 12): numbers/spaces stripped, uppercased. ✓
- `complement(1..4)` over `ACGT` → extracted `ACGT` (palindrome). ✓
- `complement(1..3)` over `AACCGGTT` → positions 1..3 = `AAC`, revcomp = `GTT`; code returns `GTT`. ✓ (confirms true reverse complement, not just complement)
- `join(1..4,9..12)` over `ACGTACGTACGT` → `ACGT`+`ACGT` = `ACGTACGT`, 2 parts. ✓
- `<1..>12` → `Is5PrimePartial=True`, `Is3PrimePartial=True`. ✓
- LOCUS `MINI 12 bp DNA linear UNK 21-JUN-1999` → Len=12, Mol=DNA, Topo=linear, Div=UNK, Date=1999-06-21. ✓

### Test quality audit
`GenBankParserTests.cs` (52 declared; 66 test cases discovered in this class): assertions use exact values — LOCUS fields, sequence equality, exact keyword list, location Start/End, partial/order flags, exact reference fields and PubMed ID, exact extracted subsequence. Edge cases covered: empty/null/whitespace content, no-LOCUS, no-features, empty `KEYWORDS` ".", minimal record, multi-record, IUPAC ambiguity codes. Tests genuinely fail on misparse (no permissive Contains-only or no-throw tautologies for core assertions).

### Findings / defects
None. Minor non-defect observations (out of scope, no spec impact): `ParseLocation` partial detection uses substring `Contains('<'/'>')` which is safe for INSDC location strings where those characters appear only as partial markers; `order()` shares the part-extraction path with `join()` but is correctly distinguished by the `IsOrder`/`IsJoin` flags.

## Verdict & follow-ups
Stage A PASS, Stage B PASS. Implementation faithfully realises the NCBI/INSDC-validated GenBank format: section parsing, 1-based-inclusive coordinates, complement (reverse-complement), join, order, partial `<`/`>`, ORIGIN number/space stripping, and `//` record termination all confirmed by source review and recomputed worked examples. No code changes required. **State: CLEAN.**

- Build: `Seqeron.Genomics.Tests` — succeeded, 0 warnings.
- Tests: `--filter FullyQualifiedName~GenBank` → 66 passed, 0 failed.

## Fix applied (2026-06-12)

Applied alongside the EMBL fix (see `PARSE-EMBL-001.md`, "Fix applied"). GenBank uses the **same** INSDC qualifier syntax as EMBL, so it shared the same defect.

- **Source:** INSDC Feature Table Definition (qualifier value format) — an embedded `"` is encoded as two consecutive quotes (`""`).
- **Defect:** `ParseFeatures` set a single-line qualifier value with `qualLine[(eqIdx + 1)..].Trim('"')`, which over-strips and does not collapse embedded `""`.
- **Fix:** added a `UnquoteQualifierValue` helper (mirrors EmblParser): removes exactly one outer quote pair and collapses `""` → `"`; unquoted values returned unchanged. Applied only to the single-line quoted-value path; locations and unquoted values untouched.
- **Tests added** (GenBankParserTests.cs): `Parse_Qualifier_EscapedDoubleQuote_Unescaped` (confirmed RED first), `Parse_Qualifier_NoEscapedQuotes_ValueUnchanged`.
- GenBank-filtered tests pass (68); full suite green at 4503, 0 failed; no regressions.
