# Validation Report: ANNOT-GFF-001 — GFF3 feature annotation / record construction + coordinate handling

- **Validated:** 2026-06-24   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.ToGff3(IEnumerable<GeneAnnotation>, string)` (export), `GenomeAnnotator.ParseGff3(IEnumerable<string>)` (parse); helpers `FormatGff3Attributes`, `EncodeGff3Value`, `ParseGff3Attributes`.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs` (records lines 40–48; `ParseGff3` 426–484; `ParseGff3Attributes` 486–502; `ToGff3` 511–524; `FormatGff3Attributes` 526–544; `EncodeGff3Value` 553–577)
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_GFF3_Tests.cs` (38 tests)
- **Stage A verdict:** PASS-WITH-NOTES (correct; documented, spec-permitted simplifications)
- **Stage B verdict:** PASS
- **State:** CLEAN

---

## Stage A — Description

### Sources opened & what they confirm
- **Sequence Ontology GFF3 Specification v1.26** — fetched live from `github.com/The-Sequence-Ontology/Specifications/.../gff3.md` and quoted directly. Confirms every format rule below (not relying on the TestSpec's citation label).
- **Wikipedia "General feature format"** and **RFC 3986** — corroborate columns, phase, and percent-encoding base standard (consistent with the SO spec).

### Format rules confirmed verbatim from the spec
1. **9 TAB-separated columns, in order:** `seqid, source, type, start, end, score, strand, phase, attributes`.
2. **1-based, inclusive coordinates:** "given in positive 1-based integer coordinates" — both ends inclusive. Feature length = `end − start + 1`.
3. **start ≤ end always:** "Start is always less than or equal to end" (exception only for circular features).
4. **Strand values:** `+`, `-`, `.`, `?` ("? … strandedness is relevant, but unknown").
5. **Phase:** "0, 1, or 2"; "REQUIRED for all CDS features". `.` used for non-CDS.
6. **Encoding:** general — tab `%09`, newline `%0A`, CR `%0D`, percent `%25`, control chars `%00–%1F,%7F`; column-9 reserved — `;`→`%3B`, `=`→`%3D`, `&`→`%26`, `,`→`%2C`; "Unescaped spaces are allowed within fields" (spaces NOT encoded).
7. **Undefined-field placeholder:** "Undefined fields are replaced with the '.' character."
8. **Version directive:** topmost line; "version number always begins with 3, the second and third numbers are optional" — so `##gff-version 3` is valid.

### Edge-case semantics
- Export: 0-based half-open internal `[Start, End)` → 1-based inclusive `[Start+1, End]` (start +1, end unchanged). Defined and sourced.
- Parse: file's 1-based start/end preserved verbatim into `GenomicFeature` (no conversion; a different record type than the export `GeneAnnotation`). Undefined score/phase (`.`) → null. Empty/invalid strand → `.`.
- Special chars in attribute values percent-encoded; spaces preserved.

### Independent cross-check — hand computation
Export `GeneAnnotation(GeneId="gene1", Start=99, End=500, Strand='+', Type="CDS", Product="test", {Note=important})`, seqId `chr1`:
Internal `[99,500)` (0-based, half-open) = 401 positions. GFF3 emits `start=100, end=500` → 1-based inclusive `[100,500]`, length 500−100+1 = **401**. Match.
Expected line: `chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=gene1;product=test;Note=important`.

### Findings / divergences
Description correct. Documented, spec-permitted simplifications (per TestSpec "Deviations"): source and score always `.`; CDS phase always `0` (valid phase; `GeneAnnotation` has no phase field); no hierarchy/FASTA/multi-value attributes; seqid written raw (caller-controlled). None violate the spec for the supported feature set.

---

## Stage B — Implementation

### Code path reviewed
- `ToGff3` (`GenomeAnnotator.cs:511–524`): emits `##gff-version 3` (515), then per feature `{seqId}\t.\t{Type}\t{Start+1}\t{End}\t.\t{Strand}\t{phase}\t{attrs}` (522); `phase = CDS ? "0" : "."` (520, OrdinalIgnoreCase).
- `FormatGff3Attributes` (526–544): emits `ID=…;product=…` then remaining attributes, skipping `translation`.
- `EncodeGff3Value` (553–577): switch encodes `\t \n \r % ; = & ,`; control chars `<0x20`/`0x7F` via `%XX`; everything else (incl. space) passthrough.
- `ParseGff3` (426–484) / `ParseGff3Attributes` (486–502): skips blank/`#`/`##` lines and `<9`-column lines; tolerant `TryParse` for start/end/score/phase (malformed → skip, no throw); `.` → null score/phase; empty strand col → `.`; decodes via `Uri.UnescapeDataString`; auto-ID `feature_{n}` when `ID` absent.

### Coordinate check (1-based inclusive) — PASS
`{ann.Start + 1}` performs the 0-based→1-based start conversion (522); end emitted unchanged, exactly correct for half-open internal `[s,e)` → 1-based inclusive `[s+1, e]`. No off-by-one; verified by hand (401 positions) and by tests M15/M17. Parser preserves file 1-based values (M2: input `1000`→`Start==1000`) — a deliberate, documented type asymmetry, not a defect.

### Cross-verification (recomputed vs code, all confirmed by passing tests)
| Check | Input | Expected | Test |
|-------|-------|----------|------|
| 1-based start | Start=99 | start=100, end=500 | M15/M17 ✓ |
| Col-9 reserved encoding | `a;b=c&d,e` | `a%3Bb%3Dc%26d%2Ce` | M19 ✓ |
| Control encoding | `a\tb\nc\r%d` | `a%09b%0Ac%0D%25d` | M22 ✓ |
| Space NOT encoded | `gene 1` | `ID=gene 1` (no `%20`) | M16 ✓ |
| CDS phase | Type=CDS | `0` | M21 ✓ |
| non-CDS phase | Type=gene | `.` | M20 ✓ |
| Parse decode | `test%20protein` | `test protein` | M11 ✓ |
| Null score/phase | `.` | null | M4/M6 ✓ |
| Strand set | `+ - . ?` | preserved | M18 ✓ |

### Variant/delegate consistency
Two record types by design: `GeneAnnotation` (0-based, export source) vs `GenomicFeature` (parse target, file 1-based). Roundtrip test S4 correctly asserts the `99 → 100` asymmetry. No `*Fast`/delegate variants.

### Test quality audit
38 tests, all green; assertions check exact sourced values (column values, encoded byte sequences, coordinates), deterministic, and cover the Stage-A edge cases (undefined `.`, malformed lines, comments/directives/blank lines, missing ID, empty input, all strands, encoding incl. spaces-not-encoded).

### Findings / defects
None. Code faithfully realises the validated description.

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES; Stage B: PASS; State: CLEAN.** No code change. Build succeeded; filtered suite `GenomeAnnotator_GFF3_Tests` 38/0. No follow-ups.
