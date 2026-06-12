# Validation Report: ANNOT-GFF-001 — GFF3 annotation output / serialization

- **Validated:** 2026-06-12   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.ToGff3(IEnumerable<GeneAnnotation>, string)` (export), `GenomeAnnotator.ParseGff3(IEnumerable<string>)` (parse), helpers `FormatGff3Attributes`, `EncodeGff3Value`, `ParseGff3Attributes`.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs` (lines 326–454)
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_GFF3_Tests.cs` (43 test cases)
- **Stage A verdict:** PASS-WITH-NOTES (correct; documented simplifications: source/score always ".", CDS phase always "0", header `##gff-version 3`)
- **Stage B verdict:** PASS
- **State:** CLEAN

---

## Stage A — Description

### Sources opened & what they confirm
- **Sequence Ontology GFF3 Specification v1.26** (https://github.com/The-Sequence-Ontology/Specifications, raw `gff3.md`) — fetched and quoted directly. Confirms every format rule below.
- **Wikipedia "General feature format"** and **RFC 3986** — corroborate columns, phase semantics, and percent-encoding base standard (cited in TestSpec; consistent with SO spec).

### Format rules confirmed from the spec (quoted)
1. **9 TAB-separated columns, in order:** `seqid, source, type, start, end, score, strand, phase, attributes`.
2. **1-based, inclusive coordinates:** "The start and end coordinates of the feature are given in positive 1-based integer coordinates" — both ends inclusive.
3. **start ≤ end always:** "Start is always less than or equal to end" — holds even on the minus strand.
4. **Strand values:** `+`, `-`, `.`, `?` ("? can be used for features whose strandedness is relevant, but unknown").
5. **Phase:** "one of the integers 0, 1, or 2" and "REQUIRED for all CDS features". (Spec does not allow "." for a CDS phase; "." is used for non-CDS.)
6. **Empty-field placeholder:** "Undefined fields are replaced with the '.' character."
7. **Encoding:** general — tab `%09`, newline `%0A`, CR `%0D`, percent `%25`, control chars `%00–%1F,%7F`; column-9 reserved — `;`→`%3B`, `=`→`%3D`, `&`→`%26`, `,`→`%2C`; "Unescaped spaces are allowed within fields" (spaces NOT encoded); "no other characters may be encoded."
8. **Version directive:** must be topmost line; major version `3` (spec example also shows `##gff-version 3.1.26`; the major-version-only form `##gff-version 3` is the universally accepted/standard header).
9. **Reserved attribute tags are capitalized:** ID, Name, Alias, Parent, Target, Gap, Derives_from, Note, Dbxref, Ontology_term, Is_circular.

### Edge-case semantics
- 0-based half-open internal `[Start, End)` → 1-based inclusive: `start = Start+1`, `end = End` (no change to end). Defined and sourced.
- Empty score/source → "." placeholder. Non-CDS phase → ".". CDS phase → integer (0/1/2).
- Special chars in attribute values percent-encoded; spaces preserved.

### Independent cross-check — hand-constructed worked example
Feature `(GeneId="gene1", Start=99, End=500, Strand='+', Type="CDS", Product="test", Attributes={Note=important})`, seqId `chr1`. Internal `[99,500)` = 0-based indices 99..499 = **1-based 100..500 inclusive**. Expected GFF3 line (TAB-separated):

```
chr1	.	CDS	100	500	.	+	0	ID=gene1;product=test;Note=important
```

### Findings / divergences
Description is correct. Documented (and spec-permitted) simplifications, per TestSpec "Deviations": source and score always ".", CDS phase always "0" (`GeneAnnotation` has no phase field — 0 is a valid phase), no hierarchy/FASTA/multi-value attributes, seqid written raw. None violate the spec for the supported feature set.

---

## Stage B — Implementation

### Code path reviewed
- `ToGff3` — `GenomeAnnotator.cs:388–401`. Emits `##gff-version 3` header (line 392), then per feature: `{seqId}\t.\t{Type}\t{Start+1}\t{End}\t.\t{Strand}\t{phase}\t{attributes}` (line 399). `phase = CDS ? "0" : "."` (line 397).
- `FormatGff3Attributes` — lines 403–421: emits `ID=…;product=…` then remaining attributes, skipping `translation`.
- `EncodeGff3Value` — lines 430–454: switch encoding `\t \n \r % ; = & ,`, control chars `<0x20`/`0x7F`, all else (incl. space) passthrough.
- `ParseGff3` / `ParseGff3Attributes` — lines 326–379 (decodes via `Uri.UnescapeDataString`).

### Coordinate check (1-based inclusive) — PASS
- `{ann.Start + 1}` performs the required 0-based→1-based `+1` conversion on start (line 399).
- End is emitted unchanged, which is exactly correct when internal coords are half-open `[Start, End)`: `[s,e)` 0-based ↦ `[s+1, e]` 1-based inclusive. Verified against ORF construction (`FindOrfsInFrame` sets `End = i+3`, exclusive). No off-by-one.
- start ≤ end holds by construction for both strands (`FindOrfs` reverse-strand adjustment yields `adjStart < adjEnd`); export does not swap, but never needs to.

### Worked example recomputed against code
For the Stage-A feature, line 399 produces exactly:
`chr1\t.\tCDS\t100\t500\t.\t+\t0\tID=gene1;product=test;Note=important` — matches the hand-derived line. Split on `\t` = 9 fields, correct column order. Confirmed by tests `ToGff3_GeneratesValidTabDelimitedOutput` (M17) and `Roundtrip_CoreDataPreserved` (S4).

### Encoding check vs spec
`EncodeGff3Value` encodes precisely the spec-required set and nothing else; spaces pass through. Verified by:
- M16 `ToGff3_EscapesSpecialCharacters`: `ID=gene 1` (space preserved), `test%3Bproduct`.
- M19 `ToGff3_EncodesAllRequiredGff3Characters`: `a%3Bb%3Dc%26d%2Ce`.
- M22 `ToGff3_EncodesControlCharacters`: `a%09b%0Ac%0D%25d`.

### Phase / strand / placeholders
- M20 `ToGff3_NonCdsFeature_PhaseIsDot` → ".", M21 `ToGff3_CdsFeature_PhaseIsZero` → "0". Matches spec NOTE 4.
- Strand emitted raw; parse round-trips +/-/./? (M18).
- Score/source columns "." (spec placeholder).

### Test quality audit
Tests assert exact sourced values (column-by-column, exact encoded strings, exact 1-based start `100`), are deterministic, and cover the Stage-A edge cases (off-by-one start, phase ".", reserved+control char encoding, spaces preserved, roundtrip). No tautological "no-throw" assertions for the core behaviour.

### Findings / defects
None. Implementation faithfully realizes the validated spec for the supported feature subset.

---

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES — spec rules confirmed verbatim; only spec-permitted simplifications (source/score "." , CDS phase "0", `##gff-version 3` header) noted.
- **Stage B:** PASS — 9-column TAB layout, **1-based-inclusive coordinate conversion correct (start = Start+1, end unchanged; no off-by-one)**, encoding matches spec exactly, phase/strand/placeholder handling correct.
- **State:** CLEAN. No code changes required.
- **Tests:** `--filter FullyQualifiedName~GFF3` → 43 passed, 0 failed. Full suite `Seqeron.Genomics.Tests` → 4461 passed, 0 failed, 1 skipped (MFE benchmark).
