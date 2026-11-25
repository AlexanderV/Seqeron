# Evidence Artifact: ANNOT-GFF-001

**Test Unit ID:** ANNOT-GFF-001
**Algorithm:** GFF3 I/O (`GenomeAnnotator.ParseGff3`, `GenomeAnnotator.ToGff3`; MCP `parse_gff3`, `to_gff3`)
**Date Collected:** 2026-11-25 (review-2026-09 B11, F8–F10)

---

## Online Sources (what was actually opened)

### Sequence Ontology GFF3 Specification v1.26 (L. Stein, 18 Aug 2020)
**URL:** https://raw.githubusercontent.com/The-Sequence-Ontology/Specifications/master/gff3.md (curl, full text, 870 lines)
**Authority rank:** 1 (the format standard).
1. "nine-column, tab-delimited"; percent-encode tab, newline, CR, `%`, control characters (%00–%1F, %7F),
   "no other characters may be encoded"; `; = & ,` "have reserved meanings in column 9"; "unescaped spaces are
   allowed within fields"; "+" for space "is no longer allowed".
2. Column 1: "must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|] … may not contain unescaped
   whitespace and must not begin with an unescaped '>'".
3. Columns 4/5: "positive 1-based integer coordinates … Start is always less than or equal to end"; zero-length
   features: "start equals end and the implied site is to the right of the indicated base".
4. Column 7: `+`, `-`, `.`, `?`. Column 8: "one of the integers 0, 1, or 2", measured from the 5' end (start on +,
   end on −); "The phase is REQUIRED for all CDS features".
5. Column 9: tag=value pairs separated by ';'; "URL escaping rules are used for tags or values containing ',=;'";
   "Multiple attributes of the same type are indicated by separating the values with the comma"
   (`Parent=AF2312,AB2812,abc-3`); Parent, Alias, Note, Dbxref, Ontology_term may be multi-valued; tags are case sensitive.
6. Directives: `##gff-version 3` topmost; `###`; `##FASTA` ends the annotation section ("Once a ##FASTA section is
   encountered no other content beyond valid FASTA sequence is allowed"); "a GFF line that begins with the
   character > creates an implied ##FASTA directive".
7. The Canonical Gene (EDEN, lines 2–24) — used verbatim (tab-delimited) as the reference dataset.

### Reference implementations (run locally)
- **GenomeTools 1.6.5** `gt gff3validator`, `gt gff3 -retainids` (Debian package `genometools`) — the reference
  validator: canonical gene → "input is valid GFF3"; per-rule error messages quoted in the tests.
- **gffutils 0.14** (`feature_from_line`) and **BCBio-GFF** (PyPI, `GFF.parse`) — attribute value lists.
- Independent Python implementation of the spec line rules (scratch oracle, `urllib.parse.unquote` decoding).

---

## Stage A verdict — 🟡 PASS-WITH-NOTES

Columns, 1-based coordinates, "." placeholder, column-9 encoding set, CDS phase rule and the `##gff-version 3`
header were correct (the phase accumulation had been fixed earlier). Gaps against the spec: multi-valued
attributes not modelled (an escaped comma and a separator comma became the same string); seqid/source dropped on
parse; no column-domain validation (strand, phase, coordinates, column count, column-9 tokens); value whitespace
trimmed; `##FASTA` not recognised; on export seqid/source/type/tags not encoded, ID/product duplicable, empty
values emitted, invalid strands and empty intervals written.

## Stage B — implementation vs references

### Defect repros (old code → reference)
| Input | Old result | Reference |
|---|---|---|
| `Note=a%2Cb` vs `Note=a,b` | both `"a,b"` | gffutils `['a,b']` vs `['a','b']` |
| `Note= x y ` | `"x y"` | gt / gffutils `" x y "` |
| strand `Z`, strand `''`, phase `3`/`-1`, start `0`, start > end, 10 columns, repeated tag | all accepted | gt gff3validator error for each |
| ToGff3 source `"my\tsrc"` | 10-column row | gt: "does not contain 9 tab separated fields" |
| ToGff3 Attributes `{"ID":"x"}` | `ID=g1;product=p;ID=x` | gt: "more than one ID attribute" |
| ToGff3 tag `a=b` | `a=b=x` | gt: "does not contain exactly one '='" |
| ToGff3 product `""` | `product=` | gt: "has no value" |
| ToGff3 strand `x` / interval [5,5) | written | gt: invalid strand / "start '6' is larger then end '5'" |
| ToGff3 seqid `chr 1`, `>chr1` | raw | spec column 1 (gt does not check seqid) |

### Oracle numbers (fixed code)
- **Parse, random corpus** (1,500 lines; seed 20261125; valid/invalid column values, escapes, multi-values,
  spaces, CR endings, 8/10 columns): C# `ParseGff3` vs independent spec oracle → 528 / 528 accepted features,
  **0 field mismatches** (FeatureId incl. `feature_n` numbering, Type, Start, End, Strand, Score, Phase,
  Attributes, Seqid, Source, AttributeValues).
- **Acceptance vs gt gff3validator** (same 1,500 lines, one file per line): 1,392 agree, 108 disagree — all
  explained by gt's C-library number parsing: 82 signed starts (`+7` read as 7; `-5` "reset to 1" with a warning),
  26 score `1,5` (strtod prefix → 1). The spec requires positive integers / a floating-point number, so these
  lines stay rejected.
- **Value lists vs gffutils 0.14:** 437 match, 58 differ; on lines without spaces in column 9: 170 / 0. All
  differences: gffutils joins a value list back into one string when a value contains spaces.
- **Value lists vs BCBio-GFF:** on lines without spaces 163 match / 7 differ — all 7 have an escaped `=` in a tag
  (BCBio decodes before splitting on `=`); with spaces BCBio strips values or switches to GFF2 parsing.
- **Canonical gene export:** 23 annotations → ToGff3 → rows identical to the spec (23/23), CDS phases 13/13,
  gt valid. Mirrored onto the − strand (x → 1497229 − x): CDS phases 13/13 equal (mirror-invariant), gt valid.
- **Export escaping corpus** (300 annotations; random seqid/source/ID/product/tag/value/multi-values over a pool
  with tab, LF, CR, `% ; = & ,`, `>`, space, quotes, U+0001, U+007F, é, 漢): line == independent spec encoder
  300/300; gt gff3validator valid 300/300; ToGff3 → ParseGff3 exact round trip 290/300 — the 10 misses are tags
  with trailing spaces, which the parser trims like gt / gffutils / BCBio.

### Choices recorded
- CDS with phase "." is accepted on parse (gt accepts it; spec says required) — export always writes the phase.
- Malformed lines are skipped, not thrown (the parser's documented discipline); gt rejects the file.
- Tags are trimmed (gt, gffutils, BCBio); values are kept verbatim.
- `Attributes` (string → string) is kept for compatibility; `AttributeValues` carries the exact lists.
