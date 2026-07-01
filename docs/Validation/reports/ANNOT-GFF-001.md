# Validation Report: ANNOT-GFF-001 — GFF3 input/output (parse + full source/score/phase export)

- **Validated:** 2026-06-25   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.ToGff3(IEnumerable<GeneAnnotation>, string)` (export), `GenomeAnnotator.ParseGff3(IEnumerable<string>)` (parse); `ComputeCdsPhases` (per-transcript cumulative CDS phase); helpers `FormatGff3Attributes`, `EncodeGff3Value`, `ParseGff3Attributes`; record fields `GeneAnnotation.Source`/`GeneAnnotation.Score`.
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs` (`GeneAnnotation` 39–48; `ParseGff3` 439–497; `ParseGff3Attributes` 499–515; `ToGff3` 538–566; `ComputeCdsPhases` 582–625; `FormatGff3Attributes` 627–645; `EncodeGff3Value` 654–678)
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_GFF3_Tests.cs` (46 tests)
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

This is a FRESH re-validation (the unit was reset to ⬜ pending after the campaign export-fidelity fix: real source/score + per-transcript cumulative CDS phase on both strands). Authoritative source — **Sequence Ontology GFF3 Specification v1.26** — was retrieved live THIS session and quoted directly (not trusting the repo's own TestSpec/Evidence/tests).

---

## Stage A — Description

### Source opened & what it confirms (retrieved this session)
Sequence Ontology GFF3 Specification v1.26, fetched from
`github.com/The-Sequence-Ontology/Specifications/.../gff3.md` (and the raw mirror). Confirmed verbatim:

1. **9 TAB-separated columns, in order:** seqid, source, type, start, end, score, strand, phase, attributes.
2. **Coordinates:** "given in positive 1-based integer coordinates, relative to the landmark given in column one" → 1-based inclusive; length = end − start + 1.
3. **source (col 2):** "a free text qualifier intended to describe the algorithm or operating procedure that generated this feature."
4. **score (col 6):** "The score of the feature, a floating point number."
5. **strand (col 7):** `+`, `-`, `.`.
6. **phase (col 8):** "For features of type 'CDS', the phase indicates where the next codon begins relative to the 5' end … the 5' end for CDS features on the plus strand is the feature's start and … on the minus strand is the feature's end." "A phase of '0' indicates that a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward) … '1' … second nucleotide … '2' … third nucleotide." "The phase is REQUIRED for all CDS features."
7. **Undefined fields** "are replaced with the '.' character."
8. **Col-9 escaping:** "tab, newline, carriage return, the percent (%) sign, and control characters must be encoded using RFC 3986 Percent-Encoding"; reserved in col 9 — `;`(%3B) `=`(%3D) `&`(%26) `,`(%2C); "unescaped spaces are allowed within fields, meaning that parsers must split on tabs, not spaces."

### Phase formula (derived from the spec definition)
The 5'-most CDS segment of a transcript opens the reading frame at its first base → phase 0. Each later
segment begins after Σ(preceding segment lengths) coding bases, so the next codon boundary is
`(3 − (Σ preceding lengths) mod 3) mod 3` bases forward. Segments are ordered 5'→3': ascending genomic
start on `+`, descending on `−` (5' end = feature end on the minus strand). This is the exact quantity the
spec's worked example tabulates.

### Independent cross-check — SO v1.26 canonical EDEN gene (ctg123), hand-recomputed
**cds00003 (+ strand), the spec's canonical multi-CDS-phase example:**

| CDS (1-based) | len = e−s+1 | Σ preceding | phase = (3 − Σ mod 3) mod 3 | Spec phase |
|---------------|-------------|-------------|------------------------------|------------|
| 3301–3902 | 602 | 0 | (3−0)%3 = **0** | 0 |
| 5000–5500 | 501 | 602 | (3−2)%3 = **1** | 1 |
| 7000–7600 | 601 | 1103 | (3−2)%3 = **1** | 1 |

→ phase sequence **0, 1, 1** — exact match to the spec.

**cds00001 (+ strand), all-multiple-of-3 segments** (1201–1500/300, 3000–3902/903, 5000–5500/501, 7000–7600/601):
Σ preceding 0, 300, 1203, 1704 (all ≡ 0 mod 3) → **0, 0, 0, 0** — exact match to the spec.

**Minus-strand worked example (hand-derived; the spec gives no minus-strand CDS example).**
Three `−` CDS segments 1-based 7000–7600 (len 601), 5000–5500 (len 501), 3301–3902 (len 602). 5'→3' order
is highest-coordinate first:

| order | seg (1-based) | len | Σ preceding | phase |
|-------|---------------|-----|-------------|-------|
| 5'-most | 7000–7600 | 601 | 0 | (3−0)%3 = **0** |
| next | 5000–5500 | 501 | 601 | (3−1)%3 = **2** |
| 3'-most | 3301–3902 | 602 | 1102 | (3−1)%3 = **2** |

Emitted in input (ascending genomic) order 3301, 5000, 7000 → phase sequence **2, 2, 0**.

### Sample emitted line (from `ToGff3`)
`chr1\tENSEMBL\tgene\t100\t500\t95.5\t+\t.\tID=gene1;product=test` — 1-based start (internal 99→100),
real source col 2, real score col 6, `.` phase for non-CDS, space in product left unescaped.

### Findings / divergences
None. Description is correct and complete against SO GFF3 v1.26. Stage A: **PASS**.

---

## Stage B — Implementation

### Code path reviewed (file:line)
- `ToGff3` (538–566): emits `##gff-version 3`; materialises input once; precomputes phases via
  `ComputeCdsPhases`; per row emits `{seqId}\t{source}\t{Type}\t{Start+1}\t{End}\t{score}\t{Strand}\t{phase}\t{attrs}`.
  `source = IsNullOrEmpty ? "." : Source`; `score = HasValue ? invariant(value) : "."`; `phase = CDS ? computed : "."`.
- `ComputeCdsPhases` (582–625): groups CDS rows by `GeneId` (Ordinal), keeps input index; orders 5'→3'
  (`−` ⇒ OrderByDescending(Start)); first segment phase 0, then `(3 − (preceding%3))%3`; accumulates
  `End − Start` (= 1-based inclusive length for half-open internal `[Start,End)`).
- `ParseGff3` (439–497): skips blank/`#`/`<9`-col lines; tolerant `TryParse` for start/end/score/phase
  (malformed → skip, no throw); `.`→null score/phase; empty strand col → `.`; `Uri.UnescapeDataString`
  decode; auto-ID `feature_{n}` when ID absent. Preserves file 1-based start/end verbatim into `GenomicFeature`.
- `EncodeGff3Value` (654–678): encodes exactly `\t \n \r % ; = & ,` + control (<0x20, 0x7F); everything
  else (incl. space) passthrough — matches the spec's "no other characters may be encoded" + spaces-allowed rule.

### Cross-verification recomputed vs code (all confirmed by passing tests)
| Check | Expected (spec) | Test |
|-------|-----------------|------|
| Plus-strand cds00003 phases | 0, 1, 1 | M28 ✓ |
| Minus-strand phases (input order) | 2, 2, 0 | M29 ✓ |
| Per-transcript independence | 0, 0 (two single-seg tx) | M30 ✓ |
| Single CDS phase | 0 | M21, M30 ✓ |
| 1-based start | internal 99 → 100 | M15/M17 ✓ |
| Real source col 2 | `ENSEMBL` | M23 ✓ |
| Source absent → `.` | `.` | M24 ✓ |
| Real score col 6 (invariant) | `95.5` | M25 ✓ |
| Score absent → `.` | `.` | M26 ✓ |
| Non-CDS phase | `.` | M20 ✓ |
| Col-9 reserved encoding | `a%3Bb%3Dc%26d%2Ce` | M19 ✓ |
| Control encoding | `a%09b%0Ac%0D%25d` | M22 ✓ |
| Space NOT encoded | `ID=gene 1` (no %20) | M16 ✓ |
| Source+score round-trip | survives parse | M27 ✓, S4 ✓ |
| Parse decode | `test protein` | M11 ✓ |
| Null score/phase parse | null | M4/M6 ✓ |
| All strands parse | + − . ? preserved | M18 ✓ |
| Malformed/short line skip | skipped | M12 ✓ |

The phase formula in code reproduces the spec's cds00003 (0,1,1) and cds00001 (0,0,0,0) tables exactly,
and the minus-strand ordering matches the spec's "5' end = feature end on the minus strand" rule.

### Variant/delegate consistency
Two record types by design: `GeneAnnotation` (0-based half-open, export source, carries Source/Score) vs
`GenomicFeature` (parse target, file 1-based). Round-trip S4 correctly asserts the 99→100 asymmetry; M27
asserts source/score survive. No `*Fast`/delegate variants.

### Numerical robustness
Score parsed/emitted with InvariantCulture (locale-safe). Phase arithmetic is `(3 − x%3)%3` on small ints —
no overflow. Parser is exception-free on malformed numerics/empty strand (tolerant TryParse).

### Test quality audit
46 tests, all green, all assert exact spec-sourced values (column values, encoded byte sequences,
coordinates, hand-derived phase sequences) — not code echoes, not "no-throw" tautologies, deterministic.
Coverage of public surface: `ParseGff3` (all 9 columns, score/phase null, strands, comments/directives/blank,
attributes incl. percent-decode, multi-feature, malformed, missing-ID, empty, only-comments), `ToGff3`
(version header, 1-based, 9-col, source present/absent, score present/absent, CDS phase plus/minus/single/
per-transcript, non-CDS `.`, encoding reserved/control/space-not-encoded, translation excluded, empty,
default seqId), round-trip both directions. Every Stage-A edge case is exercised.

Minor cosmetic note (not a defect, no fix needed): the M28 doc-comment labels its third segment "len 602"
whereas internal `(6999,7600)` is length 601; the third segment's length does not enter its own phase
(phase depends only on preceding lengths), so the assertion and result are unaffected, and 601 vs 602 both
give the same `0,1,1` because both ≡ {2 mod 3 not used here}. The asserted values match the spec.

### Findings / defects
None. Code faithfully realises the validated SO GFF3 v1.26 description on both strands. Stage B: **PASS**.

---

## Verdict & follow-ups
- **Stage A: ✅ PASS; Stage B: ✅ PASS; State: ✅ CLEAN.**
- No code or test change required this session. GFF3 fixture 46/46; full unfiltered `dotnet test Seqeron.sln`
  green (Seqeron.Genomics.Tests 18783 passed, Failed: 0; all sibling projects pass).
- Cross-check phase sequences recorded: plus-strand cds00003 = **0,1,1**; cds00001 = **0,0,0,0**;
  minus-strand (input order) = **2,2,0**. Sample line:
  `chr1\tENSEMBL\tgene\t100\t500\t95.5\t+\t.\tID=gene1;product=test`.
