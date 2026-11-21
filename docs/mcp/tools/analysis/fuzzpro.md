# fuzzpro

EMBOSS fuzzpro: protein pattern search (PROSITE-like patterns, exact).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `fuzzpro` |
| **Method ID** | `GenomicAnalyzer.FindProteinPatternHits` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 fuzzpro: every window of a protein sequence that matches a PROSITE-like pattern: amino-acid letters (case-insensitive; B = D or N, Z = E or Q, J = I or L; x = any residue), classes [ILV], exclusions {P} (plain letters only: B, Z, J, X inside [] or {} are rejected as in the binary), repeats x(3), length ranges x(2,4) (every matching window length is reported; (0,n) makes the element optional), '-' separators, the anchors '<' (N terminus) and '>' (C terminus), a terminating '.'; e.g. [DE](2)HS{P}x(2)Px(2,4)C. The search is exact (fuzzpro has no mismatch option). A start-anchored pattern with a range reports only its longest match (the binary's regular-expression search class); a single trailing stop '*' is not part of the sequence. Positions are 0-based; hits are ordered by position, then length. Verified against fuzzpro 6.6.0 on 900 random patterns (letters, classes, exclusions, repeats, ranges, anchors; sequences with X, B, Z, U, J, O, * and lower case). For nucleotide patterns see find_pattern_hits (fuzznuc), for a translated search fuzztran. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/GenomicAnalyzer.ProteinPatterns.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.ProteinPatterns.cs)
- Algorithm doc: [Known_Motif_Search.md](../../../algorithms/Motif_Analysis/Known_Motif_Search.md)
- Pattern language: [ProteinPattern.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinPattern.cs)
- Oracle: `scripts/analysis/generate_fuzzpro_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/fuzzpro_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Protein sequence (letters; a trailing * is dropped). |
| `pattern` | string | Yes | The pattern (fuzzpro -pattern). |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: fuzzpro searches the range and counts positions from the start of the original sequence; a protein is never reversed. The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. -sformat -sask -squick -sdbname -sid -iquery -ufo -fformat -fopenfile are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused. |
| `reportFormat` | string | No | EMBOSS report format (-rformat): gff or excel. The report text of the program (the whole file; #!Date is today; sequence name from sequenceName) is returned in the result field report; only for a search without sequenceQualifiers. The other formats of the program (seqtable, table, embl, genbank, motif, ...) are not written yet. Omit for no report text. |
| `sequenceName` | string | No | Sequence name printed in the report (default sequence). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pattern` | string | The normalised pattern |
| `hits` | array | Hits `{ position (0-based), length }`, ordered by position then length |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Pattern cannot be null or empty |
| 1001 | Illegal pattern '...' |

## Examples

### Example 1: Class, wildcard and class

**User Prompt:**
> Find [AG]-x(2)-[IV] in MKVLAAGIVGLLLAQPSTAEEKV.

**Expected Tool Call:**
```json
{"tool": "fuzzpro", "arguments": {"sequence": "MKVLAAGIVGLLLAQPSTAEEKV", "pattern": "[AG]-x(2)-[IV]"}}
```

**Response:**
```json
{"pattern": "[AG]x(2)[IV]", "hits": [{"position": 4, "length": 4}, {"position": 5, "length": 4}]}
```

fuzzpro prints the 1-based windows 5-8 and 6-9.

## See Also

- [find_pattern_hits](find_pattern_hits.md) - fuzznuc nucleotide patterns
- [fuzztran](fuzztran.md) - the same patterns on translated frames
- [find_protein_motifs](find_protein_motifs.md)
