# fuzztran

EMBOSS fuzztran: protein pattern search after translation of the chosen frames.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `fuzztran` |
| **Method ID** | `GenomicAnalyzer.FindTranslatedPatternHits` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 fuzztran: translates a nucleotide sequence in the chosen frames (frame 1, 2, 3, -1, -2, -3, F = forward three, R = reverse three, 6 = all; transeq conventions: reverse frame -k is the reverse complement read from its k-th base, dangling bases give X) with an EMBOSS genetic-code table (0 = standard, else the NCBI table number) and searches every translation with a fuzzpro protein pattern (see fuzzpro: letters, x, [ ], { }, repeats, ranges, anchors; exact). Each hit gives the frame, the 1-based protein start / end, the matched residues, the forward-strand nucleotide location (position 0-based, length = 3 per residue) and embossStart / embossEnd, the two numbers the program prints (exact for the forward frames; EMBOSS 6.6.0 prints negative numbers for the reverse frames - a defect of its source-feature mapping, reproduced; use position / length). Hits per frame in the order 1, 2, 3, -1, -2, -3, by protein position, then length. Verified against fuzztran 6.6.0 on 760 random DNA / pattern / frame / table cases (all frame selections, tables 0-11, N / IUPAC codes and lower case). Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Analysis/GenomicAnalyzer.ProteinPatterns.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.ProteinPatterns.cs)
- Algorithm doc: [Known_Motif_Search.md](../../../algorithms/Motif_Analysis/Known_Motif_Search.md)
- Oracle: `scripts/analysis/generate_fuzztran_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/fuzztran_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive). |
| `pattern` | string | Yes | The protein pattern (fuzztran -pattern). |
| `frame` | string | No | Frames to translate (fuzztran -frame): 1 (default), 2, 3, -1, -2, -3, F, R or 6. |
| `table` | integer | No | EMBOSS genetic code (fuzztran -table): 0 = standard (default), otherwise an NCBI table number. |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: fuzztran cuts the (reverse-complemented) text at the begin / end AS GIVEN, numbers the hits from the first position of that text and prints negative numbers for the reverse frames (a defect of the program, reproduced). The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. -sformat -sask -squick -sdbname -sid -iquery -ufo -fformat -fopenfile are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused. |
| `reportFormat` | string | No | EMBOSS report format (-rformat): any of simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq (diffseq needs 7 report tags: refused as the program does; embl / genbank are nucleotide only, pir / swiss protein only), or the report qualifiers as typed after the files, e.g. '-rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2'. The report FILE of the program (header with the command line, one block per sequence, footer; #Rundate / #!Date is now; sequence name from sequenceName) is returned in the result field report. Omit for no report text. |
| `sequenceName` | string | No | Sequence name printed in the report (default sequence). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `pattern` | string | The normalised pattern |
| `hits` | array | Hits `{ frame, proteinStart, proteinEnd, translation, position (0-based forward-strand start), length (nt), embossStart?, embossEnd? (forward frames only) }` |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1001 | Pattern cannot be null or empty |
| 1001 | frames must be 1, 2, 3, -1, -2, -3, F, R or 6 (fuzztran -frame). |

## Examples

### Example 1: Forward frame

**User Prompt:**
> Search KV in the translation of ATGAAAGTTCTGGCT.

**Expected Tool Call:**
```json
{"tool": "fuzztran", "arguments": {"sequence": "ATGAAAGTTCTGGCT", "pattern": "KV"}}
```

**Response:**
```json
{"pattern": "KV", "hits": [{"frame": 1, "proteinStart": 2, "proteinEnd": 3, "translation": "KV", "position": 3, "length": 6, "embossStart": 4, "embossEnd": 9}]}
```

fuzztran prints the feature `4 9` (frame 1, *start 2, *end 3, *translation KV).

## See Also

- [fuzzpro](fuzzpro.md) - the pattern language on a protein
- [find_pattern_hits](find_pattern_hits.md) - fuzznuc
- [translate_dna](translate_dna.md)
