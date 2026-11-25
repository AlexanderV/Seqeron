# find_pattern_hits

EMBOSS fuzznuc pattern search (classes, exclusions, repeats, ranges, anchors, mismatches, both strands).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_pattern_hits` |
| **Method ID** | `GenomicAnalyzer.FindPatternHits` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

Search a DNA sequence with an EMBOSS fuzznuc pattern (the -pattern language): IUPAC letters, classes [ACG], exclusions {AG}, repeats N(3), length ranges A(2,4) (every matching window length is reported), '-' separators, the anchors '<' (start of the sequence) and '>' (end), a terminating '.'; e.g. [CG](5)TG{A}N(1,5)C. maxMismatches = fuzznuc -pmismatch (Hamming, substitutions only; with a range the extra repeats beyond n of X(n,m) must match exactly, as in fuzznuc); bothStrands = -complement (the reverse-strand pattern is the reversed, complemented one with swapped anchors; positions are forward-strand starts). Give one pattern, or patternFile = the text of a fuzznuc pattern file (one pattern per line, or FASTA '>name' + pattern) for several named patterns. Verified against EMBOSS fuzznuc 6.6.0 on 3300 random patterns. By default a symbol other than A/C/G/T in the sequence matches nothing and costs a mismatch; with sequenceAmbiguity = true the sequence is read as the fuzznuc binary reads it (1350 more oracle rows with ambiguity codes in the sequence): a pattern IUPAC letter matches every sequence symbol whose base set lies within the letter's (pattern R matches A, G, R; sequence R matches pattern R, V, D, N but not A, M, Y), N matches every symbol, a class [ACG] matches its plain bases only, an exclusion {AG} matches every symbol not listed (ambiguity codes included), U reads as T, and the reverse strand is the IUPAC reverse complement of the sequence (the reported patterns are the plain-base normal form). Positions 0-based; hits ordered by position, length, forward before reverse. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Emboss.Motifs/EmbossMotifPrograms.PatternHits.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.Motifs/EmbossMotifPrograms.PatternHits.cs)
- Algorithm doc: [Known_Motif_Search.md](../../../algorithms/Motif_Analysis/Known_Motif_Search.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive). |
| `pattern` | string | No | One fuzznuc pattern (omit when patternFile is given). |
| `patternFile` | string | No | Text of a fuzznuc pattern file: one pattern per line, or FASTA style (>name then the pattern); several patterns are reported under their names (unnamed: pattern1, pattern2, ...). |
| `maxMismatches` | integer | No | Maximum substitutions per hit (fuzznuc -pmismatch, default 0). (default: 0) |
| `bothStrands` | boolean | No | Also search the reverse-complement strand (fuzznuc -complement, default false). (default: false) |
| `sequenceAmbiguity` | boolean | No | Read IUPAC ambiguity codes in the SEQUENCE as the fuzznuc binary does (default false: they match nothing). |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: fuzznuc searches the range of the (reverse-complemented) sequence and reports the ORIGINAL strand: positions are original, the strand is flipped for -sreverse; the sequence is read as fuzznuc reads it (ambiguity codes on). The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. -sformat selects the file format of the text; -ufo -fformat -fopenfile name the feature file of a program that takes the features of its sequences (diffseq: the files are given by the tool); -sask -squick -sdbname -sid -iquery are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused. |
| `reportFormat` | string | No | EMBOSS report format (-rformat): any of simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq (diffseq needs 7 report tags: refused as the program does; embl / genbank are nucleotide only, pir / swiss protein only), or the report qualifiers as typed after the files, e.g. '-rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2'. The report FILE of the program (header with the command line, one block per sequence, footer; #Rundate / #!Date is now; sequence name from sequenceName) is returned in the result field report; needs one pattern (pattern, not patternFile). Omit for no report text. |
| `sequenceName` | string | No | Sequence name printed in the report (default sequence). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `patterns` | object | Pattern name -> normalised pattern text |
| `hits` | array | Hits `{ pattern, position, length, strand ("forward"/"reverse"), mismatches }`, 0-based forward-strand starts, ordered by position, length, forward before reverse |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1001 | Give exactly one of pattern and patternFile |
| 1001 | Illegal pattern ... |
| 1002 | maxMismatches cannot be negative |

## Examples

### Example 1: Length range

**User Prompt:**
> Find A(1,3)C in AAAACAAAAT.

**Expected Tool Call:**
```json
{"tool": "find_pattern_hits", "arguments": {"sequence": "AAAACAAAAT", "pattern": "A(1,3)C"}}
```

**Response:**
```json
{"patterns": {"pattern1": "A(1,3)C"}, "hits": [{"pattern": "pattern1", "position": 1, "length": 4, "strand": "forward", "mismatches": 0}, {"pattern": "pattern1", "position": 2, "length": 3, "strand": "forward", "mismatches": 0}, {"pattern": "pattern1", "position": 3, "length": 2, "strand": "forward", "mismatches": 0}]}
```

fuzznuc reports every window length: 1-based 4..5, 3..5, 2..5.

### Example 2: Anchored pattern on the reverse strand

**User Prompt:**
> Find <ACG on both strands of AAACGT.

**Expected Tool Call:**
```json
{"tool": "find_pattern_hits", "arguments": {"sequence": "AAACGT", "pattern": "<ACG", "bothStrands": true}}
```

**Response:**
```json
{"patterns": {"pattern1": "<ACG"}, "hits": [{"pattern": "pattern1", "position": 3, "length": 3, "strand": "reverse", "mismatches": 0}]}
```

On the reverse strand the anchor swaps: the hit is the reverse complement CGT at the end of the sequence.

### Example 3: Pattern file

**User Prompt:**
> Search ACGT and TTGA (named p1, p2).

**Expected Tool Call:**
```json
{"tool": "find_pattern_hits", "arguments": {"sequence": "ACGTACGTTTGAACC", "patternFile": ">p1\nACGT\n>p2\nTTGA\n"}}
```

**Response:**
```json
{"patterns": {"p1": "ACGT", "p2": "TTGA"}, "hits": [{"pattern": "p1", "position": 0, "length": 4, "strand": "forward", "mismatches": 0}, {"pattern": "p1", "position": 4, "length": 4, "strand": "forward", "mismatches": 0}, {"pattern": "p2", "position": 8, "length": 4, "strand": "forward", "mismatches": 0}]}
```

## See Also

- [find_known_motifs](find_known_motifs.md) - motif sets (IUPAC, mismatches, strands)
