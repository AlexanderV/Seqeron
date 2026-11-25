# megamerger

EMBOSS megamerger: merge two large overlapping DNA sequences by word-match anchors.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `megamerger` |
| **Method ID** | `GenomicAnalyzer.MegaMerge` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 megamerger: merges two large overlapping DNA sequences by exact word matches (600 random related / unrelated pairs equal: merged sequence and report text). All exact matches of wordSize symbols (the maximal matches of seqmatchall) are ordered longest first and reduced to a minimal non-overlapping set: the longest is kept, matches inside its dead zone are dropped, those crossing it are truncated (dropped when shorter than the word), until none is left; the kept matches are ordered by position. The merged sequence = the part before the first match (from the sequence whose overlap starts at its own start; when neither does the longer prefix, upper case), the matches (first sequence, lower case), between matches the region of the sequence whose mismatch is farther from its ends (preferFirst: always the first sequence), upper case, and the part after the last match (the longer remainder in upper case when neither overlap is at an end). The report lists every overlap start / end, mismatch region and choice as the program writes it. A second sequence shorter than the word or without a match gives an empty merge; a first sequence shorter than the word (or without a word free of N) is an error (the program stops with 'No match found'). Note: this merges by anchors; the overlap merge of two short sequences by dynamic programming is the merger tool. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program). More than 4 194 304 maximal word matches of the pair are refused with an ArgumentException (SequenceAligner.WordMatchMaxMatches).

## Core Documentation Reference

- Source: [GenomicAnalyzer.MegaMerger.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.PairwiseAlignment/EmbossPairwisePrograms.MegaMerger.cs)
- Algorithm doc: [Common_Region_Detection.md](../../../algorithms/Sequence_Comparison/Common_Region_Detection.md)
- Oracle: [generate_megamerger_oracle.py](../../../../scripts/analysis/generate_megamerger_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | Yes | First DNA sequence (letters, any case). |
| `sequence2` | string | Yes | Second DNA sequence. |
| `wordSize` | integer | No | Word size (megamerger -wordsize, default 20, at least 2). (default: 20) |
| `preferFirst` | boolean | No | Use the first sequence over every mismatch region (megamerger -prefer, default false). (default: false) |
| `name1` | string | No | Name of the first sequence in the report. (default: seq1) |
| `name2` | string | No | Name of the second sequence in the report. (default: seq2) |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 and the same with 2 for the second sequence. Measured on EMBOSS 6.6.0: megamerger merges the TRIMMED sequences (positions of the report are those of the trimmed texts; -sreverse1 reverse-complements before the range is cut); -sprotein is refused (the sequences are DNA). Omit for the whole sequences. The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. -sformat selects the file format of the text; -ufo -fformat -fopenfile name the feature file of a program that takes the features of its sequences (diffseq: the files are given by the tool); -sask -squick -sdbname -sid -iquery are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `mergedSequence` | string | The merged sequence (empty when nothing was merged). |
| `report` | string | The text of the megamerger report file. |
| `matches` | array | The minimal set of anchoring matches `{ sequence1, sequence2, position1, position2, length }` (0-based starts) ordered by position. |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | The first sequence has no word of N symbols: megamerger stops with No match found |

## Examples

### Example 1: Overlap at the start, the second sequence longer

**User Prompt:**
> Merge ACGTACGTAC and ACGTACGTACGGGGGG with word size 4.

**Expected Tool Call:**
```json
{"tool": "megamerger", "arguments": {"sequence1": "ACGTACGTAC", "sequence2": "ACGTACGTACGGGGGG", "wordSize": 4, "name1": "a", "name2": "b"}}
```

**Response:**
```json
{"mergedSequence": "acgtacgtacgggggg", "report": "# Report of megamerger of: a and b\n\na overlap starts at 1\nb overlap starts at 1\n\nMatching region a 1-10 : b 1-10\nLength of match: 10\n\na overlap ends at 10\nb overlap ends at 10\n\nUsing b 11-16 as the final sequence\n\n", "matches": [{"sequence1": 0, "sequence2": 1, "position1": 0, "position2": 0, "length": 10}]}
```

megamerger writes the same merged sequence and report.

## See Also

- [merger](../alignment/merger.md) - dynamic-programming merge of two overlapping sequences
- [seqmatchall](seqmatchall.md) - all exact word matches
- [local_align_top_n](../alignment/local_align_top_n.md) - non-intersecting local alignments
