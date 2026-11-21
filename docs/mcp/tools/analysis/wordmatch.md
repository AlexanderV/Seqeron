# wordmatch

EMBOSS 6.6.0 wordmatch: every exact match of at least wordSize symbols between a set of pattern sequences (a FASTA text) and large sequences (a FASTA text), found by Rabin-Karp (260 nucleotide / protein runs: alignments in match, pair, simple, fasta, srspair, markx3, score and srs format, the statistics log and both GFF3 feature files byte for byte). The table holds every word of the pattern sequences (a pattern sequence is identified by its NAME; words around N, or X for a protein, are skipped as in seqmatchall); a window found in the table is verified and the match is extended to the right over the pattern sequence; a pattern sequence is not examined again before the end of its previous match in that scanned sequence. Returns the matches (0-based starts, per scanned sequence and pattern sequence in the program's order), the number of distinct words, the number and total length of the matches and the files of the program: the alignments (-outfile, alignFormat), the log (-logfile) and the feature files (-aoutfeat / -boutfeat). Distinct words of five or more symbols with equal hashes are not distinguished (the program picks an arbitrary one).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `wordmatch` |
| **Method ID** | `GenomicAnalyzer.Wordmatch` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 wordmatch: every exact match of at least wordSize symbols between a set of pattern sequences (a FASTA text) and large sequences (a FASTA text), found by Rabin-Karp (260 nucleotide / protein runs: alignments in match, pair, simple, fasta, srspair, markx3, score and srs format, the statistics log and both GFF3 feature files byte for byte). The table holds every word of the pattern sequences (a pattern sequence is identified by its NAME; words around N, or X for a protein, are skipped as in seqmatchall); a window found in the table is verified and the match is extended to the right over the pattern sequence; a pattern sequence is not examined again before the end of its previous match in that scanned sequence. Returns the matches (0-based starts, per scanned sequence and pattern sequence in the program's order), the number of distinct words, the number and total length of the matches and the files of the program: the alignments (-outfile, alignFormat), the log (-logfile) and the feature files (-aoutfeat / -boutfeat). Distinct words of five or more symbols with equal hashes are not distinguished (the program picks an arbitrary one).

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `patternFasta` | string | Yes | The pattern sequences as FASTA text (-asequence); every record is a pattern sequence. |
| `sequenceFasta` | string | Yes | The sequences to scan as FASTA text (-bsequence); every record is scanned. |
| `wordSize` | integer | No | Word size (wordmatch -wordsize, default 4, at least 2). (default: 4) |
| `dumpAlign` | boolean | No | Dump the matches as alignments (-dumpalign, default true). (default: true) |
| `dumpFeatures` | boolean | No | Dump the matches as GFF3 feature files (-dumpfeat, default true). (default: true) |
| `alignFormat` | string | No | The alignment format (-aformat, default match). (default: match) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `matches` | array | The matches (pattern name / start, scanned sequence name / start, length; 0-based). |
| `patternCount` | integer | Number of distinct words of the pattern set. |
| `allMatches` | integer | Number of matches. |
| `sumLength` | integer | Summed match length. |
| `alignments` | string | The -outfile text. |
| `log` | string | The -logfile statistics. |
| `featureFile1` | string | The -aoutfeat GFF3 text. |
| `featureFile2` | string | The -boutfeat GFF3 text. |

## Errors

| Code | Message |
|------|---------|
| InvalidParams | Invalid arguments (an empty or unreadable sequence, a word size below 2, a report format the program refuses for the sequence type). |

## Examples

### Example 1: Find the exact matches of a pattern set

**User Prompt:**
> Find the exact matches of a pattern set

**Expected Tool Call:**
```json
{"tool": "wordmatch", "arguments": {"patternFasta": ">p\nACGTACGTAC\n", "sequenceFasta": ">g\nTTACGTACGTACTT\n", "wordSize": 6}}
```

**Response:**
```json
{"matches": [{"pattern": "p", "patternStart": 0, "query": "g", "queryStart": 2, "length": 10}], "patternCount": 5, "allMatches": 1}
```

## See Also

- [seqmatchall](seqmatchall.md) - all-against-all word matches
- [megamerger](megamerger.md) - merge by word matches
