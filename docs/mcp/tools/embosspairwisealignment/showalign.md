# showalign

EMBOSS showalign: pretty display of a multiple alignment.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `showalign` |
| **Method ID** | `MultipleAlignmentDisplay.ShowAlign` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 showalign: the pretty text display of a multiple alignment, identical to the program (600 random alignments with random options equal). Every row except the reference is converted against the reference (the consensus unless refseq names a row or a 1-based number): show N (default) non-identities, I identities, S similarities (matrix score > 0), D dissimilarities, A all; '.' marks suppressed residues; similarCase (default true) encodes similarity in the case (non-identity / similarity / all as in the program); order I (input), A (names), S (decreasing similarity to the reference). Per block of width columns (default 60): the position numbers (every 10th, relative to the reference residues when a row is the reference, with -n / +n flanks), the tick ruler, the reference first when it is a row, the rows, the consensus (consensus = false hides it; bottom = false puts it first when it is the reference). The consensus is the EMBOSS cons algorithm (plurality in percent, setCase, identity in percent, gaps: a gap can be the consensus) on the alignment as given (case kept); a lower-case letter of the input is kept lower case. margin = width of the name column (0 = no names, -1 = longest name + 1); uppercase = 1-based column ranges put in upper case, e.g. '3-10,20-25'. The set type is that of the first row (nucleotide: A/C/G/T/U + IUPAC letters); matrix defaults to EDNAFULL / EBLOSUM62; shorter rows are padded with trailing gaps. html = true wraps the text in <pre> ... </pre> and highlight (EMBOSS range syntax, a colour after each pair) colours the ranges of every sequence and consensus line with <font color=...> (A3-13c2: 300 random texts equal). Rows are read as EMBOSS reads an alignment (A4-3): letters and '*' are kept, '.' '~' '-' are gaps, '?' '#' '+' are the unknown residue, every other character (digits, blanks, '_', punctuation) is dropped, and rows that become unequal are padded with trailing gaps. A display text above 134 217 728 characters (blocks x lines x (margin + width + 1): a huge margin, width 1 over a very long alignment) is refused with an ArgumentException (MultipleAlignmentDisplay.ShowAlignMaxOutputCharacters); a width beyond the alignment shows one block. sbegin / send (A4-9, 400 random ranges equal the program) cut the displayed columns: the consensus is that of the whole alignment and the numbers and ticks keep the positions of the whole alignment.

## Core Documentation Reference

- Source: [MultipleAlignmentDisplay.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.PairwiseAlignment/MultipleAlignmentDisplay.cs)
- Algorithm doc: [Alignment_Statistics.md](../../../algorithms/Alignment/Alignment_Statistics.md)
- Oracle: [generate_showalign_oracle.py](../../../../scripts/analysis/generate_showalign_oracle.py), [generate_showalign_range_oracle.py](../../../../scripts/analysis/generate_showalign_range_oracle.py) (-sbegin / -send)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `names` | array<string> | Yes | Sequence names (one per row; the first word is used). |
| `sequences` | array<string> | Yes | The aligned sequences (DNA or protein; gaps '-', '.', '~'; case is kept). |
| `refseq` | string | No | Reference: '0' = the consensus (default), a 1-based row number or a sequence name. |
| `show` | string | No | What to show: N non-identities (default), I identities, S similarities, D dissimilarities, A all. |
| `order` | string | No | Output order of the rows: I input (default), A alphabetical names, S similarity to the reference. |
| `matrix` | string | No | Scoring matrix name (EBLOSUM62, EDNAFULL, EPAM250, ...); default by sequence type. |
| `similarCase` | boolean | No | Case encodes similarity (default true). |
| `consensus` | boolean | No | Show the consensus line (default true). |
| `bottom` | boolean | No | Consensus at the bottom (default true) instead of the top. |
| `number` | boolean | No | Position number line (default true). |
| `ruler` | boolean | No | Tick ruler line (default true). |
| `width` | integer | No | Columns per block (default 60). |
| `margin` | integer | No | Name column width: -1 = longest name + 1 (default), 0 = no names. |
| `uppercase` | string | No | 1-based inclusive column ranges to put in upper case, e.g. '3-10,20-25' (EMBOSS range syntax; '@file' = one range per line). |
| `plurality` | number | No | Consensus plurality in percent of the total weight (default 50). |
| `setCase` | number | No | Consensus: weight at or below which the letter is lower case (default half the number of sequences). |
| `identity` | number | No | Consensus: percent of identities required at a position (default 0). |
| `gaps` | boolean | No | A gap character can be the consensus (default true). |
| `sbegin` | integer | No | -sbegin: first column shown, 1-based (0 = the first; negative = counted from the end; beyond the alignment = the last). (default: 0) |
| `send` | integer | No | -send: last column shown, 1-based (0 = the last; negative = counted from the end; before the first shown column = that column). (default: 0) |
| `html` | boolean | No | -html (default false): wrap the display in <pre> ... </pre> and colour the highlight ranges with <font color=...> tags. |
| `highlight` | string | No | -highlight: column ranges to colour with html, each followed by an HTML colour (an empty colour is red), e.g. '24-45 blue 56-78 orange'; '@file' = one range per line. Used only with html. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The display text (blocks of number line, ruler, reference, rows and consensus; lines end in a line feed). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Names and sequences are required |
| 1001 | One name per sequence is required |
| 1001 | Invalid uppercase range |
| 1001 | show must be A, I, N, S or D |
| 1001 | order must be I, A or S |
| 1001 | Reference sequence is not a sequence name or a number |

## Examples

### Example 1: Default display of a three-row nucleotide alignment

**User Prompt:**
> Show the alignment of s1 acgtac-gtaGT, s2 acgtaccgtaGA and s3 acctaccgtTGA.

**Expected Tool Call:**
```json
{"tool": "showalign", "arguments": {"names": ["s1", "s2", "s3"], "sequences": ["acgtac-gtaGT", "acgtaccgtaGA", "acctaccgtTGA"]}}
```

**Response:**
```json
{"text": "                  10  \n          ----:----|--\ns1        ......-....T\ns2        ............\ns3        ..C......T..\nConsensus acgtaccgtaGA\n\n"}
```

showalign prints the same text: identities with the consensus are dots, the position numbers are every 10th column.

## See Also

- [consambig](consambig.md) - ambiguity-code consensus of the same alignment
- [generate_emboss_consensus](generate_emboss_consensus.md) - the plurality consensus used as reference
- [multiple_alignment_statistics](multiple_alignment_statistics.md) - alignment statistics (infoalign)
