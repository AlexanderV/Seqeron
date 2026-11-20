# find_open_reading_frames

Open reading frames (ORFs) in all six reading frames.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_open_reading_frames` |
| **Method ID** | `GenomicAnalyzer.FindOpenReadingFrames` |
| **Version** | 1.2.0 |
| **Stability** | Stable |

## Description

Finds **open reading frames** in all six frames (3 forward + 3 reverse-complement).
Each ORF runs from a start codon (default `ATG`) to its first in-frame stop codon,
inclusive, and is reported only if its nucleotide length (stop included) is at least
`minLength`. The stop codons come from the NCBI genetic code `genetic_code` (default 1 =
Standard: `TAA`/`TAG`/`TGA`; e.g. 2 = vertebrate mitochondrial: `TAA`/`TAG`/`AGA`/`AGG`).
Every start is considered independently, so nested ORFs sharing a stop are all reported
(canonical Rosalind ORF semantics). The reported sequence length is divisible by 3.

Options (NCBI ORFfinder / `COrf` semantics, all additive; the defaults reproduce version 1.0):

| Option | Effect |
|--------|--------|
| `startCodons` | `AtgOnly` (default), `GeneticCodeStarts` (ATG plus the table's alternative initiation codons, e.g. TTG/CTG in table 1), `AnySenseCodon` (stop-to-stop regions) |
| `nested` | `true` (default): every in-frame start; `false`: only the most upstream start per stop (NCBI ORFfinder default, EMBOSS getorf) |
| `includeOpenEnded` | also report 3'-partial ORFs that reach the sequence end without a stop (`isTerminated` = false; end = last complete codon) |
| `minLengthExcludesStop` | compare `minLength` with the length without the stop codon (NCBI ORFfinder `-ml`, EMBOSS getorf `-minsize`) |
| `bothStrands` | `false`: forward strand only |
| `strand` | NCBI ORFfinder `-strand`: `both` \| `plus` \| `minus` (case-insensitive); alternative to `bothStrands = false` (not combinable with it) |
| `ignoreNested` | ORFfinder `-n true`: drop every ORF whose span lies completely inside another reported ORF of any frame or strand (non-strict: a shared start or end counts; identical spans: the plus-strand ORF is dropped) |
| `regionStart`, `regionEnd` | ORFfinder `-b/-e`: 1-based inclusive search region (`regionStart` 0 = 1, `regionEnd` 0 or beyond the sequence = to the end; `regionStart` > `regionEnd` or negative values are errors). The region is cut out and scanned as a sequence of its own, `start`/`end` stay in input coordinates, `position` is relative to the region |
| `orfFinderMinLength` | ORFfinder `-ml` rule: `minLength` below 30 becomes 30, then whole codons (rounded down to a multiple of 3); lifts the `minLength` ≥ 1 check |
| `allowDualCodingStopTables` | accept `genetic_code` 27, 28, 31 and let the codons the table lists as stops end ORFs (what ORFfinder does) |

NCBI ORFfinder 0.4.3 equivalent of `ORFfinder -s 1 -ml 75`: `startCodons=GeneticCodeStarts nested=false includeOpenEnded=true
minLengthExcludesStop=true orfFinderMinLength=true allowDualCodingStopTables=true minLength=75` (add `ignoreNested`, `strand`, `regionStart/End`,
`genetic_code` for `-n`, `-strand`, `-b/-e`, `-g`). These options were measured on the real ORFfinder 0.4.3 binary (1750 runs, 54 820 results, 232 162 ORFs; the
tool equals it on all of them, see the algorithm document §5.5). ORFfinder's `-c` (circular) has no effect in the binary and no counterpart here.

Without `allowDualCodingStopTables`, tables with dual-coding stop codons (27, 28, 31) are rejected: every stop codon of those tables also codes
an amino acid, so no codon unambiguously ends an ORF (as `Translator.FindOrfs`).

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs) `FindOpenReadingFrames` (scan core: [OrfScanner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/OrfScanner.cs))
- Algorithm: [Open_Reading_Frame_Detection.md](../../../algorithms/Analysis/Open_Reading_Frame_Detection.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive). An ambiguous codon ends an ORF only when every concrete codon is a stop, is never a start, and is otherwise read through. |
| `minLength` | integer | No | Minimum ORF length in nucleotides, stop codon included unless minLengthExcludesStop (default 100). |
| `genetic_code` | integer | No | NCBI translation table number whose stop codons end ORFs (default 1 = Standard; e.g. 2 = vertebrate mitochondrial, 11 = bacterial/plastid). One of the NCBI gc.prt tables (1-6, 9-16, 21-33); the dual-coding 27, 28, 31 only with allowDualCodingStopTables. |
| `startCodons` | string | No | Start-codon policy (case-insensitive): AtgOnly (default), GeneticCodeStarts (ATG and the alternative initiation codons of genetic_code, e.g. TTG/CTG in table 1), AnySenseCodon (any non-stop codon). |
| `nested` | boolean | No | true (default): report every in-frame start (nested ORFs). false: only the most upstream start per stop codon (NCBI ORFfinder default, EMBOSS getorf). |
| `includeOpenEnded` | boolean | No | true: also report ORFs that run to the sequence end without a stop codon (3'-partial, isTerminated=false). Default false. (default: false) |
| `minLengthExcludesStop` | boolean | No | true: minLength excludes the stop codon (NCBI ORFfinder -ml / EMBOSS getorf -minsize). Default false. (default: false) |
| `bothStrands` | boolean | No | true (default): scan both strands; false: forward strand only (see strand). |
| `strand` | string | No | NCBI ORFfinder -strand: both \| plus \| minus (case-insensitive). Default null = decided by bothStrands; cannot be combined with bothStrands=false. |
| `ignoreNested` | boolean | No | NCBI ORFfinder -n true: drop every ORF whose span lies completely inside another reported ORF of any frame/strand (identical spans: the plus-strand ORF is dropped). Default false. (default: false) |
| `regionStart` | integer | No | NCBI ORFfinder -b: 1-based first nucleotide of the searched region (default 1; 0 is read as 1). The region is cut out and scanned as a sequence of its own. |
| `regionEnd` | integer | No | NCBI ORFfinder -e: 1-based last nucleotide (inclusive) of the region; 0 (default) or beyond the sequence = to the end. An end below the start is an error. |
| `orfFinderMinLength` | boolean | No | NCBI ORFfinder -ml rule: minLength below 30 is changed to 30 and used in whole codons (rounded down to a multiple of 3; also lifts the minLength >= 1 check). Default false. (default: false) |
| `allowDualCodingStopTables` | boolean | No | true: genetic_code 27, 28 and 31 are accepted and the codons the table lists as stops end ORFs (what NCBI ORFfinder -g 27/28/31 does). Default false (rejected). (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | ORFs: `{ sequence, position, frame, isReverseComplement, length, codonCount, start, end, isTerminated, isStartPartial }` |

`frame` is 1–3 within the scanned strand; `position` is the 0-based start in the *scanned strand*
(the input, or its reverse complement when `isReverseComplement`); `start`/`end` are the 0-based inclusive
coordinates of the same span on the input strand (stop codon included; for a reverse-strand ORF the
reverse complement of `input[start..end]` is `sequence`); `codonCount` = length / 3. `isStartPartial` is true when the ORF begins within the first
3 nt of the scanned strand or region (ORFfinder "5'-partial", also for an ORF that begins with ATG); `isTerminated` = false is 3'-partial.

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1001 | Unknown / dual-coding genetic code, unknown `startCodons` |
| 1002 | Minimum ORF length must be at least 1 (not with `orfFinderMinLength`); negative region value or `regionStart` > `regionEnd` |
| 1001 | `strand` unknown, or combined with `bothStrands = false` |

## Examples

### Example 1: Minimal forward ORF (ATG-TAA)

**User Prompt:**
> Find ORFs of length ≥ 6 in "ATGTAA".

**Expected Tool Call:**
```json
{
  "tool": "find_open_reading_frames",
  "arguments": { "sequence": "ATGTAA", "minLength": 6 }
}
```

**Response:**
```json
{ "items": [ { "sequence": "ATGTAA", "position": 0, "frame": 1, "isReverseComplement": false, "length": 6, "codonCount": 2, "start": 0, "end": 5, "isTerminated": true } ] }
```

### Example 2: One extra codon (ATG-AAA-TAG)

**User Prompt:**
> ORFs of length ≥ 9 in "ATGAAATAG".

**Expected Tool Call:**
```json
{
  "tool": "find_open_reading_frames",
  "arguments": { "sequence": "ATGAAATAG", "minLength": 9 }
}
```

**Response:**
```json
{ "items": [ { "sequence": "ATGAAATAG", "position": 0, "frame": 1, "isReverseComplement": false, "length": 9, "codonCount": 3, "start": 0, "end": 8, "isTerminated": true } ] }
```

### Example 3: Reverse-strand ORF in input coordinates

**Expected Tool Call:**
```json
{
  "tool": "find_open_reading_frames",
  "arguments": { "sequence": "CCTTAGGGCATGG", "minLength": 9 }
}
```

**Response:** the reverse complement `CCATGCCCTAAGG` holds `ATGCCCTAA` at scanned-strand position 2 (frame 3); on the input it spans 2..10.
```json
{ "items": [ { "sequence": "ATGCCCTAA", "position": 2, "frame": 3, "isReverseComplement": true, "length": 9, "codonCount": 3, "start": 2, "end": 10, "isTerminated": true } ] }
```

### Example 4: Vertebrate mitochondrial code (TGA = Trp)

```json
{ "tool": "find_open_reading_frames", "arguments": { "sequence": "ATGTGAAAATAA", "minLength": 1, "genetic_code": 2 } }
```
Returns the 12-nt ORF `ATGTGAAAATAA` (with the Standard code the ORF is `ATGTGA`).

### Example 5: ORFfinder `-n true` (identical spans across strands)

```json
{ "tool": "find_open_reading_frames", "arguments": { "sequence": "TAATTAGCTAAAGCTGTCCACCATGCTGTCCTTCCAGGTAGCTAATTA", "minLength": 30,
  "startCodons": "AnySenseCodon", "nested": false, "includeOpenEnded": true, "minLengthExcludesStop": true, "orfFinderMinLength": true, "ignoreNested": true } }
```
Without `ignoreNested` the ORFs are `+3..44, +7..45, -3..44, -2..40, -7..45` (input spans); with it the plus ORFs are dropped (`+3..44` equals `-3..44`, `+7..45` equals `-7..45`):
`-3..44, -2..40, -7..45` (measured on ORFfinder 0.4.3, sequence `strand_equal_span_s2_12` of the oracle).

## Performance

- **Time Complexity:** O(n + output) — one pass per frame; pending starts are released at each stop (5 Mbp random DNA, all six frames: about a second in Release).
- **Space Complexity:** O(number of ORFs + pending starts of the current stop-free stretch).

## See Also

- [find_motif](find_motif.md) — locate a known motif
- [codon_frequencies](codon_frequencies.md) — codon usage in a reading frame
