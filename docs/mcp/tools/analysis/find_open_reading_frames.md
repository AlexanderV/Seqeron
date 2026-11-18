# find_open_reading_frames

Open reading frames (ORFs) in all six reading frames.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_open_reading_frames` |
| **Method ID** | `GenomicAnalyzer.FindOpenReadingFrames` |
| **Version** | 1.1.0 |
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

Tables with dual-coding stop codons (27, 28, 31) are rejected: every stop codon of those tables also codes
an amino acid, so no codon unambiguously ends an ORF (as `Translator.FindOrfs`).

## Core Documentation Reference

- Source: [GenomicAnalyzer.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs) `FindOpenReadingFrames` (scan core: [OrfScanner.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/OrfScanner.cs))
- Algorithm: [Open_Reading_Frame_Detection.md](../../../algorithms/Analysis/Open_Reading_Frame_Detection.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | DNA sequence (A/C/G/T, min length 1) |
| `minLength` | integer | No | Minimum ORF length in nucleotides, stop codon included unless `minLengthExcludesStop` (default 100, ≥ 1) |
| `genetic_code` | integer | No | NCBI translation table whose stop codons end ORFs (default 1; 27/28/31 rejected) |
| `startCodons` | string | No | `AtgOnly` (default) \| `GeneticCodeStarts` \| `AnySenseCodon` (case-insensitive) |
| `nested` | boolean | No | Report every in-frame start (default true) |
| `includeOpenEnded` | boolean | No | Also report ORFs without a stop (default false) |
| `minLengthExcludesStop` | boolean | No | NCBI ORFfinder length convention (default false) |
| `bothStrands` | boolean | No | Scan the reverse complement too (default true) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `items` | array | ORFs: `{ sequence, position, frame, isReverseComplement, length, codonCount, start, end, isTerminated }` |

`frame` is 1–3 within the scanned strand; `position` is the 0-based start in the *scanned strand*
(the input, or its reverse complement when `isReverseComplement`); `start`/`end` are the 0-based inclusive
coordinates of the same span on the input strand (stop codon included; for a reverse-strand ORF the
reverse complement of `input[start..end]` is `sequence`); `codonCount` = length / 3.

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence |
| 1001 | Unknown / dual-coding genetic code, unknown `startCodons` |
| 1002 | Minimum ORF length must be at least 1 |

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

## Performance

- **Time Complexity:** O(n + output) — one pass per frame; pending starts are released at each stop (5 Mbp random DNA, all six frames: about a second in Release).
- **Space Complexity:** O(number of ORFs + pending starts of the current stop-free stretch).

## See Also

- [find_motif](find_motif.md) — locate a known motif
- [codon_frequencies](codon_frequencies.md) — codon usage in a reading frame
