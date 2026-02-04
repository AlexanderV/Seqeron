# FASTA Parsing

## Overview

FASTA format is a text-based format for representing nucleotide or amino acid sequences. It originated from the FASTA software package (Lipman & Pearson, 1985) and has become a near-universal standard in bioinformatics.

## Format Specification

### Structure

```
>ID description (optional)
SEQUENCE_LINE_1
SEQUENCE_LINE_2
...
```

### Components

1. **Header Line (defline)**
   - Starts with `>` (greater-than symbol)
   - First word: sequence identifier
   - Remainder after whitespace: optional description
   - Single line only

2. **Sequence Data**
   - One letter per nucleotide/amino acid
   - Typically ≤80 characters per line
   - May span multiple lines
   - Invalid characters ignored (spaces, numbers, etc.)

### Multi-FASTA

Multiple sequences concatenated, each beginning with `>`:

```
>seq1 First sequence
ATGCATGC
>seq2 Second sequence  
GGCCTTAA
```

## Algorithm

### Parsing Algorithm

```
Time Complexity: O(n) where n = total characters
Space Complexity: O(m) where m = longest sequence
```

**Steps:**
1. Read lines sequentially
2. If line starts with `>`:
   - Yield previous entry (if any)
   - Extract ID (first token) and description (remainder)
   - Initialize new sequence buffer
3. Else if non-empty line:
   - Append trimmed content to sequence buffer
4. At end: yield final entry

### Output Formatting Algorithm

```
Time Complexity: O(n) where n = total sequence length
```

**Steps:**
1. Write `>` + header
2. Split sequence into chunks of `lineWidth` characters
3. Write each chunk as separate line

## Implementation

### Class: `FastaParser`

**Namespace:** `Seqeron.Genomics.IO`

### Methods

| Method | Description | Complexity |
|--------|-------------|------------|
| `Parse(string)` | Parse FASTA string | O(n) |
| `ParseFile(string)` | Parse FASTA file (sync) | O(n) |
| `ParseFileAsync(string)` | Parse FASTA file (async) | O(n) |
| `ToFasta(entries, lineWidth)` | Format to FASTA string | O(n) |
| `WriteFile(path, entries)` | Write FASTA file | O(n) |

### Return Type: `FastaEntry`

| Property | Type | Description |
|----------|------|-------------|
| `Id` | `string` | Sequence identifier |
| `Description` | `string?` | Optional description |
| `Sequence` | `DnaSequence` | The sequence data |
| `Header` | `string` | Full header (Id + Description) |

## Edge Cases

| Case | Behavior |
|------|----------|
| Empty input | Returns empty enumerable |
| Header only (no sequence) | Not yielded (requires sequence) |
| Multi-line sequence | Lines concatenated |
| Whitespace in sequence | Trimmed/ignored |
| No description | `Description` is null |
| Long sequences | Wrapped at `lineWidth` (default 80) |

## References

1. Lipman DJ, Pearson WR (1985). "Rapid and sensitive protein similarity searches". Science. 227(4693):1435–41
2. Pearson WR, Lipman DJ (1988). "Improved tools for biological sequence comparison". PNAS. 85(8):2444–8
3. Wikipedia: FASTA format - https://en.wikipedia.org/wiki/FASTA_format
4. NCBI BLAST Help - https://blast.ncbi.nlm.nih.gov/doc/blast-topics/
