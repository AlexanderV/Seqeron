# Evidence: PARSE-FASTA-001 - FASTA Parsing

## Test Unit
- **ID:** PARSE-FASTA-001
- **Area:** FileIO
- **Methods:** Parse, ParseFile, ParseFileAsync, ToFasta, WriteFile

## Authoritative Sources

### Primary Sources

1. **Wikipedia - FASTA Format**
   - URL: https://en.wikipedia.org/wiki/FASTA_format
   - Accessed: 2026-02-05
   
2. **NCBI BLAST Help - Query Input**
   - URL: https://blast.ncbi.nlm.nih.gov/doc/blast-topics/
   - Accessed: 2026-02-05

3. **NCBI FASTA Specification**
   - URL: https://www.ncbi.nlm.nih.gov/blast/fasta.shtml
   - Referenced by Wikipedia

### Original Papers
- Lipman DJ, Pearson WR (1985). "Rapid and sensitive protein similarity searches". Science. 227(4693):1435–41
- Pearson WR, Lipman DJ (1988). "Improved tools for biological sequence comparison". PNAS. 85(8):2444–8

## Format Specification (from Sources)

### Structure
1. **Description Line (defline/header)**
   - Begins with `>` (greater-than symbol)
   - First word is the sequence identifier
   - Rest (after space/tab) is optional description
   - All in a single line

2. **Sequence Lines**
   - Follow immediately after header
   - Recommended: ≤80 characters per line (Wikipedia: "typically no more than 80 characters")
   - Original format: no longer than 120 characters, usually ≤80
   - Can be multi-line (interleaved) or single-line (sequential)
   - Whitespace is ignored within sequence

### Character Sets

#### Nucleic Acid Codes (from NCBI/Wikipedia)
| Code | Meaning |
|------|---------|
| A | Adenine |
| C | Cytosine |
| G | Guanine |
| T | Thymine |
| U | Uracil (RNA) |
| N | Any (A/C/G/T/U) |
| R | Purine (A/G) |
| Y | Pyrimidine (C/T/U) |
| K | Keto (G/T/U) |
| M | Amino (A/C) |
| S | Strong (G/C) |
| W | Weak (A/T/U) |
| B | Not A (C/G/T/U) |
| D | Not C (A/G/T/U) |
| H | Not G (A/C/T/U) |
| V | Not T/U (A/C/G) |
| - | Gap |

### Multi-FASTA Format
- Multiple sequences in one file
- Each sequence starts with `>`
- Concatenation of single-FASTA entries

### Edge Cases (Documented)

1. **Empty Content**
   - Should return empty/no sequences

2. **Blank Lines**
   - "Blank lines are not allowed in the middle of FASTA input" (NCBI)
   - Common practice: skip blank lines in parsers

3. **Whitespace in Sequence**
   - "Anything other than a valid character would be ignored (including spaces, tabulators)" (Wikipedia)

4. **Case Sensitivity**
   - "lower-case letters are accepted and are mapped into upper-case" (NCBI/Wikipedia)

5. **Header without Sequence**
   - Not explicitly defined; implementation-specific

6. **Special Characters in Headers**
   - NCBI defines specific identifier formats (gi|, gb|, etc.)
   - Pipes (`|`) are part of NCBI identifier format

7. **Line Width for Output**
   - Historical: 80 characters (terminal width)
   - Modern: configurable, typically 60-80

## Test Categories from Evidence

### Must Test (Format Specification)
1. Single sequence parsing
2. Multi-sequence parsing  
3. Multi-line sequence concatenation
4. Header ID/description parsing
5. Empty input handling
6. ToFasta output formatting
7. Line width wrapping on output
8. Round-trip integrity

### Should Test (Documented Behavior)
1. Header without description
2. Whitespace/blank line handling in sequence
3. Special characters in headers (pipes, colons)
4. Very long sequences
5. File I/O operations

### Could Test (Edge Cases)
1. Unicode in descriptions
2. NCBI-format identifiers
3. Performance with large files
4. Async file reading

## Implementation Notes

Current implementation (`FastaParser.cs`):
- Uses `>` as header marker
- Splits header on space/tab to extract ID and description
- Ignores empty/whitespace-only lines in sequence
- Line width configurable (default 80)
- Returns `FastaEntry` with Id, Description, and DnaSequence
