# find_tandem_repeats_emboss

EMBOSS etandem / equicktandem tandem-repeat programs (diagonal-score models), ported exactly.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `find_tandem_repeats_emboss` |
| **Method ID** | `GenomicAnalyzer.FindTandemRepeatsEtandem` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 tandem-repeat programs ported statement by statement (diagonal-score models, not the exact-run finder find_tandem_repeats). program etandem (default): for every unit size minRepeat..maxRepeat a running per-column consensus score (+1 when the new base is the column's most frequent base, -1 otherwise; N counts as a fifth base, or as a plain -1 with mismatch) extends candidates in steps of the unit size; candidates need a best score >= threshold and, unless uniform, a consensus with two different bases; overlapping candidates lose to the highest score (ties: smaller unit). Reports score, size, count, identity = 100 * (score + n * (copies + 1)) / (2 * n * copies) and the lower-case consensus (phased as printed). program equicktandem: the quick pre-screen - the sequence compared with itself shifted by every gap 1..maxRepeat with a +1 / -1 diagonal score; reported stretches are masked for larger gaps; no consensus / identity. Hits are in the program's report order (etandem: highest score first). Non-ACGT symbols count as N (etandem) / one 'other' class (equicktandem). position is 0-based (EMBOSS Start - 1). Verified against etandem and equicktandem 6.6.0 on 1700 + 1300 random and planted sequences. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program). etandem with threshold = -2147483648 (A4-4): the program wraps threshold - 1 to 2147483647, reports its empty root and dies with SIGFPE ("Floating point exception") as soon as one repeat qualifies - refused here with InvalidOperationException; a sequence without a qualifying repeat returns no repeats as the program prints nothing; -2147483647 and above are exact.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Emboss.Repeats/EmbossRepeatPrograms.EmbossTandem.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Emboss.Repeats/EmbossRepeatPrograms.EmbossTandem.cs)
- Algorithm doc: [Tandem_Repeat_Detection.md](../../../algorithms/Repeat_Analysis/Tandem_Repeat_Detection.md)
- Oracle: `scripts/analysis/generate_emboss_tandem_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Emboss.Repeats.Tests/TestData/Analysis/emboss_tandem_oracle.tsv`

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | Yes | Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive). |
| `program` | string | No | etandem (default) \| equicktandem (case-insensitive). |
| `minRepeat` | integer | No | etandem -minrepeat: smallest unit size (>= 2, default 10 as in the program). Ignored by equicktandem. (default: 10) |
| `maxRepeat` | integer | No | Largest unit size (-maxrepeat). Default: etandem = minRepeat (as in the program), equicktandem = 600. |
| `threshold` | integer | No | Minimum score (-threshold, default 20; equicktandem needs >= 1). (default: 20) |
| `mismatch` | boolean | No | etandem -mismatch: count N as a mismatch instead of a fifth base (default false). |
| `uniform` | boolean | No | etandem -uniform: also report single-base consensus repeats (default false). |
| `sequenceQualifiers` | string | No | EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: etandem ignores -sbegin / -send (the whole sequence is searched) and numbers a -sreverse sequence on the reverse strand; equicktandem searches the range and numbers from its first position. The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. -sformat selects the file format of the text; -ufo -fformat -fopenfile name the feature file of a program that takes the features of its sequences (diffseq: the files are given by the tool); -sask -squick -sdbname -sid -iquery are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused. |
| `reportFormat` | string | No | EMBOSS report format (-rformat): any of simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq (diffseq needs 7 report tags: refused as the program does; embl / genbank are nucleotide only, pir / swiss protein only), or the report qualifiers as typed after the files, e.g. '-rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2'. The report FILE of the program (header with the command line, one block per sequence, footer; #Rundate / #!Date is now; sequence name from sequenceName) is returned in the result field report. Omit for no report text. |
| `sequenceName` | string | No | Sequence name printed in the report (default sequence). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `program` | string | etandem or equicktandem |
| `repeats` | array | Repeats `{ position (0-based), length, score, repeatSize, repeatCount, identity?, consensus? }` in the program's report order; identity / consensus only for etandem |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | Invalid DNA sequence (A/C/G/T and IUPAC codes only) |
| 1001 | program must be etandem or equicktandem |
| 1002 | minRepeat must be 2 or higher / maxRepeat must be at least minRepeat / threshold must be at least 1 (equicktandem) |

## Examples

### Example 1: etandem

**User Prompt:**
> Run etandem (unit sizes 2-8, threshold 10) on a sequence with a CAG array.

**Expected Tool Call:**
```json
{"tool": "find_tandem_repeats_emboss", "arguments": {"sequence": "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACACGTTAACGTTAACGTTAACGTTAACGTTAacgtNN", "program": "etandem", "minRepeat": 2, "maxRepeat": 8, "threshold": 10}}
```

**Response:**
```json
{"program": "etandem", "repeats": [{"position": 60, "length": 36, "score": 33, "repeatSize": 3, "repeatCount": 12, "identity": 100.0, "consensus": "cag"}, {"position": 136, "length": 36, "score": 26, "repeatSize": 6, "repeatCount": 6, "identity": 94.4, "consensus": "acgtta"}]}
```

etandem prints `33 61 96 3 12 100.0 cag` and `26 137 172 6 6 94.4 acgtta` (1-based Start / End).

### Example 2: equicktandem

**User Prompt:**
> Quick-screen (ACG)x8 for tandem repeats up to unit size 6.

**Expected Tool Call:**
```json
{"tool": "find_tandem_repeats_emboss", "arguments": {"sequence": "ACGACGACGACGACGACGACGACG", "program": "equicktandem", "maxRepeat": 6, "threshold": 10}}
```

**Response:**
```json
{"program": "equicktandem", "repeats": [{"position": 0, "length": 24, "score": 21, "repeatSize": 3, "repeatCount": 8}]}
```

equicktandem prints `21 1 24 3 8`.

## See Also

- [find_tandem_repeats](find_tandem_repeats.md) - the exact-run definition (maximal runs of a primitive unit)
- [find_approximate_tandem_repeats](find_approximate_tandem_repeats.md) - Tandem Repeats Finder (Benson 1999)
