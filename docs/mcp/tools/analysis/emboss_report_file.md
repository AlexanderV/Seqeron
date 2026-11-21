# emboss_report_file

The whole report FILE of the EMBOSS 6.6.0 programs fuzznuc, fuzzpro, fuzztran, etandem and equicktandem on a FASTA text of one or several sequences (etandem / equicktandem read the first), byte for byte as the installed program writes it: every -rformat (simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq; diffseq is refused: it needs 7 report tags; embl / genbank are nucleotide only, pir / swiss protein only), the header with the command line and the pattern / parameter text, one block per sequence, the totals footer, the report options -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq, the sequence qualifiers -sbegin -send -sreverse and the extra -origfile text of etandem / equicktandem (originalFile). reportQualifiers is the text after the files, e.g. '-rformat table -raccshow -rmaxseq 2'; the default format is seqtable (fuzznuc, fuzzpro) or table (fuzztran, etandem, equicktandem). Delegates to EmbossReportPrograms and the one EmbossReportWriter.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `emboss_report_file` |
| **Method ID** | `EmbossReportPrograms.Fuzznuc` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The whole report FILE of the EMBOSS 6.6.0 programs fuzznuc, fuzzpro, fuzztran, etandem and equicktandem on a FASTA text of one or several sequences (etandem / equicktandem read the first), byte for byte as the installed program writes it: every -rformat (simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq; diffseq is refused: it needs 7 report tags; embl / genbank are nucleotide only, pir / swiss protein only), the header with the command line and the pattern / parameter text, one block per sequence, the totals footer, the report options -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq, the sequence qualifiers -sbegin -send -sreverse and the extra -origfile text of etandem / equicktandem (originalFile). reportQualifiers is the text after the files, e.g. '-rformat table -raccshow -rmaxseq 2'; the default format is seqtable (fuzznuc, fuzzpro) or table (fuzztran, etandem, equicktandem). Delegates to EmbossReportPrograms and the one EmbossReportWriter.

## Core Documentation Reference

- Algorithm: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `program` | string | Yes | fuzznuc \| fuzzpro \| fuzztran \| etandem \| equicktandem. |
| `fasta` | string | Yes | The input file text: FASTA records (name = first word of the header, an accession-like second word is the accession, the rest the description). |
| `pattern` | string | No | The pattern (fuzznuc / fuzzpro / fuzztran -pattern). |
| `maxMismatches` | integer | No | fuzznuc -pmismatch. (default: 0) |
| `complement` | boolean | No | fuzznuc -complement. (default: false) |
| `frame` | string | No | fuzztran -frame: 1 2 3 -1 -2 -3 F R 6 (default 1). |
| `table` | integer | No | fuzztran -table (0 = standard). (default: 0) |
| `minRepeat` | integer | No | etandem -minrepeat (default 10). |
| `maxRepeat` | integer | No | etandem / equicktandem -maxrepeat (default: etandem = minRepeat, equicktandem 600). |
| `threshold` | integer | No | etandem / equicktandem -threshold (default 20). |
| `mismatch` | boolean | No | etandem -mismatch. (default: false) |
| `uniform` | boolean | No | etandem -uniform. (default: false) |
| `sequenceQualifiers` | string | No | Sequence qualifiers as typed after the files: -sbegin 3 -send -2 -sreverse -slower -supper. |
| `reportQualifiers` | string | No | Report qualifiers as typed after the files: -rformat table -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall 5 -rmaxseq 2 (negated booleans -norscoreshow). |
| `commandLine` | array<string> | No | The command-line entries printed in the report header, one per line (e.g. '[-sequence] s.fa', '-pattern ACGT'); default: synthesized from the parameters. |
| `runDate` | string | No | The run date of the header (ISO 8601); default now. |
| `fastaFile` | string | No | The input file name: sequence USAs read fasta::NAME:id (default sequence.fasta). |
| `reportFile` | string | No | The report file name printed in the header (default stdout). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `report` | string | The -outfile report text. |
| `originalFile` | string | The -origfile text of etandem / equicktandem (absent for the other programs). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Unknown report format 'name' |
| 1001 | Report format 'embl' not supported for protein sequences |
| 1001 | Report format 'diffseq' needs 7 tags, has only 2 |

## Examples

### Example 1: fuzznuc table report of two sequences

**User Prompt:**
> Write the fuzznuc table report of ACGT on two sequences.

**Expected Tool Call:**
```json
{"tool": "emboss_report_file", "arguments": {"program": "fuzznuc", "fasta": ">seq1 a test\nACGTACGTAAACGT\n>seq2\nTTACGTGG\n", "pattern": "ACGT", "reportQualifiers": "-rformat table", "commandLine": ["[-sequence] ex.fa", "-pattern ACGT", "-rformat table", "-auto", "-outfile ex.rep"], "runDate": "2023-11-14T22:13:20Z", "fastaFile": "ex.fa", "reportFile": "ex.rep"}}
```

**Response:**
```json
{"report": "########################################\n# Program: fuzznuc\n# Rundate: Tue 14 Nov 2023 22:13:20\n# Commandline: fuzznuc\n#    [-sequence] ex.fa\n#    -pattern ACGT\n#    -rformat table\n#    -auto\n#    -outfile ex.rep\n# Report_format: table\n# Report_file: ex.rep\n########################################\n\n#=======================================\n#\n# Sequence: seq1     from: 1   to: 14\n# HitCount: 3\n#\n# Pattern_name Mismatch Pattern\n# pattern             0 ACGT\n# \n# Complement: No\n#\n#=======================================\n\n  Start     End  Strand   Score Pattern      Mismatch\n      1       4       +       4 pattern:ACGT        .\n      5       8       +       4 pattern:ACGT        .\n     11      14       +       4 pattern:ACGT        .\n\n#---------------------------------------\n#---------------------------------------\n#=======================================\n#\n# Sequence: seq2     from: 1   to: 8\n# HitCount: 1\n#\n# Pattern_name Mismatch Pattern\n# pattern             0 ACGT\n# \n# Complement: No\n#\n#=======================================\n\n  Start     End  Strand   Score Pattern      Mismatch\n      3       6       +       4 pattern:ACGT        .\n\n#---------------------------------------\n#---------------------------------------\n\n#---------------------------------------\n# Total_sequences: 2\n# Total_length: 22\n# Reported_sequences: 2\n# Reported_hitcount: 4\n#---------------------------------------\n"}
```

The text is the output of the installed fuzznuc 6.6.0 (run date faked).

## See Also

- [find_pattern_hits](find_pattern_hits.md) - fuzznuc hits
- [fuzzpro](fuzzpro.md) - fuzzpro hits
