# dnadiff

MUMmer 3.23 dnadiff: the comparative analysis of two sequence sets and its report.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `dnadiff` |
| **Method ID** | `GenomicAnalyzer.Dnadiff` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUMmer 3.23 dnadiff: the comparative analysis of two sequence sets, the perl script over nucmer --maxmatch, delta-filter -1 / -m, show-coords -rclTH, show-snps -rlTHC and show-diff -rH / -qH, with its report: sequences, aligned / unaligned sequences and bases, the 1-to-1 and M-to-M alignments (count, total and average length, average identity), breakpoints, relocations, translocations, inversions, insertions, tandem insertions (count, sum, average) and the SNPs / indels with the 12 substitution classes, all also for the 'good' ones (buffer >= 20). Returns every file of the program: report, delta, 1delta, mdelta, 1coords, mcoords, snps, rdiff, qdiff and (when something is unaligned) unref / unqry. Pass delta (dnadiff -d) to analyse a precomputed nucmer delta file instead of running nucmer. The order of the per-substitution lines inside the SNP / indel blocks of the report and of the unref / unqry lines is the order of perl's hash iteration, which differs from run to run in the program (6 runs of the installed dnadiff on one input gave 6 different reports); here it is fixed (. A C G T, then other symbols). Everything else is byte for byte the program's (300 + runs, the unordered lines compared as sets). A FASTA line with white space (carriage return) is the program's 'Whitespace found in FastA' error.

## Core Documentation Reference

- Source: [GenomicAnalyzer.Dnadiff.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.Dnadiff.cs)
- Algorithm doc: [Nucmer_And_Delta_Programs.md](../../../algorithms/Alignment/Nucmer_And_Delta_Programs.md)
- Oracle: [generate_mummer_show_oracle.py](../../../../scripts/alignment/generate_mummer_show_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference` | string | Yes | Reference FASTA text. |
| `query` | string | Yes | Query FASTA text. |
| `delta` | string | No | Precomputed delta file text (dnadiff -d); null = run nucmer --maxmatch. |
| `referenceName` | string | No | Reference file name printed in the delta header (default reference.fa). |
| `queryName` | string | No | Query file name printed in the delta header (default query.fa). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `report` | string | The .report file. |
| `delta` | string | The nucmer delta file. |
| `delta1` | string | The 1-to-1 delta file (.1delta). |
| `deltaM` | string | The M-to-M delta file (.mdelta). |
| `coords1` | string | The 1-to-1 coordinates (.1coords). |
| `coordsM` | string | The M-to-M coordinates (.mcoords). |
| `snps` | string | The SNPs (.snps). |
| `referenceDiff` | string | The reference breakpoints (.rdiff). |
| `queryDiff` | string | The query breakpoints (.qdiff). |
| `unalignedReference` | string | The .unref file (null when every reference sequence aligned). |
| `unalignedQuery` | string | The .unqry file (null when every query sequence aligned). |

## Errors

| Code | Message |
|------|---------|
| 1001 | Sequence cannot be null or empty |
| 1001 | ERROR: Whitespace found in FastA, aborting. |

## Examples

### Example 1: Two contigs against one reference

**User Prompt:**
> Run dnadiff on the reference and the two query contigs.

**Expected Tool Call:**
```json
{"tool": "dnadiff", "arguments": {"reference": ">r\nTTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATATACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTATTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAAATCCCACCCGAACTAAGTTTGTCGAACCTTGGTCCAAGATCGGGACTCGGTCTCCAGGTAAGACGGGCTCATTCATAAACGTTACTAAGGGGTATAATCTTCTATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAGTCTTGTGCGC\n", "query": ">q1\nTTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGA\n>q2\nGCGCACAAGACTGCTGTAATTGGATTGCAAGTCTACTAAGTGTTCCCACCCACAAATAGAAGATTATACCCCTTAGTAACGTTTATGAATGAGCCCGTCTTACCTGGAGACCGAGTCCCGATCTTGGACCAAGGTTCGACAAACTTAGTTCGGGTGGGATTTATGCTCCGTTAGGAATACTGGGTGGGCCGCGGGGGGAATAGCTATAGCCCTGGGTTGCTAGTCAACGAAGCAAATGTATATTATTAGACGCGACCGCAGTAAGAGCGTATGTACACATTCTCCCTAGGTTGCACATGAAGAATGAGCCGGTATGCAGG\n", "delta": "ref2.fa qry2.fa\nNUCMER\n>r q1 600 320\n1 320 1 320 0 0 0\n0\n>r q2 600 320\n281 600 320 1 0 0 0\n0\n"}}
```

**Response:**
```json
{"report": "ref2.fa qry2.fa\nNUCMER\n\n                               [REF]                [QRY]\n[Sequences]\nTotalSeqs                          1                    2\nAlignedSeqs               1(100.00%)           2(100.00%)\nUnalignedSeqs               0(0.00%)             0(0.00%)\n\n[Bases]\nTotalBases                       600                  640\nAlignedBases            600(100.00%)         640(100.00%)\nUnalignedBases              0(0.00%)             0(0.00%)\n\n[Alignments]\n1-to-1                             2                    2\nTotalLength                      640                  640\nAvgLength                     320.00               320.00\nAvgIdentity                   100.00               100.00\n\nM-to-M                             2                    2\nTotalLength                      640                  640\nAvgLength                     320.00               320.00\nAvgIdentity                   100.00               100.00\n\n[Feature Estimates]\nBreakpoints                        2                    0\nRelocations                        0                    0\nTranslocations                     1                    0\nInversions                         0                    0\n\nInsertions                         0                    0\nInsertionSum                       0                    0\nInsertionAvg                    0.00                 0.00\n\nTandemIns                          0                    0\nTandemInsSum                       0                    0\nTandemInsAvg                    0.00                 0.00\n\n[SNPs]\nTotalSNPs                          0                    0\nAC                          0(0.00%)             0(0.00%)\nAT                          0(0.00%)             0(0.00%)\nAG                          0(0.00%)             0(0.00%)\nTA                          0(0.00%)             0(0.00%)\nTG                          0(0.00%)             0(0.00%)\nTC                          0(0.00%)             0(0.00%)\nGT                          0(0.00%)             0(0.00%)\nGC                          0(0.00%)             0(0.00%)\nGA                          0(0.00%)             0(0.00%)\nCT                          0(0.00%)             0(0.00%)\nCG                          0(0.00%)             0(0.00%)\nCA                          0(0.00%)             0(0.00%)\n\nTotalGSNPs                         0                    0\nTA                          0(0.00%)             0(0.00%)\nTG                          0(0.00%)             0(0.00%)\nTC                          0(0.00%)             0(0.00%)\nGT                          0(0.00%)             0(0.00%)\nGC                          0(0.00%)             0(0.00%)\nGA                          0(0.00%)             0(0.00%)\nCG                          0(0.00%)             0(0.00%)\nCT                          0(0.00%)             0(0.00%)\nCA                          0(0.00%)             0(0.00%)\nAG                          0(0.00%)             0(0.00%)\nAT                          0(0.00%)             0(0.00%)\nAC                          0(0.00%)             0(0.00%)\n\nTotalIndels                        0                    0\nA.                          0(0.00%)             0(0.00%)\nT.                          0(0.00%)             0(0.00%)\nG.                          0(0.00%)             0(0.00%)\nC.                          0(0.00%)             0(0.00%)\n.C                          0(0.00%)             0(0.00%)\n.T                          0(0.00%)             0(0.00%)\n.G                          0(0.00%)             0(0.00%)\n.A                          0(0.00%)             0(0.00%)\n\nTotalGIndels                       0                    0\nT.                          0(0.00%)             0(0.00%)\nG.                          0(0.00%)             0(0.00%)\nC.                          0(0.00%)             0(0.00%)\nA.                          0(0.00%)             0(0.00%)\n.A                          0(0.00%)             0(0.00%)\n.C                          0(0.00%)             0(0.00%)\n.T                          0(0.00%)             0(0.00%)\n.G                          0(0.00%)             0(0.00%)\n", "note": "abridged: the other files are the program's"}
```

dnadiff -d d2.delta: the report of the 1-to-1 and M-to-M alignments; the files other than the report are omitted here.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
