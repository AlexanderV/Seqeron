# est2genome

EMBOSS 6.6.0 est2genome: the spliced alignment of ESTs to a genomic sequence (exons, introns with donor / acceptor sites, both strands).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `est2genome` |
| **Method ID** | `EmbossEst2genomePrograms.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

EMBOSS 6.6.0 est2genome (Richard Mott's est_genome, embest.c): the spliced alignment of EST sequences (a FASTA text, every record is aligned) to an unspliced genomic sequence (the first record of a FASTA text). A local alignment with introns: matches +match, mismatches -mismatch, a base deleted from either sequence -gapPenalty, an intron of any length -intronPenalty, or -splicePenalty when it starts after a donor and ends on an acceptor site (usesplice, default true: GT/AG for a forward gene, CT/AC for a reversed gene). mode 'both' (default) aligns the forward and the reverse-complemented EST with forward-gene splice sites, then repeats the better one with reversed-gene splice sites and reports the best (the Note line tells which: forward est / reversed est, forward / REVERSED GENE); 'forward' and 'reverse' (case-insensitive, unique abbreviations) align one strand. An alignment is written when its score is at least minScore: the Exon lines (score, percent identity, genomic and EST coordinates, names, EST description), the +Intron / -Intron / ?Intron lines, the Span line and, after a blank line, the Segment lines; with align the alignment itself (introns as their first and last five bases and their length, > forward / < reversed splice direction) in blocks of width columns and 'Alignment Score'. reverse prints the coordinates of a reversed EST on its original strand; best = false (-nobest) also writes the other alignment (with a single strand mode the program crashes: refused with InvalidOperationException); space is the megabytes of the matrix above which the linear-space divide and conquer is used (raised to 1.5e-6 x genome length with a warning in messages); shuffle > 0 sets minScore to the best score + 1 (the program's shuffle does not change the sequence); seed has no effect; width of 0 or less with align never ends in the program (refused). Negative penalties are used as their absolute values. sequenceQualifiers: EMBOSS sequence qualifiers after the files (-sbegin1 5 -send1 90 -sreverse1 for the ESTs, -sbegin2 ... for the genomic sequence). Byte for byte equal to est2genome on 330 runs (forward / reverse genes, ESTs on either strand, every option, the linear-space recursion, the crashes).

## Core Documentation Reference

- Source: [EmbossEst2genomePrograms.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/EmbossEst2genomePrograms.cs)
- Algorithm doc: [Pairwise_Output_Formats_And_EMBOSS_Programs.md](../../../algorithms/Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md)
- Oracle: [generate_est2genome_oracle.py](../../../../scripts/alignment/generate_est2genome_oracle.py)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `estFasta` | string | Yes | The EST sequences as FASTA text (-estsequence); every record is aligned. |
| `genomeFasta` | string | Yes | The genomic sequence as FASTA text (-genomesequence); the first record is used. |
| `match` | integer | No | Score for matching two bases (-match, default 1). (default: 1) |
| `mismatch` | integer | No | Cost for mismatching two bases (-mismatch, default 1). (default: 1) |
| `gapPenalty` | integer | No | Cost for deleting a single base in either sequence, excluding introns (-gappenalty, default 2). (default: 2) |
| `intronPenalty` | integer | No | Cost for an intron, independent of length (-intronpenalty, default 40). (default: 40) |
| `splicePenalty` | integer | No | Cost for an intron between a donor and an acceptor site (-splicepenalty, default 20). (default: 20) |
| `minScore` | integer | No | Exclude alignments with scores below this threshold (-minscore, default 30). (default: 30) |
| `reverse` | boolean | No | Report the coordinates of a reversed EST on the original strand (-reverse, default false). (default: false) |
| `useSplice` | boolean | No | Use donor and acceptor splice sites (-usesplice, default true). (default: true) |
| `mode` | string | No | Comparison mode (-mode): 'both' (default), 'forward' or 'reverse'. |
| `best` | boolean | No | Print only the best alignment (-best, default true; false = -nobest). (default: true) |
| `space` | number | No | Space threshold in megabytes for the linear-space recursion (-space, default 10.0). (default: 10) |
| `shuffle` | integer | No | Shuffle count (-shuffle, default 0). (default: 0) |
| `seed` | integer | No | Random number seed (-seed, default 20825; no effect on the output). (default: 20825) |
| `align` | boolean | No | Show the alignment (-align, default false). (default: false) |
| `width` | integer | No | Alignment width (-width, default 50). (default: 50) |
| `sequenceQualifiers` | string | No | EMBOSS sequence qualifiers as typed after the files: -sbegin1 5 -send1 90 -sreverse1 -slower1 -supper1 (input 1 = the ESTs) and -sbegin2 ... (input 2 = the genomic sequence). The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. -sformat -sask -squick -sdbname -sid -iquery -ufo -fformat -fopenfile are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `text` | string | The -outfile text. |
| `messages` | string | What the program writes to standard error (the space warning). |

## Errors

| Code | Message |
|------|---------|
| 1001 | The EST sequences cannot be empty |
| 1001 | est2genome: 'x' is not a valid menu option (both, forward, reverse). |
| 1001 | est2genome: segmentation fault (-nobest with a single strand mode) |

## Examples

### Example 1: An EST spliced from three exons

**User Prompt:**
> Align this EST to the genomic sequence.

**Expected Tool Call:**
```json
{"tool": "est2genome", "arguments": {"estFasta": ">est1 an est\nGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCT\n", "genomeFasta": ">genomic some description\nACGTGGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTGTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAGTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTTTTT\n"}}
```

**Response:**
```json
{"text": "Note Best alignment is between forward est and forward genome, and splice sites imply forward gene\nExon        60 100.0     5    64 genomic          1    60 est1          an est\n+Intron    -20   0.0    65   138 genomic     \nExon        80 100.0   139   218 genomic         61   140 est1          an est\n+Intron    -20   0.0   219   262 genomic     \nExon        50 100.0   263   312 genomic        141   190 est1          an est\n\nSpan       150 100.0     5   312 genomic          1   190 est1          an est\n\nSegment     60 100.0     5    64 genomic          1    60 est1          an est\nSegment     80 100.0   139   218 genomic         61   140 est1          an est\nSegment     50 100.0   263   312 genomic        141   190 est1          an est\n", "messages": ""}
```

est2genome -estsequence e.fa -genomesequence g.fa: three exons joined by two GT/AG introns (74 and 44 bases, -20 = the splice penalty), score 150.

## See Also

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
