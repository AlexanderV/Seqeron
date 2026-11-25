# clustalw_command

ClustalW 2.1 as a function of its command line and the text of the files it names.

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `clustalw_command` |
| **Method ID** | `ClustalWProgram.Run` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

ClustalW 2.1 as a function of its command line and the text of the files it names (ClustalWProgram.Run): the arguments are those of the program (-INFILE=x.fa -TYPE=DNA|PROTEIN -ALIGN -OUTPUT=CLUSTAL|GCG|GDE|PHYLIP|PIR|NEXUS|FASTA -OUTORDER=INPUT|ALIGNED -CASE -SEQNOS -SEQNO_RANGE -RANGE=m,n -OUTFILE=name, -NEWTREE=name / -USETREE=name, -CONVERT (the input as it is in the chosen format), -PROFILE1=a -PROFILE2=b -PROFILE | -SEQUENCES with -NEWTREE1/2 / -USETREE1/2, -NOSECSTR1/2 -SECSTROUT=STRUCTURE|MASK|BOTH|NONE -HELIXGAP -STRANDGAP -LOOPGAP -TERMINALGAP -HELIXENDIN -HELIXENDOUT -STRANDENDIN -STRANDENDOUT, and every alignment parameter: -GAPOPEN -GAPEXT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE -KTUPLE -WINDOW -TOPDIAGS -PAIRGAP -SCORE -PWMATRIX -PWGAPOPEN -PWGAPEXT -MATRIX -DNAMATRIX -PWDNAMATRIX -TRANSWEIGHT -NOPGAP -NOHGAP -NOVGAP -HGAPRESIDUES -GAPDIST -CLUSTERING=NJ|UPGMA -ITERATION=TREE|ALIGNMENT|NONE -NUMITER); option names are case-insensitive, values of enumerations are matched like the program does. `files` gives the text of every file the line names (input sequences / profiles in FASTA, CLUSTAL or GDE format; guide-tree files; matrix files); the result holds the files the program writes (alignment file in the chosen format with the default name <input base>.aln / .msf / .gde / .phy / .pir / .nxs / .fasta or -OUTFILE, guide-tree files <base>.dnd or the -NEWTREE name). Guide tree files: the neighbour-joining tree written by printPhylipTree (and the UPGMA tree), the user tree read with Tree::readTree (rooted or unrooted, with or without distances; the program's refusals - not starting with '(', wrong number of leaves, unknown name, unrooted without distances, trees that send the program into undefined behaviour - are tool errors). Row order: -OUTORDER=ALIGNED (the program's default) lists the rows in the order they entered the progressive alignment. Secondary-structure masks: !SS_ / !GM_ lines of CLUSTAL profile files and "SS_ / "GM_ entries of GDE files are read as the program does and change the gap penalties of -PROFILE / -SEQUENCES. VERIFICATION: guide trees and row order on 3 300 sets, the seven alignment file formats byte for byte on 2 611 sets (301 of them -CONVERT), structure masks on 1 200 profile pairs, all against the installed clustalw 2.1. Not implemented (registered in the audit): -TREE, -BOOTSTRAP, -STATS, -PIM and the input formats PIR, MSF, RSF, EMBL, PHYLIP, NEXUS.

## Core Documentation Reference

- Source: [Seqeron.Genomics.ClustalW/ClustalWProgram.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.ClustalW/ClustalWProgram.cs), [ClustalWCommandLine.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.ClustalW/ClustalWCommandLine.cs), [ClustalWAlignmentWriter.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.ClustalW/ClustalWAlignmentWriter.cs), [SequenceAligner.ClustalWTrees.cs](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.ClustalW/ClustalWAligner.Trees.cs)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `arguments` | array<string> | Yes | The command-line arguments of clustalw, one per element (for example ["-INFILE=x.fa", "-TYPE=DNA", "-ALIGN", "-OUTPUT=GCG"]); a first argument without '-' or '/' is the input file. |
| `files` | object<string,string> | No | The text of the files the arguments name: file name -> content (input sequences or profiles, -USETREE files, matrix files). |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `files` | object<string,string> | The files the program wrote: name -> text (alignment file, guide-tree files). |
| `messages` | array<string> | What the program prints to the error stream while reading the command line (ignored values, warnings). |

## Errors

| Code | Message |
|------|---------|
| 1001 | arguments must contain at least one command-line argument; an input file named on the line is not among `files` |
| 1002 | Unknown option, missing or malformed value (the program's `Error: unknown option` / `parameter required` / bad number), output type unknown, tree file missing |
| 1003 | A guide tree refused by the program (`ClustalWTreeFileException`: not starting with '(', wrong number of leaves, unknown name, unrooted without distances, undefined behaviour of the program) |
| 1004 | A verb or input format that is not implemented (-TREE, -BOOTSTRAP; PIR / MSF / RSF / EMBL / PHYLIP / NEXUS input) |

## Examples

### Example 1: Align four DNA sequences, CLUSTAL output with residue numbers (rows in the program's ALIGNED order)

**Tool Call:**
```json
{
  "tool": "clustalw_command",
  "arguments": {
    "arguments": [
      "-INFILE=ex.fa",
      "-TYPE=DNA",
      "-ALIGN",
      "-OUTPUT=CLUSTAL",
      "-SEQNOS=ON"
    ],
    "files": {
      "ex.fa": ">human\nACGTACGTTTGACCAGTAGGCATCGATGGATTACA\n>chimp\nACGTACGATTGACCAGTAGCCATCGATGGATTACA\n>mouse\nACGTAGGATTGACCGTAGCCATCGATGGTTACAC\n>rat\nACGTAGGATTGACCGTAGCCATCGAAGGTTACAC\n"
    }
  }
}
```

**Response:**
```json
{
  "files": {
    "ex.aln": "CLUSTAL 2.1 multiple sequence alignment\n\n\nmouse           ACGTAGGATTGACC-GTAGCCATCGATGG-TTACAC 34\nrat             ACGTAGGATTGACC-GTAGCCATCGAAGG-TTACAC 34\nchimp           ACGTACGATTGACCAGTAGCCATCGATGGATTACA- 35\nhuman           ACGTACGTTTGACCAGTAGGCATCGATGGATTACA- 35\n                ***** * ****** **** ****** ** ***** \n",
    "ex.dnd": "(\n(\nhuman:0.05798,\nchimp:-0.00084)\n:0.05966,\nmouse:0.00000,\nrat:0.02941);\n"
  },
  "messages": []
}
```

`clustalw -INFILE=ex.fa -TYPE=DNA -ALIGN -OUTPUT=CLUSTAL -SEQNOS=ON` writes exactly these two files (ex.aln, ex.dnd).

### Example 2: The same sequences as a GCG (MSF) file in input order

**Tool Call:**
```json
{
  "tool": "clustalw_command",
  "arguments": {
    "arguments": [
      "-INFILE=ex.fa",
      "-TYPE=DNA",
      "-ALIGN",
      "-OUTPUT=GCG",
      "-OUTORDER=INPUT"
    ],
    "files": {
      "ex.fa": ">human\nACGTACGTTTGACCAGTAGGCATCGATGGATTACA\n>chimp\nACGTACGATTGACCAGTAGCCATCGATGGATTACA\n>mouse\nACGTAGGATTGACCGTAGCCATCGATGGTTACAC\n>rat\nACGTAGGATTGACCGTAGCCATCGAAGGTTACAC\n"
    }
  }
}
```

**Response:**
```json
{
  "files": {
    "ex.msf": "PileUp\n\n\n\n   MSF:   36  Type: N    Check:  6005   .. \n\n Name: human oo  Len:   36  Check:  6841  Weight:  35.7\n Name: chimp oo  Len:   36  Check:  6609  Weight:  0.1\n Name: mouse oo  Len:   36  Check:  6534  Weight:  21.4\n Name: rat oo  Len:   36  Check:  6021  Weight:  42.8\n\n//\n\n\n\nhuman           ACGTACGTTT GACCAGTAGG CATCGATGGA TTACA.\nchimp           ACGTACGATT GACCAGTAGC CATCGATGGA TTACA.\nmouse           ACGTAGGATT GACC.GTAGC CATCGATGG. TTACAC\nrat             ACGTAGGATT GACC.GTAGC CATCGAAGG. TTACAC\n\n",
    "ex.dnd": "(\n(\nhuman:0.05798,\nchimp:-0.00084)\n:0.05966,\nmouse:0.00000,\nrat:0.02941);\n"
  },
  "messages": []
}
```

`clustalw -INFILE=ex.fa -TYPE=DNA -ALIGN -OUTPUT=GCG -OUTORDER=INPUT` writes exactly these files; the weights 35.7 / 0.1 / 21.4 / 42.8 are the alignment's sequence weights.

## See Also

- [multiple_align_clustalw](multiple_align_clustalw.md)
- [clustalw_profile_align](clustalw_profile_align.md)
- [clustalw_sequences_to_profile](clustalw_sequences_to_profile.md)
