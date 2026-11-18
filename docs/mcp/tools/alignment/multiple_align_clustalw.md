# multiple_align_clustalw

ClustalW 2.1's multiple alignment of DNA (clustalw -TYPE=DNA -ALIGN) or PROTEIN (type='protein', clustalw -TYPE=PROTEIN -ALIGN), ported step by step and verified row-for-row against the installed clustalw 2.1 (1 000 DNA sets and 1 000 protein sets, option variants included): pairwise distances (distance='clustalw' = the full affine Smith-Waterman-Gotoh stage, 'clustalw-quick' = the fast k-tuple stage of -QUICKTREE with ktuple / window / topDiagonals / pairGap), ClustalW's neighbour-joining guide tree (rooted by Tree::reRoot), sequence weights derived from the tree branch lengths (Thompson, Higgins & Gibson 1994; noWeights = equal weights), delayed divergent sequences (maxDiv, default 30: sequences with at most that percent identity to every other sequence are aligned last), and profile-profile alignment with position-specific gap penalties (positions that already hold gaps get a reduced opening penalty and half the extension penalty, end gaps are free) and Myers-Miller linear-space alignment. DNA: the IUB weight matrix shifted to non-negative scores (negative = keep it signed) with transitions weighted by transitionWeight (clustalw -TRANSWEIGHT, default 0.5); gapOpen / gapExtend = -GAPOPEN 15 / -GAPEXT 6.66. PROTEIN: the matrix series of ClustalW's profile stage (matrix = -MATRIX gonnet (default) / blosum / pam / id; the table of the series is chosen per profile pair by its percent identity, e.g. GONNET80 / 250 / 120 / 350 / 160, BLOSUM80 / 62 / 45 / 30, PAM20 / 60 / 120 / 350), a gap opening scaled by the matrix average score and the logarithm of the profile lengths (gapOpen / gapExtend = -GAPOPEN 10 / -GAPEXT 0.2), the residue-specific gap penalties of Pascarella & Argos 1992 (noPreferredGaps = -NOPGAP), the hydrophilic-run penalties (noHydrophilicGaps = -NOHGAP; hydrophilicResidues = -HGAPRESIDUES, default GPSNDQEKR, at most 9 letters), the variable-region penalties (variableRegionGaps = clustalw's -NOVGAP, sic: it switches them ON), the gap-separation distance (gapDist = -GAPDIST, default 4; 0 = off) and the pairwise stage of F62 (pwMatrix = -PWMATRIX gonnet / blosum / pam / id, pwGapOpen 10 / pwGapExtend 0.1, ktuple 1-2). User matrices (type='protein'): matrixFile = clustalw -MATRIX=file (one matrix file, or a CLUSTAL_SERIES file whose entries 'MATRIX llimit ulimit file' are looked up in matrixFiles; the entry whose percent-identity range contains the profile pair's identity is used, with the gap-penalty scale 0.5 + (pcid - llimit) / ((ulimit - llimit) * 2)) and pwMatrixFile = clustalw -PWMATRIX=file (verified on 240 + 800 clustalw -MATRIX / -PWMATRIX runs with floating-point / integer / unlabelled / tab-separated / header-subset files; only 2 sets, both -NOVGAP, hit the gdist residual below). Known residual: with variableRegionGaps and a two-sequence profile of more than 60 % identity ClustalW reads an uninitialised C++ member (gdist), so on ~0.2 % of such runs the heap decides; the port emulates the usual value. END-GAP OPTIONS (both types): useEndGaps = UserParameters::useEndGaps and endGapPenalties = UserParameters::endGapPenalties (reachable in clustalw only through its interactive menu / a patched build; verified on 340 sets x 4 modes of the patched clustalw 2.1 binary): with endGapPenalties ClustalW 2.1 writes NO output (exit status 0) when the score of a first-pass profile alignment is negative - this tool then fails with an explicit 'abandons the alignment' error. NUCLEOTIDE MATRICES (type='dna'): dnaMatrix / pwDnaMatrix = -DNAMATRIX / -PWDNAMATRIX IUB (default) or CLUSTALW; dnaMatrixFile / pwDnaMatrixFile = -DNAMATRIX=file / -PWDNAMATRIX=file (text of a matrix file; verified on 360 sets of the installed clustalw: the command line reads and validates the file but never uses it for DNA - a defect of CommandLineParser, userDnaMatrixUse='commandline' - while the interactive menu uses it, userDnaMatrixUse='applied', the default); transitionWeight = -TRANSWEIGHT also rewrites the cells of the user matrix at the letters of header positions 0 / 4, 2 / 11 and 2 / 12. Rows are returned in input order. match / mismatch / scoreGapOpen / scoreGapExtend (DNA) or scoreMatrix / scoreGapOpen / scoreGapExtend (protein; default BLOSUM62, open -9.8, extend -0.2) define only the reported totalScore (natural affine sum of pairs, Altschul 1989).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_clustalw` |
| **Method ID** | `SequenceAligner.MultipleAlignClustalW` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

ClustalW 2.1's multiple alignment of DNA (clustalw -TYPE=DNA -ALIGN) or PROTEIN (type='protein', clustalw -TYPE=PROTEIN -ALIGN), ported step by step and verified row-for-row against the installed clustalw 2.1 (1 000 DNA sets and 1 000 protein sets, option variants included): pairwise distances (distance='clustalw' = the full affine Smith-Waterman-Gotoh stage, 'clustalw-quick' = the fast k-tuple stage of -QUICKTREE with ktuple / window / topDiagonals / pairGap), ClustalW's neighbour-joining guide tree (rooted by Tree::reRoot), sequence weights derived from the tree branch lengths (Thompson, Higgins & Gibson 1994; noWeights = equal weights), delayed divergent sequences (maxDiv, default 30: sequences with at most that percent identity to every other sequence are aligned last), and profile-profile alignment with position-specific gap penalties (positions that already hold gaps get a reduced opening penalty and half the extension penalty, end gaps are free) and Myers-Miller linear-space alignment. DNA: the IUB weight matrix shifted to non-negative scores (negative = keep it signed) with transitions weighted by transitionWeight (clustalw -TRANSWEIGHT, default 0.5); gapOpen / gapExtend = -GAPOPEN 15 / -GAPEXT 6.66. PROTEIN: the matrix series of ClustalW's profile stage (matrix = -MATRIX gonnet (default) / blosum / pam / id; the table of the series is chosen per profile pair by its percent identity, e.g. GONNET80 / 250 / 120 / 350 / 160, BLOSUM80 / 62 / 45 / 30, PAM20 / 60 / 120 / 350), a gap opening scaled by the matrix average score and the logarithm of the profile lengths (gapOpen / gapExtend = -GAPOPEN 10 / -GAPEXT 0.2), the residue-specific gap penalties of Pascarella & Argos 1992 (noPreferredGaps = -NOPGAP), the hydrophilic-run penalties (noHydrophilicGaps = -NOHGAP; hydrophilicResidues = -HGAPRESIDUES, default GPSNDQEKR, at most 9 letters), the variable-region penalties (variableRegionGaps = clustalw's -NOVGAP, sic: it switches them ON), the gap-separation distance (gapDist = -GAPDIST, default 4; 0 = off) and the pairwise stage of F62 (pwMatrix = -PWMATRIX gonnet / blosum / pam / id, pwGapOpen 10 / pwGapExtend 0.1, ktuple 1-2). User matrices (type='protein'): matrixFile = clustalw -MATRIX=file (one matrix file, or a CLUSTAL_SERIES file whose entries 'MATRIX llimit ulimit file' are looked up in matrixFiles; the entry whose percent-identity range contains the profile pair's identity is used, with the gap-penalty scale 0.5 + (pcid - llimit) / ((ulimit - llimit) * 2)) and pwMatrixFile = clustalw -PWMATRIX=file (verified on 240 + 800 clustalw -MATRIX / -PWMATRIX runs with floating-point / integer / unlabelled / tab-separated / header-subset files; only 2 sets, both -NOVGAP, hit the gdist residual below). Known residual: with variableRegionGaps and a two-sequence profile of more than 60 % identity ClustalW reads an uninitialised C++ member (gdist), so on ~0.2 % of such runs the heap decides; the port emulates the usual value. END-GAP OPTIONS (both types): useEndGaps = UserParameters::useEndGaps and endGapPenalties = UserParameters::endGapPenalties (reachable in clustalw only through its interactive menu / a patched build; verified on 340 sets x 4 modes of the patched clustalw 2.1 binary): with endGapPenalties ClustalW 2.1 writes NO output (exit status 0) when the score of a first-pass profile alignment is negative - this tool then fails with an explicit 'abandons the alignment' error. NUCLEOTIDE MATRICES (type='dna'): dnaMatrix / pwDnaMatrix = -DNAMATRIX / -PWDNAMATRIX IUB (default) or CLUSTALW; dnaMatrixFile / pwDnaMatrixFile = -DNAMATRIX=file / -PWDNAMATRIX=file (text of a matrix file; verified on 360 sets of the installed clustalw: the command line reads and validates the file but never uses it for DNA - a defect of CommandLineParser, userDnaMatrixUse='commandline' - while the interactive menu uses it, userDnaMatrixUse='applied', the default); transitionWeight = -TRANSWEIGHT also rewrites the cells of the user matrix at the letters of header positions 0 / 4, 2 / 11 and 2 / 12. Rows are returned in input order. match / mismatch / scoreGapOpen / scoreGapExtend (DNA) or scoreMatrix / scoreGapOpen / scoreGapExtend (protein; default BLOSUM62, open -9.8, extend -0.2) define only the reported totalScore (natural affine sum of pairs, Altschul 1989).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.ClustalWMultiple.cs#L69](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.ClustalWMultiple.cs#L69)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.6)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one): DNA (A, C, G, T / IUB codes) by default, amino acids (letters A-Z except J) with type='protein'. |
| `gapOpen` | number | No | ClustalW -GAPOPEN: gap-opening penalty of the profile alignment (default 15 for DNA, 10 for protein; a penalty, positive). |
| `gapExtend` | number | No | ClustalW -GAPEXT: gap-extension penalty of the profile alignment (default 6.66 for DNA, 0.2 for protein). |
| `distance` | string | No | Pairwise stage that gives the guide-tree distances: 'clustalw' (default; full DP, clustalw -ALIGN) or 'clustalw-quick' (fast k-tuple mode, clustalw -QUICKTREE). (default: clustalw) |
| `pwGapOpen` | number | No | ClustalW -PWGAPOPEN: pairwise gap-opening penalty of the full pairwise stage (default 15 for DNA, 10 for protein). |
| `pwGapExtend` | number | No | ClustalW -PWGAPEXT: pairwise gap-extension penalty of the full pairwise stage (default 6.66 for DNA, 0.1 for protein). |
| `transitionWeight` | number | No | DNA only, ClustalW -TRANSWEIGHT: weight of transitions (A-G, C-T/U) in the profile matrix (default 0.5; 0 = as mismatches, 1 = as matches). (default: 0.5) |
| `maxDiv` | integer | No | ClustalW -MAXDIV (default 30): a sequence whose highest percent identity to any other sequence is at most this value is delayed (aligned last, after the rest). (default: 30) |
| `noWeights` | boolean | No | ClustalW -NOWEIGHTS: no tree-derived sequence weights. (default: false) |
| `negative` | boolean | No | ClustalW -NEGATIVE: keep the substitution matrix signed instead of shifting it to non-negative scores (protein: one matrix BLOSUM40 / PAM120 / GONNET250 instead of the series). (default: false) |
| `ktuple` | integer | No | distance='clustalw-quick' only: k-tuple (word) size, 1-4 for DNA (clustalw -KTUPLE, default 2), 1-2 for protein (default 1). |
| `window` | integer | No | distance='clustalw-quick' only: diagonals scanned around each top diagonal (clustalw -WINDOW, default 4 for DNA, 5 for protein). |
| `topDiagonals` | integer | No | distance='clustalw-quick' only: number of best diagonals used (clustalw -TOPDIAGS, default 4 for DNA, 5 for protein). |
| `pairGap` | integer | No | distance='clustalw-quick' only: gap penalty of the fragment chain (clustalw -PAIRGAP; DNA default 5, or ktuple + 4 when ktuple is given; protein default 3, or ktuple + 3, as the clustalw command line does). |
| `match` | integer | No | Score for matching base in the reported totalScore only (DNA). (default: 1) |
| `mismatch` | integer | No | Score for mismatching base in the reported totalScore only (DNA). (default: -1) |
| `scoreGapOpen` | number | No | Gap-open score of the reported totalScore only (negative = penalty; default -2 for DNA, -9.8 for protein). |
| `scoreGapExtend` | number | No | Gap-extend score of the reported totalScore only (negative = penalty; default -1 for DNA, -0.2 for protein). |
| `type` | string | No | 'dna' (default) or 'protein' (ClustalW -TYPE). (default: dna) |
| `matrix` | string | No | type='protein': ClustalW -MATRIX, the matrix series of the profile alignment: 'gonnet' (default), 'blosum', 'pam' or 'id'. |
| `pwMatrix` | string | No | type='protein': ClustalW -PWMATRIX, the pairwise matrix: 'gonnet' (default; GONNET250), 'blosum' (BLOSUM30), 'pam' (PAM350) or 'id'. |
| `noPreferredGaps` | boolean | No | type='protein': ClustalW -NOPGAP: switch off the residue-specific (Pascarella & Argos) gap penalties. (default: false) |
| `noHydrophilicGaps` | boolean | No | type='protein': ClustalW -NOHGAP: switch off the hydrophilic-run gap penalties. (default: false) |
| `variableRegionGaps` | boolean | No | type='protein': ClustalW -NOVGAP (sic): switch the variable-region gap penalties on (clustalw's default is off). (default: false) |
| `gapDist` | integer | No | type='protein': ClustalW -GAPDIST: gap separation distance (default 4; 0 = off): columns this close to an existing gap get a higher opening penalty. (default: 4) |
| `hydrophilicResidues` | string | No | type='protein': ClustalW -HGAPRESIDUES: residues of a hydrophilic stretch (default GPSNDQEKR; at most the first 9 letters are used). (default: GPSNDQEKR) |
| `scoreMatrix` | string | No | type='protein': substitution matrix of the reported totalScore only: BLOSUM62 (default) or PAM250. |
| `matrixFile` | string | No | type='protein': ClustalW -MATRIX=file, the TEXT of a user matrix file for the profile alignment (NCBI / EMBOSS layout: '#' comments, a header line of residue letters, one row per residue with an optional leading label, one extra '*' row, integer or floating-point values (a row whose first value contains '.' is scaled by 10); only the lower triangle is read) or of a CLUSTAL_SERIES file (lines 'MATRIX llimit ulimit file', file names resolved through matrixFiles). Replaces 'matrix'. |
| `matrixFiles` | object | No | type='protein': the matrix files named by a CLUSTAL_SERIES matrixFile: file name -> text of the matrix file. |
| `pwMatrixFile` | string | No | type='protein': ClustalW -PWMATRIX=file, the TEXT of a user matrix file for the pairwise (distance) stage (same layout as matrixFile, a single matrix). Replaces 'pwMatrix'. |
| `useEndGaps` | boolean | No | End-gap option of the profile stage, UserParameters::useEndGaps (the interactive menu item 'Toggle End Gap Separation'; the command-line flag -ENDGAPS only ever writes false): gaps at the sequence ends count as gaps in the position-specific penalties. Default false. (default: false) |
| `endGapPenalties` | boolean | No | End-gap option of the profile stage, UserParameters::endGapPenalties: the end gaps of a profile alignment are penalised instead of free. Default false. With it ClustalW 2.1 abandons the alignment silently (no output file) when the score of a first-pass profile alignment is negative; this tool then fails with a 'ClustalW 2.1 abandons the alignment ...' error (ClustalWAbortedException). (default: false) |
| `dnaMatrix` | string | No | type='dna': ClustalW -DNAMATRIX built-in: 'iub' (default) or 'clustalw' (ClustalW 1.6's identity matrix, gap penalties scaled by 0.66). Ignored when dnaMatrixFile is given. |
| `pwDnaMatrix` | string | No | type='dna': ClustalW -PWDNAMATRIX built-in for the pairwise stage: 'iub' (default) or 'clustalw' (opening x 0.6667, extension x 0.751). Ignored when pwDnaMatrixFile is given. |
| `dnaMatrixFile` | string | No | type='dna': ClustalW -DNAMATRIX=file, the TEXT of a nucleotide matrix file for the profile stage (same layout as matrixFile, a single matrix; the transition weights of getProfileAlignMatrix index the matrix with the header letters at positions 0, 2, 4, 11 and 12 - give at least 13 letters). How it is used: see userDnaMatrixUse. |
| `pwDnaMatrixFile` | string | No | type='dna': ClustalW -PWDNAMATRIX=file, the TEXT of a nucleotide matrix file for the full pairwise stage (a single matrix, x100, never shifted). See userDnaMatrixUse. |
| `userDnaMatrixUse` | string | No | type='dna', with dnaMatrixFile / pwDnaMatrixFile: 'applied' (default; what ClustalW's interactive menu does - the matrix is used) or 'commandline' (what clustalw -DNAMATRIX=file -PWDNAMATRIX=file does: the file is read and validated but registered under the PROTEIN slot - a defect of CommandLineParser - so the DNA stages keep the IUB / built-in matrix). (default: applied) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA (input order). |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a nucleotide beats a gap on ties). |
| `totalScore` | integer | Natural affine sum-of-pairs score of the alignment under match / mismatch / scoreGapOpen / scoreGapExtend (ClustalW itself reports no such score). |

## Errors

| Code | Message |
|------|---------|
| 1010 | Invalid DNA at sequences[i] (type=dna); residue not in A-Z except J (type=protein) |
| 1002 | gapOpen, gapExtend, pwGapOpen, pwGapExtend, transitionWeight must be non-negative; dnaMatrix / pwDnaMatrix must be iub or clustalw; dnaMatrixFile / pwDnaMatrixFile must be a single ClustalW matrix file; userDnaMatrixUse must be applied or commandline; matrixFile / pwMatrixFile must be a ClustalW matrix file (or series with its matrixFiles); matrixFile / pwMatrixFile must be a ClustalW matrix file (or series with its matrixFiles); distance must be clustalw or clustalw-quick; type must be dna or protein; matrix / pwMatrix / scoreMatrix must be a known name |

## Examples

### Example 1: Five related 30-mers

**Tool Call:**
```json
{
  "tool": "multiple_align_clustalw",
  "arguments": {
    "sequences": [
      "ACGTACGTTTGACCAGTAGGCATCGATGGA",
      "ACGTACGATTGACCAGTAGCCATCGATGGA",
      "TTGACCAGTAGGCATCGATCCCAAGT",
      "ACGTAGGATTGACCGTAGCCATCGATGGA",
      "GGGTTTACCAGTAGGCATCGATGGAAACCT"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTACGTTTGACCAGTAGGCATCGATGGA-----",
    "ACGTACGATTGACCAGTAGCCATCGATGGA-----",
    "--------TTGACCAGTAGGCATCGATCCCAAGT-",
    "ACGTAGGATTGACC-GTAGCCATCGATGGA-----",
    "-----GGGTTTACCAGTAGGCATCGATGGAAACCT"
  ],
  "consensus": "ACGTACGATTGACCAGTAGGCATCGATGGA-----",
  "totalScore": 77
}
```

### Example 2: Four protein sequences (type='protein')

**Tool Call:**
```json
{
  "tool": "multiple_align_clustalw",
  "arguments": {
    "sequences": [
      "LFFFFFFFMFFFFFDKFFFFFFFFFFFFFFFF",
      "FLAFLVKTYAHTGLKDKEGF",
      "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF",
      "FFFRNPFFPFGFFFFFFFFFFFFSFFHFFFFFF"
    ],
    "type": "protein"
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "LFFF---FFFFMFFFFFDKFFFFFFFFFFFFFFFF",
    "------------FLAFLVKTYAHTGLKDKEGF---",
    "-FFF---FFFFFFFFFFFFFFFFFFFFFFFFFFFF",
    "-FFFRNPFFPFGFFFFFFFFFFFFSFFHFFFFFF-"
  ],
  "consensus": "-FFF---FFFFFFFFFFFFFFFFFFFFFFFFFFFF",
  "totalScore": 290
}
```

## Worked Example

The rows are exactly what `clustalw -TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` writes for these five sequences (ClustalW 2.1, run here): the five sequences are aligned along ClustalW's neighbour-joining guide tree with tree-derived weights and position-specific gap penalties; sequences with at most 30 % identity to every other would be delayed and added last.

## Protein stage (type='protein')

The protein example is set 20 of `clustalw_protein_msa_oracle.tsv`: the rows are exactly what `clustalw -TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` writes (GONNET series, `-GAPOPEN 10 -GAPEXT 0.2`, residue-specific / hydrophilic gap penalties, `-GAPDIST 4`); `totalScore` is the natural affine sum of pairs under BLOSUM62 (open -9.8, extend -0.2, EMBOSS convention). The port is verified on 1 000 protein sets (option variants incl. matrix series, `-NOPGAP/-NOHGAP/-NOVGAP/-GAPDIST/-HGAPRESIDUES`); the only known residual is the undefined `gdist` member of ClustalW with `-NOVGAP` (see the algorithm doc §4.6).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.ClustalWMultiple.cs#L69](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.ClustalWMultiple.cs#L69)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
