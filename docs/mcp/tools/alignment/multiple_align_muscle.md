# multiple_align_muscle

MUSCLE 3.8.1551's multiple alignment (Edgar 2004, NAR 32:1792; muscle -in x.fa -out y.fa), ported step by step from the program's source: letters upper-cased, gaps / white space dropped, letters outside the alphabet replaced by the wildcard (X protein, N nucleotide), a leading M of a protein aligned as X and restored (methionine hack); first guide tree from a k-mer distance (distance1: kmer6_6 over the six MAFFT residue groups for proteins, kmer4_6 for nucleotides, kmer20_3, kmer20_4, kbit20_3, pwkimura, pwscoredist) clustered by UPGMB (0.1 x average + 0.9 x minimum linkage; cluster1: upgma, upgmamin, upgmamax, neighborjoining) and rooted (root1 pseudo / midlongestspan / minavgleafdist), ClustalW tree weights (weight1: none, henikoff, henikoffpb, gsc, clustalw, threeway), profile-profile dynamic programming with position-specific gap open / close penalties and half-open terminal gaps (termGaps half / full / ext), the log-expectation score (VTML200 likelihood ratios, center -0.52, gap open -2.9) for proteins or MUSCLE's nucleotide sum-of-pairs score (BLASTZ matrix, gap open -400) for DNA / RNA (scoreModel sp = PAM200, sv = VTML240 for proteins, a matrix file replaces the matrix of the SP score); maxIters >= 2 rebuilds the tree from the distances of the draft (distance2 pctidkimura / pctidlog / edit / scoredist, cluster2, root2) and realigns only the subtrees whose leaf sets changed (maxTrees); maxIters >= 3 adds maxIters - 2 rounds of restricted-partition refinement (every guide-tree edge splits the sequences in two, the two profiles are re-aligned and kept when the objective score objScore sp / dp / xp / ps / spf / spm improves, inside the blocks between anchor columns unless anchors = false). The program's own default is maxIters 8 (params.cpp g_uMaxIters = 8; verified with the installed muscle3, whose -help text 'default 16' is stale: a plain run writes the rows of -maxiters 8, not 16); the default here is 2 (draft + tree refinement), pass maxIters = 8 for the program's behaviour. Options accepted but without effect, because the 3.8.1551 source compiles them out: diags / diags1 / diags2 / diagLength / diagBreak / diagMargin, hydro / hydroFactor / tomHydro, gapOpen2 / gapExtend2 / gapAmbig; maxMb / maxHours (resource limits). useTree = a rooted binary Newick tree with edge lengths whose labels are the names (skips the tree refinement); clusterOnly / computeWeights stop early. Rows are returned in input order (treeOrder = the order muscle writes) with the files the program can write (tree1, tree2, distance matrices, weights). totalScore is the natural affine sum of pairs (BLOSUM62 10 / 0.2 for proteins, EDNAFULL 10 / 0.5 for nucleotides), not MUSCLE's objective. VERIFICATION: maxIters 1 and 2 with the option set of the first Docker round (le / sp / sv, gapOpen / gapExtend / center, distance1 kmer6_6, cluster1 / cluster2 upgma*, sueff, maxTrees 0) equal the MUSCLE 3.8.1551 binary built from the same source on 420 / 420 user-run sets (rows and output order). Everything else (maxIters >= 3, objScore, weights, the other distances, neighborjoining and roots, anchors, termGaps, brenner, useTree, the tree / distance / weight files) is verified against the same binary on 630 further user-run sets (docs/Validation/review-2026-09/handoff/A1-10c-muscle-r2): rows in the written order, files and the program's refusals equal; the user matrix (matrix) is verified for proteins with a '*' column (the file without it and nucleotides: round 3).

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_muscle` |
| **Method ID** | `SequenceAligner.MultipleAlignMuscle` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

MUSCLE 3.8.1551's multiple alignment (Edgar 2004, NAR 32:1792; muscle -in x.fa -out y.fa), ported step by step from the program's source: letters upper-cased, gaps / white space dropped, letters outside the alphabet replaced by the wildcard (X protein, N nucleotide), a leading M of a protein aligned as X and restored (methionine hack); first guide tree from a k-mer distance (distance1: kmer6_6 over the six MAFFT residue groups for proteins, kmer4_6 for nucleotides, kmer20_3, kmer20_4, kbit20_3, pwkimura, pwscoredist) clustered by UPGMB (0.1 x average + 0.9 x minimum linkage; cluster1: upgma, upgmamin, upgmamax, neighborjoining) and rooted (root1 pseudo / midlongestspan / minavgleafdist), ClustalW tree weights (weight1: none, henikoff, henikoffpb, gsc, clustalw, threeway), profile-profile dynamic programming with position-specific gap open / close penalties and half-open terminal gaps (termGaps half / full / ext), the log-expectation score (VTML200 likelihood ratios, center -0.52, gap open -2.9) for proteins or MUSCLE's nucleotide sum-of-pairs score (BLASTZ matrix, gap open -400) for DNA / RNA (scoreModel sp = PAM200, sv = VTML240 for proteins, a matrix file replaces the matrix of the SP score); maxIters >= 2 rebuilds the tree from the distances of the draft (distance2 pctidkimura / pctidlog / edit / scoredist, cluster2, root2) and realigns only the subtrees whose leaf sets changed (maxTrees); maxIters >= 3 adds maxIters - 2 rounds of restricted-partition refinement (every guide-tree edge splits the sequences in two, the two profiles are re-aligned and kept when the objective score objScore sp / dp / xp / ps / spf / spm improves, inside the blocks between anchor columns unless anchors = false). The program's own default is maxIters 8 (params.cpp g_uMaxIters = 8; verified with the installed muscle3, whose -help text 'default 16' is stale: a plain run writes the rows of -maxiters 8, not 16); the default here is 2 (draft + tree refinement), pass maxIters = 8 for the program's behaviour. Options accepted but without effect, because the 3.8.1551 source compiles them out: diags / diags1 / diags2 / diagLength / diagBreak / diagMargin, hydro / hydroFactor / tomHydro, gapOpen2 / gapExtend2 / gapAmbig; maxMb / maxHours (resource limits). useTree = a rooted binary Newick tree with edge lengths whose labels are the names (skips the tree refinement); clusterOnly / computeWeights stop early. Rows are returned in input order (treeOrder = the order muscle writes) with the files the program can write (tree1, tree2, distance matrices, weights). totalScore is the natural affine sum of pairs (BLOSUM62 10 / 0.2 for proteins, EDNAFULL 10 / 0.5 for nucleotides), not MUSCLE's objective. VERIFICATION: maxIters 1 and 2 with the option set of the first Docker round (le / sp / sv, gapOpen / gapExtend / center, distance1 kmer6_6, cluster1 / cluster2 upgma*, sueff, maxTrees 0) equal the MUSCLE 3.8.1551 binary built from the same source on 420 / 420 user-run sets (rows and output order). Everything else (maxIters >= 3, objScore, weights, the other distances, neighborjoining and roots, anchors, termGaps, brenner, useTree, the tree / distance / weight files) is verified against the same binary on 630 further user-run sets (docs/Validation/review-2026-09/handoff/A1-10c-muscle-r2): rows in the written order, files and the program's refusals equal; the user matrix (matrix) is verified for proteins with a '*' column (the file without it and nucleotides: round 3).

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.Muscle.cs#L67](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Muscle.cs#L67)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.9)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one, none empty): letters, '-', '.' and white space only. |
| `type` | string | No | 'auto' (default; MUSCLE's guess: nucleotide when >= 95 % of the first 100 residues are A, C, G, T, U or N), 'protein', 'dna' or 'rna' (MUSCLE -seqtype). (default: auto) |
| `maxIters` | integer | No | MUSCLE -maxiters (at least 1): 1 = progressive draft, 2 (default) = plus tree refinement, 3 and above = plus restricted-partition refinement iterations (the program's default is 8: pass 8 for a plain 'muscle -in x.fa'; the -help text 'default 16' of 3.8.1551 is stale). (default: 2) |
| `scoreModel` | string | No | Profile-profile score (-le / -sp / -sv / -spn): 'default' (LE for proteins, SPN for nucleotides), 'le', 'sp' (PAM200), 'sv' (VTML240), 'spn'. A score of the wrong alphabet is replaced like the program does (SPN for nucleotides, LE for proteins). (default: default) |
| `gapOpen` | number | No | MUSCLE -gapopen (a negative score; default -2.9 LE, -1439 SP, -300 SV, -400 SPN DNA, -420 RNA). |
| `gapExtend` | number | No | MUSCLE -gapextend (default 0). |
| `center` | number | No | MUSCLE -center (default -0.52 LE, 0 otherwise, -300 RNA). |
| `distance1` | string | No | MUSCLE -distance1: 'default' (kmer6_6 protein, kmer4_6 nucleotide), 'kmer6_6', 'kmer4_6' (nucleotides only), 'kmer20_3', 'kmer20_4', 'kbit20_3', 'pwkimura', 'pwscoredist' (protein only). kmer20_3 / kmer20_4 / kbit20_3 need sequences of at least 4 residues (the program reads outside its buffer for shorter ones). (default: default) |
| `cluster1` | string | No | MUSCLE -cluster1 (first guide tree): 'upgmb' (default), 'upgma', 'upgmamin', 'upgmamax' or 'neighborjoining'. (default: upgmb) |
| `cluster2` | string | No | MUSCLE -cluster2 (tree of iteration 2): 'upgmb' (default), 'upgma', 'upgmamin', 'upgmamax' or 'neighborjoining'. (default: upgmb) |
| `sueff` | number | No | MUSCLE -sueff (default 0.1): weight of the average linkage in UPGMB. |
| `maxTrees` | integer | No | MUSCLE -maxtrees (default 1): 0 switches the tree refinement of iteration 2 off; more than 1 repeats it only for a weighting other than ClustalW. (default: 1) |
| `treeOrder` | boolean | No | true: rows in the order muscle writes them (depth-first leaf order of the final guide tree, -group) instead of the input order. (default: false) |
| `distance2` | string | No | MUSCLE -distance2 (distance of the aligned pairs of iteration 2): 'default' (pctidkimura), 'pctidkimura', 'pctidlog', 'edit', 'scoredist' (protein only). (default: default) |
| `root1` | string | No | MUSCLE -root1: 'pseudo' (default: the root of the clustering), 'midlongestspan' or 'minavgleafdist'. (default: pseudo) |
| `root2` | string | No | MUSCLE -root2 (tree of iteration 2): 'pseudo', 'midlongestspan' or 'minavgleafdist'. (default: pseudo) |
| `weight1` | string | No | MUSCLE -weight1 (sequence weights of the progressive stage): 'clustalw' (default), 'none', 'henikoff', 'henikoffpb', 'gsc', 'threeway'. (default: clustalw) |
| `weight2` | string | No | MUSCLE -weight2 (sequence weights of the refinement): same values as weight1. (default: clustalw) |
| `objScore` | string | No | MUSCLE -objscore (refinement objective): 'spm' (default: xp up to 100 sequences, spf above), 'sp', 'dp', 'xp', 'ps' (LE score only), 'spf'. (default: spm) |
| `termGaps` | string | No | MUSCLE -termgaps: 'half' (default), 'full', 'ext'. (default: half) |
| `anchors` | boolean | No | MUSCLE -anchors (true, default) / -noanchors (false): refine inside the blocks between anchor columns. (default: true) |
| `anchorSpacing` | integer | No | MUSCLE -anchorspacing (default 32). (default: 32) |
| `smoothWindow` | integer | No | MUSCLE -smoothwindow (odd, default 7). (default: 7) |
| `smoothScoreCeil` | number | No | MUSCLE -smoothscoreceil (default 3 LE, 200 SP, 90 SV, 999 SPN). |
| `minBestColScore` | number | No | MUSCLE -minbestcolscore (default 2 LE, 300 SP, 130 SV, 90 SPN). |
| `minSmoothScore` | number | No | MUSCLE -minsmoothscore (default 1 LE, 125 SP, 40 SV, 90 SPN). |
| `brenner` | boolean | No | MUSCLE -brenner: build the root alignment by aligning every sequence to the root profile. (default: false) |
| `diags1` | boolean | No | MUSCLE -diags1 (no effect: compiled out of 3.8.1551; -diags sets diags1 and diags2). (default: false) |
| `diags2` | boolean | No | MUSCLE -diags2 (no effect). (default: false) |
| `diagLength` | integer | No | MUSCLE -diaglength (default 24; no effect). (default: 24) |
| `diagBreak` | integer | No | MUSCLE -diagbreak (default 1; no effect). (default: 1) |
| `diagMargin` | integer | No | MUSCLE -diagmargin (default 5; no effect). (default: 5) |
| `hydro` | integer | No | MUSCLE -hydro (default 5; no effect). (default: 5) |
| `hydroFactor` | number | No | MUSCLE -hydrofactor (default 1.2; no effect). (default: 1.2) |
| `tomHydro` | boolean | No | MUSCLE -tomhydro (no effect). (default: false) |
| `gapOpen2` | number | No | MUSCLE -gapopen2 (no effect). |
| `gapExtend2` | number | No | MUSCLE -gapextend2 (no effect). |
| `gapAmbig` | number | No | MUSCLE -gapambig (no effect). |
| `maxMb` | integer | No | MUSCLE -maxmb (no effect on the result). |
| `maxHours` | number | No | MUSCLE -maxhours (no effect on the result). |
| `matrix` | string | No | MUSCLE -matrix: the text of a score matrix file (NCBI format) replacing the matrix of the SP score. For DNA/RNA the program reads the file under the nucleotide alphabet (the A/C/G/T entries of the matrix) and then runs with the amino-acid SP score. |
| `useTree` | string | No | MUSCLE -usetree: a rooted binary Newick tree with edge lengths, leaf labels = names; skips the first tree and the tree refinement. |
| `names` | array<string> | No | FASTA labels of the sequences (default s0, s1, ...): the leaf names of useTree and of the returned trees. |
| `clusterOnly` | boolean | No | MUSCLE -clusteronly with -tree1: return only the first guide tree (tree1), no alignment. (default: false) |
| `computeWeights` | boolean | No | MUSCLE -computeweights: stop after the progressive alignment and return the weights file instead of an alignment. (default: false) |
| `stable` | boolean | No | MUSCLE -stable: the 3.8.1551 program refuses it ("-stable not supported in this version of muscle"), so true is rejected; treeOrder = false (default) already returns the rows in input order. (default: false) |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA (input order; the order MUSCLE writes with treeOrder). Empty with computeWeights (the program then writes only the weights). |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a residue beats a gap on ties). |
| `totalScore` | integer | Natural affine sum-of-pairs score of the alignment (BLOSUM62 10 / 0.2 for proteins, EDNAFULL 10 / 0.5 for nucleotides; MUSCLE itself reports no such score), rounded. |
| `exactTotalScore` | number | The same score unrounded. |
| `names` | array<string> | FASTA labels of the rows (input order; the names parameter, default s0, s1, ...). |
| `outputOrder` | array<integer> | Indices of the input rows in the order the program writes them (treeOrder). |
| `tree1` | string | The first guide tree in the program's -tree1 text format, or null. |
| `tree2` | string | The tree of iteration 2 (-tree2), or null. |
| `distanceMatrix1` | string | The first distance matrix (-distmx1), or null. |
| `distanceMatrix2` | string | The second distance matrix (-distmx2), or null. |
| `weightsFile` | string | The sequence weights (-computeweights), or null. |
| `referenceBehaviourUndefined` | boolean | True when MUSCLE 3.8.1551's own behaviour on this input is undefined (-matrix with DNA/RNA reads never-written profile memory; the reference may crash); the alignment is the zero-memory result. |

## Errors

| Code | Message |
|------|---------|
| 1010 | sequence is empty or contains a character other than a letter, '-', '.' or white space |
| 1002 | type must be auto, protein, dna or rna; maxIters must be at least 1; scoreModel / distance1 / distance2 / cluster1 / cluster2 / root1 / root2 / weight1 / weight2 / objScore / termGaps must be a known name; scoreModel and distance1 must suit the alphabet; options the program refuses (stable, ...) are reported as errors |

## Examples

### Example 1: Five related protein sequences (maxIters 2)

**Tool Call:**
```json
{
  "tool": "multiple_align_muscle",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQPWTTQAEAKP",
      "MKVLAGIVGLLAQPWTQAEAKPLL",
      "MRVLAAGIIGLLLQPWTTQAEAK",
      "KVLAAGIVALLLAQPWTTAEAKPS",
      "MKVAAGIVGLLLAQPWTTQEAKP"
    ],
    "type": "protein"
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "MKVLAAGIVGLLLAQPWTTQAEAKP--",
    "MKV-LAGIVG-LLAQPWTQA-EAKPLL",
    "MRVLAAGIIG-LLLQPWTTQAEAK---",
    "-KVLAAGIVALLLAQPWTTA-EAKPS-",
    "MKV-AAGIVGLLLAQPWTTQ-EAKP--"
  ],
  "consensus": "MKVLAAGIVGLLLAQPWTTQ-EAKP--",
  "totalScore": 656
}
```

### Example 2: Four DNA sequences (progressive draft, maxIters 1)

**Tool Call:**
```json
{
  "tool": "multiple_align_muscle",
  "arguments": {
    "sequences": [
      "ACGTACGTTAGCTAGCTAGGATCGATCGGATC",
      "ACGTACGTAGCTAGCTAGGATCGATCGATC",
      "ACGTTCGTTAGCTAGCTAGGATCGGATCGGATC",
      "ACGTACGTTAGCTTAGCTAGGATCGATCGGATC"
    ],
    "type": "dna",
    "maxIters": 1
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTACGTTAGC-TAGCTAGGATC-GATCGGATC",
    "ACGTACG-TAGC-TAGCTAGGATC-GATC-GATC",
    "ACGTTCGTTAGC-TAGCTAGGATCGGATCGGATC",
    "ACGTACGTTAGCTTAGCTAGGATC-GATCGGATC"
  ],
  "consensus": "ACGTACGTTAGC-TAGCTAGGATC-GATCGGATC",
  "totalScore": 783
}
```

## Verification status

The installed `muscle` of the audit machine is 5.x (the PPP / probcons-like pipeline), a different program, so it cannot be the oracle. The oracle is the MUSCLE 3.8.1551 binary built from the same Debian source by the Docker hand-off package `docs/Validation/review-2026-09/handoff/A1-10c-muscle` (420 sets): **awaiting its run**. Until it is integrated the alignments are not verified against muscle 3.8; the score formulas (log-expectation / SPN, VTML200 / BLASTZ tables), the ClustalW weights example of the source, the UPGMA test matrix of the source, the Kimura / Dayhoff distance, the k-mer distance definition, the wildcard rules and the optimality of the profile dynamic program (exhaustive enumeration, 480 random profile pairs) are tested.

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.Muscle.cs#L67](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Muscle.cs#L67)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
- Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research 32(5): 1792-1797.
