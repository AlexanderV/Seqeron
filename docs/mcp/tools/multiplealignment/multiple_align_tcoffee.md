# multiple_align_tcoffee

The multiple alignment of the T-Coffee 13.41 program (t_coffee; Notredame, Higgins & Heringa 2000) reproduced end to end: its own primary library built by the program's pair-wise methods - 'proba_pair' (the program's default: ProbCons pair-HMM posteriors > 0.01, at most 4 x min(length) pairs per pair of sequences), 'slow_pair' (Myers-Miller global alignment, dna_idmat 10 / 0 with gap open -10 / extend -1 for DNA, BLOSUM62 for proteins) and 'lalign_id_pair' (Huang-Miller SIM, the lalignTopAlignments best non-intersecting local alignments), several methods merged the way the program's forked library merge does (which is not their union) - then library relaxation, the program's guide tree (k-tuple distances, neighbour joining with the trifurcating root turned into a binary tree; fast neighbour joining from 100 sequences) and the sparse linked-pair-wise progressive stage. Rows equal t_coffee -seq s.fa [-method m1,m2] -output fasta_aln on 1 500+ committed DNA and protein sets (libraries pair by pair on 800+ more), for any -n_core. The type (DNA when more than 98 % of the letters are A G C T U N R Y M K S W H B V D, otherwise protein) selects the matrices, penalties and models. RNA input (a U) uses the program's separate RNA pair-HMM (5 matrix types, 2 insert state pairs, get_tot_prob). Further command-line options of the program are accepted in 'options' (a t_coffee flag string, e.g. '-tree_mode upgma -distance_matrix_mode fast -seq_weight t_coffee -relax_lib 2 -filter_lib 20 -gapopen -50 -gapext -2 -nomatch 5 -outorder aligned'): -relax_lib N, -filter_lib T, -tree_mode nj|upgma (cwph, kmeans and unknown names are the program's refusals: TCoffeeAbortedException), -distance_matrix_mode ktup|fast|idscore|no|slow|cscore|clustalw|geometric|very_fast|random (random = rand() % 100 per pair from the glibc generator seeded by randomSeed, default 10 = the program's DEBUG_SRAND run; the five before it are computed from the RELAXED library / the FASTA DP between the relaxation and the re-weighting, as the program does; -quicktree = very_fast; aln, ktup2 and unknown names are refused as the program does; random is not reproducible: the program's rand is seeded differently on every run), -distance_matrix_sim_mode (idmat_sim1/2/3, gap1, logid), -seq_weight t_coffee|no|@w1,w2,... (weights in input order), -weight no, -gapopen, -gapext, -cosmetic_penalty, -nomatch, -do_normalise, -inorder aligned|input, -outorder input|aligned, -lalign_n_top, -method; flags without influence on the rows (-evaluate_mode, -shrink_lib, -n_core, ...) are accepted and ignored; any other flag is rejected by name. The progressive stage runs other dynamic programs with -dp_mode gotoh_pair_wise (ssec_pair_wise, gotoh_pair_wise_test) | myers_miller_pair_wise (default) | fasta_pair_wise | cfasta_pair_wise | very_fast_pair_wise on the library column scores of -extend_mode very_fast_triplet (default) | fast_triplet | slow_triplet | triplet | relative_triplet | rna0, with -max_n_pair, -tg_mode, -ndiag, -diag_threshold; the dp modes and extend modes the program refuses or dies on raise TCoffeeAbortedException, the ones that only run in the program but are not ported (viterbiL / viterbiD, co_pair_wise = Clustal Omega, linked_pair_wise_collapse, sim_pair_wise_lalign; the rna contact libraries) raise an error naming them. Every option is verified row for row against the installed t_coffee on 150+ committed sets (4 300+ sets in all). 'names' = the sequence names: the program sorts the sequences by name (strcmp order) before building the library and the tree (-inorder aligned, its default) and restores the input order on output, which is the 'name-prefix quirk' (s1 / s10 trees differ); without names the names are taken to sort in input order. -matrix NAME (38 matrices of the program: blosum30..80mt, pam120..350mt, md_*, idmat, est_idmat, ...). Further pair-wise methods: 'fast_pair' (FASTA heuristic: banded Gotoh DP on the k-tuple diagonals), 'exon_pair' / 'exon2_pair' / 'exon3_pair' (Myers-Miller alignment of the genomic _G_ templates given in genomicTemplates, matrix exon2mt, expressed in the residues of the sequences), 'biphasic_pair' (biphasic_pair_wise: two affine gap models, opening x 2 / extension / 2 and opening / 2 / extension x 2), 'lalign_rs_s_pair' (lalign_id_pair on every pair INCLUDING each sequence against itself), 'promo_pair' / 'blastr_pair' (slow_pair on the dinucleotide-extended alphabet, matrices promoter_tf1 / blosumR, gap open -60 / -20; the program's shared sequence array is changed lossily by the extension, which is reproduced); -distance_matrix_sim_mode also accepts the group names of make_group_aa (simple, mafft, clustalw, polarity, vasiliky, make_all, an unknown name = no group) and idscore (a gap in the first row of an aligned pair under -distance_matrix_mode fast makes the program segfault, exit status 11). A set on which the program itself dies - lalign_id_pair with a very small gap open (e.g. 'lalign_id_pair@EP@GOP@-3': SIM's edit-script array overflows its heap chunk), the methods clean_slow_pair / hash_pair / subop1_pair / subop2_pair / cfast_pair, the idscore segfault - raises TCoffeeAbortedException (ExitStatus, Reason). Methods that run an external program (T-Coffee's side of the interface is ported, the external program is NOT reimplemented): 'muscle_pair' (MUSCLE 3.8.1551 on every pair, the library of its alignment), 'TMalign_pair' / 'mustang_pair' / 'sap_pair' and -mode 3dcoffee (structure alignment of the PDB structures given as structureTemplates: the alignment TM-align / MUSTANG wrote is given in structureAlignerOutputs; sap_pair = TM-align because SAP is not installed; a pair without two usable structures falls back to proba_pair), 'blastp_msa' / 'tblastx_msa' (the HSPs of the all-against-all BLAST / tblastx search given as blastXml become the library, weight = truncated HSP percent identity); 'dali_pair' and 'fugue_pair' stop as in the installed program (DaliLite / fugueali are not installed, TCoffeeAbortedException), 'mafft_pair' is the MAFFT 7.505 port of the pair stage (multiple_align_mafft; 460 / 460 sets of the program equal) and -mode expresso is reported as not available (expresso downloads structures from the RCSB). More dynamic programs: -dp_mode gotoh_pair_wise_lgp (position specific gap opening / closing penalties from the gaps of the two profiles), ktup_pair_wise (k-tuple score only: the rows stay as they are, unpadded), precomputed_pair_wise (stops unless the first rows of two groups have the same length), gotoh_pair_wise_sw (Smith-Waterman form: every merge keeps only the columns of the best local alignment, rows of different length; 10 % of the committed sets are undefined behaviour of the program) test_pair_wise (a debugging routine: rows padded with gaps; a second group shorter than n columns is a segmentation fault, longer rows than the longest sequence overflow the program's row buffers), procoffee_pair_wise (sparse linked DP, older diagonal cap, groups swapped) and viterbi_pair_wise (Viterbi path of the protein pair-HMM on the residue frequency profiles); extend modes pc, mixt, g_coffee, g_coffee_quadruplets, quadruplet, test and rna / rna1-4 (no contact library = slow_triplet). totalScore is the natural affine sum of pairs (EDNAFULL / BLOSUM62), not the program's objective (0 and an empty consensus when the rows have different lengths).

## Overview

| Property | Value |
|----------|-------|
| **Server** | MultipleAlignment |
| **Tool Name** | `multiple_align_tcoffee` |
| **Method ID** | `SequenceAligner.MultipleAlignTCoffeeProgram` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The multiple alignment of the T-Coffee 13.41 program (t_coffee; Notredame, Higgins & Heringa 2000) reproduced end to end: its own primary library built by the program's pair-wise methods - 'proba_pair' (the program's default: ProbCons pair-HMM posteriors > 0.01, at most 4 x min(length) pairs per pair of sequences), 'slow_pair' (Myers-Miller global alignment, dna_idmat 10 / 0 with gap open -10 / extend -1 for DNA, BLOSUM62 for proteins) and 'lalign_id_pair' (Huang-Miller SIM, the lalignTopAlignments best non-intersecting local alignments), several methods merged the way the program's forked library merge does (which is not their union) - then library relaxation, the program's guide tree (k-tuple distances, neighbour joining with the trifurcating root turned into a binary tree; fast neighbour joining from 100 sequences) and the sparse linked-pair-wise progressive stage. Rows equal t_coffee -seq s.fa [-method m1,m2] -output fasta_aln on 1 500+ committed DNA and protein sets (libraries pair by pair on 800+ more), for any -n_core. The type (DNA when more than 98 % of the letters are A G C T U N R Y M K S W H B V D, otherwise protein) selects the matrices, penalties and models. RNA input (a U) uses the program's separate RNA pair-HMM (5 matrix types, 2 insert state pairs, get_tot_prob). Further command-line options of the program are accepted in 'options' (a t_coffee flag string, e.g. '-tree_mode upgma -distance_matrix_mode fast -seq_weight t_coffee -relax_lib 2 -filter_lib 20 -gapopen -50 -gapext -2 -nomatch 5 -outorder aligned'): -relax_lib N, -filter_lib T, -tree_mode nj|upgma (cwph, kmeans and unknown names are the program's refusals: TCoffeeAbortedException), -distance_matrix_mode ktup|fast|idscore|no|slow|cscore|clustalw|geometric|very_fast|random (random = rand() % 100 per pair from the glibc generator seeded by randomSeed, default 10 = the program's DEBUG_SRAND run; the five before it are computed from the RELAXED library / the FASTA DP between the relaxation and the re-weighting, as the program does; -quicktree = very_fast; aln, ktup2 and unknown names are refused as the program does; random is not reproducible: the program's rand is seeded differently on every run), -distance_matrix_sim_mode (idmat_sim1/2/3, gap1, logid), -seq_weight t_coffee|no|@w1,w2,... (weights in input order), -weight no, -gapopen, -gapext, -cosmetic_penalty, -nomatch, -do_normalise, -inorder aligned|input, -outorder input|aligned, -lalign_n_top, -method; flags without influence on the rows (-evaluate_mode, -shrink_lib, -n_core, ...) are accepted and ignored; any other flag is rejected by name. The progressive stage runs other dynamic programs with -dp_mode gotoh_pair_wise (ssec_pair_wise, gotoh_pair_wise_test) | myers_miller_pair_wise (default) | fasta_pair_wise | cfasta_pair_wise | very_fast_pair_wise on the library column scores of -extend_mode very_fast_triplet (default) | fast_triplet | slow_triplet | triplet | relative_triplet | rna0, with -max_n_pair, -tg_mode, -ndiag, -diag_threshold; the dp modes and extend modes the program refuses or dies on raise TCoffeeAbortedException, the ones that only run in the program but are not ported (viterbiL / viterbiD, co_pair_wise = Clustal Omega, linked_pair_wise_collapse, sim_pair_wise_lalign; the rna contact libraries) raise an error naming them. Every option is verified row for row against the installed t_coffee on 150+ committed sets (4 300+ sets in all). 'names' = the sequence names: the program sorts the sequences by name (strcmp order) before building the library and the tree (-inorder aligned, its default) and restores the input order on output, which is the 'name-prefix quirk' (s1 / s10 trees differ); without names the names are taken to sort in input order. -matrix NAME (38 matrices of the program: blosum30..80mt, pam120..350mt, md_*, idmat, est_idmat, ...). Further pair-wise methods: 'fast_pair' (FASTA heuristic: banded Gotoh DP on the k-tuple diagonals), 'exon_pair' / 'exon2_pair' / 'exon3_pair' (Myers-Miller alignment of the genomic _G_ templates given in genomicTemplates, matrix exon2mt, expressed in the residues of the sequences), 'biphasic_pair' (biphasic_pair_wise: two affine gap models, opening x 2 / extension / 2 and opening / 2 / extension x 2), 'lalign_rs_s_pair' (lalign_id_pair on every pair INCLUDING each sequence against itself), 'promo_pair' / 'blastr_pair' (slow_pair on the dinucleotide-extended alphabet, matrices promoter_tf1 / blosumR, gap open -60 / -20; the program's shared sequence array is changed lossily by the extension, which is reproduced); -distance_matrix_sim_mode also accepts the group names of make_group_aa (simple, mafft, clustalw, polarity, vasiliky, make_all, an unknown name = no group) and idscore (a gap in the first row of an aligned pair under -distance_matrix_mode fast makes the program segfault, exit status 11). A set on which the program itself dies - lalign_id_pair with a very small gap open (e.g. 'lalign_id_pair@EP@GOP@-3': SIM's edit-script array overflows its heap chunk), the methods clean_slow_pair / hash_pair / subop1_pair / subop2_pair / cfast_pair, the idscore segfault - raises TCoffeeAbortedException (ExitStatus, Reason). Methods that run an external program (T-Coffee's side of the interface is ported, the external program is NOT reimplemented): 'muscle_pair' (MUSCLE 3.8.1551 on every pair, the library of its alignment), 'TMalign_pair' / 'mustang_pair' / 'sap_pair' and -mode 3dcoffee (structure alignment of the PDB structures given as structureTemplates: the alignment TM-align / MUSTANG wrote is given in structureAlignerOutputs; sap_pair = TM-align because SAP is not installed; a pair without two usable structures falls back to proba_pair), 'blastp_msa' / 'tblastx_msa' (the HSPs of the all-against-all BLAST / tblastx search given as blastXml become the library, weight = truncated HSP percent identity); 'dali_pair' and 'fugue_pair' stop as in the installed program (DaliLite / fugueali are not installed, TCoffeeAbortedException), 'mafft_pair' is the MAFFT 7.505 port of the pair stage (multiple_align_mafft; 460 / 460 sets of the program equal) and -mode expresso is reported as not available (expresso downloads structures from the RCSB). More dynamic programs: -dp_mode gotoh_pair_wise_lgp (position specific gap opening / closing penalties from the gaps of the two profiles), ktup_pair_wise (k-tuple score only: the rows stay as they are, unpadded), precomputed_pair_wise (stops unless the first rows of two groups have the same length), gotoh_pair_wise_sw (Smith-Waterman form: every merge keeps only the columns of the best local alignment, rows of different length; 10 % of the committed sets are undefined behaviour of the program) test_pair_wise (a debugging routine: rows padded with gaps; a second group shorter than n columns is a segmentation fault, longer rows than the longest sequence overflow the program's row buffers), procoffee_pair_wise (sparse linked DP, older diagonal cap, groups swapped) and viterbi_pair_wise (Viterbi path of the protein pair-HMM on the residue frequency profiles); extend modes pc, mixt, g_coffee, g_coffee_quadruplets, quadruplet, test and rna / rna1-4 (no contact library = slow_triplet). totalScore is the natural affine sum of pairs (EDNAFULL / BLOSUM62), not the program's objective (0 and an empty consensus when the rows have different lengths).

## Core Documentation Reference

- Source: [Seqeron.Genomics.TCoffee/TCoffeeAligner.Library.cs#L656](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.TCoffee/TCoffeeAligner.Library.cs#L446)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.10)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one, none empty): letters only (the program removes gaps; here they are rejected). |
| `methods` | array<string> | No | Primary-library methods in command-line order (t_coffee -method m1,m2): 'proba_pair' (default; = plain t_coffee seqs.fa), 'slow_pair', 'lalign_id_pair', 'clustalw_pair' (the ClustalW 2.1 pair alignments), 'biphasic_pair', 'lalign_rs_s_pair', 'promo_pair', 'blastr_pair', 'fast_pair', 'exon_pair', 'exon2_pair', 'exon3_pair', each optionally with the program's parameter overrides, e.g. 'slow_pair@EP@GOP@-30@GEP@-3@MATRIX@blosum50mt'; 'clean_slow_pair', 'hash_pair', 'subop1_pair', 'subop2_pair', 'cfast_pair' are refused by the program itself and raise TCoffeeAbortedException. |
| `lalignTopAlignments` | integer | No | t_coffee -lalign_n_top (default 10): alignments reported by lalign_id_pair for every pair of sequences. |
| `relaxLibrary` | boolean | No | t_coffee -relax_lib: true = 1 (default, library relaxation), false = 0. |
| `type` | string | No | 'auto' (default; the program's decision), 'dna' or 'protein' (forces the matrices, penalties and models). |
| `options` | string | No | Further t_coffee flags, e.g. '-tree_mode upgma -relax_lib 2 -seq_weight t_coffee -gapext -2' (see the tool description for the supported list); applied after the other parameters. |
| `names` | array<string> | No | Sequence names (one per sequence, in input order): the program sorts the sequences by name before computing anything (-inorder aligned); omitted = names that sort in input order. |
| `genomicTemplates` | array<string> | No | The genomic (_G_) template sequence of every input sequence for the methods exon_pair / exon2_pair / exon3_pair (t_coffee -template_file: one entry per sequence in input order, null / empty = the sequence has no template; only pairs whose two sequences have one are aligned). Omitted = no templates. |
| `structureTemplates` | array<string> | No | The sequence of the PDB structure (a _P_ template) of every input sequence for the structure methods TMalign_pair / mustang_pair / sap_pair (one entry per sequence in input order, null / empty = no structure). A structure is attached only when it is alike enough (identity > 35 %, coverage > 50 %); a pair without two structures is aligned by proba_pair, as the program does. |
| `structureAlignerOutputs` | array<object> | No | What the EXTERNAL structure aligner (TM-align / MUSTANG, not reimplemented) wrote for pairs of structures: method 'TMalign_pair', 'mustang_pair' or 'sap_pair', the input indices of the two sequences in the order the program calls it (names sorted), and the output (TMalign: its standard output, the last four lines matter; mustang: the content of results.afasta; sap_pair: the content of the file <template1>_<template2>.sap_results = the raw output of SAP 1.1.3). A missing pair = the aligner failed (no library entries). Without any 'sap_pair' item the method sap_pair uses the TMalign entries (SAP not installed: the program switches to TM-align); with at least one 'sap_pair' item SAP is installed and the program's own sap_pair glue reads the SAP files (a pair of two usable structures without a file stops the program, as it did in the oracle). |
| `blastXml` | string | No | The XML (blastall -p blastp -m7 / NCBI BLAST+ -outfmt 5) of the search of every sequence against a local BLAST database of all of them, one Iteration per sequence in order, for the methods blastp_msa and tblastx_msa (the BLAST search itself is the external program and is not reimplemented). The hits are matched by the names (see names; default q000, q001, ...). |
| `format` | string | No | Optional output format; the alignment written in that format is returned in the result field 'formatted'. Biopython Bio.Align names: fasta, a2m, clustal, phylip, nexus, stockholm (byte for byte as Bio.Align.write; 'emboss:' + name selects the EMBOSS writer of the same name). EMBOSS -aformat names: pair, simple, srspair, srs, markx0, markx1, markx2, markx3, markx10, match, score, sam, tcoffee (two rows) and msf, mega, meganon, nexusnon, paup, paupnon, phylipnon, selex, treecon (any number of rows; EMBOSS prints the run date: now, UTC). Omit for no formatted text. |
| `randomSeed` | integer | No | The srand seed of -distance_matrix_mode random (rand() % 100 per pair, glibc rand). The program seeds from /dev/urandom (two runs differ) unless the environment variable DEBUG_SRAND is set, which is srand(10): the default 10 is that deterministic run (verified on 310 sets). Ignored by every other mode. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA (input order). |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a residue beats a gap on ties). |
| `totalScore` | integer | Natural affine sum-of-pairs score of the alignment (EDNAFULL 10 / 0.5 for DNA, BLOSUM62 10 / 0.2 for proteins; the program's own objective is not reported). |
| `formatted` | string | No | The alignment written in the requested format (null when format is omitted). |

## Errors

| Code | Message |
|------|---------|
| 1010 | a sequence is empty or contains a character other than a letter |
| 1002 | methods must be proba_pair, slow_pair or lalign_id_pair; type must be auto, dna or protein; lalignTopAlignments must be at least 1; options contains a t_coffee flag or value that is not supported (named in the message); names does not have one entry per sequence |

## Examples

### Example 1: Four DNA sequences, `-method slow_pair,lalign_id_pair`

**Tool Call:**
```json
{
  "tool": "multiple_align_tcoffee",
  "arguments": {
    "sequences": [
      "ACGTTGCATGCCGATAGC",
      "ACGTGCATGCGATAGGC",
      "ACTTGCATCCGATTGC",
      "CGTTGCATGGCCGATAGC"
    ],
    "methods": [
      "slow_pair",
      "lalign_id_pair"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "ACGTTGCAT-GCCGATAG-C",
    "ACGT-GCAT-GC-GATAGGC",
    "AC-TTGCAT--CCGATTG-C",
    "-CGTTGCATGGCCGATAG-C"
  ],
  "consensus": "ACGTTGCAT-GCCGATAG-C",
  "totalScore": 238
}
```

### Example 2: Five related protein sequences, the program's default (`t_coffee seqs.fa` = proba_pair)

**Tool Call:**
```json
{
  "tool": "multiple_align_tcoffee",
  "arguments": {
    "sequences": [
      "MKVLAAGIVGLLLAQPWTTQAEAKP",
      "MKVLAGIVGLLAQPWTQAEAKPLL",
      "MRVLAAGIIGLLLQPWTTQAEAK",
      "KVLAAGIVALLLAQPWTTAEAKPS",
      "MKVAAGIVGLLLAQPWTTQEAKP"
    ]
  }
}
```

**Response:**
```json
{
  "alignedSequences": [
    "MKVLAAGIVGLLLAQPWTTQAEAK-P-",
    "MK-VLAGIVGLL-AQPWTQ-AEAKPLL",
    "MRVLAAGIIGLL-LQPWTTQAEA--K-",
    "KV-LAAGIVALLLAQPWTT-AEAKPS-",
    "MK-VAAGIVGLLLAQPWTT-QEAK-P-"
  ],
  "consensus": "MK-LAAGIVGLLLAQPWTT-AEAK-P-",
  "totalScore": 559
}
```

## Verification status

Rows equal the installed `t_coffee` 13.41 (Debian t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1, `-n_core 1`) on 1 480 committed DNA and protein sets of 2-8 sequences (the plain default, `slow_pair,lalign_id_pair`, the single methods, the reverse order and a three-method merge: guide tree and every row), on two sets of 100 / 101 sequences (fast neighbour joining) and the primary libraries pair by pair, weight by weight on 820 further sets; `scripts/alignment/generate_tcoffee_program_oracle.py`. Not available: `proba_pair` of RNA, the other `-method`s of the program (the methods the tool supports are listed in its description above).

## References

- Algorithm source: [Seqeron.Genomics.TCoffee/TCoffeeAligner.Library.cs#L656](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.TCoffee/TCoffeeAligner.Library.cs#L446)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
- Notredame, C., Higgins, D. G., Heringa, J. (2000). T-Coffee: A novel method for fast and accurate multiple sequence alignment. J Mol Biol 302(1): 205-217.
