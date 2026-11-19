# multiple_align_tcoffee

The multiple alignment of the T-Coffee 13.41 program (t_coffee; Notredame, Higgins & Heringa 2000) reproduced end to end: its own primary library built by the program's pair-wise methods - 'proba_pair' (the program's default: ProbCons pair-HMM posteriors > 0.01, at most 4 x min(length) pairs per pair of sequences), 'slow_pair' (Myers-Miller global alignment, dna_idmat 10 / 0 with gap open -10 / extend -1 for DNA, BLOSUM62 for proteins) and 'lalign_id_pair' (Huang-Miller SIM, the lalignTopAlignments best non-intersecting local alignments), several methods merged the way the program's forked library merge does (which is not their union) - then library relaxation, the program's guide tree (k-tuple distances, neighbour joining with the trifurcating root turned into a binary tree; fast neighbour joining from 100 sequences) and the sparse linked-pair-wise progressive stage. Rows equal t_coffee -seq s.fa [-method m1,m2] -output fasta_aln on 1 500+ committed DNA and protein sets (libraries pair by pair on 800+ more), for any -n_core. The type (DNA when more than 98 % of the letters are A G C T U N R Y M K S W H B V D, otherwise protein) selects the matrices, penalties and models. RNA input (a U) uses the program's separate RNA pair-HMM (5 matrix types, 2 insert state pairs, get_tot_prob). Further command-line options of the program are accepted in 'options' (a t_coffee flag string, e.g. '-tree_mode upgma -distance_matrix_mode fast -seq_weight t_coffee -relax_lib 2 -filter_lib 20 -gapopen -50 -gapext -2 -nomatch 5 -outorder aligned'): -relax_lib N, -filter_lib T, -tree_mode nj|upgma (cwph, kmeans and unknown names are the program's refusals: TCoffeeAbortedException), -distance_matrix_mode ktup|fast|idscore|no (aln, ktup2 and unknown names are refused as the program does; slow / cscore / clustalw / geometric are not ported yet), -distance_matrix_sim_mode (idmat_sim1/2/3, gap1, logid), -seq_weight t_coffee|no|@w1,w2,... (weights in input order), -weight no, -gapopen, -gapext, -cosmetic_penalty, -nomatch, -do_normalise, -inorder aligned|input, -outorder input|aligned, -lalign_n_top, -method; flags without influence on the rows (-extend_mode, -evaluate_mode, -shrink_lib, -tg_mode, -ktuple, -n_core, ...) are accepted and ignored; any other flag is rejected by name. Every option is verified row for row against the installed t_coffee on 150+ committed sets (4 300+ sets in all). 'names' = the sequence names: the program sorts the sequences by name (strcmp order) before building the library and the tree (-inorder aligned, its default) and restores the input order on output, which is the 'name-prefix quirk' (s1 / s10 trees differ); without names the names are taken to sort in input order. -matrix NAME (38 matrices of the program: blosum30..80mt, pam120..350mt, md_*, idmat, est_idmat, ...). Further pair-wise methods (audit A2-3 / A2-4): 'biphasic_pair' (biphasic_pair_wise: two affine gap models, opening x 2 / extension / 2 and opening / 2 / extension x 2), 'lalign_rs_s_pair' (lalign_id_pair on every pair INCLUDING each sequence against itself), 'promo_pair' / 'blastr_pair' (slow_pair on the dinucleotide-extended alphabet, matrices promoter_tf1 / blosumR, gap open -60 / -20; the program's shared sequence array is changed lossily by the extension, which is reproduced); -distance_matrix_sim_mode also accepts the group names of make_group_aa (simple, mafft, clustalw, polarity, vasiliky, make_all, an unknown name = no group) and idscore (a gap in the first row of an aligned pair under -distance_matrix_mode fast makes the program segfault, exit status 11). A set on which the program itself dies - lalign_id_pair with a very small gap open (e.g. 'lalign_id_pair@EP@GOP@-3': SIM's edit-script array overflows its heap chunk), the methods clean_slow_pair / hash_pair / subop1_pair / subop2_pair / cfast_pair, the idscore segfault - raises TCoffeeAbortedException (ExitStatus, Reason). totalScore is the natural affine sum of pairs (EDNAFULL / BLOSUM62), not the program's objective.

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `multiple_align_tcoffee` |
| **Method ID** | `SequenceAligner.MultipleAlignTCoffeeProgram` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

The multiple alignment of the T-Coffee 13.41 program (t_coffee; Notredame, Higgins & Heringa 2000) reproduced end to end: its own primary library built by the program's pair-wise methods - 'proba_pair' (the program's default: ProbCons pair-HMM posteriors > 0.01, at most 4 x min(length) pairs per pair of sequences), 'slow_pair' (Myers-Miller global alignment, dna_idmat 10 / 0 with gap open -10 / extend -1 for DNA, BLOSUM62 for proteins) and 'lalign_id_pair' (Huang-Miller SIM, the lalignTopAlignments best non-intersecting local alignments), several methods merged the way the program's forked library merge does (which is not their union) - then library relaxation, the program's guide tree (k-tuple distances, neighbour joining with the trifurcating root turned into a binary tree; fast neighbour joining from 100 sequences) and the sparse linked-pair-wise progressive stage. Rows equal t_coffee -seq s.fa [-method m1,m2] -output fasta_aln on 1 500+ committed DNA and protein sets (libraries pair by pair on 800+ more), for any -n_core. The type (DNA when more than 98 % of the letters are A G C T U N R Y M K S W H B V D, otherwise protein) selects the matrices, penalties and models. RNA input (a U) uses the program's separate RNA pair-HMM (5 matrix types, 2 insert state pairs, get_tot_prob). Further command-line options of the program are accepted in 'options' (a t_coffee flag string, e.g. '-tree_mode upgma -distance_matrix_mode fast -seq_weight t_coffee -relax_lib 2 -filter_lib 20 -gapopen -50 -gapext -2 -nomatch 5 -outorder aligned'): -relax_lib N, -filter_lib T, -tree_mode nj|upgma (cwph, kmeans and unknown names are the program's refusals: TCoffeeAbortedException), -distance_matrix_mode ktup|fast|idscore|no (aln, ktup2 and unknown names are refused as the program does; slow / cscore / clustalw / geometric are not ported yet), -distance_matrix_sim_mode (idmat_sim1/2/3, gap1, logid), -seq_weight t_coffee|no|@w1,w2,... (weights in input order), -weight no, -gapopen, -gapext, -cosmetic_penalty, -nomatch, -do_normalise, -inorder aligned|input, -outorder input|aligned, -lalign_n_top, -method; flags without influence on the rows (-extend_mode, -evaluate_mode, -shrink_lib, -tg_mode, -ktuple, -n_core, ...) are accepted and ignored; any other flag is rejected by name. Every option is verified row for row against the installed t_coffee on 150+ committed sets (4 300+ sets in all). 'names' = the sequence names: the program sorts the sequences by name (strcmp order) before building the library and the tree (-inorder aligned, its default) and restores the input order on output, which is the 'name-prefix quirk' (s1 / s10 trees differ); without names the names are taken to sort in input order. -matrix NAME (38 matrices of the program: blosum30..80mt, pam120..350mt, md_*, idmat, est_idmat, ...). Further pair-wise methods (audit A2-3 / A2-4): 'biphasic_pair' (biphasic_pair_wise: two affine gap models, opening x 2 / extension / 2 and opening / 2 / extension x 2), 'lalign_rs_s_pair' (lalign_id_pair on every pair INCLUDING each sequence against itself), 'promo_pair' / 'blastr_pair' (slow_pair on the dinucleotide-extended alphabet, matrices promoter_tf1 / blosumR, gap open -60 / -20; the program's shared sequence array is changed lossily by the extension, which is reproduced); -distance_matrix_sim_mode also accepts the group names of make_group_aa (simple, mafft, clustalw, polarity, vasiliky, make_all, an unknown name = no group) and idscore (a gap in the first row of an aligned pair under -distance_matrix_mode fast makes the program segfault, exit status 11). A set on which the program itself dies - lalign_id_pair with a very small gap open (e.g. 'lalign_id_pair@EP@GOP@-3': SIM's edit-script array overflows its heap chunk), the methods clean_slow_pair / hash_pair / subop1_pair / subop2_pair / cfast_pair, the idscore segfault - raises TCoffeeAbortedException (ExitStatus, Reason). totalScore is the natural affine sum of pairs (EDNAFULL / BLOSUM62), not the program's objective.

## Core Documentation Reference

- Source: [Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeLibrary.cs#L446](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeLibrary.cs#L446)
- Algorithm doc: [Multiple_Sequence_Alignment.md](../../../algorithms/Alignment/Multiple_Sequence_Alignment.md) (section 4.10)

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequences` | array<string> | Yes | Sequences to align (at least one, none empty): letters only (the program removes gaps; here they are rejected). |
| `methods` | array<string> | No | Primary-library methods in command-line order (t_coffee -method m1,m2): 'proba_pair' (default; = plain t_coffee seqs.fa), 'slow_pair', 'lalign_id_pair', 'clustalw_pair' (the ClustalW 2.1 pair alignments), 'biphasic_pair', 'lalign_rs_s_pair', 'promo_pair', 'blastr_pair', each optionally with the program's parameter overrides, e.g. 'slow_pair@EP@GOP@-30@GEP@-3@MATRIX@blosum50mt'; 'clean_slow_pair', 'hash_pair', 'subop1_pair', 'subop2_pair', 'cfast_pair' are refused by the program itself and raise TCoffeeAbortedException. |
| `lalignTopAlignments` | integer | No | t_coffee -lalign_n_top (default 10): alignments reported by lalign_id_pair for every pair of sequences. (default: 10) |
| `relaxLibrary` | boolean | No | t_coffee -relax_lib: true = 1 (default, library relaxation), false = 0. (default: true) |
| `type` | string | No | 'auto' (default; the program's decision), 'dna' or 'protein' (forces the matrices, penalties and models). (default: auto) |
| `options` | string | No | Further t_coffee flags, e.g. '-tree_mode upgma -relax_lib 2 -seq_weight t_coffee -gapext -2' (see the tool description for the supported list); applied after the other parameters. |
| `names` | array<string> | No | Sequence names (one per sequence, in input order): the program sorts the sequences by name before computing anything (-inorder aligned); omitted = names that sort in input order. |

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
| `alignedSequences` | array<string> | The aligned rows of the MSA (input order). |
| `consensus` | string | Majority-vote consensus of the aligned rows (gaps vote; a residue beats a gap on ties). |
| `totalScore` | integer | Natural affine sum-of-pairs score of the alignment (EDNAFULL 10 / 0.5 for DNA, BLOSUM62 10 / 0.2 for proteins; the program's own objective is not reported). |

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

Rows equal the installed `t_coffee` 13.41 (Debian t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1, `-n_core 1`) on 1 480 committed DNA and protein sets of 2-8 sequences (the plain default, `slow_pair,lalign_id_pair`, the single methods, the reverse order and a three-method merge: guide tree and every row), on two sets of 100 / 101 sequences (fast neighbour joining) and the primary libraries pair by pair, weight by weight on 820 further sets; `scripts/alignment/generate_tcoffee_program_oracle.py`. Not available: `proba_pair` of RNA, the other `-method`s of the program (open items A1-13c / A1-13d).

## References

- Algorithm source: [Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeLibrary.cs#L446](../../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.TCoffeeLibrary.cs#L446)
- Binding: [AlignmentTools.cs](../../../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment/Tools/AlignmentTools.cs)
- Notredame, C., Higgins, D. G., Heringa, J. (2000). T-Coffee: A novel method for fast and accurate multiple sequence alignment. J Mol Biol 302(1): 205-217.
