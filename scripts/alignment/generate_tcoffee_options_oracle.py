#!/usr/bin/env python3
"""Oracle generator for the command-line options of SequenceAligner.MultipleAlignTCoffeeProgram (audit A1-13c / A1-13d, B09 F71+).

Runs the installed T-Coffee program (Debian t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1, `t_coffee`) on random sets with a FLAG string and writes
  tcoffee_program_options_oracle.tsv : feature <TAB> type <TAB> flags <TAB> seq0|seq1|... <TAB> row0|row1|... <TAB> names
                                       rows in INPUT order of the sequences (in the order of the program's output file when the flags contain `-outorder aligned`);
                                       FAIL:<rc> when the program fails deterministically; names = the sequence names (empty = q000, q001, ...).
A flag value `@random` after `-seq_weight` is replaced by random weights: the program reads them from a SINGLE_SEQ_WEIGHT_FORMAT_01 file, the recorded flags carry the
in-memory form `@w1,w2,...`. Every run uses `-n_core 1 -no_warning` and a private HOME.

usage: generate_tcoffee_options_oracle.py OUT_DIR [FEATURE ...] [SCALE=1.0 via env SCALE] [SEED via env SEED]   (needs `t_coffee` on PATH)
"""
import os, random, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_tcoffee_program_oracle import make_set, DNA, PROTEIN

RNA = 'ACGU'
ALPHA = {'dna': DNA, 'protein': PROTEIN, 'rna': RNA}
MIXED_ALPHABETS = ['ACGTN', 'ACGTRYMKSWN', PROTEIN, 'ACGT', 'ACDEFGHIKLMNPQRSTVWY']

# feature -> list of (type, flags, sets, (kmin, kmax), (lmin, lmax), name mode)   name mode: '' = q000.. | 'shuffled' | 'seq' (s0, s1, ... s12: sorted as strings) | 'prefix'
PLAN = {
    'rna': [('rna', '', 260, (2, 8), (3, 36), ''),
            ('rna', '-method proba_pair', 60, (2, 7), (3, 30), ''),
            ('rna', '-method slow_pair,lalign_id_pair', 40, (2, 7), (3, 30), ''),
            ('rna', '-method proba_pair,slow_pair,lalign_id_pair', 40, (2, 6), (3, 30), ''),
            ('rna', '-method lalign_id_pair,proba_pair', 30, (2, 6), (3, 30), '')],
    'relax': [(t, f, n, (2, 8), (3, 36), '') for t in ('dna', 'protein') for f, n in (('-relax_lib 0', 40), ('-relax_lib 2', 50), ('-relax_lib 3', 30), ('-filter_lib 20', 40),
                                                                                   ('-filter_lib 100', 30), ('-relax_lib 2 -filter_lib 50', 30),
                                                                                   ('-method slow_pair,lalign_id_pair -relax_lib 2', 30))],
    'gap': [(t, f, n, (2, 8), (3, 36), '') for t in ('dna', 'protein') for f, n in (('-gapopen -50', 30), ('-gapext -2', 40), ('-gapext 1', 20), ('-cosmetic_penalty 5', 20),
                                                                                  ('-gapopen -30 -gapext -1', 30), ('-nomatch 5', 30), ('-nomatch 50 -gapext -1', 30),
                                                                                  ('-nomatch 200 -do_normalise 0', 20), ('-gapext -20 -gapopen -100', 30),
                                                                                  ('-method slow_pair,lalign_id_pair -gapext -2', 30))],
    'tree': [(t, f, n, (3, 9), (3, 36), '') for t in ('dna', 'protein') for f, n in (('-tree_mode upgma', 50), ('-distance_matrix_mode fast', 40), ('-distance_matrix_mode idscore', 30),
                                                                                  ('-distance_matrix_mode no', 20), ('-distance_matrix_mode fast -tree_mode upgma', 40),
                                                                                  ('-distance_matrix_mode idscore -tree_mode upgma', 30),
                                                                                  ('-distance_matrix_mode fast -distance_matrix_sim_mode idmat_sim2', 25),
                                                                                  ('-distance_matrix_mode fast -distance_matrix_sim_mode idmat_sim3', 25),
                                                                                  ('-distance_matrix_mode fast -distance_matrix_sim_mode idmat_gap1', 25),
                                                                                  ('-distance_matrix_mode fast -distance_matrix_sim_mode idmat_logid', 25),
                                                                                  ('-distance_matrix_mode idscore -distance_matrix_sim_mode idmat_sim3', 25),
                                                                                  ('-method slow_pair,lalign_id_pair -distance_matrix_mode fast -tree_mode upgma', 30))],
    'weight': [(t, f, n, (3, 9), (3, 36), '') for t in ('dna', 'protein') for f, n in (('-seq_weight t_coffee', 60), ('-seq_weight @random', 60),
                                                                                    ('-method slow_pair -seq_weight t_coffee', 40),
                                                                                    ('-seq_weight t_coffee -nomatch 20', 30), ('-seq_weight t_coffee -weight no', 15),
                                                                                    ('-seq_weight t_coffee -distance_matrix_mode idscore -tree_mode upgma', 30))],
    'order': [(t, f, n, (2, 8), (3, 30), nm) for t in ('dna', 'protein') for f, n, nm in (('', 40, 'shuffled'), ('-outorder aligned', 40, 'shuffled'), ('-inorder input', 40, 'shuffled'),
                                                                                       ('-inorder input -outorder aligned', 30, 'shuffled'),
                                                                                       ('', 30, 'prefix'), ('-tree_mode upgma', 20, 'shuffled'),
                                                                                       ('-method slow_pair,lalign_id_pair', 30, 'shuffled'))] +
             [(t, f, n, (11, 13), (4, 20), 'seq') for t in ('dna', 'protein') for f, n in (('', 30), ('-tree_mode upgma', 20))],
    'matrix': [(t, f, n, (2, 7), (3, 30), '') for t in ('dna', 'protein', 'rna') for f, n in (('-matrix blosum50mt', 25), ('-matrix pam250mt', 25), ('-matrix idmat', 20), ('-matrix blosum80mt', 20),
                                                                                         ('-matrix md_250mt', 20), ('-matrix blosum30mt', 20), ('-matrix est_idmat', 15),
                                                                                         ('-method slow_pair -matrix blosum50mt', 25), ('-method slow_pair -matrix pam120mt', 20),
                                                                                         ('-method slow_pair,lalign_id_pair,proba_pair -matrix blosum50mt', 20),
                                                                                         ('-matrix blosum62mt -tree_mode upgma -seq_weight t_coffee', 15))],
    'method': [(t, f, n, (2, 7), (3, 30), '') for t in ('dna', 'protein') for f, n in (('-method slow_pair@EP@GOP@-30', 25), ('-method slow_pair@EP@GEP@-3', 25),
                                                                                    ('-method slow_pair@EP@MATRIX@blosum50mt', 25),
                                                                                    ('-method slow_pair@EP@GOP@-5@GEP@-2@MATRIX@pam250mt', 25),
                                                                                    ('-method lalign_id_pair@EP@GOP@-5', 25), ('-method lalign_id_pair@EP@MATRIX@blosum62mt@GEP@-2', 25),
                                                                                    ('-method slow_pair@EP@GOP@-20,lalign_id_pair@EP@GEP@-2', 25),
                                                                                    ('-method lalign_id_pair@EP@GOP@-20@GEP@-6 -lalign_n_top 4', 20))],
    # A2-7: clustalw_pair decides the residue type of EVERY pair file itself (the type guess of the program on the SECOND sequence), not T-Coffee's type of the set:
    # sets whose sequences come from different alphabets (nucleotide-only letters like YAGR / NSMYC are protein for clustalw; a protein first sequence with a nucleotide second one is DNA)
    # A2-7: the float sums of the pair list depend on the order of the group's sequences (ascending index): large sets with a nomatch penalty that makes near-ties visible
    'floatorder': [('protein', f, n, (6, 10), (6, 30), '') for f, n in (('-nomatch 40 -gapext -1', 100), ('-nomatch 52 -gapext -1', 100), ('-nomatch 15 -gapext -2 -seq_weight t_coffee', 60))],
    # A2-3 / A2-4: further pair-wise methods of the program (biphasic_pair_wise, SIM on every pair incl. the self pair, slow_pair on the dinucleotide alphabet)
    'biphasic': [(t, f, n, (2, 7), (3, 34), '') for t in ('dna', 'protein', 'rna') for f, n in (('-method biphasic_pair', 70), ('-method biphasic_pair,slow_pair', 25), ('-method biphasic_pair@EP@GOP@-30@GEP@-2', 20),
                                                                                              ('-method proba_pair,biphasic_pair -tree_mode upgma', 15))],
    'lalign_rs_s': [(t, f, n, (2, 7), (3, 34), '') for t in ('dna', 'protein', 'rna') for f, n in (('-method lalign_rs_s_pair', 70), ('-method lalign_rs_s_pair,slow_pair', 25), ('-method lalign_rs_s_pair -lalign_n_top 3', 20),
                                                                                                  ('-method lalign_id_pair,lalign_rs_s_pair', 15))],
    'promo_blastr': [(t, f, n, (2, 7), (3, 34), '') for t in ('dna', 'protein', 'rna') for f, n in (('-method promo_pair', 45), ('-method blastr_pair', 45), ('-method promo_pair,blastr_pair', 15), ('-method slow_pair,promo_pair', 15),
                                                                                                   ('-method blastr_pair@EP@GOP@-40@GEP@-1', 15))],
    'refused': [(t, f, n, (2, 5), (3, 20), '') for t in ('dna', 'protein') for f, n in (('-method clean_slow_pair', 16), ('-method hash_pair', 16), ('-method subop1_pair', 12), ('-method subop2_pair', 12), ('-method cfast_pair', 12),
                                                                                       ('-method slow_pair,hash_pair', 12))],
    # A2-6: -distance_matrix_sim_mode containing idscore: `fast` hands the ALIGNED rows to idscore_pairseq; a gap in the FIRST row of a pair segfaults the program (exit status 11), sets without such a gap are scored
    'simidscore': [(t, f, n, (2, 6), (5, 30), '') for t in ('dna', 'protein') for f, n in (('-distance_matrix_mode fast -distance_matrix_sim_mode idscore', 60), ('-distance_matrix_mode fast -distance_matrix_sim_mode idscore_sim3', 40),
                                                                                         ('-distance_matrix_mode ktup -distance_matrix_sim_mode idscore', 15), ('-distance_matrix_mode idscore -distance_matrix_sim_mode idscore_sim2', 15))],
    'simidscore_gapless': [(t, f, n, (3, 8), (6, 30), '') for t in ('dna', 'protein') for f, n in (('-distance_matrix_mode fast -distance_matrix_sim_mode idscore', 60), ('-distance_matrix_mode fast -distance_matrix_sim_mode idscore_sim2 -tree_mode upgma', 40),
                                                                                                  ('-distance_matrix_mode fast -distance_matrix_sim_mode idscore -seq_weight t_coffee', 40))],
    # A2-6: the residue groups of -distance_matrix_sim_mode (make_group_aa: first token before '_'): simple, mafft, clustalw, polarity, vasiliky, make_all, unknown names (no group), with and without a library method
    # that already built the groups of the matrix idmat in the process
    'simgroups': [(t, f'-distance_matrix_mode fast -distance_matrix_sim_mode {m}' + extra, n, (4, 9), (6, 30), '')
                  for t, n in (('protein', 14), ('dna', 6)) for m in ('simple', 'mafft', 'clustalw', 'polarity', 'vasiliky', 'make_all', 'bogus', 'simple_sim2', 'mafft_gap1', 'vasiliky_logid', 'cov', 'blosum62mt_sim1')
                  for extra in ('',)] + [('protein', '-method slow_pair -distance_matrix_mode fast -distance_matrix_sim_mode vasiliky', 25, (4, 9), (6, 30), '')],
    # A2-6: tree / distance modes the program refuses or dies on (cwph and unknown tree modes: core dump for 3+ sequences; kmeans: no alignment written; aln: segfault; ktup2 / unknown distance mode: exit status 1)
    'treemodes': [(t, f, 30, (2, 6), (4, 24), '') for t in ('dna', 'protein') for f in ('-tree_mode cwph', '-tree_mode kmeans', '-tree_mode bogus')],
    'distrefused': [(t, f, 30, (2, 6), (4, 24), '') for t in ('dna', 'protein') for f in ('-distance_matrix_mode aln', '-distance_matrix_mode ktup2', '-distance_matrix_mode bogus')],
    # A2-2: fast_pair = fasta_gotoh_pair_wise (k-tuple diagonals, z-score threshold / -ndiag N best diagonals, banded Gotoh DP); longer sequences than the other features so that the band matters
    'fast_pair': [(t, f, n, (2, 7), lr, '') for t in ('dna', 'protein', 'rna') for lr in ((8, 70),) for f, n in (('-method fast_pair', 50), ('-method fast_pair,slow_pair', 15), ('-method fast_pair@EP@GOP@-20@GEP@-2', 15),
                                                                                       ('-method fast_pair -ndiag 3', 20), ('-method fast_pair -ndiag 10', 15), ('-method fast_pair -nomatch 5', 15),
                                                                                       ('-method proba_pair,fast_pair -tree_mode upgma', 10))] +
                 [(t, '-method fast_pair', 30, (2, 6), (3, 14), '') for t in ('dna', 'protein', 'rna')],
    # A2-6: distance modes that need the RELAXED library (slow / cscore / clustalw: the linked pair-wise DP of every pair; very_fast / -quicktree: very_fast_gotoh_pair_wise with library-scored diagonals) and geometric (fasta_pair_wise, pam250mt)
    'dist_slow': [(t, f, n, (2, 8), (6, 40), '') for t in ('dna', 'protein') for f, n in (('-distance_matrix_mode slow', 50), ('-distance_matrix_mode slow -seq_weight t_coffee', 30), ('-distance_matrix_mode slow -tree_mode upgma', 25),
                                                                                         ('-distance_matrix_mode slow -relax_lib 0', 15), ('-distance_matrix_mode slow -do_normalise 0', 15), ('-distance_matrix_mode slow -do_normalise 0 -tree_mode upgma', 15),
                                                                                         ('-distance_matrix_mode slow -distance_matrix_sim_mode idmat_sim2', 15), ('-distance_matrix_mode slow -method slow_pair,lalign_id_pair', 20),
                                                                                         ('-distance_matrix_mode slow -nomatch 20', 15), ('-distance_matrix_mode slow -gapopen -20 -gapext -5', 15))],
    'dist_cscore': [(t, f, n, (2, 8), (6, 40), '') for t in ('dna', 'protein') for f, n in (('-distance_matrix_mode cscore', 50), ('-distance_matrix_mode cscore -seq_weight t_coffee', 30), ('-distance_matrix_mode cscore -tree_mode upgma', 25),
                                                                                          ('-distance_matrix_mode cscore -relax_lib 0', 15), ('-distance_matrix_mode cscore -filter_lib 600', 15), ('-distance_matrix_mode cscore -filter_lib 2000', 8),
                                                                                          ('-distance_matrix_mode cscore -method slow_pair,lalign_id_pair', 20), ('-distance_matrix_mode cscore -nomatch 20', 15))],
    'dist_clustalw': [(t, f, n, (2, 8), (6, 40), '') for t in ('dna', 'protein') for f, n in (('-distance_matrix_mode clustalw', 50), ('-distance_matrix_mode clustalw -seq_weight t_coffee', 30), ('-distance_matrix_mode clustalw -tree_mode upgma', 25),
                                                                                            ('-distance_matrix_mode clustalw -relax_lib 0', 15), ('-distance_matrix_mode clustalw -do_normalise 0', 15),
                                                                                            ('-distance_matrix_mode clustalw -method slow_pair,lalign_id_pair', 20), ('-distance_matrix_mode clustalw -nomatch 20', 15))],
    'dist_veryfast': [(t, f, n, (2, 8), (6, 50), '') for t in ('dna', 'protein') for f, n in (('-distance_matrix_mode very_fast', 50), ('-quicktree', 40), ('-quicktree -tree_mode upgma', 25), ('-quicktree -seq_weight t_coffee', 30),
                                                                                             ('-quicktree -relax_lib 0', 15), ('-quicktree -tg_mode 0', 15), ('-quicktree -tg_mode 2', 15), ('-quicktree -nomatch 20', 15),
                                                                                             ('-distance_matrix_mode very_fast -distance_matrix_sim_mode idmat_sim2', 15), ('-quicktree -method slow_pair,lalign_id_pair', 15))],
    # B09 A3-6: `-distance_matrix_mode random` draws rand() % 100 per pair; t_coffee seeds with vsrand(0) = 4 bytes of /dev/urandom (two runs differ) unless the environment variable DEBUG_SRAND is set (srand(10), util.c): run_tc
    # sets it for these runs, so every set is a deterministic oracle row of the port's default RandomSeed = 10.
    'dist_random': [(t, f, n, (2, 8), (6, 40), '') for t in ('dna', 'protein') for f, n in (('-distance_matrix_mode random', 60), ('-distance_matrix_mode random -seq_weight t_coffee', 40), ('-distance_matrix_mode random -tree_mode upgma', 30),
                                                                                           ('-distance_matrix_mode random -method slow_pair,lalign_id_pair', 15), ('-distance_matrix_mode random -relax_lib 0', 10))],
    'dist_geometric': [(t, f, n, (2, 8), (6, 50), '') for t in ('dna', 'protein') for f, n in (('-distance_matrix_mode geometric', 60), ('-distance_matrix_mode geometric -seq_weight t_coffee', 40), ('-distance_matrix_mode geometric -tree_mode upgma', 30),
                                                                                            ('-distance_matrix_mode geometric -method slow_pair,lalign_id_pair', 15), ('-distance_matrix_mode geometric -relax_lib 0', 10))],
    # A2-4: exon_pair / exon2_pair / exon3_pair = fast_pair on the genomic (_G_) templates of `-template_file` (matrix exon2mt); `@templates` is replaced by generated templates (recorded as `-template_file @t0,t1,...`)
    'exon': [(t, f, n, (2, 6), (8, 40), '') for t in ('dna', 'protein') for f, n in (('-method exon_pair -template_file @templates', 40), ('-method exon2_pair -template_file @templates', 50), ('-method exon3_pair -template_file @templates', 50),
                                                                                     ('-method exon_pair,slow_pair -template_file @templates', 20), ('-method exon2_pair,exon3_pair -template_file @templates', 20),
                                                                                     ('-method exon2_pair@EP@GOP@-30@GEP@-2 -template_file @templates', 15), ('-method proba_pair,exon3_pair -tree_mode upgma -template_file @templates', 15))],
    'exon_sparse': [(t, f, n, (3, 7), (8, 40), '') for t in ('dna', 'protein') for f, n in (('-method exon_pair -template_file @templates', 40), ('-method exon3_pair -template_file @templates', 40), ('-method exon_pair,slow_pair -template_file @templates', 30))],
    # A2-5: -dp_mode of the progressive stage other than linked_pair_wise, with the library column scores of -extend_mode (very_fast_triplet default / slow_triplet / triplet / fast_triplet), -max_n_pair, gap penalties and -tg_mode
    'dp_gotoh': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode gotoh_pair_wise', 60), ('-dp_mode ssec_pair_wise', 10), ('-dp_mode gotoh_pair_wise_test', 30), ('-dp_mode gotoh_pair_wise -gapopen -20 -gapext -2', 30),
                                                                                         ('-dp_mode gotoh_pair_wise -gapopen -50 -gapext -5 -tg_mode 0', 20), ('-dp_mode gotoh_pair_wise -gapext -3 -tg_mode 2', 20), ('-dp_mode gotoh_pair_wise -nomatch 20', 20),
                                                                                         ('-dp_mode gotoh_pair_wise -extend_mode slow_triplet', 20), ('-dp_mode gotoh_pair_wise -extend_mode triplet -nomatch 10', 20), ('-dp_mode gotoh_pair_wise -extend_mode fast_triplet', 25),
                                                                                         ('-dp_mode gotoh_pair_wise -max_n_pair 2', 20), ('-dp_mode gotoh_pair_wise -max_n_pair 0', 20), ('-dp_mode gotoh_pair_wise -seq_weight t_coffee', 25),
                                                                                         ('-dp_mode gotoh_pair_wise -relax_lib 0', 15), ('-dp_mode gotoh_pair_wise -tree_mode upgma -distance_matrix_mode fast', 25), ('-dp_mode gotoh_pair_wise -method slow_pair,lalign_id_pair', 20),
                                                                                         ('-dp_mode gotoh_pair_wise -do_normalise 0', 10), ('-dp_mode gotoh_pair_wise -outorder aligned', 10))],
    'dp_mm': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode myers_miller_pair_wise', 60), ('-dp_mode default', 10), ('-dp_mode myers_miller_pair_wise -gapopen -20 -gapext -2', 30),
                                                                                    ('-dp_mode myers_miller_pair_wise -gapopen -50 -gapext -5', 20), ('-dp_mode myers_miller_pair_wise -nomatch 20', 20), ('-dp_mode myers_miller_pair_wise -extend_mode slow_triplet', 20),
                                                                                    ('-dp_mode myers_miller_pair_wise -extend_mode fast_triplet', 20), ('-dp_mode myers_miller_pair_wise -max_n_pair 2 -seq_weight t_coffee', 25),
                                                                                    ('-dp_mode myers_miller_pair_wise -tree_mode upgma -method slow_pair,lalign_id_pair', 20), ('-dp_mode myers_miller_pair_wise -relax_lib 0', 15))],
    'dp_fasta': [(t, f, n, (2, 8), (6, 50), '') for t in ('dna', 'protein') for f, n in (('-dp_mode fasta_pair_wise', 60), ('-dp_mode fasta_pair_wise -ndiag 5', 25), ('-dp_mode fasta_pair_wise -ndiag 20', 20), ('-dp_mode fasta_pair_wise -diag_threshold 1', 20),
                                                                                       ('-dp_mode fasta_pair_wise -gapopen -20 -gapext -2', 25), ('-dp_mode fasta_pair_wise -tg_mode 0 -gapopen -30 -gapext -3', 15), ('-dp_mode fasta_pair_wise -extend_mode slow_triplet', 15),
                                                                                       ('-dp_mode cfasta_pair_wise', 60), ('-dp_mode cfasta_pair_wise -ndiag 4', 25), ('-dp_mode cfasta_pair_wise -gapopen -20 -gapext -2', 25),
                                                                                       ('-dp_mode very_fast_pair_wise', 60), ('-dp_mode very_fast_pair_wise -gapopen -20 -gapext -2', 25), ('-dp_mode very_fast_pair_wise -seq_weight t_coffee -tree_mode upgma', 20))],
    # A2-5: biphasic_pair_wise (two affine gap models) and glocal_pair_wise (match / deletion / insertion / unmatched pair) as the dynamic program of the progressive stage
    'dp_other': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode biphasic_pair_wise', 60), ('-dp_mode biphasic_pair_wise -gapopen -20 -gapext -2', 40), ('-dp_mode biphasic_pair_wise -extend_mode slow_triplet -nomatch 10', 25),
                                                                                       ('-dp_mode glocal_pair_wise', 60), ('-dp_mode glocal_pair_wise -gapopen -20 -gapext -2', 40), ('-dp_mode glocal_pair_wise -extend_mode slow_triplet -seq_weight t_coffee', 25),
                                                                                       ('-dp_mode glocal_pair_wise -gapext -5 -tree_mode upgma', 25))],
    # A2-5 (WP30): gotoh_pair_wise_lgp = the Gotoh DP with the position specific gap opening / closing penalties of aln2local_penalties
    'dp_lgp': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode gotoh_pair_wise_lgp', 40), ('-dp_mode gotoh_pair_wise_lgp -gapopen -20 -gapext -2', 60), ('-dp_mode gotoh_pair_wise_lgp -gapopen -50 -gapext -5', 40),
                                                                                    ('-dp_mode gotoh_pair_wise_lgp -gapext -3', 20), ('-dp_mode gotoh_pair_wise_lgp -gapopen -30 -nomatch 20', 20), ('-dp_mode gotoh_pair_wise_lgp -gapopen -20 -gapext -2 -extend_mode slow_triplet', 25),
                                                                                    ('-dp_mode gotoh_pair_wise_lgp -gapopen -20 -gapext -2 -extend_mode fast_triplet', 15), ('-dp_mode gotoh_pair_wise_lgp -gapopen -20 -gapext -2 -extend_mode triplet', 15),
                                                                                    ('-dp_mode gotoh_pair_wise_lgp -gapopen -40 -gapext -4 -seq_weight t_coffee', 25), ('-dp_mode gotoh_pair_wise_lgp -gapopen -20 -gapext -2 -relax_lib 0', 15),
                                                                                    ('-dp_mode gotoh_pair_wise_lgp -gapopen -20 -gapext -2 -tree_mode upgma -distance_matrix_mode fast', 20), ('-dp_mode gotoh_pair_wise_lgp -gapopen -20 -gapext -2 -method slow_pair,lalign_id_pair', 15))],
    # A2-5 (WP30): gotoh_pair_wise_sw (LOCAL: the rows keep only the best local alignment of every merge), ktup_pair_wise (only a score, the rows stay as they are), precomputed_pair_wise (stops unless the first rows of two groups have the same length)
    'dp_sw': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode gotoh_pair_wise_sw', 60), ('-dp_mode gotoh_pair_wise_sw -gapopen -20 -gapext -2', 40), ('-dp_mode gotoh_pair_wise_sw -gapopen -50 -gapext -5', 30),
                                                                                  ('-dp_mode gotoh_pair_wise_sw -gapopen -10 -gapext -1 -nomatch 5', 25), ('-dp_mode gotoh_pair_wise_sw -gapopen -20 -gapext -2 -extend_mode slow_triplet', 25),
                                                                                  ('-dp_mode gotoh_pair_wise_sw -gapopen -20 -gapext -2 -extend_mode fast_triplet', 15), ('-dp_mode gotoh_pair_wise_sw -gapopen -20 -gapext -2 -seq_weight t_coffee', 20),
                                                                                  ('-dp_mode gotoh_pair_wise_sw -gapopen -20 -gapext -2 -relax_lib 0', 15), ('-dp_mode gotoh_pair_wise_sw -gapopen -20 -gapext -2 -tree_mode upgma', 20))],
    'dp_ktup': [(t, f, n, (2, 8), (4, 40), '') for t in ('dna', 'protein') for f, n in (('-dp_mode ktup_pair_wise', 50), ('-dp_mode ktup_pair_wise -ktuple 2', 20), ('-dp_mode ktup_pair_wise -ktuple 3', 15), ('-dp_mode ktup_pair_wise -extend_mode slow_triplet', 15),
                                                                                 ('-dp_mode ktup_pair_wise -seq_weight t_coffee -tree_mode upgma', 15), ('-dp_mode ktup_pair_wise -gapopen -20 -gapext -2', 10))],
    'dp_precomputed_equal': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode precomputed_pair_wise', 60), ('-dp_mode precomputed_pair_wise -tree_mode upgma', 30), ('-dp_mode precomputed_pair_wise -seq_weight t_coffee', 20))],
    'dp_precomputed': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode precomputed_pair_wise', 40), ('-dp_mode precomputed_pair_wise -tree_mode upgma', 20))],
    # A2-5 (WP30): test_pair_wise (no alignment: every row padded with gaps to l0 + l1 - n; a heap overflow of the row buffers kills part of the sets, a negative generate_null size the rest of the short ones)
    'dp_test': [(t, f, n, kr, lr, '') for t in ('dna', 'protein') for f, n, kr, lr in (('-dp_mode test_pair_wise', 110, (2, 8), (1, 36)), ('-dp_mode test_pair_wise', 60, (2, 4), (1, 6)), ('-dp_mode test_pair_wise -tree_mode upgma', 40, (2, 8), (3, 30)))],
    # A2-5 (WP30): the column scores pc / mixt / g_coffee / g_coffee_quadruplets of -extend_mode (slow_get_dp_cost with another evaluate_residue_pair) on generic dynamic programs
    'dp_extend2': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in tuple((f'-dp_mode {dp} -extend_mode {em}{extra}', n) for em in ('pc', 'mixt', 'g_coffee', 'g_coffee_quadruplets')
                                                                                              for dp, extra, n in (('gotoh_pair_wise', '', 30), ('gotoh_pair_wise', ' -gapopen -20 -gapext -2', 20), ('myers_miller_pair_wise', '', 15), ('gotoh_pair_wise', ' -nomatch 10', 10),
                                                                                                                   ('gotoh_pair_wise', ' -seq_weight t_coffee -tree_mode upgma', 10), ('fasta_pair_wise', '', 10)))],
    # A2-5 (WP30): procoffee_pair_wise (the sparse linked DP with the older cl2diag_cap390, the two groups swapped)
    'dp_procoffee': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode procoffee_pair_wise', 80), ('-dp_mode procoffee_pair_wise -gapopen -20 -gapext -2', 40), ('-dp_mode procoffee_pair_wise -relax_lib 0', 30),
                                                                                       ('-dp_mode procoffee_pair_wise -seq_weight t_coffee -tree_mode upgma', 30), ('-dp_mode procoffee_pair_wise -nomatch 10', 20), ('-dp_mode procoffee_pair_wise -method slow_pair,lalign_id_pair', 30))],
    # A2-5 (WP30): viterbi_pair_wise (the Viterbi path of proba_pair's pair-HMM on the frequency profiles of the two groups)
    'dp_viterbi': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode viterbi_pair_wise', 100), ('-dp_mode viterbi_pair_wise -relax_lib 0', 30), ('-dp_mode viterbi_pair_wise -seq_weight t_coffee -tree_mode upgma', 30),
                                                                                    ('-dp_mode viterbi_pair_wise -method slow_pair,lalign_id_pair', 30), ('-dp_mode viterbi_pair_wise -gapopen -20 -gapext -2', 20))],
    # A2-5 (WP30): -extend_mode test (slow_get_dp_cost_test over residue_pair_extended_list_raw: the raw sum of the smaller weights x 10 / (ns1 x ns2)) and rna / rna1-4 (no contact library given)
    'dp_extend3': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in tuple([(f'-dp_mode {dp} -extend_mode test{extra}', n) for dp, extra, n in (('gotoh_pair_wise', '', 40), ('gotoh_pair_wise', ' -gapopen -20 -gapext -2', 25), ('myers_miller_pair_wise', '', 25), ('fasta_pair_wise', '', 10), ('gotoh_pair_wise', ' -nomatch 10', 10))]
                                                                                              + [(f'-dp_mode gotoh_pair_wise -extend_mode {em}', 25) for em in ('rna', 'rna1', 'rna2', 'rna3', 'rna4')])],
    # A2-5 (WP30): -extend_mode quadruplet (get_dp_cost_quadruplet: slow_get_dp_cost for a single sequence in a group, else fast_get_dp_cost_quadruplet over residue_pair_extended_list_quadruplet)
    'dp_quad': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode gotoh_pair_wise -extend_mode quadruplet', 70), ('-dp_mode gotoh_pair_wise -extend_mode quadruplet -gapopen -20 -gapext -2', 30), ('-dp_mode myers_miller_pair_wise -extend_mode quadruplet', 30),
                                                                                  ('-dp_mode gotoh_pair_wise -extend_mode quadruplet -nomatch 10', 15), ('-dp_mode gotoh_pair_wise -extend_mode quadruplet -seq_weight t_coffee -tree_mode upgma', 20), ('-dp_mode fasta_pair_wise -extend_mode quadruplet', 15))],
    # A2-5: -extend_mode relative_triplet (residue_pair_relative_extended_list reads the list of the first residue at the positions of the second residue's list: past the end of the shorter list = heap), rna0 (= slow_triplet) and the modes with no effect on linked_pair_wise
    'dp_extend': [(t, f, n, (2, 8), (6, 36), '') for t in ('dna', 'protein') for f, n in (('-dp_mode gotoh_pair_wise -extend_mode relative_triplet', 60), ('-dp_mode myers_miller_pair_wise -extend_mode relative_triplet', 30),
                                                                                       ('-dp_mode gotoh_pair_wise -extend_mode rna0', 30), ('-extend_mode relative_triplet', 20), ('-extend_mode pc', 15), ('-extend_mode mixt -extend_mode very_fast_triplet', 15))],
    # A2-5: dynamic programs / extend modes the program refuses or dies on (exit statuses 1, 6, 11), measured; unknown names stop every run, ssp / tm / curvature as soon as there are two sequences
    'dp_refused': [(t, f, 20, (2, 5), (5, 20), '') for t in ('dna', 'protein') for f in ('-dp_mode idscore_pair_wise', '-dp_mode subop1_pair_wise', '-dp_mode subop2_pair_wise', '-dp_mode proba_pair_wise', '-dp_mode pavie_pair_wise',
                                                                                  '-dp_mode domain_pair_wise', '-dp_mode fasta_cdna_pair_wise', '-dp_mode cfasta_cdna_pair_wise', '-dp_mode gotoh_pair_wise_lalign', '-dp_mode seq_viterbi_pair_wise',
                                                                                  '-dp_mode hh_pair_wise', '-dp_mode cfasta_sw_pair_wise', '-dp_mode gotoh_pair_wise_lgp_sticky', '-dp_mode bogus_pair_wise', '-extend_mode bogus', '-extend_mode ssp', '-extend_mode tm',
                                                                                  '-extend_mode curvature', '-extend_mode matrix')],
    # A2-11: muscle_pair = the installed muscle3 (3.8.1551) on every pair (`muscle -in -out -quiet -maxmb 0`); needs the private plugin directory of run_tc
    'muscle_pair': [(t, f, n, (2, 6), (4, 36), '') for t in ('dna', 'protein') for f, n in (('-method muscle_pair', 90), ('-method muscle_pair,slow_pair', 30), ('-method proba_pair,muscle_pair', 20),
                                                                                   ('-method muscle_pair -tree_mode upgma', 20), ('-method muscle_pair -relax_lib 0', 20))],
    'muscle_mixed': [('mixed', f, n, (2, 5), (3, 24), '') for f, n in (('-method muscle_pair', 80), ('-method muscle_pair,slow_pair', 20))],
    # A2-11: mafft_pair = the installed mafft 7.505 on every pair (`mafft --anysymbol in > out`; the program sets MAFFT_BINARIES to the private plugin directory, see run_tc); the port of MAFFT is A2-19
    'mafft_pair': [(t, f, n, (2, 6), (4, 36), '') for t in ('dna', 'protein') for f, n in (('-method mafft_pair', 90), ('-method mafft_pair,slow_pair', 30), ('-method proba_pair,mafft_pair', 20),
                                                                                  ('-method mafft_pair -tree_mode upgma', 20), ('-method mafft_pair -relax_lib 0', 20))],
    'mafft_mixed': [('mixed', f, n, (2, 5), (3, 24), '') for f, n in (('-method mafft_pair', 80), ('-method mafft_pair,slow_pair', 20))],
    'clustalw_mixed': [('mixed', f, n, (2, 5), (3, 24), '') for f, n in (('-method clustalw_pair', 160), ('-method clustalw_pair,slow_pair', 40), ('-method proba_pair,clustalw_pair', 30))],
    'clustalw': [(t, f, n, (2, 6), (4, 28), '') for t in ('dna', 'protein') for f, n in (('-method clustalw_pair', 100), ('-method clustalw_pair,slow_pair', 40), ('-method proba_pair,clustalw_pair', 30),
                                                                                    ('-method clustalw_pair -tree_mode upgma', 30))],
}


def make_templates(rng, seqs, typ, none_prob):
    """A2-4: a genomic (_G_) template for every sequence: the sequence with substitutions, lower-cased stretches and inserted runs of the exon letters o b j (the program's exon2mt matrix scores them 100 against themselves);
    sets with several letters of the program's quirky set "subset_objOBJ-" so that exon2_pair / exon3_pair get library entries; a sequence has no template with probability none_prob."""
    alpha = ALPHA[typ].lower()
    out = []
    for s in seqs:
        if rng.random() < none_prob:
            out.append('')
            continue
        if rng.random() < 0.2:                # a template equal to the sequence: fix_seq_seq needs no reconciling alignment (so no aln2sim call decides the residue groups of the process)
            out.append(s)
            continue
        t = []
        for ch in s:
            r = rng.random()
            if r < 0.12: ch = rng.choice(alpha + 'objsetu')
            elif r < 0.2: continue
            if rng.random() < 0.7: ch = ch.lower()
            t.append(ch)
            if rng.random() < 0.12: t.extend(rng.choice('objOBJsetu') for _ in range(rng.randint(1, 3)))
        out.append(''.join(t) or 'o')
    return out


def run_tc(seqs, flags, names, weights, templates=None):
    with tempfile.TemporaryDirectory() as d:
        with open(d + '/s.fa', 'w') as f:
            for n, s in zip(names, seqs): f.write(f'>{n}\n{s}\n')
        if templates is not None:
            with open(d + '/tpl.fa', 'w') as f:
                for i, (n, t) in enumerate(zip(names, templates)):
                    if not t: continue
                    with open(d + f'/g{i}.fa', 'w') as g: g.write(f'>{n}\n{t}\n')
                    f.write(f'>{n} _G_ {d}/g{i}.fa\n')
            flags = flags.replace('@templates', d + '/tpl.fa')
        if weights is not None:
            with open(d + '/w.txt', 'w') as f:
                f.write('# SINGLE_SEQ_WEIGHT_FORMAT_01\n' + ''.join(f'{n} {w}\n' for n, w in zip(names, weights)))
            flags = flags.replace('@random', d + '/w.txt')
        env = dict(os.environ, HOME=d)
        if '-distance_matrix_mode random' in flags: env['DEBUG_SRAND'] = '1'      # srand(10): see 'dist_random'
        plug = d + '/.t_coffee/plugins/linux'   # the program points MAFFT_BINARIES / PATH at this private directory: link the installed helpers (mafft: /usr/lib/mafft/bin + lib/mafft; muscle = muscle3) so that mafft_pair / muscle_pair run (A2-11)
        os.makedirs(plug, exist_ok=True)
        for src in ('/usr/lib/mafft/bin', '/usr/lib/mafft/lib/mafft'):
            if os.path.isdir(src):
                for fn in os.listdir(src):
                    if not os.path.exists(plug + '/' + fn): os.symlink(src + '/' + fn, plug + '/' + fn)
        if os.path.exists('/usr/bin/muscle3') and not os.path.exists(plug + '/muscle'): os.symlink('/usr/bin/muscle3', plug + '/muscle')
        if os.path.exists('/usr/bin/legacy_blast') and not os.path.exists(plug + '/legacy_blast.pl'): os.symlink('/usr/bin/legacy_blast', plug + '/legacy_blast.pl')   # ncbi-blast+-legacy names it legacy_blast, T-Coffee looks for legacy_blast.pl (A2-9)
        for tool in ('TMalign', 'mustang'):
            if os.path.exists('/usr/bin/' + tool) and not os.path.exists(plug + '/' + tool): os.symlink('/usr/bin/' + tool, plug + '/' + tool)
        cmd = ['t_coffee', '-seq', d + '/s.fa', '-n_core', '1', '-no_warning', '-outfile', d + '/o.fa', '-output', 'fasta_aln'] + flags.split()
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=d, env=env)
        rows, cur = {}, None
        try:
            for l in open(d + '/o.fa'):
                l = l.strip()
                if l.startswith('>'): cur = l[1:].split()[0]; rows[cur] = ''
                elif cur: rows[cur] += l
            if '-outorder aligned' in flags: return '|'.join(rows.values())
            return '|'.join(rows[n] for n in names)
        except (FileNotFoundError, KeyError):
            return f'FAIL:{r.returncode}'


def make_names(rng, k, mode):
    if mode == 'shuffled': return ['n' + ''.join(rng.choice('abcdefgh') for _ in range(3)) + str(i) for i in range(k)]
    if mode == 'seq': return [f's{i}' for i in range(k)]
    if mode == 'prefix':
        names = list(dict.fromkeys(f's{i}' + ('0' if rng.random() < 0.3 else '') for i in range(k)))
        return names + [f'z{i}' for i in range(k - len(names))]
    return [f'q{i:03d}' for i in range(k)]


def job(a):
    feature, typ, flags, idx, kr, lr, nm, seed = a
    rng = random.Random(seed * 1000003 + idx * 7919 + sum(map(ord, feature + typ + flags + nm)) * 104729)
    if feature in ('simidscore_gapless', 'dp_precomputed_equal'):
        # substitution-only variants of one base sequence: equal lengths, so the Myers-Miller alignment of every pair has no gap and the program does not crash
        alpha = ALPHA[typ]; base = ''.join(rng.choice(alpha) for _ in range(rng.randint(*lr)))
        seqs = [''.join(c if rng.random() > 0.2 else rng.choice(alpha) for c in base) for _ in range(rng.randint(*kr))]
    elif typ == 'mixed':
        seqs = [''.join(rng.choice(rng.choice(MIXED_ALPHABETS)) for _ in range(rng.randint(max(2, lr[0] - 1), lr[1]))) for _ in range(rng.randint(*kr))]
    else:
        seqs = make_set(rng, ALPHA[typ], kr, lr)
    names = make_names(rng, len(seqs), nm)
    weights = [round(rng.uniform(0.05, 2.2), 3) for _ in seqs] if '@random' in flags else None
    rec_flags = flags.replace('@random', '@' + ','.join(map(str, weights))) if weights else flags
    templates = None
    if '@templates' in flags:
        templates = make_templates(rng, seqs, typ, 0.15 if 'sparse' in feature else 0.0)
        rec_flags = rec_flags.replace('@templates', '@' + ','.join(templates))
    return '\t'.join([feature, typ, rec_flags, '|'.join(seqs), run_tc(seqs, flags, names, weights, templates), '|'.join(names) if nm else ''])


def main():
    out = sys.argv[1]
    feats = sys.argv[2:] or list(PLAN)
    scale = float(os.environ.get('SCALE', '1')); seed = int(os.environ.get('SEED', '7'))
    os.makedirs(out, exist_ok=True)
    jobs = [(f, t, fl, i, kr, lr, nm, seed) for f in feats for t, fl, n, kr, lr, nm in PLAN[f] for i in range(max(1, int(n * scale)))]
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as ex:
        res = list(ex.map(job, jobs))
    path = os.path.join(out, 'tcoffee_program_options_oracle.tsv')
    new = not os.path.exists(path)
    with open(path, 'a') as f:
        if new: f.write('feature\ttype\tflags\tseqs\trows\tnames\n')
        for r in res: f.write(r + '\n')


if __name__ == '__main__':
    main()
