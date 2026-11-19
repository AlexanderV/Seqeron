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
    'clustalw_mixed': [('mixed', f, n, (2, 5), (3, 24), '') for f, n in (('-method clustalw_pair', 160), ('-method clustalw_pair,slow_pair', 40), ('-method proba_pair,clustalw_pair', 30))],
    'clustalw': [(t, f, n, (2, 6), (4, 28), '') for t in ('dna', 'protein') for f, n in (('-method clustalw_pair', 100), ('-method clustalw_pair,slow_pair', 40), ('-method proba_pair,clustalw_pair', 30),
                                                                                    ('-method clustalw_pair -tree_mode upgma', 30))],
}


def run_tc(seqs, flags, names, weights):
    with tempfile.TemporaryDirectory() as d:
        with open(d + '/s.fa', 'w') as f:
            for n, s in zip(names, seqs): f.write(f'>{n}\n{s}\n')
        if weights is not None:
            with open(d + '/w.txt', 'w') as f:
                f.write('# SINGLE_SEQ_WEIGHT_FORMAT_01\n' + ''.join(f'{n} {w}\n' for n, w in zip(names, weights)))
            flags = flags.replace('@random', d + '/w.txt')
        env = dict(os.environ, HOME=d)
        cmd = ['t_coffee', '-seq', d + '/s.fa', '-n_core', '1', '-no_warning', '-outfile', d + '/o.fa', '-output', 'fasta_aln'] + flags.split()
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=d, env=env)
        rows, cur = {}, None
        try:
            for l in open(d + '/o.fa'):
                l = l.strip()
                if l.startswith('>'): cur = l[1:]; rows[cur] = ''
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
    if feature == 'simidscore_gapless':
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
    return '\t'.join([feature, typ, rec_flags, '|'.join(seqs), run_tc(seqs, flags, names, weights), '|'.join(names) if nm else ''])


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
