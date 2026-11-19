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
