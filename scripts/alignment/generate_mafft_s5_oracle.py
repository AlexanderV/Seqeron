#!/usr/bin/env python3
"""Builds the oracles of the remaining options of the script mafft (B09 WP72, A2-19c-s5): --skipiterate, --core / --coreext, ... One group per option family; every case is one run of the installed
script `mafft` (Debian 7.505-1), repeated twice (a case whose runs differ is dropped); the oracle holds return code, standard output, standard error and every file the run wrote.

It reuses the infrastructure of WP70 (generate_mafft_command_oracle.py: sets, option runner; mafft_command_oracle_lib.py).

Usage: python3 scripts/alignment/generate_mafft_s5_oracle.py [group ...]   (no group = all)   env MC_SCALE (default 1), MC_OUT (JSONL path), MC_TSV (plain TSV), MC_THREADS
"""
import hashlib, json, os, random, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L
import generate_mafft_command_oracle as C

SCALE = float(os.environ.get('MC_SCALE', '1'))
OUT = os.environ.get('MC_OUT', '/tmp/mafft_s5_oracle.jsonl')
THREADS = int(os.environ.get('MC_THREADS', '8'))


def build(rng, want):
    out = []
    basic = [('small', 14), ('tiny', 3), ('pair', 2), ('gapped', 3), ('odd', 4), ('mid', 1)]
    n = int(210 * SCALE)

    def add(gid, opts, mix=basic, count=None, anysymbol=None):
        if want and gid not in want:
            return
        out.extend(C.cases(gid, opts, mix, count or n, rng, anysymbol))

    thr = ['0', '0.05', '0.1', '0.2', '0.3', '0.5', '0.7', '0.9', '0.99', '1', '1.5', '-0.5', '2e-1']
    iters = [['--maxiterate', '1'], ['--maxiterate', '2'], ['--maxiterate', '1000'], ['--maxiterate', '5', '--retree', '1'], ['--maxiterate', '2', '--nofft'], ['--maxiterate', '3', '--treeout'],
             ['--localpair', '--maxiterate', '2'], ['--globalpair', '--maxiterate', '1000'], ['--genafpair', '--maxiterate', '2'], ['--maxiterate', '2', '--anysymbol'], ['--maxiterate', '2', '--weighti', '0', '--localpair'], ['--maxiterate', '2', '--reorder'], ['--maxiterate', '2', '--distout'], ['--maxiterate', '3', '--clustalout'], ['--maxiterate', '2', '--quiet']]
    add('skipiterate', [['--skipiterate', t] + i for t in thr for i in iters] + [['--skipiterate', t] for t in thr[:6]],
        [('small', 14), ('tiny', 2), ('pair', 2), ('gapped', 3), ('odd', 2), ('mid', 4)])
    # --core / --coreext / --corethr / --corewin: setcore reads the file hat2 (the iterative methods, --distout)
    core_bases = [['--core', '--maxiterate', '2'], ['--core', '--maxiterate', '1000'], ['--core', '--maxiterate', '3', '--retree', '1'], ['--core', '--localpair', '--maxiterate', '2'], ['--core', '--globalpair', '--maxiterate', '2'],
                  ['--core', '--genafpair', '--maxiterate', '2'], ['--core', '--distout'], ['--core', '--localpair', '--distout'], ['--core', '--maxiterate', '2', '--nofft'], ['--core'], ['--core', '--localpair'], ['--core', '--retree', '1'],
                  ['--core', '--parttree'], ['--core', '--maxiterate', '2', '--weighti', '0', '--localpair']]
    core_extra = [[], [], [], ['--coreext'], ['--corewin', '5'], ['--corewin', '10', '--coreext'], ['--corewin', '20', '--corethr', '0.3'], ['--corethr', '0.9'], ['--corethr', '0.1', '--corewin', '8'], ['--corewin', '1'],
                  ['--corewin', '0'], ['--corethr', '0'], ['--corethr', '1'], ['--clustalout'], ['--phylipout'], ['--reorder'], ['--anysymbol'], ['--outnum'], ['--quiet'], ['--corewin', '7', '--clustalout', '--coreext'],
                  ['--namelength', '5', '--clustalout'], ['--nuc'], ['--amino'], ['--inputorder'], ['--corewin', '3', '--reorder'], ['--corewin', 'x'], ['--corethr', 'abc'], ['--corewin', '12', '--phylipout', '--coreext']]
    add('core', [b + e for b in core_bases for e in core_extra], [('small', 10), ('tiny', 2), ('pair', 2), ('gapped', 3), ('odd', 2), ('mid', 6)])
    # --aamatrix: a user-defined 20 x 20 amino-acid matrix (the lower half, the letters in the first line in any order, optionally "frequency" or "norescale")
    if not want or 'aamatrix' in want:
        letters = 'ARNDCQEGHILKMFPSTWYV'

        def matrix_text(rng):
            order = list(letters)
            if rng.random() < 0.7:
                rng.shuffle(order)
            style = rng.random()
            m = [[0] * 20 for _ in range(20)]
            for i in range(20):
                for j in range(i + 1):
                    m[i][j] = rng.randint(-4, 3) if i != j else rng.randint(4, 11)
            lines = []
            if rng.random() < 0.5:
                lines.append('# user-defined matrix')
            hdr = ' '.join(order) if style < 0.6 else ('  '.join(order) if style < 0.85 else ''.join(order))
            if rng.random() < 0.15:
                hdr = hdr.lower()
            lines.append('   ' + hdr)
            for i in range(20):
                vals = []
                for j in range(i + 1):
                    v = m[i][j]
                    vals.append(str(v) if rng.random() < 0.75 else ('%.1f' % v))
                lines.append(order[i] + ' ' + ' '.join(vals))
                if rng.random() < 0.05:
                    lines.append('# comment line')
            tail = rng.random()
            if tail < 0.25:
                fr = [rng.randint(1, 9) / 100.0 for _ in range(20)]
                lines.append('frequency ' + ' '.join('%.2f' % x for x in fr))
            elif tail < 0.4:
                lines.append('norescale')
            text = '\n'.join(lines) + ('\n' if rng.random() < 0.9 else '')
            if rng.random() < 0.1:
                text = text.replace('\n', '\r\n')
            return text, order

        am_opts = [[], [], ['--maxiterate', '2'], ['--localpair'], ['--globalpair', '--maxiterate', '2'], ['--genafpair'], ['--fmodel'], ['--op', '2', '--ep', '0.1'], ['--retree', '1'], ['--treeout'], ['--anysymbol'],
                   ['--clustalout'], ['--quiet'], ['--parttree'], ['--retree', '0', '--treeout'], ['--maxiterate', '1000', '--nofft'], ['--bl', '62'], ['--6merpair'], ['--unweight'], ['--reorder'], ['--distout'], ['--auto'], ['--large']]
        mixa = [('small', 12), ('tiny', 2), ('pair', 2), ('gapped', 3), ('odd', 3), ('mid', 3)]
        for k in range(n):
            text, order = matrix_text(rng)
            kind = rng.random()
            if kind < 0.06:
                order2 = order[:]
                text = text.replace(order2[3], '', 1)  # a letter is missing in the first line
            elif kind < 0.1:
                text = '\n'.join(text.split('\n')[:rng.randint(2, 12)]) + '\n'  # too few rows
            sname = rng.choices([m[0] for m in mixa], [m[1] for m in mixa])[0]
            seqs, itext = C.make_input(rng, C.SETS[sname], rng.choice([None, None, 30, 61]))
            v = ['--aamatrix', 'mtx.txt'] + list(rng.choice(am_opts))
            if sname == 'odd' and '--anysymbol' not in v:
                v = ['--anysymbol'] + v
            out.append(('aamatrix-%04d' % k, 'aamatrix', v, itext, {'mtx.txt': text}))
    # --add / --addprofile with --maxiterate (dvtditr -K n: the branches that cut the original alignment are not refined)
    if not want or 'additer' in want:
        variants = [['--maxiterate', '2'], ['--maxiterate', '1000'], ['--maxiterate', '1', '--retree', '1'], ['--maxiterate', '3', '--nofft'], ['--maxiterate', '2', '--anysymbol'], ['--maxiterate', '2', '--treeout'], ['--maxiterate', '2', '--reorder'],
                    ['--maxiterate', '2', '--clustalout'], ['--maxiterate', '2', '--distout'], ['--maxiterate', '2', '--op', '2', '--ep', '0.1'], ['--maxiterate', '5', '--skipiterate', '0.3'], ['--maxiterate', '2', '--bunkatsunashi'],
                    ['--maxiterate', '2', '--keeplength'], ['--maxiterate', '2', '--mapout'], ['--maxiterate', '2', '--6merpair'], ['--maxiterate', '2', '--unweight'], ['--maxiterate', '2', '--quiet']]
        for k in range(n):
            v = list(rng.choice(variants))
            aligned, new = C.aligned_family(rng, rng.randint(2, 9))
            profile = rng.random() < 0.15
            if '--anysymbol' in v and rng.random() < 0.5:
                new[0] = new[0][:3] + 'X' + new[0][3:]
            opt = '--addprofile' if profile else '--add'
            files = {'new.fa': L.fasta(new, ['n%d' % i for i in range(len(new))])}
            if profile:
                al2, nw2 = C.aligned_family(rng, rng.randint(2, 5))
                files = {'new.fa': al2}
            out.append(('additer-%04d' % k, 'additer', [opt, 'new.fa'] + v, aligned, files))
    # --addprofile (disttbfast -I n): the sequences of new.fa are an alignment that must form a group of the guide tree
    if not want or 'addprofile' in want:
        variants = [[], [], ['--retree', '1'], ['--retree', '3'], ['--6merpair'], ['--nofft'], ['--op', '2', '--ep', '0.1'], ['--anysymbol'], ['--reorder'], ['--clustalout'], ['--treeout', '--distout'], ['--keeplength'], ['--mapout'],
                    ['--bl', '80'], ['--unweight'], ['--quiet'], ['--outnum'], ['--inputorder']]
        for k in range(n):
            v = list(rng.choice(variants))
            aligned, new = C.aligned_family(rng, rng.randint(2, 9))
            if len(new) > 1 and rng.random() < 0.6:
                rc, o, e, f2 = L.run_once(['--quiet'], L.fasta(new, ['n%d' % i for i in range(len(new))]))
                prof = o
            elif rng.random() < 0.5:
                al2, nw2 = C.aligned_family(rng, rng.randint(2, 5))
                prof = al2
            else:
                prof = L.fasta(new, ['n%d' % i for i in range(len(new))])  # not aligned: the program stops
            if '--anysymbol' in v and rng.random() < 0.5:
                prof = prof.replace('>n0\n', '>n0\n', 1)
            out.append(('addprofile-%04d' % k, 'addprofile', ['--addprofile', 'new.fa'] + v, aligned, {'new.fa': prof}))
    return out


def main():
    rng = random.Random(7250816)
    want = set(sys.argv[1:])
    allcases = build(rng, want)
    n = len(allcases)
    dropped = []
    done = 0
    tsv = os.environ.get('MC_TSV')
    ft = open(tsv, 'w', encoding='utf-8', newline='\n') if tsv else None
    full = 0
    with ThreadPoolExecutor(THREADS) as ex, open(OUT, 'w') as fo:
        for res, bad in ex.map(C.run_case, allcases):
            done += 1
            if res is None:
                dropped.append(bad)
                continue
            fo.write(json.dumps(res) + '\n')
            if ft:
                sha = '-' if res['nostderr'] else hashlib.sha256(res['err'].encode('utf-8')).hexdigest()
                keep = (not res['nostderr']) and len(res['err']) < 200000  # the text of every standard error (the tests count the program's messages)
                full += 1
                cols = [res['id'], res['group'], json.dumps(res['args']), res['input'], json.dumps(res['files']), str(res['rc']), res['out'], sha,
                        res['err'] if keep else '', json.dumps(res['outfiles'])]
                ft.write('\t'.join(L.esc(c) for c in cols) + '\n')
            if done % 200 == 0:
                print(done, '/', n, file=sys.stderr)
    if ft:
        ft.close()
    print('cases', n, 'written', n - len(dropped), 'non-deterministic (dropped)', len(dropped), dropped[:10], file=sys.stderr)


if __name__ == '__main__':
    main()
