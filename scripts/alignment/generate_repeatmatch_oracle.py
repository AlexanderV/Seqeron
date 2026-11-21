#!/usr/bin/env python3
"""A3-9e3 oracle: `repeat-match` and the script `exact-tandems` of MUMmer 3.23 (Debian 3.23+dfsg-8).
Rows (JSON lines):
  {"kind":"rm","fasta","args":[...],"out","err","rc"}                 repeat-match <args> g.fa (-n N, -f, -t, -E, -V k): standard output / error (latin-1), exit status
  {"kind":"rm","ub":true,...}                                         a row whose result is undefined in the program (an empty file / a file without a '>' line: the sequence buffer is never written; valgrind: use of
                                                                      uninitialised values; the output changes under MALLOC_PERTURB_) - recorded, not compared
  {"kind":"et","fasta","n","out","err","rc","left":[names]}          the INSTALLED script exact-tandems g.fa N (broken: `awk -f $` and the file name on the next line): standard output / error and the files left behind
  {"kind":"etx","fasta","n","out","err","rc"}                         the INTENDED pipeline `repeat-match -t -n N g.fa | tail +3 | sort -k1n -k2n | awk -f /usr/lib/mummer/tandem-repeat.awk`, run by hand
Sequences: random, planted tandem repeats (unit lengths 1 - 60, copies 2 - 12, with and without point errors), dispersed repeats, reverse-complement repeats, palindromes, homopolymers, N / IUPAC characters,
upper case, several records (only the first is used), characters that are not letters, line widths.
Usage: python3 scripts/alignment/generate_repeatmatch_oracle.py [seed=20261205] [n_rm=330] [n_et=215] | gzip > repeatmatch_oracle.jsonl.gz"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09rm_%d' % os.getpid()
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261205
    n_rm = int(sys.argv[2]) if len(sys.argv) > 2 else 330
    n_et = int(sys.argv[3]) if len(sys.argv) > 3 else 215
    rnd = random.Random(seed)

    def rseq(k, alpha='acgt'): return ''.join(rnd.choice(alpha) for _ in range(k))
    def rc(s): return s.translate(COMP)[::-1]
    def mutate(s, rate):
        o = []
        for c in s:
            r = rnd.random()
            if r < rate * .6: o.append(rnd.choice('acgt'))
            elif r < rate * .75: pass
            elif r < rate: o.append(c + rnd.choice('acgt'))
            else: o.append(c)
        return ''.join(o)
    def fasta(recs, width=60, upper=False):
        out = []
        for nm, sq in recs:
            out.append('>' + nm + '\n')
            if upper: sq = sq.upper()
            for i in range(0, len(sq), width): out.append(sq[i:i + width] + '\n')
        return ''.join(out)
    def reset():
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
    def lat(b): return b.decode('latin-1')
    def run(cmd, **kw): return subprocess.run(cmd, cwd=WORK, capture_output=True, **kw)

    def genome(kind):
        """one sequence with planted structure"""
        n = rnd.choice([40, 80, 150, 300, 600, 1200, 2500])
        s = rseq(n)
        parts = [s]
        for _ in range(rnd.randint(0, 4)):
            u = rseq(rnd.randint(1, 60))
            t = u * rnd.randint(2, 12)
            if rnd.random() < .4: t = mutate(t, .03)
            pos = rnd.randint(0, len(s))
            s = s[:pos] + t + s[pos:]
        for _ in range(rnd.randint(0, 3)):                       # dispersed repeats (forward and reverse complement)
            seg = s[rnd.randint(0, max(0, len(s) - 100)):][:rnd.randint(15, 120)]
            if rnd.random() < .5: seg = rc(seg)
            if rnd.random() < .3: seg = mutate(seg, .05)
            pos = rnd.randint(0, len(s))
            s = s[:pos] + seg + s[pos:]
        r = rnd.random()
        if kind % 17 == 3: s = 'a' * rnd.randint(10, 90)                        # homopolymer
        elif kind % 17 == 5: s = 'ac' * rnd.randint(8, 60)
        elif kind % 17 == 7: x = rseq(rnd.randint(8, 80)); s = x + rc(x)      # a reverse-complement palindrome: the program stops
        elif kind % 17 == 9: x = rseq(rnd.randint(8, 80)); s = x + rseq(rnd.randint(0, 5)) + rc(x)
        elif kind % 17 == 11: s = ''.join(rnd.choice('acgtnnryk') for _ in range(rnd.randint(30, 300)))
        elif kind % 17 == 13: s = rseq(rnd.randint(1, 6))
        return s

    rows = []

    def rm_row(fa, args):
        reset()
        open(WORK + '/g.fa', 'w').write(fa)
        p = run(['repeat-match'] + args + ['g.fa'])
        if len(p.stdout) > 30000: return None                                   # keep the oracle small
        row = {'kind': 'rm', 'fasta': fa, 'args': args, 'out': lat(p.stdout), 'err': lat(p.stderr), 'rc': p.returncode}
        if p.returncode == 0 or p.returncode == 255:
            q = run(['repeat-match'] + args + ['g.fa'], env=dict(os.environ, MALLOC_PERTURB_='165'))
            if q.stdout != p.stdout or q.stderr != p.stderr: row['ub'] = True
        return row

    # ---------- repeat-match
    for k in range(n_rm):
        s = genome(k)
        recs = [('g%d description text' % k, s)]
        if k % 9 == 4: recs.append(('second', rseq(rnd.randint(5, 80))))
        if k % 12 == 6: recs[0] = (recs[0][0], s.upper())
        fa = fasta(recs, rnd.choice([50, 60, 70, 1000]))
        if k % 23 == 8:
            fa = fa.replace('\n', '\r\n')
        if k % 29 == 10:                                                           # characters that are not letters: x with a warning
            fa = fa.replace('a', '1', 2).replace('c', '-', 1)
        args = []
        r = rnd.random()
        if r < .6: args += ['-n', str(rnd.choice([3, 5, 8, 10, 12, 15, 20, 25, 30, 40]))]
        elif r < .65 and len(s) <= 200: args += ['-n', '1']
        if rnd.random() < .25: args.append('-f')
        if rnd.random() < .2: args.append('-t')
        if len(s) <= 700 and rnd.random() < .3: args.append('-E')
        if len(s) <= 150 and rnd.random() < .25: args += ['-V', str(rnd.choice([1, 2]))]
        elif rnd.random() < .05: args += ['-V', '1']
        row = rm_row(fa, args)
        if row is not None: rows.append(row)
    # large genomes (the tree of a long repeat is deep: the traversals must not recurse) and long tandem arrays
    for k in range(6):
        u = rseq(rnd.randint(30, 400))
        g = rseq(rnd.randint(3000, 8000)) + u * rnd.randint(3, 10) + rseq(rnd.randint(2000, 9000)) + rc(u) + rseq(rnd.randint(500, 3000)) + u[:rnd.randint(10, 20)] * 6
        row = rm_row(fasta([('big%d' % k, g)], 70), ['-n', str(rnd.choice([25, 40, 60])), rnd.choice(['-f', '-t', '-f'])])
        if row is not None: rows.append(row)
    for g, args in [('a' * 4000, ['-n', '3900']), ('a' * 3000, ['-t', '-n', '1400']), ('ac' * 2500, ['-n', '4900']), ('acg' * 1500, ['-f', '-n', '4400']), ('a' * 1200, ['-n', '1190', '-V', '1'])]:
        row = rm_row(fasta([('homopolymer', g)], 60), args)
        if row is not None: rows.append(row)
    # the special inputs
    for fa, args in [('', []), ('acgtacgtacgtacgt\n', []), ('>x\n', []), ('>x\nacgtacgtacgtacgtacgtacgtac\n', ['-n', '0']), ('>x\nacgtacgtacgtacgtacgtacgtac\n', ['-n', 'x']),
                     ('>' + 'x' * 600 + '\nacgtacgtacgtacgtacgtacgtac\n', ['-n', '5']), ('>' + 'x' * 499 + '\nacgtacgtacgtacgtacgtacgtacgtac\n', ['-n', '5']), ('>' + 'x' * 500 + '\nacgtacgtacgtacgtacgtacgtac\n', ['-n', '5']),
                     ('>x\nacgtacgtacgtacgtacgtacgtac', ['-n', '5']), ('>x', []), ('>x\nA\n', ['-n', '1']), ('>x\nacgt\n', ['-n', '1']), ('>x\nacgtaa\n', ['-n', '2', '-V', '2']),
                     ('>x\n\n>y\nacgtacgtacgtac\n', ['-n', '4']), ('text before\n>x\nacgtacgtacgtacgtacgtac\n', ['-n', '4']), ('>x\nacgt acgt\tacgtac\n>y\nggg\n', ['-n', '3'])]:
        reset()
        open(WORK + '/g.fa', 'w').write(fa)
        if args[-2:] == ['-n', 'x']:
            continue
        p = run(['repeat-match'] + args + ['g.fa'])
        row = {'kind': 'rm', 'fasta': fa, 'args': args, 'out': lat(p.stdout), 'err': lat(p.stderr), 'rc': p.returncode}
        q = run(['repeat-match'] + args + ['g.fa'], env=dict(os.environ, MALLOC_PERTURB_='165'))
        if q.stdout != p.stdout or q.stderr != p.stderr or (not fa.lstrip().startswith('>') and '>' not in fa):
            row['ub'] = True
            vg = run(['valgrind', '-q', '--error-exitcode=9', 'repeat-match'] + args + ['g.fa'])
            row['valgrind'] = vg.returncode == 9 or 'uninitialised' in lat(vg.stderr)
        rows.append(row)

    # ---------- exact-tandems: the installed script and the intended pipeline
    for k in range(n_et):
        r = rnd.random()
        u = rseq(rnd.randint(1, 50))
        s = rseq(rnd.randint(0, 300)) + u * rnd.randint(2, 14)
        if rnd.random() < .5: s = mutate(s, .01)
        s += rseq(rnd.randint(0, 200)) + (rseq(rnd.randint(1, 30)) * rnd.randint(2, 6))
        if rnd.random() < .3: s += rseq(rnd.randint(10, 60))
        if k % 19 == 7: s = 'a' * rnd.randint(20, 80)
        if k % 23 == 9: x = rseq(rnd.randint(8, 60)); s = x + rc(x)
        n = rnd.choice([3, 5, 8, 10, 12, 15, 20, 30, 40])
        fa = fasta([('t%d' % k, s)], 60, upper=(k % 13 == 5))
        reset()
        open(WORK + '/g.fa', 'w').write(fa)
        p = run(['exact-tandems', 'g.fa', str(n)])
        if len(run(['repeat-match', '-t', '-n', str(n), 'g.fa']).stdout) > 20000: continue
        left = sorted(f for f in os.listdir(WORK) if f != 'g.fa')
        row = {'kind': 'et', 'fasta': fa, 'n': n, 'out': lat(p.stdout), 'err': lat(p.stderr), 'rc': p.returncode, 'left': [f.split('.', 1)[1] if f[0].isdigit() else f for f in left]}
        rows.append(row)
        reset()
        open(WORK + '/g.fa', 'w').write(fa)
        p = run('repeat-match -t -n %d g.fa | tail +3 > m.tmp; sort -k1n -k2n m.tmp | awk -f /usr/lib/mummer/tandem-repeat.awk' % n, shell=True)
        rows.append({'kind': 'etx', 'fasta': fa, 'n': n, 'out': lat(p.stdout), 'err': lat(p.stderr), 'rc': p.returncode})
    # the usage errors of the script
    for args in [[], ['g.fa'], ['g.fa', '10', 'x'], ['nofile.fa', '10']]:
        reset()
        open(WORK + '/g.fa', 'w').write('>x\nacgtacgtacgtacgtacgtac\n')
        p = run(['exact-tandems'] + args)
        rows.append({'kind': 'et', 'fasta': '>x\nacgtacgtacgtacgtacgtac\n', 'n': 0, 'args': args, 'out': lat(p.stdout), 'err': lat(p.stderr), 'rc': p.returncode, 'left': []})
    shutil.rmtree(WORK, ignore_errors=True)
    for r in rows: print(json.dumps(r))
    kinds = {}
    for r in rows:
        key = r['kind'] + (':rc%d' % r['rc'] if r['rc'] else '') + (':ub' if r.get('ub') else '')
        kinds[key] = kinds.get(key, 0) + 1
    sys.stderr.write('rows=%d %s\n' % (len(rows), kinds))


main()
