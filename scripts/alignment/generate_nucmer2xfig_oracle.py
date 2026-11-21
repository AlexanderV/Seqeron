#!/usr/bin/env python3
"""A3-9e oracle: `nucmer2xfig` (perl script of MUMmer 3.23) on coords files of `show-coords -c -l` form: synthetic alignments of many contigs against references of 0.02 - 3 Mb (several alignments per contig, contigs
split by a gap of more than twice their length, identities below 50, equal and adjacent starts, a repeated contig id later in the file), real coords files of nucmer runs, and malformed input (short lines, blank lines,
text instead of numbers, an empty file). Rows: {"coords","out","rc"}. Usage: python3 generate_nucmer2xfig_oracle.py [n=260] [seed=20261205] | gzip > nucmer2xfig_oracle.jsonl.gz"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09xfig'

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 260
    rnd = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 20261205)
    shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
    rows = 0
    for case in range(n):
        kind = rnd.random()
        lines = []
        if kind < .75:
            rlen = rnd.choice([20000, 600000, 1250000, 3000000, rnd.randint(1000, 4000000)])
            rid = rnd.choice(['2R', 'chr1', 'NC_000913', 'ref'])
            for c in range(rnd.randint(1, 14)):
                qlen = rnd.randint(300, 80000); qid = 'contig%d' % (c + 1 if rnd.random() < .9 else rnd.randint(1, 5))
                pos = rnd.randint(1, max(2, rlen - qlen))
                for _ in range(rnd.randint(1, 5)):
                    l = rnd.randint(50, min(qlen, 30000)); s1 = pos if rnd.random() < .8 else max(1, pos - rnd.randint(0, 200)); e1 = min(rlen, s1 + l)
                    pos = e1 + rnd.choice([0, 10, qlen // 2, qlen * 3])
                    pct = rnd.choice([rnd.uniform(40, 100), rnd.uniform(90, 100), 100.0, 49.0, 50.0])
                    s2 = rnd.randint(1, 1000); e2 = s2 + l
                    cov1 = 100.0 * l / rlen; cov2 = 100.0 * l / qlen
                    lines.append('%d  %d | %d %d | %d  %d | %.2f  | %d  %d | %.2f  %.2f | %s %s' % (s1, e1, s2, e2, e1 - s1 + 1, l, pct, rlen, qlen, cov1, cov2, rid, qid))
            if rnd.random() < .6: lines.sort(key=lambda l: (l.split()[18], int(l.split()[3])))
        elif kind < .9:
            ref = ''.join(rnd.choice('acgt') for _ in range(rnd.randint(2000, 6000)))
            qs = []
            for _ in range(rnd.randint(1, 5)):
                a = rnd.randint(0, len(ref) - 800); piece = list(ref[a:a + rnd.randint(300, 800)])
                for i in range(len(piece)):
                    if rnd.random() < .03: piece[i] = rnd.choice('acgt')
                qs.append(''.join(piece))
            open(WORK + '/r.fa', 'w').write('>chrA\n' + ref + '\n'); open(WORK + '/q.fa', 'w').write(''.join('>c%d\n%s\n' % (i, q) for i, q in enumerate(qs)))
            if subprocess.run(['nucmer', '-p', 'o', '-l', '12', '-c', '20', 'r.fa', 'q.fa'], cwd=WORK, capture_output=True).returncode != 0: continue
            lines = subprocess.run(['show-coords', '-c', '-l', '-H', '-r', 'o.delta'], cwd=WORK, capture_output=True, text=True).stdout.splitlines()
        else:
            lines = rnd.choice([[], [''], ['abc def'], ['1 2 | 3 4 | 5 6 | x | 7 8 | 9 10 | a b'], ['10  50 | 1 2 | 3  4 | 70.5  | 1000000  400 | 0.1  0.2 | R Q', '', '60 90 | 1 2 | 3 4 | 95.00 | 1000000 400 | 0.1 0.2 | R Q', 'short line'], ['  5  9 | 1 2 | 3  4 | 1e2  | 12abc  400 | 0.1  0.2 | R Q']])
        text = '\n'.join(lines) + ('\n' if lines and rnd.random() < .9 else '')
        open(WORK + '/in.coords', 'w').write(text)
        p = subprocess.run(['nucmer2xfig', 'in.coords'], cwd=WORK, capture_output=True)
        print(json.dumps({'coords': text, 'out': p.stdout.decode('latin-1'), 'rc': p.returncode})); rows += 1
    sys.stderr.write('rows=%d\n' % rows)

main()
