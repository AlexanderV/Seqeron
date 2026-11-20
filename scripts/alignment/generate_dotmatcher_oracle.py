#!/usr/bin/env python3
"""ALIGN-LOCAL-001 oracle (audit A2-28d): EMBOSS dotmatcher 6.6.0 (windowed diagonal score dot plot) vs DotMatcher.Compute.

`dotmatcher a b -windowsize W -threshold T -matrixfile M -graph data` writes the plot as a text file: after the title lines the `Line x1 y1 x2 y2 colour 0` entries are the diagonal segments in the
order the program draws them (coordinates 1-based, x on sequence a, y on sequence b), followed by the frame and axes (`Rectangle ...`). This script records the segments before the frame.
Cases: related DNA pairs (EDNAFULL, EDNASIMPLE) and protein pairs (EBLOSUM62, EPAM250, EBLOSUM30), 12-200 residues (always longer than the window), W 3-20, T chosen around the score of an identical window.
TSV: kind  matrix  W  T  seq1  seq2  segments ('x1:y1:x2:y2' separated by spaces, '-' = none)
Never aborts: a failing call is skipped.
Usage: python3 scripts/alignment/generate_dotmatcher_oracle.py [n=500] > dotmatcher_oracle.tsv
"""
import os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261131)
AA = 'ACDEFGHIKLMNPQRSTVWY'


def mutate(s, rate, alpha):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.6: out.append(rnd.choice(alpha))
        elif r < rate * 0.8: pass
        elif r < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    rows = 0
    with tempfile.TemporaryDirectory() as d:
        for _ in range(n):
            protein = rnd.random() < 0.4
            alpha = AA if protein else 'ACGT'
            base = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(12, 200)))
            if rnd.random() < 0.7:
                a = mutate(base, rnd.choice([0.05, 0.15, 0.3]), alpha) + ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, 20)))
                b = mutate(base, rnd.choice([0.05, 0.15, 0.3]), alpha)
            else:
                a = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(12, 120))); b = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(12, 120)))
            w = rnd.randint(3, min(20, len(a) - 1, len(b) - 1))
            mat = rnd.choice(['EBLOSUM62', 'EBLOSUM62', 'EPAM250', 'EBLOSUM30']) if protein else rnd.choice(['EDNAFULL', 'EDNAFULL', 'EDNASIMPLE'])
            unit = 5 if not protein else 5
            t = rnd.randint(max(1, w * unit // 3), w * unit)
            open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
            for f in os.listdir(d):
                if f.startswith('dm'): os.remove(os.path.join(d, f))
            r = subprocess.run(['dotmatcher', d + '/a.fa', d + '/b.fa', '-windowsize', str(w), '-threshold', str(t), '-matrixfile', mat, '-graph', 'data', '-goutfile', d + '/dm', '-auto'], capture_output=True, text=True)
            files = [f for f in os.listdir(d) if f.startswith('dm') and f.endswith('.dat')]
            if r.returncode != 0 or not files: continue
            segs = []
            for line in open(os.path.join(d, files[0])):
                if line.startswith('Rectangle'): break
                m = re.match(r'Line x1 (\S+) y1 (\S+) x2 (\S+) y2 (\S+) colour', line)
                if m: segs.append(':'.join(str(int(float(g))) for g in m.groups()))
            print('\t'.join(['prot' if protein else 'dna', mat, str(w), str(t), a, b, ' '.join(segs) or '-']))
            rows += 1
    sys.stderr.write('rows=%d\n' % rows)


main()
