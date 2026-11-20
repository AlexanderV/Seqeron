#!/usr/bin/env python3
"""ALIGN-SEMI-001 oracle (A3-13): EMBOSS merger 6.6.0 on random overlapping read pairs.

Pairs of reads of a random genome fragment (120-450 nt): overlapping (dovetail either way, containment), with point errors, 1-3 base indels, N and IUPAC ambiguity codes
(R Y K M S W), lower / mixed case, and a few unrelated pairs; -gapopen / -gapextend from {50/5 (default), 10/0.5, 20/2, 100/10, 5/1, 0/0}, default EDNAFULL.
`merger -asequence a.fa -bsequence b.fa -outfile m.txt -outseq m.fa`: merged sequence (outseq, case as written), the `Conflicts` lines and the warning lines of the report, the score.
TSV (tests/.../TestData/Alignment/merger_oracle.tsv):  seqA \t seqB \t gapopen \t gapextend \t merged \t conflicts ('pos:baseA:pos:baseB:using' joined by ';' or '-') \t warnings (lines joined by '|' or '-') \t score
Usage: python3 scripts/alignment/generate_merger_oracle.py [rows] > merger_oracle.tsv
"""
import os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261125)


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate:
            out.append(rnd.choice('ACGT'))
        elif r < rate * 1.15:
            continue
        elif r < rate * 1.3:
            out.append(c); out.append(rnd.choice('ACGT'))
        else:
            out.append(c)
    return ''.join(out)


def ambiguity(s, rate):
    return ''.join((rnd.choice('NNNRYKMSW') if rnd.random() < rate else c) for c in s)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 420
    done = 0
    with tempfile.TemporaryDirectory() as tmp:
        for _ in range(n):
            g = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(120, 450)))
            L = len(g)
            kind = rnd.random()
            if kind < 0.6:   # dovetail overlap
                x2 = rnd.randint(L // 4, L // 2); y1 = rnd.randint(x2 + 15, min(L, x2 + L // 2)); a, b = g[0:y1], g[x2:L]
            elif kind < 0.78:  # containment
                x = rnd.randint(1, L // 3); y = rnd.randint(L // 2, L - 1); a, b = g, g[x:y]
            elif kind < 0.90:  # unaligned flanks on both sequences (warnings)
                core = g[rnd.randint(0, L // 3):][:rnd.randint(40, 120)]
                junk = lambda: ''.join(rnd.choice('ACGT') for _ in range(rnd.choice([0, rnd.randint(1, 25), 12])))
                a, b = junk() + core + junk(), junk() + core + junk()
            elif kind < 0.95:  # same start / same end overhangs
                x = rnd.randint(10, 40); a, b = g[x:], g[x:L - rnd.randint(5, 30)]
            else:
                a, b = g[:rnd.randint(30, 100)], ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(30, 100)))
            a = mutate(a, rnd.choice([0.0, 0.01, 0.03])); b = mutate(b, rnd.choice([0.0, 0.01, 0.03]))
            a = ambiguity(a, rnd.choice([0, 0, 0.01, 0.03])); b = ambiguity(b, rnd.choice([0, 0, 0.01, 0.03]))
            if rnd.random() < 0.15: a = a.lower()
            if rnd.random() < 0.15: b = ''.join(c.lower() if rnd.random() < 0.5 else c for c in b)
            if rnd.random() < 0.3: a, b = b, a
            go, ge = rnd.choice([(50, 5), (50, 5), (10, 0.5), (20, 2), (100, 10), (5, 1), (0, 0)])
            fa, fb = os.path.join(tmp, 'a.fa'), os.path.join(tmp, 'b.fa')
            open(fa, 'w').write('>a\n' + a + '\n'); open(fb, 'w').write('>b\n' + b + '\n')
            out, outseq = os.path.join(tmp, 'm.txt'), os.path.join(tmp, 'm.fa')
            r = subprocess.run(['merger', '-asequence', fa, '-bsequence', fb, '-outfile', out, '-outseq', outseq, '-gapopen', str(go), '-gapextend', str(ge), '-auto'], capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(outseq):
                sys.stderr.write('FAIL %s %s: %s\n' % (a, b, r.stderr.strip()[:200])); continue
            merged = ''.join(l.strip() for l in open(outseq) if not l.startswith('>'))
            txt = open(out).read()
            conflicts = ';'.join('%s:%s:%s:%s:%s' % m for m in re.findall(r"^#\s+(\d+)\s+'(.)'\s+(\d+)\s+'(.)'\s+'(.)'", txt, re.M)) or '-'
            warns = []
            for line in txt.split('\n'):
                l = line[1:].strip() if line.startswith('#') else line.strip()
                if l.startswith(('The region of alignment', 'Only the sequence of', 'There is an equal amount', 'Sequence a is being', 'Sequence b is being')):
                    warns.append(re.sub(r'sequence ([ab])\b', lambda m: 'sequence ' + m.group(1).upper(), re.sub(r'Sequence ([ab]) is', lambda m: 'Sequence ' + m.group(1).upper() + ' is', l)))
            score = re.search(r'# Score: (\S+)', txt).group(1)
            print('\t'.join([a, b, str(go), str(ge), merged, conflicts, '|'.join(warns) or '-', score]))
            done += 1
    sys.stderr.write('rows: %d\n' % done)


if __name__ == '__main__':
    main()
