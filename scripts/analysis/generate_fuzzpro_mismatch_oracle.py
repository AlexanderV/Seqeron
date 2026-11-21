#!/usr/bin/env python3
"""A5-8c2 (B09 F195) oracle: EMBOSS fuzzpro / fuzztran 6.6.0 protein patterns WITH MISMATCHES (-pmismatch n: embPatGetTypeII picks BYP / TUB / brute force for a protein pattern, embPatBruteForce with the protein alphabet).

Same random syntax-tree patterns as generate_fuzzpro_oracle.py (letters, x, [classes], {exclusions}, repeats, ranges, anchors), -pmismatch 1 - 3.
`fuzzpro -sequence s.fa -pattern P -pmismatch n -rformat excel` -> (start, end, score, mismatch) per hit (the excel format has the Mismatch column; empty = 0); fuzztran the same on a nucleotide sequence with -frame 6 (hits as (frame, start, end, mismatch)).
TSV: kind \t sequence \t pattern \t mismatch \t hits ('start0:length:mm' 0-based, sorted, space separated; '-' = none; 'ERR' = the binary rejects the pattern or aborts).  fuzztran rows: 'start:end:mm@frame' (signed frame) with the program's own 1-based protein start / end.
Usage: python3 scripts/analysis/generate_fuzzpro_mismatch_oracle.py [rows=800] > fuzzpro_mismatch_oracle.tsv     (needs EMBOSS 6.6.0)
"""
import os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import generate_fuzzpro_oracle as G  # noqa: E402

rnd = random.Random(20261128)
G.rnd = rnd
AA = G.AA


def run(prog, args, tmp):
    out = os.path.join(tmp, 'o.xls')
    if os.path.exists(out): os.remove(out)
    r = subprocess.run([prog] + args + ['-outfile', out, '-rformat', 'excel', '-auto'], capture_output=True, text=True)
    text = r.stdout + r.stderr
    if r.returncode != 0 or 'Illegal pattern' in text or 'Fatal' in text or 'Error' in text:
        return None
    hits = []
    if os.path.exists(out):
        for line in open(out):
            f = line.rstrip('\n').split('\t')
            if line.startswith('SeqName') or len(f) < 7: continue
            hits.append(f)
    return hits


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 800
    done = 0
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, 's.fa')
        for i in range(n):
            tran = i % 4 == 3
            if tran:
                L = rnd.randint(30, 300); s = ''.join(rnd.choice('ACGT') for _ in range(L))
                pat, sample = G.gen_pattern(s)
            else:
                L = rnd.randint(12, 140)
                s = ''.join(rnd.choice(AA) for _ in range(L))
                pat, sample = G.gen_pattern(s)
                if sample and rnd.random() < 0.6:
                    p = rnd.randrange(len(s)); s = s[:p] + sample + s[p:]
                if rnd.random() < 0.3:
                    s = ''.join(rnd.choice('XBZUJO*') if rnd.random() < 0.05 else c for c in s)
                # plant near-misses so that mismatch hits exist
                if sample and rnd.random() < 0.6:
                    q = list(s); p = rnd.randrange(len(q)); c = list(sample)
                    for _ in range(rnd.randint(1, 2)):
                        if c: c[rnd.randrange(len(c))] = rnd.choice(AA)
                    s = ''.join(q[:p]) + ''.join(c) + ''.join(q[p:])
                if rnd.random() < 0.2:
                    s = s.lower()
            mm = rnd.choice([1, 1, 2, 3])
            open(fa, 'w').write('>s\n' + s + '\n')
            if tran:
                hits = run('fuzztran', ['-sequence', fa, '-pattern', pat, '-pmismatch', str(mm), '-frame', '6'], tmp)
            else:
                hits = run('fuzzpro', ['-sequence', fa, '-pattern', pat, '-pmismatch', str(mm)], tmp)
            if hits is None:
                print('%s\t%s\t%s\t%d\tERR' % ('tran' if tran else 'pro', s, pat, mm)); done += 1; continue
            if tran:
                # excel columns of fuzztran: SeqName Start End Score Strand Pattern Mismatch Frame PStart PEnd Translation
                rows = sorted('%s:%s:%s@%s%s' % (h[8], h[9], h[6] if h[6] != '.' else '0', '-' if h[4] == '-' else '', h[7]) for h in hits)
            else:
                rows = sorted('%d:%d:%s' % (int(h[1]) - 1, int(h[2]) - int(h[1]) + 1, h[6] if h[6] != '.' else '0') for h in hits)
            print('%s\t%s\t%s\t%d\t%s' % ('tran' if tran else 'pro', s, pat, mm, ' '.join(rows) or '-'))
            done += 1
    sys.stderr.write('rows: %d\n' % done)


if __name__ == '__main__':
    main()
