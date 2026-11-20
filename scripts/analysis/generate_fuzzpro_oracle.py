#!/usr/bin/env python3
"""GENOMIC-MOTIFS-001 oracle (A3-12): EMBOSS fuzzpro 6.6.0 (protein pattern search) on random syntax-tree patterns.

fuzzpro has no -pmismatch (6.6.0 help): the search is exact. Protein sequences 12-140 aa over the 20 amino acids with occasional X, B, Z, U, J, O, * and lower case; patterns drawn
from the documented language: letters, `x` (any), classes [ILV], exclusions {P}, repeats X(3), ranges X(2,4), `-` separators, the anchors `<` / `>`, a trailing period; planted instances.
`fuzzpro -sequence s.fa -pattern P -rformat gff` -> (start, end) 1-based per hit (the binary lists every matching window). The binary reads a trailing stop `*` (one) off the sequence.
TSV (tests/.../TestData/Analysis/fuzzpro_oracle.tsv):  sequence \t pattern \t hits ('start0:length' 0-based, sorted, space separated; '-' = none; 'ERR' = the binary rejects the pattern)
Usage: python3 scripts/analysis/generate_fuzzpro_oracle.py [rows] > fuzzpro_oracle.tsv
"""
import os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261127)
AA = 'ACDEFGHIKLMNPQRSTVWY'


def gen_pattern(seq):
    n = rnd.randint(1, 6)
    parts, sample = [], ''
    for _ in range(n):
        k = rnd.random()
        if k < 0.4:
            c = rnd.choice(AA + 'x' * 2 + 'XBZ'); atom = c; s = rnd.choice(AA)
        elif k < 0.7:
            L = rnd.sample(AA + ('XBZ' if rnd.random() < 0.1 else ''), rnd.randint(1, 4)); atom = '[' + ''.join(L) + ']'; s = rnd.choice([a for a in L if a in AA] or ['A'])
        else:
            L = rnd.sample(AA + ('XBZ' if rnd.random() < 0.1 else ''), rnd.randint(1, 3)); atom = '{' + ''.join(L) + '}'; s = rnd.choice([a for a in AA if a not in L] or ['A'])
        r = rnd.random(); reps = 1
        if r < 0.2:
            reps = rnd.randint(2, 3); atom += '(%d)' % reps
        elif r < 0.4:
            lo = rnd.randint(0, 2); hi = lo + rnd.randint(1, 2); atom += '(%d,%d)' % (lo, hi); reps = lo
        parts.append(atom); sample += s * reps
    sep = '-' if rnd.random() < 0.4 else ''
    pat = sep.join(parts)
    if all(p.endswith(',0)') or '(0,' in p for p in parts) and all('(0,' in p for p in parts):
        return gen_pattern(seq)
    if rnd.random() < 0.15: pat = '<' + pat
    if rnd.random() < 0.15: pat = pat + '>'
    if rnd.random() < 0.3: pat += '.'
    return pat, sample


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 700
    done = 0
    with tempfile.TemporaryDirectory() as tmp:
        for _ in range(n):
            L = rnd.randint(12, 140)
            s = ''.join(rnd.choice(AA) for _ in range(L))
            pat, sample = gen_pattern(s)
            if sample and rnd.random() < 0.6:
                p = rnd.randrange(len(s)); s = s[:p] + sample + s[p:]
            if rnd.random() < 0.45:
                s = ''.join(rnd.choice('XBZUJO*') if rnd.random() < 0.06 else c for c in s)
            if rnd.random() < 0.25:
                s = s.lower() if rnd.random() < 0.5 else ''.join(c.lower() if rnd.random() < 0.5 else c for c in s)
            fa, out = os.path.join(tmp, 's.fa'), os.path.join(tmp, 'o.gff')
            open(fa, 'w').write('>s\n' + s + '\n')
            if os.path.exists(out): os.remove(out)
            r = subprocess.run(['fuzzpro', '-sequence', fa, '-pattern', pat, '-outfile', out, '-rformat', 'gff', '-auto'], capture_output=True, text=True)
            if r.returncode != 0:
                sys.stderr.write('FAIL %s %s: %s\n' % (s, pat, r.stderr.strip()[:200])); continue
            if 'Illegal pattern' in (r.stdout + r.stderr):
                print('%s\t%s\tERR' % (s, pat))   # the binary rejects the pattern (B / Z / X inside [] or {} : "Nesting '[' in [] not allowed" / "Non alpha character")
                done += 1
                continue
            hits = []
            if os.path.exists(out):
                for line in open(out):
                    if line.startswith('#') or not line.strip(): continue
                    f = line.split('\t')
                    hits.append((int(f[3]) - 1, int(f[4]) - int(f[3]) + 1))
            hits.sort()
            print('%s\t%s\t%s' % (s, pat, ' '.join('%d:%d' % h for h in hits) or '-'))
            done += 1
    sys.stderr.write('rows: %d\n' % done)


if __name__ == '__main__':
    main()
