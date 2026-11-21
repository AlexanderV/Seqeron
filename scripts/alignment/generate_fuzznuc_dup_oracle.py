#!/usr/bin/env python3
"""A5-8c2b (B09 F195) oracle: fuzznuc 6.6.0 reports for patterns with ranges and -pmismatch (embPatBruteForce -> patOUBrute pushes one hit per successful unrolling, so a window that
several unrollings match is listed several times), both strands (-complement), excel format (Start End Score Strand Pattern Mismatch).
JSONL.GZ rows: id, seq (one FASTA record text), pattern, mismatch, complement, report (the excel text, whole file).
Usage: python3 scripts/alignment/generate_fuzznuc_dup_oracle.py OUT.jsonl.gz [n=300] [seed]     (needs EMBOSS 6.6.0)
"""
import gzip, json, os, random, subprocess, sys, tempfile


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            s = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(20, 90)))
            parts = []
            for _ in range(rnd.randint(2, 4)):
                a = rnd.choice(['A', 'C', 'G', 'T', 'N', 'R', 'Y', '[AC]', '[GT]', '{A}', '{CG}', '[ACG]'])
                r = rnd.random()
                if r < 0.55:
                    lo = rnd.randint(0, 2); hi = lo + rnd.randint(1, 2); a += '(%d,%d)' % (lo, hi)
                elif r < 0.7:
                    a += '(%d)' % rnd.randint(1, 3)
                parts.append(a)
            pat = ''.join(parts)
            if all('(0,' in a for a in parts): continue          # every element optional: the program also lists empty hits (not ported: registered)
            if rnd.random() < 0.1: pat = '<' + pat
            if rnd.random() < 0.1: pat = pat + '>'
            mm = rnd.choice([1, 1, 2]); comp = rnd.random() < 0.5
            open(os.path.join(d, 's.fa'), 'w').write('>s\n' + s + '\n')
            o = os.path.join(d, 'o.xls')
            if os.path.exists(o): os.remove(o)
            p = subprocess.run(['fuzznuc', 's.fa', '-pattern', pat, '-pmismatch', str(mm), '-complement', 'Y' if comp else 'N', '-rformat', 'excel', '-outfile', 'o.xls', '-auto'], cwd=d, capture_output=True, text=True, timeout=20)
            if p.returncode != 0 or 'Illegal' in p.stderr + p.stdout or 'Fatal' in p.stderr:
                continue
            rep = open(o, encoding='latin-1').read() if os.path.exists(o) else ''
            rows.append({'id': len(rows), 'seq': '>s\n' + s + '\n', 'pattern': pat, 'mismatch': mm, 'complement': comp, 'report': rep})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 20261130)
