#!/usr/bin/env python3
"""A5-8c3b (B09 F217) oracle: fuzznuc, fuzzpro and fuzztran 6.6.0 reports for patterns whose elements are ALL optional ranges (C(0,2)C(0,1), [AC](0,3), ...), with and without -pmismatch, both strands, anchors < and >:
the recursive matcher (embPatBruteForce -> patOUBrute) also pushes EMPTY hits (end = start - 1, score 0) and hits that start past the end of the sequence. excel format (Start End Score Strand Pattern Mismatch), whole file.
JSONL.GZ rows as generate_fuzznuc_dup_oracle.py plus `program` (fuzznuc | fuzzpro | fuzztran; only fuzznuc has `complement`, a fuzztran row has `frame`). Usage: python3 scripts/alignment/generate_fuzznuc_optional_oracle.py OUT.jsonl.gz [n=300] [seed]     (needs EMBOSS 6.6.0)
"""
import gzip, json, os, random, subprocess, sys, tempfile


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            tran = len(rows) >= n * 3 // 4
            prot = len(rows) >= n // 2
            alpha = 'ACDEFGHIKLMNPQRSTVWY' if prot and not tran else 'ACGT'
            s = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(3 if tran else 1, 60 if tran else 40)))
            frame = rnd.choice(['1', '2', '3', 'F', '6', 'R', '-1', '-2']) if tran else None
            parts = []
            for _ in range(rnd.randint(1, 4)):
                a = rnd.choice(['A', 'C', 'D', 'X', 'L', '[AC]', '[DE]', '{A}', '{CG}', '[ACG]']) if prot else rnd.choice(['A', 'C', 'G', 'T', 'N', 'R', 'Y', '[AC]', '[GT]', '{A}', '{CG}', '[ACG]'])
                lo = 0; hi = rnd.randint(1, 3)
                parts.append(a + '(%d,%d)' % (lo, hi))
            if rnd.random() < 0.25: parts[rnd.randrange(len(parts))] = rnd.choice(['A', 'C', 'D', 'X', 'L'] if prot else ['A', 'C', 'G', 'T', 'N']) + rnd.choice(['(0,1)', '(0,2)'])
            pat = ''.join(parts)
            if rnd.random() < 0.15: pat = '<' + pat
            if rnd.random() < 0.15: pat = pat + '>'
            mm = rnd.choice([0, 0, 0, 1, 2]); comp = rnd.random() < 0.5
            open(os.path.join(d, 's.fa'), 'w').write('>s\n' + s + '\n')
            o = os.path.join(d, 'o.xls')
            if os.path.exists(o): os.remove(o)
            args = ['fuzztran', 's.fa', '-pattern', pat, '-frame', frame, '-rformat', 'excel', '-outfile', 'o.xls', '-auto'] if tran else ['fuzzpro', 's.fa', '-pattern', pat, '-rformat', 'excel', '-outfile', 'o.xls', '-auto'] if prot else ['fuzznuc', 's.fa', '-pattern', pat, '-complement', 'Y' if comp else 'N', '-rformat', 'excel', '-outfile', 'o.xls', '-auto']
            if mm: args += ['-pmismatch', str(mm)]
            try:
                p = subprocess.run(args, cwd=d, capture_output=True, text=True, timeout=20)
            except subprocess.TimeoutExpired:
                continue
            if p.returncode != 0 or 'Illegal' in p.stderr + p.stdout or 'Fatal' in p.stderr:
                continue
            rep = open(o, encoding='latin-1').read() if os.path.exists(o) else ''
            rows.append({'id': len(rows), 'program': 'fuzztran' if tran else 'fuzzpro' if prot else 'fuzznuc', 'frame': frame, 'seq': '>s\n' + s + '\n', 'pattern': pat, 'mismatch': mm, 'complement': comp and not prot, 'report': rep})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 20261223)
