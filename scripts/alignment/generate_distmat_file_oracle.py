#!/usr/bin/env python3
"""A5-8b (B09 F184) oracle: the whole `-outfile` text of EMBOSS 6.6.0 `distmat` (title, method line, `Using ambiguity codes`, `Using base positions`, `Gap weighting is`, the matrix with names), compared as TEXT by
MultipleAlignmentStatistics.DistanceMatrixFile on random gapped alignments with every -nucmethod / -protmethod / -gapweight / -ambiguous / -position (the alignments of generate_distmat_oracle.py).

JSONL.GZ rows: rows (list), names, kind (dna|protein), method, gapweight, position, ambiguous, calculatea, parametera, report (the file text).
Usage: python3 scripts/alignment/generate_distmat_file_oracle.py OUT.jsonl.gz [n=300] [seed]   (needs distmat on PATH; rows whose run fails or times out are skipped)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_distmat_oracle as g  # noqa: E402


def main(out, n, seed):
    g.rnd = random.Random(seed)
    rnd = g.rnd
    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            kind = rnd.choice(['dna', 'protein'])
            msa = g.make_msa(kind, gaps=rnd.random() < 0.5, ambig=rnd.random() < 0.5)
            method = rnd.choice([0, 1, 2] if kind == 'protein' else [0, 1, 2, 3, 4, 5])
            gapw = rnd.choice([0, 0.5, 1, 0.25]); pos = rnd.choice([123, 123, 12, 1, 2, 3]) if kind == 'dna' else rnd.choice([123, 12])
            if method in (2, 4) and kind == 'protein' and pos != 123: continue    # protein Kimura / Tajima-Nei read past their buffers for another position
            if method == 4 and pos != 123: continue
            amb = rnd.random() < 0.3 and (kind == 'protein' or not any('-' in r for r in msa)); calca = method == 5 and rnd.random() < 0.5; parama = rnd.choice([1.0, 0.5, 2.0])
            names = [rnd.choice(['s%d', 'seq_%d', 'X%05d.1', 'a-b%d']) % (i + 1) if True else '' for i in range(len(msa))]
            names = [n_ % (i + 1) if '%' in n_ else n_ for i, n_ in enumerate(names)]
            open(d + '/a.fa', 'w').write(''.join(f'>{nm}\n{r}\n' for nm, r in zip(names, msa)))
            cmd = ['distmat', '-sequence', d + '/a.fa', '-outfile', d + '/o.txt', '-auto', '-gapweight', str(gapw), '-position', str(pos), '-nucmethod', str(method), '-protmethod', str(method if method <= 2 else 0), '-parametera', str(parama)]
            if amb: cmd.append('-ambiguous')
            if calca: cmd.append('-calculatea')
            try:
                p = subprocess.run(cmd, capture_output=True, timeout=4)
            except subprocess.TimeoutExpired:
                continue
            if p.returncode != 0 or not os.path.exists(d + '/o.txt'): continue
            rows.append({'rows': msa, 'names': names, 'kind': kind, 'method': method, 'gapweight': gapw, 'position': pos, 'ambiguous': amb, 'calculatea': calca, 'parametera': parama,
                         'report': open(d + '/o.txt', encoding='latin-1').read()})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n') as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 20261222)
