#!/usr/bin/env python3
"""Runs ClustalW 2.1 (the -O0 build) with -CLUSTERING=UPGMA on every set of inputs/sets.tsv and writes out/upgma_iteration_oracle.tsv. Never aborts: a failing run is recorded (rc, empty rows).

Columns: id  mode  type  profile1(;)  profile2(;)  options  rc  rows(;)  stdout(first 300 chars)
`align-*`: profile1 = the sequences (unaligned), profile2 empty; command `-ALIGN -CLUSTERING=UPGMA -ITERATION=alignment|tree [-NUMITER]`. `profile` / `sequences`: profile 1 and 2 are the aligned rows produced
by the program itself (`-ALIGN`, NJ, same options minus the command-specific ones) / the unaligned sequences of profile 2 (`sequences`, unaligned); command `-PROFILE1 -PROFILE2 -PROFILE|-SEQUENCES -CLUSTERING=UPGMA`.
"""
import os, subprocess, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.environ.get('CW_BIN', '/opt/cw/clustalw_O0')
CLUSTERING = os.environ.get('CW_CLUSTERING', 'UPGMA')   # only for a plumbing test with a binary that cannot run UPGMA
LIMIT = int(os.environ.get('CW_LIMIT', '0'))
FLAGS = {'gapopen': 'GAPOPEN', 'gapext': 'GAPEXT', 'maxdiv': 'MAXDIV', 'noweights': 'NOWEIGHTS', 'matrix': 'MATRIX', 'quicktree': 'QUICKTREE', 'numiter': 'NUMITER'}


def read_fasta(path):
    rows, cur = [], None
    for l in open(path, errors='replace'):
        l = l.rstrip('\n')
        if l.startswith('>'): rows.append(''); cur = len(rows) - 1
        elif cur is not None: rows[cur] += l.strip()
    return rows


def write_fasta(path, rows, first=0):
    open(path, 'w').write(''.join(f'>s{first + j:03d}\n{s}\n' for j, s in enumerate(rows)))


def opt_args(opts):
    return ['-' + FLAGS[k] + ('=' + v if v else '') for k, _, v in (kv.partition('=') for kv in opts) if k in FLAGS]


def run(cmd, d):
    p = subprocess.run(cmd, capture_output=True, timeout=300, cwd=d)
    return p.returncode, (p.stdout + p.stderr).decode('utf-8', errors='backslashreplace')[:300].replace('\t', ' ').replace('\n', '|')


def main():
    os.makedirs('out', exist_ok=True)
    ok = bad = 0
    with open(os.path.join(HERE, 'inputs', 'sets.tsv')) as f, open('out/upgma_iteration_oracle.tsv', 'w') as o, tempfile.TemporaryDirectory() as d:
        f.readline()
        o.write('id\tmode\ttype\tprofile1\tprofile2\toptions\trc\trows\tstdout\n')
        for n, line in enumerate(f):
            if LIMIT and n >= LIMIT: break
            sid, mode, t, seqs, options = line.rstrip('\n').split('\t')
            seqs = seqs.split(';')
            opts = options.split(',') if options else []
            kv = {k: v for k, _, v in (x.partition('=') for x in opts)}
            for fn in os.listdir(d): os.remove(os.path.join(d, fn))
            rc, rows, p1, p2, msg = -999, [], [], [], ''
            try:
                base = [BIN, f'-TYPE={t.upper()}', '-OUTPUT=FASTA', '-OUTORDER=INPUT', '-QUIET']
                if mode.startswith('align'):
                    write_fasta(d + '/x.fa', seqs)
                    it = 'alignment' if mode == 'align-alignment' else 'tree'
                    rc, msg = run(base + [f'-INFILE={d}/x.fa', '-ALIGN', '-CLUSTERING=' + CLUSTERING, f'-ITERATION={it}', f'-OUTFILE={d}/o.fa'] + opt_args([x for x in opts if not x.startswith('n1=')]), d)
                    rows = read_fasta(d + '/o.fa') if os.path.exists(d + '/o.fa') else []
                    p1 = seqs
                else:
                    n1 = int(kv['n1'])
                    plain = [x for x in opts if not x.startswith('n1=') and not x.startswith('numiter')]
                    for part, name, first in ((seqs[:n1], 'a', 0), (seqs[n1:], 'b', n1)):
                        write_fasta(f'{d}/{name}.raw.fa', part, first)
                        if len(part) > 1 and (name == 'a' or mode == 'profile'):
                            subprocess.run(base + [f'-INFILE={d}/{name}.raw.fa', '-ALIGN', f'-OUTFILE={d}/{name}.fa'] + opt_args(plain), capture_output=True, timeout=300, cwd=d)
                            rows_part = read_fasta(f'{d}/{name}.fa') if os.path.exists(f'{d}/{name}.fa') else part
                        else:
                            rows_part = part
                        write_fasta(f'{d}/{name}.fa', rows_part, first)
                        if name == 'a': p1 = rows_part
                        else: p2 = rows_part
                    rc, msg = run(base + [f'-PROFILE1={d}/a.fa', f'-PROFILE2={d}/b.fa', '-PROFILE' if mode == 'profile' else '-SEQUENCES', '-CLUSTERING=' + CLUSTERING, f'-OUTFILE={d}/o.fa'] + opt_args(plain), d)
                    rows = read_fasta(d + '/o.fa') if os.path.exists(d + '/o.fa') else []
            except Exception as e:
                rc, msg = -999, str(e)[:200]
            good = rc == 0 and bool(rows)
            ok += good; bad += not good
            o.write('\t'.join([sid, mode, t, ';'.join(p1), ';'.join(p2), options, str(rc), ';'.join(rows), msg]) + '\n'); o.flush()
    print(f'{ok} sets with rows, {bad} without')


if __name__ == '__main__':
    main()
