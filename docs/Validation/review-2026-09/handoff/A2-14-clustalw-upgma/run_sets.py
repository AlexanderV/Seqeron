#!/usr/bin/env python3
"""Runs ClustalW 2.1 `-CLUSTERING=UPGMA -ALIGN` (the program's own UPGMA guide tree) on every set of inputs/sets.tsv with the O0 build (falls back to O2 when O0 fails)
and writes out/clustalw_upgma_oracle.tsv. Never aborts: a failing run is recorded (rc, 'NOOUT').

Columns: id  type  seqs  options  binary  rc  rows(';')  tree(the .dnd text, newlines as '|')  stdout(first 300 chars)
"""
import os, subprocess, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BINS = [b for b in ('/opt/cw/clustalw_O0', '/opt/cw/clustalw_O2') if os.path.exists(b)]


def fasta_rows(path):
    rows, cur = [], None
    for l in open(path, errors='replace'):
        l = l.rstrip('\n')
        if l.startswith('>'):
            rows.append(''); cur = len(rows) - 1
        elif cur is not None:
            rows[cur] += l.strip()
    return rows


def run(binary, d, t, seqs, opts):
    open(d + '/x.fa', 'w').write(''.join(f'>s{j:03d}\n{s}\n' for j, s in enumerate(seqs)))
    for f in ('o.fa', 'o.dnd', 'x.dnd'):
        if os.path.exists(d + '/' + f):
            os.remove(d + '/' + f)
    cmd = [binary, '-INFILE=' + d + '/x.fa', '-TYPE=' + t.upper(), '-ALIGN', '-CLUSTERING=UPGMA', '-OUTPUT=FASTA', '-OUTFILE=' + d + '/o.fa', '-NEWTREE=' + d + '/o.dnd', '-OUTORDER=INPUT', '-QUIET']
    for kv in opts.split(',') if opts else []:
        k, _, v = kv.partition('=')
        cmd.append('-' + k.upper() + ('=' + v if v else ''))
    p = subprocess.run(cmd, capture_output=True, timeout=300, cwd=d)
    rows = fasta_rows(d + '/o.fa') if os.path.exists(d + '/o.fa') and os.path.getsize(d + '/o.fa') > 0 else []
    tp = d + '/o.dnd' if os.path.exists(d + '/o.dnd') and os.path.getsize(d + '/o.dnd') > 0 else d + '/x.dnd'
    tree = open(tp, errors='replace').read().replace('\n', '|') if os.path.exists(tp) else ''
    return p.returncode, rows, tree, (p.stdout + p.stderr).decode('utf-8', errors='backslashreplace')[:300].replace('\t', ' ').replace('\n', '|')


def main():
    os.makedirs('out', exist_ok=True)
    ok = failed = 0
    with open(os.path.join(HERE, 'inputs', 'sets.tsv')) as f, open('out/clustalw_upgma_oracle.tsv', 'w') as o, tempfile.TemporaryDirectory() as d:
        f.readline()
        o.write('id\ttype\tseqs\toptions\tbinary\trc\trows\ttree\tstdout\n')
        for line in f:
            sid, t, seqs, opts = line.rstrip('\n').split('\t')
            result = None
            for b in BINS:
                try:
                    rc, rows, tree, msg = run(b, d, t, seqs.split(';'), opts)
                except Exception as e:
                    rc, rows, tree, msg = -999, [], '', str(e)[:200]
                result = (os.path.basename(b), rc, rows, tree, msg)
                if rc == 0 and rows:
                    break
            b, rc, rows, tree, msg = result if result else ('none', -1, [], '', 'no binary built')
            ok += rc == 0 and bool(rows); failed += not (rc == 0 and rows)
            o.write('\t'.join([sid, t, seqs, opts, b, str(rc), ';'.join(rows), tree, msg]) + '\n'); o.flush()
    print(f'{ok} sets with rows, {failed} without')


if __name__ == '__main__':
    main()
