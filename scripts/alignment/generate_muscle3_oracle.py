#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/muscle3_oracle.tsv (B09 audit A1-10c) from the Docker package output (MUSCLE v3.8.1551, user-run).

Inputs : docs/Validation/review-2026-09/handoff/A1-10c-muscle/inputs/sets.tsv                     (id, type, seqs, options; round 1, 420 sets)
         docs/Validation/review-2026-09/handoff/A1-10c-muscle/returned/A1-10c-muscle-out.zip       (out/muscle3_oracle.tsv: id, type, options, order, rows; out/ENV.txt; out/run.log)
         [round 2, when present] docs/Validation/review-2026-09/handoff/A1-10c-muscle-r2/ (same layout; its zip is A1-10c-muscle-r2-out.zip)
Output : id <TAB> type <TAB> options <TAB> seqs <TAB> order <TAB> rows
         (seqs = ';'-joined input sequences, order = the ids of the rows in the order MUSCLE wrote them, rows = ';'-joined aligned rows in INPUT order)
The files are parsed as data only (never executed); the binary is MUSCLE v3.8.1551 built unmodified from the Debian muscle3 3.8.1551-2 source (ENV.txt).
"""
import os, sys, zipfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
base = os.path.join(root, 'docs', 'Validation', 'review-2026-09', 'handoff')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'muscle3_oracle.tsv')

PACKAGES = [('A1-10c-muscle', 'A1-10c-muscle-out.zip', 'M'), ('A1-10c-muscle-r2', 'A1-10c-muscle-r2-out.zip', 'R')]
rows_out = []
for pkg, zipname, _ in PACKAGES:
    zpath = os.path.join(base, pkg, 'returned', zipname)
    if not os.path.exists(zpath):
        print('skip (not returned yet):', pkg)
        continue
    sets = {}
    with open(os.path.join(base, pkg, 'inputs', 'sets.tsv')) as f:
        next(f)
        for line in f:
            p = line.rstrip('\n').split('\t')
            sets[p[0]] = (p[1], p[2], p[3] if len(p) > 3 else '')
    with zipfile.ZipFile(zpath) as z:
        names = z.namelist()
        env = z.read(next(n for n in names if n.endswith('ENV.txt'))).decode()
        assert 'MUSCLE v3.8.1551' in env, env
        tsv = z.read(next(n for n in names if n.endswith('muscle3_oracle.tsv'))).decode().splitlines()
    n = 0
    for line in tsv[1:]:
        f = line.split('\t')
        sid, typ, opts, order, rows = f[0], f[1], f[2], f[3], f[4]
        t, seqs, o = sets[sid]
        assert t == typ and o == opts, sid
        if rows.startswith('FAILED'):
            rows = 'FAILED'
        rows_out.append((sid, typ, opts, seqs, order, rows))
        n += 1
    assert n == len(sets), (pkg, n, len(sets))
    print(pkg, n, 'sets')

with open(out_path, 'w', newline='\n') as f:
    f.write('id\ttype\toptions\tseqs\torder\trows\n')
    for r in rows_out:
        f.write('\t'.join(r) + '\n')
print(len(rows_out), 'rows written;', sum(1 for r in rows_out if r[5] == 'FAILED'), 'FAILED')
