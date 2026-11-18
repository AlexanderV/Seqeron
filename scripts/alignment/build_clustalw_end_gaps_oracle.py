#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/clustalw_end_gaps_oracle.tsv (B09 audit A1-10d part 2) from the Docker package output.

Inputs : docs/Validation/review-2026-09/handoff/A1-10d/inputs/sets.tsv           (id, type, seqs, options, control_rows of the INSTALLED clustalw 2.1)
         docs/Validation/review-2026-09/handoff/A1-10d/returned/A1-10d-out.zip   (out/end_gaps_oracle.tsv: id, type, mode, rows; out/control_check.txt)
Output : id <TAB> type <TAB> mode <TAB> seqs <TAB> options <TAB> rows       (mode in control / useendgaps / endgappen / both;
         rows = ';'-joined aligned rows in input order, or ABORT when the binary wrote no alignment)

The package binary is clustalw 2.1+lgpl-7 (Debian source) with ONE edit: the `-ENDGAPS` handler reads CW_USEENDGAPS / CW_ENDGAPPEN and sets
UserParameters::useEndGaps / endGapPenalties; its control runs equal the installed binary on 340 / 340 sets (control_check.txt).
The user's run reported FAILED (no output file, exit status 0) for 34 sets in the modes endgappen and both: ClustalW's first-pass profile score was negative
(MSA::multiSeqAlign `score < 0` -> return -1); those rows are written as ABORT.
"""
import csv, io, os, sys, zipfile

root = os.path.join(os.path.dirname(__file__), '..', '..')
pkg = os.path.join(root, 'docs', 'Validation', 'review-2026-09', 'handoff', 'A1-10d')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'clustalw_end_gaps_oracle.tsv')

sets = {}
with open(os.path.join(pkg, 'inputs', 'sets.tsv')) as f:
    next(f)
    for line in f:
        p = line.rstrip('\n').split('\t')
        sets[p[0]] = (p[1], p[2], p[3] if len(p) > 3 else '', p[4] if len(p) > 4 else '')

with zipfile.ZipFile(os.path.join(pkg, 'returned', 'A1-10d-out.zip')) as z:
    names = z.namelist()
    tsv = next(n for n in names if n.endswith('end_gaps_oracle.tsv'))
    check = z.read(next(n for n in names if n.endswith('control_check.txt'))).decode()
    assert 'control rows equal the committed binary rows: 340 of 340' in check, check
    lines = z.read(tsv).decode().splitlines()

rows_out = []
for line in lines[1:]:
    id_, typ, mode, rows = line.split('\t')
    t, seqs, opts, control = sets[id_]
    assert t == typ
    if mode == 'control':
        assert rows == control, id_
    rows_out.append((id_, typ, mode, seqs, opts, 'ABORT' if rows == 'FAILED' else rows))

with open(out_path, 'w', newline='\n') as f:
    f.write('id\ttype\tmode\tseqs\toptions\trows\n')
    for r in rows_out:
        f.write('\t'.join(r) + '\n')
print(len(rows_out), 'rows;', sum(1 for r in rows_out if r[5] == 'ABORT'), 'ABORT')
