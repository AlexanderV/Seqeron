#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/muscle3_r2_oracle.tsv (B09 audit A1-10c part a2) from the ROUND 2 Docker package output.

Inputs : docs/Validation/review-2026-09/handoff/A1-10c-muscle-r2/inputs/sets.tsv                      (id, kind, type, data, options)
         docs/Validation/review-2026-09/handoff/A1-10c-muscle-r2/inputs/matrices/*.txt                  (the matrix files named by @blosum62 / @blosum62nostar / @pam250)
         docs/Validation/review-2026-09/handoff/A1-10c-muscle-r2/returned/A1-10c-muscle-r2-out.zip      (out/muscle3_r2_oracle.tsv, out/ENV.txt: MUSCLE v3.8.1551)
Output : id kind type options data rc names rows tree1 tree2 dist1 dist2 weights usedtree stdout
         = the oracle row of the package with the input data (data) and the matrix file text substituted into the options (matrix=@name -> matrix=<TEXT:name>, the text is appended as
         the last column `matrixtext`), so the test needs no other file. The test (to be written when the archive is returned) maps every kind to the port:
           align     SequenceAligner.MuscleRun(seqs, options) -> Alignment rows in the written order (names / OutputOrder), Tree1 / Tree2 / DistanceMatrix1 / 2 / WeightsFile equal the
                     files, rc != 0 -> the port must throw the same way (FAILED rows are the program's refusals, e.g. -stable)
           profile   MuscleProfileAlign        refine  MuscleRefine           spscore  MuscleScoreAlignment ('SP=%.4g' in stdout)   ppscore  MuscleProfileScore ('Score=%.4g')
           maketree  MuscleGuideTree           refinew MuscleRefineWindows    sw       MuscleLocalAlign                              profdb   (the program's behaviour is recorded)
         Protein ppscore / sw / refinew and profdb are expected to FAIL in the program (no score matrix / alphabet); the port rejects them with NotSupportedException.
"""
import os, re, sys, zipfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
pkg = os.path.join(root, 'docs', 'Validation', 'review-2026-09', 'handoff', 'A1-10c-muscle-r2')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'muscle3_r2_oracle.tsv')
zpath = os.path.join(pkg, 'returned', 'A1-10c-muscle-r2-out.zip')
if not os.path.exists(zpath):
    sys.exit('not returned yet: ' + zpath)

sets = {}
with open(os.path.join(pkg, 'inputs', 'sets.tsv')) as f:
    next(f)
    for line in f:
        sid, kind, typ, data, opts = line.rstrip('\n').split('\t')
        sets[sid] = (kind, typ, data, opts)
matrices = {'@' + n[:-4]: open(os.path.join(pkg, 'inputs', 'matrices', n)).read() for n in os.listdir(os.path.join(pkg, 'inputs', 'matrices'))}

with zipfile.ZipFile(zpath) as z:
    names = z.namelist()
    env = z.read(next(n for n in names if n.endswith('ENV.txt'))).decode()
    assert 'MUSCLE v3.8.1551' in env, env
    lines = z.read(next(n for n in names if n.endswith('muscle3_r2_oracle.tsv'))).decode().splitlines()

esc = lambda s: s.replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n')
rows = []
for line in lines[1:]:
    f = line.split('\t')
    assert len(f) == 14, (len(f), line[:80])
    sid, kind, typ, opts = f[0], f[1], f[2], f[3]
    k, t, data, o = sets[sid]
    assert (k, t, o) == (kind, typ, opts), sid
    mtext = ''
    for key, text in matrices.items():
        if re.search('=' + re.escape(key) + '(,|$)', opts):
            mtext = text
    rows.append([sid, kind, typ, opts, data] + f[4:] + [esc(mtext)])
with open(out_path, 'w', newline='\n') as fh:
    fh.write('id\tkind\ttype\toptions\tdata\trc\tnames\trows\ttree1\ttree2\tdist1\tdist2\tweights\tusedtree\tstdout\tmatrixtext\n')
    for r in rows:
        fh.write('\t'.join(r) + '\n')
print(len(rows), 'rows;', sum(1 for r in rows if not r[5].isdigit()), 'failed/no-output')
