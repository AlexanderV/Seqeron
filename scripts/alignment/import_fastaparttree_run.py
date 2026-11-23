#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/mafft_fastaparttree_oracle.tsv.gz from a returned Docker run of `mafft --fastaparttree` (docs/Validation/review-2026-09/handoff/A2-19c-fastaparttree/returned/*.zip, run 3: MAFFT 7.505-1 +
fasta34 built from fasta-34.26.5; the archive is untrusted data: it is read, never executed). Columns: id, options (with --quiet, as run), seqs ('|' separated), rc, rows ('|' separated).
Usage: python3 -I scripts/alignment/import_fastaparttree_run.py <returned.zip> [<more.zip> ...]   (every archive is unpacked into its own temporary directory)"""
import gzip, os, sys, tempfile, zipfile

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'mafft_fastaparttree_oracle.tsv.gz')
rows = []
for archive in sys.argv[1:]:
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(archive) as z:
            for n in z.namelist():
                assert not n.startswith('/') and '..' not in n, n
            z.extractall(td)
        path = os.path.join(td, 'out', 'results.tsv')
        for line in open(path, encoding='latin-1'):
            f = line.rstrip('\n').split('\t')
            if len(f) == 7:        # run 3: id, options, seqs, rc, rows, sha256 of stderr, first 200 bytes (all runs used --quiet)
                rows.append([f[0], f[1] + ' --quiet', f[2], f[3], f[4]])
            else:
                raise SystemExit('unknown format of ' + path)
with gzip.open(out_path, 'wt', encoding='utf-8') as f:
    f.write('\n'.join('\t'.join(r) for r in rows) + '\n')
print('wrote', len(rows), 'rows to', out_path)
