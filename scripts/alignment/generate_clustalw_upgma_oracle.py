#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A2-14, B09 F101) oracle: ClustalW 2.1 `-CLUSTERING=UPGMA -ALIGN` (the program's own UPGMA guide tree, `src/tree/UPGMA`).

The installed Debian `clustalw` (2.1+lgpl-7) segmentation-faults on every UPGMA run, and so does the unmodified source built at -O2 (undefined behaviour of the program); the same unmodified
source built at -O0 runs to the end. The user ran docker package `docs/Validation/review-2026-09/handoff/A2-14-clustalw-upgma/` (600 sets; build of the Debian source, O0 and O2) and returned
`returned/A2-14-clustalw-upgma-out.zip`. This script only converts that archive's `out/clustalw_upgma_oracle.tsv` into the committed test data (nothing models the binary):

  id  type(dna|protein)  seqs(;)  options(,)  rows(;)  dnd(the guide-tree text, newlines as '|')

Usage: python3 -I scripts/alignment/generate_clustalw_upgma_oracle.py RETURNED_OUT_TSV OUTFILE
       (RETURNED_OUT_TSV = out/clustalw_upgma_oracle.tsv of the extracted archive; data, never executed)
"""
import csv, sys


def main(src, dst):
    n = 0
    with open(src, newline='') as f, open(dst, 'w', newline='\n') as o:
        r = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        header = next(r)
        assert header[:9] == ['id', 'type', 'seqs', 'options', 'binary', 'rc', 'rows', 'tree', 'stdout'], header
        o.write('id\ttype\tseqs\toptions\trows\tdnd\n')
        for row in r:
            sid, t, seqs, opts, binary, rc, rows, tree = row[:8]
            if rc != '0' or not rows:
                continue          # never happens in the returned archive (600 sets with rows, 0 without)
            o.write('\t'.join([sid, t, seqs, opts, rows, tree]) + '\n')
            n += 1
    print(n, 'sets written')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
