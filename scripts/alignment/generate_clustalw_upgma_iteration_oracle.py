#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A2-14d, B09 F297) oracle: ClustalW 2.1 `-CLUSTERING=UPGMA` together with `-ITERATION=alignment|tree`, `-PROFILE` and `-SEQUENCES`
(`TreeInterface::getWeightsForProfileAlignUPGMA`, `Iteration::removeFirstIterate` / `iterationOnTreeNode` with UPGMA tree weights).

The installed Debian `clustalw` (2.1+lgpl-7) and an -O2 build of the same source segmentation-fault on every UPGMA run (undefined behaviour of the program); the unmodified source built at -O0
runs to the end. The user ran the docker package `docs/Validation/review-2026-09/handoff/A2-14-upgma-iteration/` (800 sets, `-O0` build `clustalw_O0`, sha256 1d8af3c9...d0aa0 in ENV.txt) and returned
`returned/A2-14-upgma-iteration-out.zip`. This script only converts that archive's `out/upgma_iteration_oracle.tsv` into the committed test data (nothing models the binary):

  id  mode(align-alignment|align-tree|profile|sequences)  type(dna|protein)  profile1(;)  profile2(;)  options(,)  rows(;)
  (align-*: profile1 = the unaligned sequences, profile2 empty; profile / sequences: the two parts as given to the program, n1 = the number of sequences of the first part)

Usage: python3 -I scripts/alignment/generate_clustalw_upgma_iteration_oracle.py RETURNED_OUT_TSV OUTFILE
       (RETURNED_OUT_TSV = out/upgma_iteration_oracle.tsv of the extracted archive; data, never executed)
"""
import csv, sys


def main(src, dst):
    n = 0
    with open(src, newline='') as f, open(dst, 'w', newline='\n') as o:
        r = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        header = next(r)
        assert header[:9] == ['id', 'mode', 'type', 'profile1', 'profile2', 'options', 'rc', 'rows', 'stdout'], header
        o.write('id\tmode\ttype\tprofile1\tprofile2\toptions\trows\n')
        for row in r:
            sid, mode, t, p1, p2, opts, rc, rows = row[:8]
            if rc != '0' or not rows:
                continue          # never happens in the returned archive (800 sets with rows, 0 without)
            o.write('\t'.join([sid, mode, t, p1, p2, opts, rows]) + '\n')
            n += 1
    print(n, 'sets written')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
