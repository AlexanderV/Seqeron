#!/usr/bin/env python3
"""Maintainer tool of the hand-off package (work package 1b.9f): collects the VERIFIED containers of a data directory (e8-verify exit code 0, the container hashed) into reference_hashes.tsv, the table that
the checker compares the owner's containers with. The builder's output depends on the input and on PSI only, not on the machine, the thread count or Resume, so a container that the owner builds from the
same text must be byte-identical to the one that was built and verified on the reference machine: that is the cross-architecture check of the arm64 run (x86-64 reference, arm64 result).
usage: s1b_refs.py <data directory> [--source TEXT] [--table reference_hashes.tsv]      (rows with the same kind / n / psi are replaced)"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    data = argv[0]
    opt = dict(zip(argv[1::2], argv[2::2]))
    table_path = opt.get('--table', os.path.join(HERE, 'reference_hashes.tsv'))
    head = ['kind', 'n', 'psi', 'container_sha256', 'container_bytes', 'input_sha256', 'source']
    rows = {}
    if os.path.exists(table_path):
        with open(table_path) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) == len(head):
                    rows[(f[0], f[1], f[2])] = f
    added = 0
    for d in sorted(glob.glob(os.path.join(data, 'runs', '*'))):
        try:
            meta = json.load(open(os.path.join(d, 'run.json')))
            sha, size = open(os.path.join(d, 'container.sha256')).read().split()[:2]
            vrc = int(open(os.path.join(d, 'verify.rc')).read().strip())
            sums = dict(line.split('=', 1) for line in open(os.path.join(d, 'input.sums')).read().splitlines() if '=' in line)
        except (OSError, ValueError, IndexError):
            continue
        if vrc != 0:
            continue
        key = (meta['kind'], str(meta['n']), 'psi' if meta.get('psi', True) else 'nopsi')
        row = [key[0], key[1], key[2], sha, size, sums.get('sha256', ''), opt.get('--source', 'the reference machine')]
        if key in rows and rows[key][3] != sha:
            print(f'WARNING: {key} had container {rows[key][3][:16]}..., this run has {sha[:16]}...: the builder changed, or the earlier row was wrong')
        rows[key] = row
        added += 1
    with open(table_path, 'w') as fh:
        fh.write('\t'.join(head) + '\n')
        for key in sorted(rows, key=lambda k: (k[0], int(k[1]), k[2])):
            fh.write('\t'.join(rows[key]) + '\n')
    print(f'{added} verified container(s) written to {table_path} ({len(rows)} rows)')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
