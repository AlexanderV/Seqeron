#!/usr/bin/env python3
"""GENOMIC-SIMILARITY-001 oracle (A3-11): `mash dist` 2.3 (Debian mash 2.3+dfsg-6) vs GenomicAnalyzer.CompareSimilarity(sketchSize) on random related sequences.

For every case: two sequences (200-6000 nt, the second a point-mutated / indel copy of the first or an unrelated one, optional N runs), k in 3..31, sketch size s in 10..1000,
canonical k-mers (mash default) or `-n`; `mash sketch -k K -s S [-n] -o X seq.fa` for both, `mash dist a.msh b.msh` -> distance, p-value, shared/denominator.
TSV (tests/.../TestData/Analysis/mash_dist_oracle.tsv):  seq1 \t seq2 \t k \t s \t canonical(1|0) \t distance \t p-value \t shared \t denominator
(mash prints distance / p-value with 5 significant digits; the test compares at that precision).
Never aborts: a failing mash call is reported on stderr and the row skipped.
Install route of the reference program: `apt-get install mash` (Debian 2.3+dfsg-6); no source tarball is committed to the repository.
Usage: python3 scripts/analysis/generate_mash_dist_oracle.py [rows] > mash_dist_oracle.tsv
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261122)


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate:
            out.append(rnd.choice('ACGT'))
        elif r < rate * 1.1:
            continue
        elif r < rate * 1.2:
            out.append(c)
            out.append(rnd.choice('ACGT'))
        else:
            out.append(c)
    return ''.join(out)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    done = 0
    with tempfile.TemporaryDirectory() as tmp:
        for i in range(n):
            L = rnd.choice([200, 500, 1000, 3000, 6000])
            a = ''.join(rnd.choice('ACGT') for _ in range(L))
            if rnd.random() < 0.15:
                b = ''.join(rnd.choice('ACGT') for _ in range(rnd.choice([L, L // 2, L * 2])))
            else:
                b = mutate(a, rnd.choice([0.0, 0.005, 0.01, 0.03, 0.08, 0.15, 0.3]))
            if rnd.random() < 0.2:
                p = rnd.randrange(len(a)); a = a[:p] + 'N' * rnd.randint(1, 20) + a[p:]
            k = rnd.choice([3, 5, 7, 9, 11, 13, 15, 16, 17, 19, 21, 25, 31])
            s = rnd.choice([10, 50, 100, 200, 500, 1000])
            canon = rnd.random() < 0.7
            ok = True
            for name, seq in (('a', a), ('b', b)):
                with open(os.path.join(tmp, name + '.fa'), 'w') as f:
                    f.write('>' + name + '\n' + seq + '\n')
                cmd = ['mash', 'sketch', '-k', str(k), '-s', str(s)] + ([] if canon else ['-n']) + ['-o', os.path.join(tmp, name), os.path.join(tmp, name + '.fa')]
                r = subprocess.run(cmd, capture_output=True, text=True)
                if r.returncode != 0:
                    sys.stderr.write('FAIL sketch %s: %s\n' % (cmd, r.stderr.strip()[:200])); ok = False
            if not ok:
                continue
            r = subprocess.run(['mash', 'dist', os.path.join(tmp, 'a.msh'), os.path.join(tmp, 'b.msh')], capture_output=True, text=True)
            t = r.stdout.split()
            if r.returncode != 0 or len(t) < 5:
                sys.stderr.write('FAIL dist k=%d s=%d: %s\n' % (k, s, r.stderr.strip()[:200])); continue
            print('\t'.join([a, b, str(k), str(s), '1' if canon else '0', t[2], t[3], t[4].split('/')[0], t[4].split('/')[1]]))
            done += 1
    sys.stderr.write('rows: %d\n' % done)


if __name__ == '__main__':
    main()
