#!/usr/bin/env python3
"""B09 WP68 (F310): converts the returned run 2 of the Docker package docs/Validation/review-2026-09/handoff/A2-29b-ksw2 (the real ksw2 / minimap2 functions of lh3/ksw2 master and minimap2 v2.26 on the 2 000 cases of cases.txt)
into the test data tests/.../TestData/Alignment/ksw2_real_oracle.tsv.gz (one row per case: the inputs of the harness line and the printed result of the real program).
Inputs (all read as data, nothing is executed): CASES (cases.txt of the package), OUT (the returned out/ directory: ksw2_sse41.txt, ksw2_sse2.txt, ksw2_scalar.txt, ksw2_scalar.txt.crashes, ksw2_oracle.txt, versions.txt).
Row columns: id variant flag w zdrop end_bonus q e q2 e2 m mat query target status max zdropped max_q max_t mqe mqe_t mte mte_q score n_cigar reach_end cigar sse2
  status = ok | crash:<exit code of the aborted scalar harness>;  result fields are '-' for a crashed case;  sse2 = '=' when the SSE2-only build (KSW_SSE2_ONLY) printed the same line as the SSE4.1 build, '-' for the scalar variants.
Usage: python3 -I convert_ksw2_returned.py CASES OUT_DIR OUTFILE.tsv.gz
"""
import gzip, sys

RESULT = 13


def read_lines(path):
    with open(path) as f:
        return [l.rstrip('\n') for l in f if l.strip()]


def main(cases_path, out_dir, outfile):
    cases = {}
    for line in read_lines(cases_path):
        t = line.split(' ')
        cid, variant, flag, w, zdrop, eb, q, e, q2, e2, m = t[:11]
        m = int(m)
        mat = t[11:11 + m * m]
        query, target = t[11 + m * m:11 + m * m + 2]
        cases[int(cid)] = [cid, variant, flag, w, zdrop, eb, q, e, q2, e2, str(m), ','.join(mat), query, target]
    sse41 = {int(l.split(' ')[0]): l for l in read_lines(out_dir + '/ksw2_sse41.txt')}
    sse2 = {int(l.split(' ')[0]): l for l in read_lines(out_dir + '/ksw2_sse2.txt')}
    scalar = {int(l.split(' ')[0]): l for l in read_lines(out_dir + '/ksw2_scalar.txt')}
    crashes = {int(l.split(' ')[0]): l.split(' ')[1] for l in read_lines(out_dir + '/ksw2_scalar.txt.crashes')}
    merged = {int(l.split(' ')[0]): l for l in read_lines(out_dir + '/ksw2_oracle.txt')}
    assert not read_lines(out_dir + '/ksw2_sse41.txt.crashes') and not read_lines(out_dir + '/ksw2_sse2.txt.crashes'), 'the SSE harnesses did not crash'
    assert len(cases) == 2000
    rows = []
    for cid in sorted(cases):
        c = cases[cid]
        variant = c[1]
        simd = variant.endswith('2')
        if simd:
            line = sse41[cid]
            assert merged[cid] == line
            same = '=' if sse2.get(cid) == line else sse2.get(cid, 'missing')
            status = 'ok'
        elif cid in crashes:
            assert cid not in scalar, cid
            line, same, status = None, '-', 'crash:' + crashes[cid]
        else:
            line = scalar[cid]
            assert merged[cid] == line
            same, status = '-', 'ok'
        if line is None:
            res = ['-'] * 12
        else:
            f = line.split(' ')
            assert f[0] == c[0] and f[1] == variant, (f[:2], c[:2])
            res = f[2:]
            assert len(res) == 12 and res[-1] == (res[-1] if res[-1] else '*')
        rows.append('\t'.join(c + [status] + res + [same]))
    prov = [
        '# Seqeron review B09 WP68: output of the REAL ksw2 programs (Docker package A2-29b-ksw2, run 2, 2026-11-23, x86_64, gcc 13.3.0 -O2): ksw_extz2_sse / ksw_extd2_sse of minimap2 v2.26 (SSE4.1 build; the SSE2-only build printed the same line, column sse2 = "=") and ksw_extz / ksw_extd of lh3/ksw2 master.',
        '# Sources: ksw2_extz2_sse.c, ksw2_extd2_sse.c sha256 as recorded for the port; mm2/ksw2.h of v2.26 (ce9dae2a...) differs from the recorded minimap2-master header (e6499bda...) only in KSW_EZ_SPLICE_SCORE, KSW_SPSC_OFFSET and the prototype of ksw_exts2_sse (not used here).',
        '# status crash:N = the scalar harness aborted (exit N) on this case: ksw2_extz.c:116 / ksw2_extd.c writes eh[j] with j > qlen (heap buffer overflow, ASan + valgrind); the result columns are "-".',
        '# columns: id variant flag w zdrop end_bonus q e q2 e2 m mat query target status max zdropped max_q max_t mqe mqe_t mte mte_q score n_cigar reach_end cigar sse2   (max is the unsigned 31-bit field ksw_extz_t.max; sequences are digit strings, "-" = empty)',
    ]
    data = ('\n'.join(prov) + '\n' + '\n'.join(rows) + '\n').encode()
    with open(outfile, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', compresslevel=9, mtime=0, filename='') as f:
        f.write(data)
    print(len(rows), 'rows,', len(crashes), 'crashed')


if __name__ == '__main__':
    main(*sys.argv[1:4])
