#!/usr/bin/env python3
"""A3-9c2 (B09 WP45) oracle of the ksw2 kernels that minimap2 reaches only through align.c: ksw_ll_i16 (the 16-bit striped local alignment of `mm_test_zdrop` / `mm_seed_ext_score` / `mm_align1_inv`) and ksw_exts2_sse (the spliced extension of `-x splice`):
the C functions of minimap2 v2.26 (https://raw.githubusercontent.com/lh3/minimap2/v2.26/ ksw2.h ksw2_ll_sse.c ksw2_exts2_sse.c) compiled with ksw2_oracle_driver.c (-msse4.1) and run on random inputs.
Output (TSV, stdout): Z2 <query> <target> <q> <e> <w> <zdrop> <end_bonus> <flag> <fields as S2>   D2 <query> <target> <q> <e> <q2> <e2> <w> <zdrop> <end_bonus> <flag> <fields>   LL <query> <target> <gapo> <gape> <score> <qe> <te>   and   S2 <query> <target> <q> <e> <q2> <noncan> <zdrop> <junc_bonus> <flag> <junc|-> <max> <zdropped> <max_q> <max_t> <mqe> <mqe_t> <mte> <mte_q> <score> <reach_end> <n_cigar> <cigar...>
Usage: python3 generate_ksw2_oracle.py [n_ll=500 [n_s2=500 [n_z2=500 per variant]]] > ksw2_oracle.tsv     (needs gcc and network access to raw.githubusercontent.com for the three sources, or a directory SRC with them in the environment variable KSW2_SRC)
"""
import os, random, subprocess, sys, tempfile, urllib.request

rnd = random.Random(20261224)
FILES = ['ksw2.h', 'ksw2_ll_sse.c', 'ksw2_exts2_sse.c', 'ksw2_extz2_sse.c', 'ksw2_extd2_sse.c']


def build():
    d = tempfile.mkdtemp()
    src = os.environ.get('KSW2_SRC')
    for f in FILES:
        if src:
            open(os.path.join(d, f), 'wb').write(open(os.path.join(src, f), 'rb').read())
        else:
            req = urllib.request.Request('https://raw.githubusercontent.com/lh3/minimap2/v2.26/' + f, headers={'User-Agent': 'Seqeron-review/1.0'})
            open(os.path.join(d, f), 'wb').write(urllib.request.urlopen(req).read())
    drv = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ksw2_oracle_driver.c')
    exe = os.path.join(d, 'drv')
    subprocess.check_call(['gcc', '-O1', '-msse4.1', '-I', d, '-o', exe, drv, os.path.join(d, 'ksw2_ll_sse.c'), os.path.join(d, 'ksw2_exts2_sse.c'), os.path.join(d, 'ksw2_extz2_sse.c'), os.path.join(d, 'ksw2_extd2_sse.c')])
    return exe


def rseq(n):
    return ''.join(rnd.choice('ACGT') for _ in range(n))


def mutate(s, rate):
    out = []
    for c in s:
        x = rnd.random()
        if x < rate / 2:
            out.append(rnd.choice('ACGTN' if rnd.random() < 0.05 else 'ACGT'))
        elif x < rate * 0.75:
            continue
        elif x < rate:
            out.append(c)
            out.extend(rnd.choice('ACGT') for _ in range(rnd.randint(1, 4)))
        else:
            out.append(c)
    return ''.join(out)


def main():
    n_ll = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    n_s2 = int(sys.argv[2]) if len(sys.argv) > 2 else 500
    n_z2 = int(sys.argv[3]) if len(sys.argv) > 3 else 500
    exe = build()
    lines = []
    keys = []
    for _ in range(n_ll):
        t = rseq(rnd.choice([1, 2, 5, 17, 40, 90, 200, 500]))
        q = mutate(rnd.choice([t, t[rnd.randint(0, len(t) // 2):], rseq(rnd.randint(1, 120))]), rnd.choice([0.0, 0.1, 0.3])) or 'A'
        if rnd.random() < 0.3:
            q = q[::-1]
        gapo, gape = rnd.choice([(4, 2), (6, 2), (2, 1), (12, 2), (39, 3)])
        lines.append('LL %s %s %d %d' % (q, t, gapo, gape))
        keys.append('LL\t%s\t%s\t%d\t%d' % (q, t, gapo, gape))
    for _ in range(n_s2):
        exons = [rseq(rnd.randint(20, 120)) for _ in range(rnd.randint(1, 4))]
        kind = rnd.choice(['GT..AG', 'CT..AC', 'GC..AG', 'AT..AC', 'rand'])
        target = exons[0]
        for e in exons[1:]:
            ln = rnd.randint(30, 400)
            body = rseq(ln)
            intr = {'GT..AG': 'GT' + body[2:-2] + 'AG', 'CT..AC': 'CT' + body[2:-2] + 'AC', 'GC..AG': 'GC' + body[2:-2] + 'AG', 'AT..AC': 'AT' + body[2:-2] + 'AC', 'rand': body}[kind]
            target += intr + e
        q = mutate(''.join(exons), rnd.choice([0.0, 0.02, 0.06]))
        if rnd.random() < 0.3:
            q = rseq(rnd.randint(0, 20)) + q
        target += rseq(rnd.randint(0, 40))
        if not q:
            q = 'A'
        if rnd.random() < 0.3:
            q, target = q[::-1], target[::-1]
        q1, e1, q2 = rnd.choice([(4, 2, 24), (2, 1, 32), (6, 2, 25), (12, 2, 40)])
        noncan = rnd.choice([0, 3, 9, 20])
        zdrop = rnd.choice([-1, 100, 200, 1000])
        jb = rnd.choice([0, 5, 9])
        flag = 0
        flag |= rnd.choice([0, 0x100, 0x200, 0x100, 0x200])
        flag |= rnd.choice([0, 0x400])
        flag |= rnd.choice([0, 0x800, 0x800])
        flag |= rnd.choice([0, 0x40, 0x40 | 0x80, 0x80, 0x1])        # EXTZ_ONLY, EXTZ_ONLY | REV_CIGAR, REV_CIGAR, SCORE_ONLY-like bits
        flag |= rnd.choice([0, 0x02])                                  # RIGHT
        flag |= rnd.choice([0, 0x08])                                  # APPROX_MAX
        if rnd.random() < 0.4:
            junc = ''.join(rnd.choice('0000000001248') if rnd.random() < 0.05 else '0' for _ in range(len(target)))
        else:
            junc = '-'
        lines.append('S2 %s %s %d %d %d %d %d %d %d %s' % (q, target, q1, e1, q2, noncan, zdrop, jb, flag, junc))
        keys.append('S2\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s' % (q, target, q1, e1, q2, noncan, zdrop, jb, flag, junc))
    for kind in ('Z2', 'D2'):
        for _ in range(n_z2):
            base = rseq(rnd.choice([5, 20, 60, 150, 400, 900]))
            q = mutate(base, rnd.choice([0.0, 0.05, 0.15, 0.3])) or 'A'
            t = mutate(base, rnd.choice([0.0, 0.05, 0.15, 0.3])) or 'C'
            if rnd.random() < 0.2:
                t = t + rseq(rnd.randint(0, 200))
            if rnd.random() < 0.15:
                q = q + rseq(rnd.randint(10, 300))
            flag = 0
            for bit, pr in ((0x1, 0.15), (0x2, 0.4), (0x8, 0.3), (0x10, 0.3), (0x40, 0.5), (0x80, 0.3), (0x4, 0.1)):
                if rnd.random() < pr:
                    flag |= bit
            w = rnd.choice([-1, -1, 1, 5, 30, 200])
            zdrop = rnd.choice([-1, 15, 50, 120, 400])
            eb = rnd.choice([-1, 0, 5, 10, 30])
            q1, e1 = rnd.choice([(4, 2), (6, 2), (2, 1), (12, 3), (0, 1)])
            if kind == 'Z2':
                lines.append('Z2 %s %s %d %d %d %d %d %d' % (q, t, q1, e1, w, zdrop, eb, flag))
                keys.append('Z2\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d' % (q, t, q1, e1, w, zdrop, eb, flag))
            else:
                q2, e2 = rnd.choice([(24, 1), (26, 1), (16, 0), (40, 1), (81, 1)])
                if not (e1 > e2 and q1 + e1 < q2 + e2):
                    q1, e1 = 4, 2
                lines.append('D2 %s %s %d %d %d %d %d %d %d %d' % (q, t, q1, e1, q2, e2, w, zdrop, eb, flag))
                keys.append('D2\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % (q, t, q1, e1, q2, e2, w, zdrop, eb, flag))
    res = subprocess.run([exe], input='\n'.join(lines) + '\n', capture_output=True, text=True, check=True).stdout.strip().split('\n')
    assert len(res) == len(keys)
    for k, r in zip(keys, res):
        print(k + '\t' + r.replace(' ', '\t', 10) if not k.startswith('LL') else k + '\t' + r.replace(' ', '\t'))


if __name__ == '__main__':
    main()
