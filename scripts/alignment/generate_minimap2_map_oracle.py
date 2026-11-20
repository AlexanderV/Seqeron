#!/usr/bin/env python3
"""ALIGN-MULTI-001 oracle (audit A4-7): minimap2 2.26 up to the end of the chaining stage - the primary / secondary marking (`mm_set_parent`, `mm_select_sub`), the long-join re-chaining of `-r b,B` (`mg_lchain_rmq`),
`--rmq` chaining, `-N` / `-p` / `--mask-level` / `--mask-len` / `--hard-mask-level`, the divergence estimate (`mm_est_err`) and the mapping quality (`mm_set_mapq`) - vs Minimap2Chaining.Map.

Every case runs `minimap2 [options] ref.fa query.fa` (PAF without -c) and records every line: query start / end, strand, target name, target start / end, mlen, blen, mapq and the tags tp (P / S), cm, s1, s2 (primary
hits only), dv and rl. JSON lines:
  {"group": "...", "opts": {"k","w","g","r": [b, B], "n","m","N","p","mask_level","mask_len","hard_mask","rmq","chain_gap_scale","chain_skip_scale","max_chain_skip","max_chain_iter","no_long_join"? ...}, "names": [...], "qname": "...",
   "targets": [...], "query": "...", "hits": [[qs, qe, "+|-", tname, ts, te, mlen, blen, mapq, "P|S", cm, s1, s2|null, dv|null, rl], ...]}
Groups (>= 300 cases each): qname (audit A5-1: query names whose first byte is >= 0x80 - UTF-8 multi-byte names - or that contain U+0000, over references with 2 - 4 identical copies of the query so that the hash of the name decides the order of
equal-score chains: khash.h __ac_X31_hash_string runs over a signed char, and a C string ends at NUL), sec (secondary-chain marking: repeats on both strands, -N -p --mask-level --mask-len --hard-mask-level, -r b,b), lj (long-join, -r b,B with B > b on queries made of distant pieces), rmq (--rmq), misc (-g -n -m
--chain-gap-scale --chain-skip-scale --max-chain-skip --max-chain-iter on related pairs, -r b,b).
Usage: python3 scripts/alignment/generate_minimap2_map_oracle.py [per_group=330] > minimap2_map_oracle.jsonl
"""
import json, random, subprocess, sys, tempfile

rnd = random.Random(20261206)
QNAME_SEED = 20261207  # the qname group (A5-1) re-seeds the generator, so the four earlier groups are reproduced exactly
COMP = str.maketrans('ACGT', 'TGCA')


def rseq(n): return ''.join(rnd.choice('ACGT') for _ in range(n))
def revcomp(s): return s.translate(COMP)[::-1]


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice('ACGT'))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice('ACGT'))
        else: out.append(c)
    return ''.join(out)


def base():
    return {'k': 15, 'w': 10, 'g': None, 'r': None, 'n': None, 'm': None, 'N': None, 'p': None, 'mask_level': None, 'mask_len': None, 'hard_mask': False, 'rmq': False,
            'chain_gap_scale': None, 'chain_skip_scale': None, 'max_chain_skip': None, 'max_chain_iter': None}


def command(o):
    c = ['minimap2', '-k', str(o['k']), '-w', str(o['w'])]
    for key, flag in (('g', '-g'), ('n', '-n'), ('m', '-m'), ('N', '-N'), ('p', '-p')):
        if o[key] is not None: c += [flag, str(o[key])]
    if o['r'] is not None: c += ['-r', '%d,%d' % tuple(o['r'])]
    for key, flag in (('mask_level', '--mask-level'), ('mask_len', '--mask-len'), ('chain_gap_scale', '--chain-gap-scale'), ('chain_skip_scale', '--chain-skip-scale'), ('max_chain_skip', '--max-chain-skip'), ('max_chain_iter', '--max-chain-iter')):
        if o[key] is not None: c += [flag + '=' + str(o[key])]
    if o['hard_mask']: c.append('--hard-mask-level')
    if o['rmq']: c.append('--rmq=yes')
    return c


def run(o, names, targets, qname, query, d):
    open(d + '/r.fa', 'w').write(''.join('>%s\n%s\n' % (n, t) for n, t in zip(names, targets)))
    open(d + '/q.fa', 'wb').write(('>%s\n%s\n' % (qname, query)).encode('utf-8'))
    r = subprocess.run(command(o) + [d + '/r.fa', d + '/q.fa'], capture_output=True, text=True, encoding='utf-8', errors='replace')
    if r.returncode != 0: return None
    hits = []
    for line in r.stdout.split('\n'):
        f = line.split('\t')
        if len(f) < 12: continue
        tags = {t.split(':')[0]: t.split(':', 2)[2] for t in f[12:]}
        hits.append([int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10]), int(f[11]), tags['tp'], int(tags['cm']), int(tags['s1']),
                     int(tags['s2']) if 's2' in tags else None, tags.get('dv'), int(tags.get('rl', -1))])
    return hits


def random_qname():
    # code points outside ASCII (2-, 3- and 4-byte UTF-8), optionally with U+0000; no white space (kseq ends the name there), no surrogates
    pools = [(0xA1, 0x24F), (0x370, 0x3FF), (0x400, 0x4FF), (0x4E00, 0x9FFF), (0x3040, 0x30FF), (0x1F300, 0x1F5FF), (0x80 + 0x20, 0x7FF)]
    def ch():
        lo, hi = rnd.choice(pools)
        while True:
            c = chr(rnd.randint(lo, hi))
            if not c.isspace() and c.isprintable(): return c
    kind = rnd.choice(['lead', 'lead', 'lead', 'inner', 'nul-first', 'nul-inner', 'ascii'])
    ascii_part = ''.join(rnd.choice('abcdefgXYZ_/.019') for _ in range(rnd.randint(1, 5)))
    wide = ''.join(ch() for _ in range(rnd.randint(1, 4)))
    if kind == 'lead': return wide + ascii_part
    if kind == 'inner': return ascii_part + wide
    if kind == 'nul-first': return '\x00' + ascii_part + wide
    if kind == 'nul-inner': return wide + '\x00' + ascii_part
    return ascii_part


def make(group):
    o = base(); names = ['r0']; qname = rnd.choice(['q', 'read1', 'seq_7', 'query/1'])
    if group == 'qname':
        qname = random_qname()
        unit = rseq(rnd.randint(300, 600))
        copy_count = rnd.choice([2, 3, 3, 4])
        ref = ''.join(rseq(rnd.randint(50, 400)) + unit for _ in range(copy_count)) + rseq(rnd.randint(50, 300))
        query = unit if rnd.random() < 0.8 else revcomp(unit)
        if rnd.random() < 0.3: o['N'] = rnd.choice([1, 2, 5, 10])
        if rnd.random() < 0.3: o['p'] = rnd.choice([0.5, 0.8, 1.0])
        return o, ['r0'], [ref], qname, query
    o['k'] = rnd.choice([11, 13, 15, 15, 15, 17, 19]); o['w'] = rnd.choice([3, 5, 10, 10, 10])
    if group == 'sec':
        o['r'] = [500, 500] if rnd.random() < 0.7 else [100, 100]
        unit = rseq(rnd.randint(150, 700))
        copies = [mutate(unit, rnd.choice([0, 0.01, 0.04, 0.1])) for _ in range(rnd.randint(2, 5))]
        ref = ''.join(rseq(rnd.randint(30, 400)) + (revcomp(c) if rnd.random() < 0.35 else c) for c in copies) + rseq(rnd.randint(30, 300))
        targets = [ref] + ([rseq(rnd.randint(300, 1500)) + mutate(unit, 0.03)] if rnd.random() < 0.3 else [])
        names = ['r%d' % i for i in range(len(targets))]
        query = rseq(rnd.randint(0, 80)) + mutate(unit, rnd.choice([0.01, 0.05, 0.1])) + rseq(rnd.randint(0, 80))
        if rnd.random() < 0.25: query = revcomp(query)
        if rnd.random() < 0.5: o['N'] = rnd.choice([0, 1, 2, 3, 5, 10, 50])
        if rnd.random() < 0.5: o['p'] = rnd.choice([0, 0.2, 0.5, 0.8, 0.95, 1.0])
        if rnd.random() < 0.4: o['mask_level'] = rnd.choice([0.0, 0.1, 0.3, 0.5, 0.8, 1.0])
        if rnd.random() < 0.2: o['mask_len'] = rnd.choice([0, 20, 100, 500])
        if rnd.random() < 0.15: o['hard_mask'] = True
        if rnd.random() < 0.3: o['m'] = rnd.choice([20, 30, 40, 60]); o['n'] = rnd.choice([2, 3, 5])
    elif group == 'lj':
        # distant pieces of one long region on the target: re-chain with -r b,B (B > b)
        o['r'] = [rnd.choice([100, 200, 500]), rnd.choice([1000, 5000, 20000])]
        segs = [rseq(rnd.randint(150, 500)) for _ in range(rnd.randint(2, 5))]
        gaps = [rseq(rnd.randint(200, 3000)) for _ in range(len(segs))]
        ref = ''.join(g + s for g, s in zip(gaps, segs)) + rseq(rnd.randint(50, 500))
        err = rnd.choice([0.01, 0.04, 0.1])
        q = ''
        for s in segs: q += mutate(s, err) + rseq(rnd.randint(0, rnd.choice([30, 300, 1500])))
        query = q if rnd.random() < 0.7 else revcomp(q)
        targets = [ref]
        o['k'] = rnd.choice([11, 13, 15, 15]); 
        if rnd.random() < 0.3: o['N'] = rnd.choice([0, 2, 5])
        if rnd.random() < 0.3: o['m'] = rnd.choice([30, 40, 80])
    elif group == 'rmq':
        o['rmq'] = True
        o['r'] = rnd.choice([[500, 500], [1000, 100000], [200, 20000], [500, 20000]])
        segs = [rseq(rnd.randint(150, 600)) for _ in range(rnd.randint(2, 5))]
        ref = ''.join(rseq(rnd.randint(100, 2500)) + s for s in segs) + rseq(rnd.randint(50, 500))
        err = rnd.choice([0.01, 0.04, 0.1])
        q = ''.join(mutate(s, err) + rseq(rnd.randint(0, rnd.choice([30, 300, 1500]))) for s in segs)
        if rnd.random() < 0.3:  # an inversion
            i = rnd.randrange(len(q)); j = min(len(q), i + rnd.randint(100, 500)); q = q[:i] + revcomp(q[i:j]) + q[j:]
        query = q
        targets = [ref] + ([mutate(ref[:2000], 0.05)] if rnd.random() < 0.3 else [])
        names = ['r%d' % i for i in range(len(targets))]
        if rnd.random() < 0.3: o['N'] = rnd.choice([0, 2, 5, 50])
        if rnd.random() < 0.2: o['max_chain_skip'] = rnd.choice([1, 5, 25, 100])
    else:  # misc
        o['r'] = [rnd.choice([30, 100, 500]), 0]; o['r'][1] = o['r'][0]
        o['g'] = rnd.choice([100, 400, 1000, 5000, 5000])
        o['m'] = rnd.choice([20, 30, 40, 40, 80]); o['n'] = rnd.choice([2, 3, 3, 5])
        if rnd.random() < 0.5: o['chain_gap_scale'] = rnd.choice([0.0, 0.3, 0.8, 0.8, 1.5, 3.0])
        if rnd.random() < 0.4: o['chain_skip_scale'] = rnd.choice([0.0, 0.2, 0.5, 1.0])
        if rnd.random() < 0.4: o['max_chain_skip'] = rnd.choice([1, 3, 10, 25, 25])
        if rnd.random() < 0.3: o['max_chain_iter'] = rnd.choice([1, 5, 50, 5000])
        segs = [rseq(rnd.randint(60, 400)) for _ in range(rnd.randint(2, 4))]
        ref = ''.join(s + rseq(rnd.randint(10, 150)) for s in segs)
        err = rnd.choice([0.01, 0.03, 0.08, 0.15])
        q = ''.join((lambda p: revcomp(p) if rnd.random() < 0.25 else p)(mutate(s, err)) + rseq(rnd.randint(0, 120)) for s in segs)
        query = q; targets = [ref]
        if rnd.random() < 0.3: targets.append(rseq(rnd.randint(200, 800)) + mutate(segs[0], 0.05))
        names = ['r%d' % i for i in range(len(targets))]
    return o, names, targets, qname, query


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 330
    counts = {}
    with tempfile.TemporaryDirectory() as d:
        for group in ['sec', 'lj', 'rmq', 'misc', 'qname']:
            if group == 'qname': rnd.seed(QNAME_SEED)
            made = attempts = 0
            while made < per and attempts < per * 4:
                attempts += 1
                o, names, targets, qname, query = make(group)
                if len(query) < 40: continue
                hits = run(o, names, targets, qname, query, d)
                if hits is None or not hits: continue
                print(json.dumps({'group': group, 'opts': o, 'names': names, 'qname': qname, 'targets': targets, 'query': query, 'hits': hits}))
                made += 1
            counts[group] = made
    sys.stderr.write('rows %s\n' % counts)


main()
