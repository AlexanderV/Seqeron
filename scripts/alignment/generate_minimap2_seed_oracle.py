#!/usr/bin/env python3
"""ALIGN-MULTI-001 oracle (audit A4-7): the SEED stage of minimap2 2.26 beyond the defaults - frequent-minimizer handling (`mm_seed_select`, `mid_occ`, `max_max_occ`, `occ_dist`, `mm_seed_mz_flt`), `-f` / `-U`, `-H`,
the hidden `-T` (sdust masking), `-X` (skip self / dual hits) and `--for-only` / `--rev-only` - vs Minimap2Chaining.CollectSeeds.

For every case `minimap2 [options] --print-seeds ref.fa query.fa` is run: stderr gives the anchor array (SD lines: target name, end on the target, strand, end on the query, span), the repeat length (RS line, the PAF `rl` tag)
and the occurrence limit (`mid_occ = N` of the log). JSON lines:
  {"group": "freq|fu|hpc|dust|ava|strand", "opts": {"k":..,"w":..,"H":bool,"T":int,"f":float|null,"fmax":int|null,"e":int|null,"U":[min,max]|null,"q_occ_frac":float|null,"X":bool,"for_only":bool,"rev_only":bool},
   "names": [...], "qname": "...", "targets": [...], "query": "...", "mid_occ": N, "rep_len": R, "seeds": [[tid, strand(0|1), tpos, qpos, span], ...]}
Groups (>= 300 cases each): freq = defaults on references with many repeats (minimizers above mid_occ, runs of them: mm_seed_select), fu = random -f / -U / --q-occ-frac, hpc = -H on homopolymer-rich sequences,
dust = -T 5..40 on low-complexity queries, ava = -X with targets of the query's name (same / different length) and names that sort before / after it, strand = --for-only / --rev-only.
Cases with more than 700 seeds are skipped. Never aborts.
Usage: python3 scripts/alignment/generate_minimap2_seed_oracle.py [per_group=330] > minimap2_seed_oracle.jsonl
"""
import json, os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261205)
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


def homopolymer_rich(n):
    out = []
    while len(out) < n:
        c = rnd.choice('ACGT'); out.extend(c * rnd.choice([1, 1, 2, 3, 4, 6, 9]))
    return ''.join(out[:n])


def low_complexity(n):
    parts = []
    while sum(map(len, parts)) < n:
        kind = rnd.random()
        if kind < 0.4: parts.append(rnd.choice(['AT', 'CA', 'GA', 'CAG', 'AAC', 'TTAGGG']) * rnd.randint(8, 40))
        elif kind < 0.55: parts.append(rnd.choice('ACGT') * rnd.randint(12, 50))
        else: parts.append(rseq(rnd.randint(30, 160)))
    return ''.join(parts)[:n]


def repeat_rich_reference():
    # unique stretches plus several repeat families copied 8-60 times with 0-6 % divergence
    fams = [rseq(rnd.randint(25, 140)) for _ in range(rnd.randint(1, 3))]
    parts = []
    for _ in range(rnd.randint(2, 4)): parts.append(rseq(rnd.randint(150, 700)))
    for f in fams:
        for _ in range(rnd.randint(8, 60)):
            parts.append(mutate(f, rnd.choice([0, 0.02, 0.06])) + rseq(rnd.randint(0, 40)))
    rnd.shuffle(parts)
    ref = ''.join(parts)
    return ref, fams


def query_from(ref, fams, repeat_bias=0.5):
    q = []
    for _ in range(rnd.randint(2, 5)):
        r = rnd.random()
        if r < repeat_bias and fams: q.append(mutate(rnd.choice(fams) * rnd.randint(1, 4), 0.03))
        elif len(ref) > 300:
            i = rnd.randrange(0, len(ref) - 200); q.append(mutate(ref[i:i + rnd.randint(80, min(500, len(ref) - i))], 0.04))
        q.append(rseq(rnd.randint(0, 60)))
    s = ''.join(q)
    return revcomp(s) if rnd.random() < 0.3 else s


def run(opts, names, targets, qname, query, d):
    open(d + '/r.fa', 'w').write(''.join('>%s\n%s\n' % (n, t) for n, t in zip(names, targets)))
    open(d + '/q.fa', 'w').write('>%s\n%s\n' % (qname, query))
    cmd = ['minimap2', '-k', str(opts['k']), '-w', str(opts['w'])]
    if opts['H']: cmd.append('-H')
    if opts['T']: cmd += ['-T', str(opts['T'])]
    if opts['f'] is not None: cmd += ['-f', repr(opts['f']) + (',%d' % opts['fmax'] if opts['fmax'] is not None else '')]
    if opts['e'] is not None: cmd += ['-e', str(opts['e'])]
    if opts['U'] is not None: cmd += ['-U', '%d,%d' % tuple(opts['U'])]
    if opts['q_occ_frac'] is not None: cmd += ['--q-occ-frac', repr(opts['q_occ_frac'])]
    if opts['X']: cmd.append('-X')
    if opts['for_only']: cmd.append('--for-only')
    if opts['rev_only']: cmd.append('--rev-only')
    r = subprocess.run(cmd + ['--print-seeds', d + '/r.fa', d + '/q.fa'], capture_output=True, text=True)
    if r.returncode != 0: return None
    mid = None; rep = 0; seeds = []
    tid = {n: i for i, n in enumerate(names)}
    for line in r.stderr.split('\n'):
        m = re.search(r'mid_occ = (\d+)', line)
        if m: mid = int(m.group(1))
        p = line.split('\t')
        if p[0] == 'RS': rep = int(p[1])
        elif p[0] == 'SD': seeds.append([tid[p[1]], 1 if p[3] == '-' else 0, int(p[2]), int(p[4]), int(p[5])])
    if mid is None: return None
    return mid, rep, seeds


def base_opts():
    return {'k': 15, 'w': 10, 'H': False, 'T': 0, 'f': None, 'fmax': None, 'U': None, 'e': None, 'q_occ_frac': None, 'X': False, 'for_only': False, 'rev_only': False}


def make(group):
    opts = base_opts()
    names = ['r0']; qname = 'q'
    if group == 'freq':
        # default options: mid_occ is max(10, the second largest occurrence count + 1) once the index holds 5001 - 10000 keys (w = 2: a reference of 8 - 9.5 kb), so a key is frequent only when it clearly dominates: homopolymer runs
        ref = rseq(rnd.randint(8000, 9500)); runs = []
        for n_run in range(rnd.choice([1, 1, 1, 2])):
            c = rnd.choice('ACGT'); run = c * (rnd.randint(40, 400) if n_run == 0 else rnd.randint(14, 24)); runs.append(run)
            i = rnd.randrange(0, len(ref)); ref = ref[:i] + run + ref[i:]
        targets = [ref]
        qp = []
        opts_k = rnd.choice([11, 13, 15, 15, 17])
        for _ in range(rnd.randint(2, 4)):
            if rnd.random() < 0.8: qp.append(runs[0][0] * (opts_k + rnd.randint(0, 9)) if rnd.random() < 0.85 else rnd.choice('ACGT') * rnd.randint(14, 40))
            i = rnd.randrange(0, len(ref) - 300); qp.append(mutate(ref[i:i + rnd.randint(80, 300)], 0.03))
            qp.append(rseq(rnd.randint(0, 50)))
        query = ''.join(qp)
        if rnd.random() < 0.3: query = revcomp(query)
        opts['k'] = opts_k; opts['w'] = 2
        names = ['r0']
    elif group in ('fu', 'strand'):
        ref, fams = repeat_rich_reference(); targets = [ref]
        if rnd.random() < 0.3: targets.append(rseq(rnd.randint(200, 900)))
        query = query_from(ref, fams)
        opts['k'] = rnd.choice([9, 11, 13, 15, 15, 15, 17, 19]); opts['w'] = rnd.choice([3, 5, 10, 10, 10, 15])
        if group == 'fu':
            if rnd.random() < 0.8:
                opts['f'] = rnd.choice([0.0, 0.00005, 0.0001, 0.0002, 0.001, 0.005, 0.02, 0.1, 0.5, 0.9, 1.0, 2.0, 5.0, 12.0, 40.0, 200.0])
                if rnd.random() < 0.2: opts['fmax'] = rnd.choice([20, 100, 1000])
            if rnd.random() < 0.3: opts['e'] = rnd.choice([0, 50, 100, 250, 500, 2000])
            if rnd.random() < 0.7: lo = rnd.choice([1, 2, 3, 5, 10, 20, 50]); opts['U'] = [lo, rnd.choice([lo, lo + 5, 30, 100, 500, 1000000])]
            if rnd.random() < 0.5: opts['q_occ_frac'] = rnd.choice([0.0, 0.001, 0.01, 0.05, 0.3, 1.0])
        if group == 'strand':
            opts['for_only' if rnd.random() < 0.5 else 'rev_only'] = True
        names = ['r%d' % i for i in range(len(targets))]
    elif group == 'hpc':
        base = homopolymer_rich(rnd.randint(300, 1500)); ref = rseq(rnd.randint(50, 300)) + base + rseq(rnd.randint(50, 300))
        targets = [ref]
        q = mutate(base, rnd.choice([0.01, 0.05, 0.1])); query = rseq(rnd.randint(0, 40)) + (revcomp(q) if rnd.random() < 0.3 else q)
        opts['H'] = True; opts['k'] = rnd.choice([9, 11, 15, 19, 19, 21]); opts['w'] = rnd.choice([5, 10, 10, 19])
        if opts['k'] > 28: opts['k'] = 19
        names = ['r0']
    elif group == 'dust':
        lc = low_complexity(rnd.randint(300, 900)); ref = rseq(rnd.randint(100, 400)) + lc + rseq(rnd.randint(100, 400))
        targets = [ref]; query = rseq(rnd.randint(0, 80)) + mutate(lc, 0.02) + rseq(rnd.randint(0, 80))
        opts['T'] = rnd.choice([2, 5, 10, 15, 20, 30, 40, 60]); opts['k'] = rnd.choice([9, 11, 13, 15, 15]); opts['w'] = rnd.choice([3, 5, 10, 10])
    elif group == 'ava':
        core = rseq(rnd.randint(300, 900)); other = mutate(core, 0.05) + rseq(rnd.randint(0, 100))
        query = core
        names_pool = {'before': 'a_seq', 'after': 'z_seq', 'same': 'q'}
        targets = [core, other, core + rseq(rnd.randint(1, 50)), revcomp(core)]
        names = [rnd.choice(['q', 'a_seq', 'z_seq', 'q', 'a_seq']) for _ in targets]
        names = ['%s%d' % (n, i) if n != 'q' or rnd.random() < 0.3 else n for i, n in enumerate(names)]
        if 'q' not in names: names[0] = 'q'
        if len(set(names)) != len(names): names = ['%s_%d' % (n, i) for i, n in enumerate(names)]; names[rnd.randrange(len(names))] = 'q'
        # make every name unique again
        seen = set(); fixed = []
        for i, n in enumerate(names):
            fixed.append(n if n not in seen else '%s_dup%d' % (n, i)); seen.add(n)
        names = fixed
        opts['X'] = True; opts['k'] = rnd.choice([11, 13, 15, 15]); opts['w'] = rnd.choice([5, 10, 10])
    return opts, names, targets, qname, query


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 330
    counts = {}
    with tempfile.TemporaryDirectory() as d:
        for group in ['freq', 'fu', 'hpc', 'dust', 'ava', 'strand']:
            made = 0; attempts = 0
            while made < per and attempts < per * 5:
                attempts += 1
                opts, names, targets, qname, query = make(group)
                if len(query) < 30: continue
                res = run(opts, names, targets, qname, query, d)
                if res is None: continue
                mid, rep, seeds = res
                if len(seeds) > 700: continue
                print(json.dumps({'group': group, 'opts': opts, 'names': names, 'qname': qname, 'targets': targets, 'query': query, 'mid_occ': mid, 'rep_len': rep, 'seeds': seeds}))
                made += 1
            counts[group] = made
    sys.stderr.write('rows %s\n' % counts)


main()
