#!/usr/bin/env python3
"""The table of gate G3 "PH5 time must not exceed the PH4 time on the same input": medians of the PH4 reference runs (e5_ph4ref.tsv, the E4 base mode, scale-down) and of the PH5 runs (e5_adv.tsv, e5_shm.tsv, e5_par.tsv).
usage: python3 -I ratio_e5.py  (in this directory)"""
import statistics, collections

def load(path, metric):
    d = collections.defaultdict(list)
    for line in open(path):
        c = line.rstrip('\n').split('\t')
        if len(c) >= 7 and c[3].isdigit() and c[4] == metric:
            try:
                d[c[2]].append(float(c[5]))
            except ValueError:
                pass
    return {k: statistics.median(v) for k, v in d.items()}

ph4 = load('e5_ph4ref.tsv', 'ph4_s')
adv = load('e5_adv.tsv', 'ph5_s')
advfwd = load('e5_adv.tsv', 'forward_s')
shm = load('e5_shm.tsv', 'ph5_s')
par = load('e5_par.tsv', 'ph5_s')
names = ['a10', 'a50', 'nrun30', 'fib10', 'fib30', 'rep10', 'dna10', 'dna_100m', 'genome_100m']
print('| input | PH4 s | A | A, spill on tmpfs | F, NSV after | F, NSV concurrent | A + PSV + NSV (both), NSV after | both, concurrent |')
print('|---|---|---|---|---|---|---|---|')
def cell(v, base):
    if v is None:
        return '-'
    r = v / base
    return ('**%.2f (%.2f)**' if r > 1 else '%.2f (%.2f)') % (v, r)
for n in names:
    b = ph4.get(n + '.ph4')
    if b is None:
        continue
    print('| %s | %.2f | %s | %s | %s | %s | %s | %s |' % (n, b, cell(adv.get(n + '.a'), b), cell(shm.get(n + '.a.shm'), b), cell(adv.get(n + '.f'), b), cell(par.get(n + '.f.par'), b), cell(adv.get(n + '.both'), b), cell(par.get(n + '.both.par'), b)))
