#!/usr/bin/env python3
"""Round 2 / round 3 lists of the mutation audit of 1b.7a: takes the result files of the previous round, drops the mutants that were killed, drops the mutants that are argued equivalent (EQUIVALENT, with the
   argument; they are reported, not re-run) and writes 1b7a_<round>_<file>_listed.json for the rest (survivors and mutants that did not compile), to be re-run after the tests were strengthened.
   usage: make_1b7a_round2.py <round name, e.g. r2> <previous prefix, '' for round 1 or r2> [extra json with replacement mutants]"""
import glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
rnd = sys.argv[1]
prev = sys.argv[2] if len(sys.argv) > 2 else ''
EQUIVALENT = {}  # id -> argument (filled in by the audit, see the ADR subsection)
exec(open(os.path.join(HERE, '1b7a_equivalent.py')).read()) if os.path.exists(os.path.join(HERE, '1b7a_equivalent.py')) else None
pre = '1b7a_%s_' % prev if prev else '1b7a_'
total = 0
for lst in sorted(glob.glob(os.path.join(HERE, pre + '*_listed.json'))):
    base = os.path.basename(lst)[len(pre):-len('_listed.json')]
    if not prev and base.startswith(('r2_', 'r3_')):
        continue
    res = os.path.join(HERE, 'results_1b7a_%s%s.tsv' % (prev + '_' if prev else '', base))
    if not os.path.exists(res):
        continue
    status = {}
    for line in open(res):
        p = line.rstrip('\n').split('\t')
        status[p[0]] = p[2:]
    keep = []
    for m in json.load(open(lst)):
        st = status.get(m['id'])
        if st is None:
            keep.append(m)
            continue
        first = st[0]
        second = st[1] if len(st) > 1 else ''
        killed = first.startswith('killed') or second.startswith('killed')
        if killed or m['id'] in EQUIVALENT:
            continue
        keep.append(m)
    if keep:
        out = os.path.join(HERE, '1b7a_%s_%s_listed.json' % (rnd, base))
        json.dump(keep, open(out, 'w'), indent=1)
        total += len(keep)
        print('%-32s %3d to re-run' % (base, len(keep)))
# replacement mutants for the ones that did not compile (round 2 only)
repl_path = os.path.join(HERE, '1b7a_replacements.json')
if rnd == 'r2' and os.path.exists(repl_path):
    repl = json.load(open(repl_path))
    for base, muts in repl.items():
        out = os.path.join(HERE, '1b7a_r2_%s_listed.json' % base)
        cur = json.load(open(out)) if os.path.exists(out) else []
        cur = [m for m in cur if m['id'] not in {r['id'] for r in muts}] + muts
        json.dump(cur, open(out, 'w'), indent=1)
        print('%-32s +%d replacement(s)' % (base, len(muts)))
        total += len(muts)
print('total', total)
