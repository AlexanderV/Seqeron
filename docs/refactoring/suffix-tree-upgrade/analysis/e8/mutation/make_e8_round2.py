#!/usr/bin/env python3
"""Round 2 of the verifier audit: the mutants that survived round 1 (and the three that did not compile, replaced by others that do) are run again against the strengthened tests.
   usage: make_e8_round2.py <repo root>   (reads e8_<file>_listed.json and results_round1_<file>.tsv, writes e8_round2_<file>_listed.json)"""
import json, os, sys
ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
REPLACE = {   # the three that gave analyzer errors: another site of the same intent
    'B02': ('the unused bits of the last data byte are not checked', 'if (usedBits != 0 && (file.Pointer[dataBytes - 1]', 'if (usedBits > 8 && (file.Pointer[dataBytes - 1]'),
    'C23': ('a pair left over in a block is not reported', 'if (lcp.HasEscape && lcp.PairCount > 0 && pair != lcp.LowerBound(end)) totals.Lcp.Add(', 'if (lcp.HasEscape && lcp.PairCount > 0 && pair < 0) totals.Lcp.Add('),
    'E01': ('a statistics file without a final line feed passes', "if (!text.EndsWith('\\n'))", 'if (text.Length == 0)'),
}
DROP = {'D14'}   # vacuous (string concatenation with the empty string): not a mutant
files = {'outputverifier': 'OutputVerifier.cs', 'metafiles': 'MetaFiles.cs', 'lcpview': 'LcpView.cs', 'intervalstack': 'IntervalStack.cs', 'e8': 'E8.cs'}
total = 0
for key in files:
    listed = json.load(open('e8_%s_listed.json' % key))
    res = {}
    for line in open('results_round1_%s.tsv' % key):
        p = line.rstrip('\n').split('\t')
        res[p[0]] = (p[2], p[3] if len(p) > 3 else '')
    out = []
    for m in listed:
        r1, r2 = res.get(m['id'], ('', ''))
        killed = r1.startswith('killed') or r2.startswith('killed')
        if m['id'] in DROP or killed:
            continue
        if m['id'] in REPLACE:
            m = dict(m, desc=REPLACE[m['id']][0], old=REPLACE[m['id']][1], new=REPLACE[m['id']][2])
        out.append(m)
    json.dump(out, open('e8_round2_%s_listed.json' % key, 'w'), indent=1)
    total += len(out)
    print(key, len(out), [m['id'] for m in out])
print(total, 'mutants for round 2')
