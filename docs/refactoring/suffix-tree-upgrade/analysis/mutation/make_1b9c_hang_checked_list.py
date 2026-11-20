#!/usr/bin/env python3
"""Work package 1b.9c: the list for the second pass over the kills that the quiet re-run (results_1b9c_hang.tsv) classified as a crash or a hang of the host: the same mutants run with CHECKED=1 (rerun_1b9c_hang.py: the checked build
   -p:EsaBounds=true and ESA_EXPECT_BOUNDS=true, where a wrong cursor is an exception of the bounds checks and not a read outside the mapping). Writes 1b9c_h_crash_listed.json. usage: make_1b9c_hang_checked_list.py"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
kinds = {}
for line in open(os.path.join(HERE, 'results_1b9c_hang.tsv')):
    p = line.rstrip('\n').split('\t')
    if len(p) >= 3:
        kinds[p[0]] = p[2]
out = [m for m in json.load(open(os.path.join(HERE, '1b9c_h_all_listed.json'))) if kinds.get(m['package'] + ':' + m['id'], '') in ('crash', 'hang', 'run-limit', 'no-result')]
json.dump(out, open(os.path.join(HERE, '1b9c_h_crash_listed.json'), 'w'), indent=1)
print(len(out), [m['package'] + ':' + m['id'] for m in out])
