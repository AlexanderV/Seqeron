#!/usr/bin/env python3
"""Splits the 561 listed mutants of work package 1b.7d (1b7d_<file>_listed.json, written by make_1b7d_mutants.py) into two halves per file: round 1a = the mutants at even positions of the list, round 1b = those at odd positions (the lists
   follow the structure of each source file, so every region of the code is in both halves). The whole list does not fit the budget of the package (one mutant is a build and a test run of about a minute, a survivor of the first stage costs
   a further build of EmProto.Tests), so the halves are run one after the other and the report says which were run.
   usage: split_1b7d_halves.py   (in this directory; writes 1b7d_r1a_<file>_listed.json and 1b7d_r1b_<file>_listed.json)"""
import glob, json, os

here = os.path.dirname(os.path.abspath(__file__))
total = 0
for f in sorted(glob.glob(os.path.join(here, '1b7d_*_listed.json'))):
    base = os.path.basename(f)
    if base.startswith('1b7d_r'):
        continue
    rest = base[len('1b7d_'):]
    data = json.load(open(f))
    for tag, part in (('r1a', data[0::2]), ('r1b', data[1::2])):
        json.dump(part, open(os.path.join(here, '1b7d_%s_%s' % (tag, rest)), 'w'), indent=1)
    total += len(data)
print('split', total, 'mutants')
