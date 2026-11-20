#!/usr/bin/env python3
"""Work package 1b.9c: the lists for stage 2 of the audit, the whole default set of SuffixTree.Esa.Tests, of every mutant that survived stage 1 (the tests of the list and the classes that name the file) in results_1b9c_<d|e>_<file>.tsv and
   has no line yet in results_1b9c_s2_<d|e>_<file>.tsv. Writes 1b9c_s2_<d|e>_<file>_listed.json (the same mutants as the lists they come from)."""
import glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '../../../../..'))
# Stage 1 of the audit ran at the head of the package until the first gap tests were committed and at that commit (the head plus the first batch of gap tests and the seams) after it, so a kill whose failing tests are ALL
# tests of this package may be a mutant that the tests of the head do not kill: it goes to stage 2 (at the head) like a survivor.
NEW = []
for f in ('SaGap1b9cTests.cs', 'SesaGap1b9cTests.cs'):
    NEW += re.findall(r'public void (\w+)', open(os.path.join(ROOT, 'tests/SuffixTree/SuffixTree.Esa.Tests/External', f)).read())


def only_new(verdict):
    m = re.match(r'killed\((\d+) failed: ([^)]*)\)', verdict)
    return bool(m) and int(m.group(1)) <= 3 and all(any(n.startswith(x) or x.startswith(n) for n in NEW) for x in m.group(2).split(','))


n = 0
for f in sorted(glob.glob(os.path.join(HERE, '1b9c_[de]_*_listed.json'))):
    name = os.path.basename(f)[len('1b9c_'):-len('_listed.json')]
    res = os.path.join(HERE, 'results_1b9c_%s.tsv' % name)
    if not os.path.exists(res):
        continue
    surv = set()
    for line in open(res):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3 and (p[2] == 'survived' or p[2].startswith('NO-RESULT') or only_new(p[2])):
            surv.add(p[0])
    out = [m for m in json.load(open(f)) if m['id'] in surv]
    if out:
        json.dump(out, open(os.path.join(HERE, '1b9c_s2_%s_listed.json' % name), 'w'), indent=1)
        n += len(out)
print('stage 2 lists:', n)
