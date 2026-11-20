#!/usr/bin/env python3
"""Final status of every listed mutant of 1b.7d (561) and 1b.7e (600) after work package 1b.9c, over the earlier audits (summary_1b7d.tsv, summary_1b7e.tsv) and the runs of this package:
     results_1b9c_<d|e>_<file>.tsv      stage 1: the 409 mutants that the earlier audits did not run (make_1b9c_lists.py): the filter of the list; for the first part of the run also the classes that name the file
     results_1b9c_s2_<d|e>_<file>.tsv   stage 2: the survivors of stage 1 (and the kills that only a test of this package made), at the head of the audit: the classes that name the file, then the whole default set of SuffixTree.Esa.Tests
     results_1b9c_hang.tsv              the one-by-one re-run of every kill that only a watchdog or a crash of the host counted (17 of the earlier audits, those of this package)
     results_1b9c_c_<d|e>_<file>.tsv    the mutants that a new test closes, run against the final code with the new test classes only: the kill must reproduce
   Statuses, in the order of the rules:
     killed (earlier audit)       an assertion failed in a run of 1b.7d / 1b.7e
     killed (stage 1 / stage 2)   an assertion failed in a run of this package (stage 2: a survivor of the listed filter that the tests of the head kill)
     killed, assertion (re-run)   a kill that a watchdog or a crash had counted, and that the re-run alone shows to be a failing assertion
     killed, hang / crash         a kill that the re-run alone confirms as a hang or a crash of the host (the mutant does not terminate, or ends the process)
     closed by a test             survived the whole default set at the head; a new test closes it and the re-run against the final tests was killed by a test of this package
     equivalent                   argued in 1b7d_equivalent.py / 1b7e_equivalent.py (EQUIVALENT earlier, EQUIVALENT_1B9C in this package)
     OPEN                         anything else (the gate is zero)
   Writes summary_1b9c.tsv (package, id, status, description) and summary_1b9c_files.tsv (the arithmetic per file); prints both. usage: summarize_1b9c.py"""
import collections, glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, '1b9c_closed.py')).read())  # CLOSED
ns = {'d': {}, 'e': {}}
exec(open(os.path.join(HERE, '1b7d_equivalent.py')).read(), ns['d'])
exec(open(os.path.join(HERE, '1b7e_equivalent.py')).read(), ns['e'])
EQ = {k: dict(ns[k]['EQUIVALENT'], **ns[k]['EQUIVALENT_1B9C']) for k in 'de'}

earlier = {}
for l in open(os.path.join(HERE, 'summary_1b7d.tsv')):
    b, i, st, d = l.rstrip('\n').split('\t')
    earlier[('d', i)] = st
for l in open(os.path.join(HERE, 'summary_1b7e.tsv')):
    i, st, d = l.rstrip('\n').split('\t')
    earlier[('e', i)] = st
listed = collections.OrderedDict()
for f in sorted(glob.glob(os.path.join(HERE, '1b7d_*_listed.json'))):
    b = os.path.basename(f)
    if not re.match(r'1b7d_r\d', b):
        for m in json.load(open(f)):
            listed.setdefault('1b7d/' + b[len('1b7d_'):-len('_listed.json')], []).append(('d', m['id'], m['desc']))
tmp = collections.OrderedDict()
for f in sorted(glob.glob(os.path.join(HERE, '1b7e_r[12]_*_listed.json'))):
    for m in json.load(open(f)):
        tmp.setdefault('1b7e/' + os.path.basename(m['file']), []).append(('e', m['id'], m['desc']))
listed.update(tmp)


def read(pattern, rx):
    out = collections.defaultdict(list)
    for f in sorted(glob.glob(os.path.join(HERE, pattern))):
        m = re.match(rx, os.path.basename(f))
        if m:
            for line in open(f):
                p = line.rstrip('\n').split('\t')
                if len(p) >= 3:
                    out[(m.group(1), p[0])].append(p[2:])
    return out


s1 = read('results_1b9c_[de]_*.tsv', r'results_1b9c_([de])_.*\.tsv')
s2 = read('results_1b9c_s2_[de]_*.tsv', r'results_1b9c_s2_([de])_.*\.tsv')
cl = read('results_1b9c_c_[de]_*.tsv', r'results_1b9c_c_([de])_.*\.tsv')
hang = {}
if os.path.exists(os.path.join(HERE, 'results_1b9c_hang.tsv')):
    for line in open(os.path.join(HERE, 'results_1b9c_hang.tsv')):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3:
            k, i = p[0].split(':')
            hang[(k, i)] = p[2]


# tests that fail under load (several mutants at once) and pass alone: DistTests.ScratchIsReleasedChunkByChunk_... polls the files of the scratch directory while a pass runs, and the reservation test of the system source reads
# the blocks that the file system booked; a kill whose failing tests are all of these is not a kill (the mutant is run again alone)
FLAKY = ('ScratchIsReleasedChunkByChunk', 'TheSystemSource_OnLinux_ReservesTheBlocks')


def failed(rs):
    for r in rs:
        m = re.search(r'killed\((\d+) failed: ([^)]*)', ' '.join(r))
        if m:
            names = [x for x in m.group(2).split(',') if x]
            if int(m.group(1)) > len(names) or any(not x.startswith(FLAKY) for x in names):
                return True
    return False


def dogged(rs):
    return any('killed(hang' in ' '.join(r) or 'killed(host crash' in ' '.join(r) for r in rs)


# the tests of this package: a kill by them alone at stage 1 does not count (stage 1 ran at the head only for the first part of the audit)
NEW = []
for f in sorted(glob.glob(os.path.join(HERE, '../../../../..', 'tests/SuffixTree/SuffixTree.Esa.Tests/External', '*Gap1b9cTests.cs'))):
    NEW += re.findall(r'public void (\w+)', open(f).read())


def only_new(rs):
    for r in rs:
        m = re.match(r'killed\((\d+) failed: ([^)]*)\)', r[0])
        if m and int(m.group(1)) <= 3 and all(any(n.startswith(x) or x.startswith(n) for n in NEW) for x in m.group(2).split(',')):
            return True
    return False


STATUSES = ['killed (earlier audit)', 'killed (stage 1)', 'killed (stage 2)', 'killed, assertion (re-run)', 'killed, hang / crash', 'closed by a test', 'equivalent', 'OPEN']
out = open(os.path.join(HERE, 'summary_1b9c.tsv'), 'w')
files = open(os.path.join(HERE, 'summary_1b9c_files.tsv'), 'w')
files.write('file\tlisted\trun_earlier\tunrun\t' + '\t'.join(STATUSES) + '\n')
total = collections.Counter()
grand = collections.Counter()
open_ids = []
for fname, ms in listed.items():
    c = collections.Counter()
    for k, i, d in ms:
        key = (k, i)
        st0 = earlier.get(key, 'not run')
        h = hang.get(key)
        status = None
        if h is not None:  # a kill that a watchdog or a crash had counted (or that this package's stage 1 counted that way): the verdict of the re-run alone
            if h == 'assertion':
                status = 'killed, assertion (re-run)'
            elif h in ('crash', 'hang', 'run-limit'):
                status = 'killed, hang / crash'
            # "none": the mutant survived the quiet re-run, so it is a survivor that a test or an argument must close (below)
        elif st0 == 'killed':
            status = 'killed (earlier audit)'
        elif st0 == 'not run':
            r1, r2 = s1.get(key, []), s2.get(key, [])
            if failed(r1) and not only_new(r1):
                status = 'killed (stage 1)'
            elif failed(r2):
                status = 'killed (stage 2)'
        if status is None:
            if key in CLOSED and failed(cl.get(key, [])):
                status = 'closed by a test'
            elif i in EQ[k] and st0 != 'killed':
                status = 'equivalent'
            elif st0 == 'equivalent':
                status = 'equivalent'
        if status is None:
            status = 'OPEN'
        c[status] += 1
        out.write('%s\t%s\t%s\t%s\n' % (k, i, status, d))
        if status == 'OPEN':
            open_ids.append((k, i))
    n_unrun = sum(1 for k, i, d in ms if earlier.get((k, i), 'not run') == 'not run')
    files.write('%s\t%d\t%d\t%d\t%s\n' % (fname, len(ms), len(ms) - n_unrun, n_unrun, '\t'.join(str(c[s]) for s in STATUSES)))
    grand.update(c)
    total['listed'] += len(ms)
    total['unrun'] += n_unrun
out.close()
files.write('TOTAL\t%d\t%d\t%d\t%s\n' % (total['listed'], total['listed'] - total['unrun'], total['unrun'], '\t'.join(str(grand[s]) for s in STATUSES)))
files.close()
print(open(os.path.join(HERE, 'summary_1b9c_files.tsv')).read())
print('open:', len(open_ids), open_ids[:12])
