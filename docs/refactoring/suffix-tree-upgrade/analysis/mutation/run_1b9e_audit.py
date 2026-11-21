#!/usr/bin/env python3
"""Work package 1b.9e, the mutation audit of the code that the package adds: the lists 1b9e_*_listed.json (make_1b9e_mutants.py) are run by run_1b9c_mutants.py in N git worktrees at the same time, Release, watchdog,
analyzers off, no NUnit.StopOnError (see the docstring of run_1b9c_mutants.py).
  stage1   the lists of the format, the lock and the reopening of the container (sesaresumestate, sesabuildlock, mappedoutputfile) with the quick set SesaResumeFormatTests (a fraction of a second of tests); the lists of
           the recording, the builder and PH1-PH3 (sesaresume, sesabuilder, sabuilder) with the whole set of the package (SesaResumeFormatTests, SesaResume1b9eTests); results_1b9e_<list>.tsv
  stage1b  the survivors of the quick set run again with the whole set of the package: 1b9e_b_<list>_listed.json, results_1b9e_b_<list>.tsv
  closed   the survivors run again at the head that holds the tests written for them (a survivor is a gap of the tests unless it is equivalent): 1b9e_c_<list>_listed.json, results_1b9e_c_<list>.tsv
  stage2   the survivors of the earlier stages against the whole default set of SuffixTree.Esa.Tests (every test that is not [Explicit]): 1b9e_s2_<list>_listed.json, results_1b9e_s2_<list>.tsv
A tmpfs that a mutant left mounted is unmounted after every list. The runner appends to a result file after every mutant and skips what the file holds, so this can be stopped and started again.
usage: run_1b9e_audit.py <worktree,worktree,...> <commit> [stage1,stage1b,closed,stage2]"""
import glob, json, os, queue, subprocess, sys, threading
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = subprocess.check_output(['git', '-C', HERE, 'rev-parse', '--show-toplevel'], text=True).strip()
trees, commit = sys.argv[1].split(','), sys.argv[2]
n = len(trees)
what = (sys.argv[3] if len(sys.argv) > 3 else 'stage1,stage2').split(',')
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', SEQERON_BENCH_DIR='/tmp/stbench')
env.pop('DOTNET_GCHeapHardLimit', None)
QUICK = 'FullyQualifiedName~SesaResumeFormatTests'
NEW = 'FullyQualifiedName~SesaResume1b9eTests|FullyQualifiedName~SesaResumeFormatTests|FullyQualifiedName~SaCheckpoint1b9eTests'
QUICK_LISTS = ('sesaresumestate', 'sesabuildlock', 'mappedoutputfile')
SLICE = 13  # the largest list is cut in slices of this many mutants so that the worktrees share it
for w in trees:
    if not os.path.isdir(w):
        subprocess.run(['git', '-C', ROOT, 'worktree', 'add', '--detach', w, commit], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(['git', '-C', w, 'checkout', '--detach', '-q', commit], check=True)
    subprocess.run(['git', '-C', w, 'checkout', '--', '.'], check=True)


def unmount_leftovers(tmp):
    out = subprocess.run(['findmnt', '-rn', '-o', 'TARGET'], capture_output=True, text=True).stdout.split('\n')
    for t in out:
        if t.startswith(tmp):
            subprocess.run(['umount', '-l', t], capture_output=True)


def todo(listfile, resname):
    res = os.path.join(HERE, resname)
    done = set(l.split('\t')[0] for l in open(res)) if os.path.exists(res) else set()
    return [m for m in json.load(open(listfile)) if m['id'] not in done]


def survivors(name, first=False):
    """The survivors of the last stage that this list went through (a quick list goes through stage1 and stage1b); with first=True those of stage1."""
    res = os.path.join(HERE, ('results_1b9e_b_%s.tsv' if name in QUICK_LISTS and not first else 'results_1b9e_%s.tsv') % name)
    if not os.path.exists(res):
        return set()
    return set(l.split('\t')[0] for l in open(res) if len(l.split('\t')) > 2 and l.split('\t')[2].startswith('survived'))


def run_phase(jobs):
    q = queue.Queue()
    for job in jobs:
        q.put(job)

    def worker(k):
        tmp = '/dev/shm/1b9e-aud%d' % k
        while True:
            try:
                listfile, resname, extra, sl = q.get_nowait()
            except queue.Empty:
                return
            e = dict(env, MUTANT_TMPDIR=tmp, **extra)
            cmd = ['python3', '-I', os.path.join(HERE, 'run_1b9c_mutants.py'), trees[k - 1], listfile, os.path.join(HERE, resname)] + ([sl] if sl else [])
            subprocess.run(cmd, env=e, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            unmount_leftovers(tmp)

    ts = [threading.Thread(target=worker, args=(k,)) for k in range(1, n + 1)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()


def slices(listfile, count):
    return ['%d:%d' % (i, min(i + SLICE, count)) for i in range(0, count, SLICE)] if count > SLICE else ['']


if 'stage1' in what:
    jobs = []
    for f in sorted(glob.glob(os.path.join(HERE, '1b9e_*_listed.json'))):
        name = os.path.basename(f)[len('1b9e_'):-len('_listed.json')]
        if name.startswith(('s2_', 'b_', 'c_')):
            continue
        if not todo(f, 'results_1b9e_%s.tsv' % name):
            continue
        filt = QUICK if name in QUICK_LISTS else NEW
        for sl in slices(f, len(json.load(open(f)))):
            jobs.append((f, 'results_1b9e_%s.tsv' % name, {'MUTANT_FILTER_ONLY': filt, 'MUTANT_TEST_LIMIT': '900'}, sl))
    run_phase(jobs)
if 'stage1b' in what:
    jobs = []
    for f in sorted(glob.glob(os.path.join(HERE, '1b9e_*_listed.json'))):
        name = os.path.basename(f)[len('1b9e_'):-len('_listed.json')]
        if name.startswith(('s2_', 'b_', 'c_')) or name not in QUICK_LISTS:
            continue
        surv = survivors(name, first=True)
        sel = [m for m in json.load(open(f)) if m['id'] in surv]
        if sel:
            f2 = os.path.join(HERE, '1b9e_b_%s_listed.json' % name)
            json.dump(sel, open(f2, 'w'), indent=1)
            if todo(f2, 'results_1b9e_b_%s.tsv' % name):
                jobs.append((f2, 'results_1b9e_b_%s.tsv' % name, {'MUTANT_FILTER_ONLY': NEW, 'MUTANT_TEST_LIMIT': '900'}, ''))
    run_phase(jobs)
if 'closed' in what:
    jobs = []
    for f in sorted(glob.glob(os.path.join(HERE, '1b9e_*_listed.json'))):
        name = os.path.basename(f)[len('1b9e_'):-len('_listed.json')]
        if name.startswith(('s2_', 'b_', 'c_')):
            continue
        surv = survivors(name)
        sel = [m for m in json.load(open(f)) if m['id'] in surv]
        if sel:
            f2 = os.path.join(HERE, '1b9e_c_%s_listed.json' % name)
            json.dump(sel, open(f2, 'w'), indent=1)
            if todo(f2, 'results_1b9e_c_%s.tsv' % name):
                jobs.append((f2, 'results_1b9e_c_%s.tsv' % name, {'MUTANT_FILTER_ONLY': NEW, 'MUTANT_TEST_LIMIT': '900'}, ''))
    run_phase(jobs)
if 'stage2' in what:
    jobs = []
    for f in sorted(glob.glob(os.path.join(HERE, '1b9e_*_listed.json'))):
        name = os.path.basename(f)[len('1b9e_'):-len('_listed.json')]
        if name.startswith(('s2_', 'b_', 'c_')):
            continue
        surv = survivors(name)
        resc = os.path.join(HERE, 'results_1b9e_c_%s.tsv' % name)
        if os.path.exists(resc):
            surv &= set(l.split('\t')[0] for l in open(resc) if len(l.split('\t')) > 2 and l.split('\t')[2].startswith('survived'))
        sel = [m for m in json.load(open(f)) if m['id'] in surv]
        if sel:
            f2 = os.path.join(HERE, '1b9e_s2_%s_listed.json' % name)
            json.dump(sel, open(f2, 'w'), indent=1)
            if todo(f2, 'results_1b9e_s2_%s.tsv' % name):
                jobs.append((f2, 'results_1b9e_s2_%s.tsv' % name, {'MUTANT_STAGES': '2', 'MUTANT_TEST_LIMIT': '1500'}, ''))
    run_phase(jobs)
open(os.path.join(HERE, '.audit_1b9e.done'), 'w').write('done\n')
