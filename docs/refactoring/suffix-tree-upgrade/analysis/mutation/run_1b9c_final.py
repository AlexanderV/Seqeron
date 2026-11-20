#!/usr/bin/env python3
"""Work package 1b.9c, the last runs of the audit, at the final head, in N git worktrees at the same time (one source file per run, so each list is dealt whole to a worktree; the largest first):
     closed   1b9c_c_*_listed.json   the mutants that a new test closes, run with the new test classes only ("the kill must reproduce"); results_1b9c_c_*.tsv
     fresh 1  1b9c_f_*_listed.json   the fresh mutants, the new test classes only; results_1b9c_f_*.tsv
     eq       1b9c_q_*_listed.json   the argued equivalents, the whole default set (an equivalent mutant must survive it); results_1b9c_q_*.tsv
     fresh 2  1b9c_f2_*_listed.json  the survivors of fresh 1, the whole default set; results_1b9c_f2_*.tsv
   The runner (run_1b9c_mutants.py) appends to a result file after every mutant and skips what the file holds, so this can be started again.      s2       1b9c_s2_*_listed.json  the survivors of stage 1 that no test of this package closes and no argument covers, against the whole default set at the head of the audit; results_1b9c_s2_*.tsv
   usage: run_1b9c_final.py <worktree,worktree,...> <commit> [s2,closed,fresh,eq]"""
import glob, json, os, queue, subprocess, sys, threading
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = subprocess.check_output(['git', '-C', HERE, 'rev-parse', '--show-toplevel'], text=True).strip()
trees, commit = sys.argv[1].split(','), sys.argv[2]
n = len(trees)
what = (sys.argv[3] if len(sys.argv) > 3 else 'closed,fresh,eq').split(',')
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', SEQERON_BENCH_DIR='/tmp/stbench')
env.pop('DOTNET_GCHeapHardLimit', None)
NEW = 'FullyQualifiedName~Gap1b9cTests'
for w in trees:
    if not os.path.isdir(w):
        subprocess.run(['git', '-C', ROOT, 'worktree', 'add', '--detach', w, commit], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(['git', '-C', w, 'checkout', '--detach', '-q', commit], check=True)
    subprocess.run(['git', '-C', w, 'checkout', '--', '.'], check=True)


def todo(listfile, resname):
    res = os.path.join(HERE, resname)
    done = set(l.split('\t')[0] for l in open(res)) if os.path.exists(res) else set()
    return [m for m in json.load(open(listfile)) if m['id'] not in done]


def run_phase(jobs):
    q = queue.Queue()
    for size, job in sorted(jobs, key=lambda x: -x[0]):
        q.put(job)

    def worker(k):
        while True:
            try:
                listfile, resname, extra = q.get_nowait()
            except queue.Empty:
                return
            e = dict(env, MUTANT_TMPDIR='/dev/shm/1b9c-fin%d' % k, **extra)
            subprocess.run(['python3', '-I', os.path.join(HERE, 'run_1b9c_mutants.py'), trees[k - 1], listfile, os.path.join(HERE, resname)], env=e, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    ts = [threading.Thread(target=worker, args=(k,)) for k in range(1, n + 1)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()


jobs = []
if 's2' in what:
    for f in sorted(glob.glob(os.path.join(HERE, '1b9c_s2_*_listed.json'))):
        name = os.path.basename(f)[:-len('_listed.json')]
        t = todo(f, 'results_%s.tsv' % name)
        if t:
            jobs.append((3 * len(t), (f, 'results_%s.tsv' % name, {'MUTANT_STAGES': '2', 'MUTANT_TEST_LIMIT': '1500'})))
if 'closed' in what:
    for f in sorted(glob.glob(os.path.join(HERE, '1b9c_c_*_listed.json'))):
        name = os.path.basename(f)[:-len('_listed.json')]
        t = todo(f, 'results_%s.tsv' % name)
        if t:
            jobs.append((len(t), (f, 'results_%s.tsv' % name, {'MUTANT_FILTER_ONLY': NEW, 'MUTANT_TEST_LIMIT': '600'})))
if 'fresh' in what:
    for f in sorted(glob.glob(os.path.join(HERE, '1b9c_f_*_listed.json'))):
        name = os.path.basename(f)[:-len('_listed.json')]
        t = todo(f, 'results_%s.tsv' % name)
        if t:
            jobs.append((len(t), (f, 'results_%s.tsv' % name, {'MUTANT_FILTER_ONLY': NEW, 'MUTANT_TEST_LIMIT': '600'})))
if 'eq' in what:
    for f in sorted(glob.glob(os.path.join(HERE, '1b9c_q_*_listed.json'))):
        name = os.path.basename(f)[:-len('_listed.json')]
        t = todo(f, 'results_%s.tsv' % name)
        if t:
            jobs.append((3 * len(t), (f, 'results_%s.tsv' % name, {'MUTANT_STAGES': '2', 'MUTANT_TEST_LIMIT': '1500'})))
run_phase(jobs)
if 'fresh' in what:  # pass 2 of the fresh mutants: the survivors of pass 1 against the whole default set
    jobs = []
    for f in sorted(glob.glob(os.path.join(HERE, '1b9c_f_*_listed.json'))):
        name = os.path.basename(f)[:-len('_listed.json')]
        res = os.path.join(HERE, 'results_%s.tsv' % name)
        surv = set(l.split('\t')[0] for l in open(res) if len(l.split('\t')) > 2 and l.split('\t')[2] == 'survived') if os.path.exists(res) else set()
        sel = [m for m in json.load(open(f)) if m['id'] in surv]
        if sel:
            f2 = os.path.join(HERE, '1b9c_f2_%s_listed.json' % name[len('1b9c_f_'):])
            json.dump(sel, open(f2, 'w'), indent=1)
            t = todo(f2, 'results_1b9c_f2_%s.tsv' % name[len('1b9c_f_'):])
            if t:
                jobs.append((3 * len(t), (f2, 'results_1b9c_f2_%s.tsv' % name[len('1b9c_f_'):], {'MUTANT_STAGES': '2', 'MUTANT_TEST_LIMIT': '1500'})))
    run_phase(jobs)
open(os.path.join(HERE, '.final_1b9c.done'), 'w').write('done\n')
