#!/usr/bin/env python3
"""Single-site mutants of the production engine against the BUDGETED HARNESS of work package 1b.8 (the explicit class FullDifferentialFuzzTests of SuffixTree.Esa.Tests, a short budget per mutant). Every build and every
   test run has a WATCHDOG (a mutant can make a walk endless or crash the test host: the whole process group is killed on expiry, a hang or a host crash counts as killed); the runs do NOT use NUnit.StopOnError; analyzers are off
   for the mutant builds; the memory of a test host is capped (DOTNET_GCHeapHardLimit) because a mutant can make a list grow without end. The source file is restored after every mutant (and at the end, also on an exception).
   The mutants run in a separate git worktree, never in the working checkout. A result line: id, description, outcome, and for a kill the first failure of the harness (the operation: engine, step, member).
   usage: run_1b8_mutants.py <worktree> <1b8_<file>_listed.json> <result.tsv> [<start>:<end> slice]
   environment: MUTANT_BUILD_LIMIT (s, default 400), MUTANT_TEST_LIMIT (s, default 360), FUZZ_BUDGET_S (default 60), FUZZ_THREADS (default 2), FUZZ_HEAP_LIMIT (default 0xC0000000),
   FUZZ_FILTER (a test filter instead of the harness: used to show that the existing unit tests kill a mutant that the harness does not)."""
import json, os, re, signal, subprocess, sys

ROOT = sys.argv[1]
BUILD_LIMIT = int(os.environ.get('MUTANT_BUILD_LIMIT', '400'))
TEST_LIMIT = int(os.environ.get('MUTANT_TEST_LIMIT', '360'))
BUDGET = os.environ.get('FUZZ_BUDGET_S', '60')
FILTER = os.environ.get('FUZZ_FILTER', 'FullyQualifiedName~FullDifferentialFuzzTests')
THREADS = os.environ.get('FUZZ_THREADS', '2')
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', SEQERON_BENCH_DIR='/tmp/stbench')
mutants_file, out = sys.argv[2], sys.argv[3]
sl = sys.argv[4] if len(sys.argv) > 4 else ''
PROJ = ROOT + '/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj'
mutants = json.load(open(mutants_file))
if sl:
    a, b = sl.split(':')
    mutants = mutants[int(a or 0):int(b) if b else None]
path = os.path.join(ROOT, mutants[0]['file'])
orig = open(path).read()
FLAGS = ['-p:RunAnalyzers=false', '-p:EnforceCodeStyleInBuild=false', '-p:TreatWarningsAsErrors=false']  # a mutant that leaves a field unassigned or a local unused is still a valid mutant


def sh(args, limit, extra_env=None):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=dict(env, **(extra_env or {})), start_new_session=True)
    try:
        txt, _ = p.communicate(timeout=limit)
        return txt, False
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        txt, _ = p.communicate()
        return txt, True


def build():
    txt, hung = sh(['dotnet', 'build', PROJ, '-c', 'Release', '--nologo', '-v', 'q'] + FLAGS, BUILD_LIMIT)
    if hung:
        return 'BUILD-TIMEOUT'
    if re.search(r'error (CS|IDE|CA|S)\d+', txt) or 'Build FAILED' in txt:
        mm = re.search(r'error ((CS|IDE|CA|S)\d+)', txt)
        return 'BUILD-ERROR ' + (mm.group(1) if mm else '')
    return ''


def test(mid):
    rep = '%s.%s.rep' % (out, mid)
    for f in (rep, rep + '.log'):
        if os.path.exists(f):
            os.remove(f)
    xenv = {'ESA_FUZZ_BUDGET_S': BUDGET, 'ESA_FUZZ_THREADS': THREADS, 'ESA_FUZZ_MAXFAIL': '1', 'ESA_FUZZ_MINIMIZE_S': '4', 'ESA_FUZZ_HANG_S': '90',
            'DOTNET_GCHeapHardLimit': os.environ.get('FUZZ_HEAP_LIMIT', '0xC0000000'), 'ESA_FUZZ_REPORT': rep}
    txt, hung = sh(['dotnet', 'test', PROJ, '-c', 'Release', '--nologo', '--no-build', '--filter', FILTER], TEST_LIMIT, xenv)
    first = ''
    if os.path.exists(rep + '.log'):
        for line in open(rep + '.log'):
            if line.startswith('FAIL '):
                first = re.sub(r'\s+', ' ', line.strip())[:230]
                break
    if hung:
        return 'killed(hang)'
    if 'HANG:' in txt:
        return 'killed(watchdog of the harness)'
    if 'Test host process crashed' in txt or 'test run was aborted' in txt:
        mm = re.search(r'(System\.\w+Exception|Fatal error[^\n]*)', txt)
        return 'killed(host crash: %s)' % (mm.group(1) if mm else '?')
    if 'Failed!' in txt:
        if first:
            return 'killed(%s)' % first
        failed = re.findall(r'^\s+Failed (\S+)', txt, re.M)
        if failed and 'FullDifferentialFuzzTests' not in FILTER:
            return 'killed(by the unit tests: %d failed, first %s)' % (len(failed), failed[0])
        mm = re.search(r'Error Message:\s*\n\s*(.{0,200})', txt)
        return 'killed(test failed: %s)' % (mm.group(1).strip() if mm else '?')
    if 'Passed!' in txt:
        return 'survived'
    return 'NO-RESULT'


results = []
try:
    for m in mutants:
        if os.path.join(ROOT, m['file']) != path:
            raise SystemExit('one source file per run')
        if orig.count(m['old']) != 1:
            results.append((m['id'], m['desc'], 'PATTERN-NOT-UNIQUE'))
            continue
        open(path, 'w').write(orig.replace(m['old'], m['new']))
        b = build()
        results.append((m['id'], m['desc'], b if b else test(m['id'])))
        open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
finally:
    open(path, 'w').write(orig)
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', mutants[0]['file']], env=env)
