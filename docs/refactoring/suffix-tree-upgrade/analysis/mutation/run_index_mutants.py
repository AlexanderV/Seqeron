#!/usr/bin/env python3
"""Single-site mutants of the index-contract implementation (work package 1b.5): each mutant is built and tested in the Release configuration and, if it survives, again in the checked configuration
   (-p:EsaBounds=true), exactly as run_kernel_mutants.py does, with a WATCHDOG: a wrong index can make a loop non-terminating, so every build and every test run has a time limit, the whole process group is killed
   on expiry (the test host included) and the mutant counts as killed(hang). The source file is restored after every mutant.
   usage: ESA_ROOT=<git worktree> run_index_mutants.py <file relative to the repo root> <mutants.json> <result.tsv> [<start>:<end> slice of the list] [nunit filter, default FullyQualifiedName~Esa.Tests.Contract]
   environment: MUTANT_BUILD_LIMIT (s, default 300), MUTANT_TEST_LIMIT (s, default 150)."""
import json, os, re, signal, subprocess, sys

ROOT = os.environ.get('ESA_ROOT', '/home/user/Seqeron')
BUILD_LIMIT = int(os.environ.get('MUTANT_BUILD_LIMIT', '300'))
TEST_LIMIT = int(os.environ.get('MUTANT_TEST_LIMIT', '150'))
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', SEQERON_BENCH_DIR='/tmp/stbench')
kernel, mutants_file, out = sys.argv[1], sys.argv[2], sys.argv[3]
sl = sys.argv[4] if len(sys.argv) > 4 else ''
filt = sys.argv[5] if len(sys.argv) > 5 else 'FullyQualifiedName~Esa.Tests.Contract'
path = os.path.join(ROOT, kernel)
proj = ROOT + '/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj'
orig = open(path).read()
mutants = json.load(open(mutants_file))
if sl:
    a, b = sl.split(':')
    mutants = mutants[int(a or 0):int(b) if b else None]


def sh(args, limit, e):
    """Runs a command in its own process group; returns (text, timed_out)."""
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=e, start_new_session=True)
    try:
        txt, _ = p.communicate(timeout=limit)
        return txt, False
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        txt, _ = p.communicate()
        return txt, True


def run(bounds):
    e = dict(env)
    if bounds:
        e['ESA_EXPECT_BOUNDS'] = 'true'
    flag = '-p:EsaBounds=' + ('true' if bounds else 'false')
    txt, hung = sh(['dotnet', 'build', proj, '-c', 'Release', '--nologo', '-v', 'q', flag], BUILD_LIMIT, e)
    if hung:
        return 'BUILD-TIMEOUT'
    if re.search(r'error (CS|IDE|CA|S)\d+', txt) or 'Build FAILED' in txt:
        return 'BUILD-ERROR'
    txt, hung = sh(['dotnet', 'test', proj, '-c', 'Release', '--nologo', '--no-build', flag, '--filter', filt], TEST_LIMIT, e)
    if hung:
        return 'killed(hang)'
    if 'Test host process crashed' in txt or 'test run was aborted' in txt:
        return 'killed(host crash)'
    if 'Failed!' in txt:
        m = re.search(r'Failed!.*?Failed:\s+(\d+)', txt)
        return 'killed(%s)' % (m.group(1) if m else '?')
    if 'Passed!' in txt:
        return 'survived'
    return 'NO-RESULT'


results = []
try:
    for m in mutants:
        src = orig
        pairs = [(m['old'], m['new'])] + [tuple(x) for x in m.get('also', [])]
        if any(src.count(o) != 1 for o, _ in pairs):
            results.append((m['id'], m['desc'], 'PATTERN-NOT-UNIQUE', ''))
            continue
        for o, nw in pairs:
            src = src.replace(o, nw)
        open(path, 'w').write(src)
        r1 = run(False)
        r2 = run(True) if r1 == 'survived' else ''
        results.append((m['id'], m['desc'], r1, r2))
        open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
finally:
    open(path, 'w').write(orig)
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', kernel], env=env)
