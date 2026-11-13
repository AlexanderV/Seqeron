#!/usr/bin/env python3
"""Single-site mutants of the new code of the driver of the end-to-end pipeline (prototype E6: Pipeline/EndToEnd.cs, Prepare.cs, the scale rules, Experiments/E6.cs): each mutant is built and tested in the Release configuration with a WATCHDOG: a wrong scan can make a
   loop non-terminating (a window that never advances, a stack that never retreats), so every build and every test run has a time limit, the whole process group (the test host included) is killed on
   expiry, and the mutant counts as killed(hang). Stage 1 runs the WholePipeline tests without the slow ones (the exhaustive enumeration, the families at 10^5 and 10^6, the seeded fuzz sample, the thread and K identity); a mutant that survives it is run again against stage 2 (all of them but the exhaustive enumeration). The source file is restored after every mutant.
   usage: E6_ROOT=<git worktree> run_e6_mutants.py <file relative to the repo root> <mutants.json> <result.tsv> [<start>:<end> slice] [stage-2 filter]
   environment: MUTANT_BUILD_LIMIT (s, default 300), MUTANT_TEST_LIMIT (s, default 240)."""
import json, os, re, signal, subprocess, sys

ROOT = os.environ.get('E6_ROOT', '/home/user/Seqeron')
BUILD_LIMIT = int(os.environ.get('MUTANT_BUILD_LIMIT', '300'))
TEST_LIMIT = int(os.environ.get('MUTANT_TEST_LIMIT', '240'))
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', SEQERON_BENCH_DIR='/tmp/stbench')
kernel, mutants_file, out = sys.argv[1], sys.argv[2], sys.argv[3]
sl = sys.argv[4] if len(sys.argv) > 4 else ''
stage2 = sys.argv[5] if len(sys.argv) > 5 else 'FullyQualifiedName~EmProto.Tests.WholePipeline&FullyQualifiedName!~Exhaustive_'
stage1 = 'FullyQualifiedName~EmProto.Tests.WholePipeline&FullyQualifiedName!~Exhaustive_&FullyQualifiedName!~TheDeepFamilies&FullyQualifiedName!~At100000&FullyQualifiedName!~ASeeded&FullyQualifiedName!~OutputIsIdentical'
path = os.path.join(ROOT, kernel)
proj = ROOT + '/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj'
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


def build():
    txt, hung = sh(['dotnet', 'build', proj, '-c', 'Release', '--nologo', '-v', 'q'], BUILD_LIMIT, env)
    if hung:
        return 'BUILD-TIMEOUT'
    if re.search(r'error (CS|IDE|CA|S)\d+', txt) or 'Build FAILED' in txt:
        return 'BUILD-ERROR ' + (re.search(r'error ((CS|IDE|CA|S)\d+)', txt).group(1) if re.search(r'error ((CS|IDE|CA|S)\d+)', txt) else '')
    return ''


def test(filt):
    txt, hung = sh(['dotnet', 'test', proj, '-c', 'Release', '--nologo', '--no-build', '--filter', filt], TEST_LIMIT, env)
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
        b = build()
        if b:
            results.append((m['id'], m['desc'], b, ''))
        else:
            r1 = test(stage1)
            r2 = test(stage2) if r1 == 'survived' else ''
            results.append((m['id'], m['desc'], r1, r2))
        open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
finally:
    open(path, 'w').write(orig)
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', kernel], env=env)
