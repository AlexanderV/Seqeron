#!/usr/bin/env python3
"""Single-site mutants of the new code of E9 (apps/SuffixTree.EmProto/Ph6), each built and tested in the Release configuration with a WATCHDOG (a wrong cursor can make a scan non-terminating or crash the test host: the whole
   process group is killed on expiry, a hang or a host crash counts as killed). Stage 1 runs the Ph6 tests without the slow exhaustive / family ones (about 25 s of tests); a mutant that survives it is run again against stage 2
   (all Ph6 tests, about 90 s). The source file is restored after every mutant (and at the end, also on an exception).
   usage: E9_ROOT=<checkout> run_e9_mutants.py <e9_<file>_listed.json> <result.tsv> [<start>:<end> slice]
   environment: MUTANT_BUILD_LIMIT (s, default 300), MUTANT_TEST_LIMIT (s, default 240)."""
import json, os, re, signal, subprocess, sys

ROOT = os.environ.get('E9_ROOT', '/home/user/Seqeron')
BUILD_LIMIT = int(os.environ.get('MUTANT_BUILD_LIMIT', '300'))
TEST_LIMIT = int(os.environ.get('MUTANT_TEST_LIMIT', '240'))
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', SEQERON_BENCH_DIR='/tmp/stbench')
mutants_file, out = sys.argv[1], sys.argv[2]
sl = sys.argv[3] if len(sys.argv) > 3 else ''
stage2 = 'FullyQualifiedName~EmProto.Tests.Ph6&FullyQualifiedName!~ExhaustiveSmallStrings'
stage1 = 'FullyQualifiedName~EmProto.Tests.Ph6&FullyQualifiedName!~ExhaustiveSmallStrings&FullyQualifiedName!~ExhaustiveTinyTexts&FullyQualifiedName!~OneText_EveryWindow&FullyQualifiedName!~RandomTexts_TablesAgree&FullyQualifiedName!~AdversarialFamilies&FullyQualifiedName!~LongPeriodic'
proj = ROOT + '/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj'
mutants = json.load(open(mutants_file))
if sl:
    a, b = sl.split(':')
    mutants = mutants[int(a or 0):int(b) if b else None]
path = os.path.join(ROOT, mutants[0]['file'])
orig = open(path).read()


def sh(args, limit, e):
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
        mm = re.search(r'error ((CS|IDE|CA|S)\d+)', txt)
        return 'BUILD-ERROR ' + (mm.group(1) if mm else '')
    return ''


def test(filt):
    txt, hung = sh(['dotnet', 'test', proj, '-c', 'Release', '--nologo', '--no-build', '--filter', filt], TEST_LIMIT, env)
    if hung:
        return 'killed(hang)'
    if 'Test host process crashed' in txt or 'test run was aborted' in txt:
        return 'killed(host crash)'
    if 'Failed!' in txt:
        mm = re.search(r'Failed!.*?Failed:\s+(\d+)', txt)
        return 'killed(%s)' % (mm.group(1) if mm else '?')
    if 'Passed!' in txt:
        return 'survived'
    return 'NO-RESULT'


results = []
try:
    for m in mutants:
        if os.path.join(ROOT, m['file']) != path:
            raise SystemExit('one source file per run')
        if orig.count(m['old']) != 1:
            results.append((m['id'], m['desc'], 'PATTERN-NOT-UNIQUE', ''))
            continue
        open(path, 'w').write(orig.replace(m['old'], m['new']))
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
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', m['file']], env=env)
