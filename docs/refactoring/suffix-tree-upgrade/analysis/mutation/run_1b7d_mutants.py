#!/usr/bin/env python3
"""Single-site mutants of the production code of work package 1b.7d (src/SuffixTree/Algorithms/SuffixTree.Esa/External), each built and tested in the Release configuration with a WATCHDOG (a wrong cursor can make a
   loop non-terminating or crash the test host: the whole process group is killed on expiry, a hang or a host crash counts as killed). Stage 1 builds the slim copy of SuffixTree.Esa.Tests (the new External tests of 1b.7d only) and runs them (about 15 s of tests);
   a mutant that survives it is built again with SuffixTree.EmProto.Tests and run against the quick differential tests of that project (stage 2; the data-file and 10^6-family tests are left out). The source file is
   restored after every mutant (and at the end, also on an exception). The mutants run in a separate git worktree, never in the working checkout. Analyzers are off for the mutant builds (a mutant that an analyzer
   rejects is still a valid test of the tests).
   usage: run_1b7d_mutants.py <worktree> <1b7d_<file>_listed.json> <result.tsv> [<start>:<end> slice]
   environment: MUTANT_BUILD_LIMIT (s, default 400), MUTANT_TEST_LIMIT (s, default 200), MUTANT_STAGES ("12" default; "1" = stage 1 only, "2" = stage 2 only, for the survivors of a run of stage 1)."""
import json, os, re, signal, subprocess, sys

ROOT = sys.argv[1]
BUILD_LIMIT = int(os.environ.get('MUTANT_BUILD_LIMIT', '400'))
TEST_LIMIT = int(os.environ.get('MUTANT_TEST_LIMIT', '200'))
STAGES = os.environ.get('MUTANT_STAGES', '12')
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', SEQERON_BENCH_DIR='/tmp/stbench')
mutants_file, out = sys.argv[2], sys.argv[3]
sl = sys.argv[4] if len(sys.argv) > 4 else ''
stage1_default = '(FullyQualifiedName~ExternalDcsTests|FullyQualifiedName~BucketSorterTests|FullyQualifiedName~SaBuilderTests)'  # the listed mutants carry the filter of their own file
stage2 = '(FullyQualifiedName~SmallTexts_WithRotatingTinyParameters|FullyQualifiedName~TheRankFile_IsTheE3bFile|FullyQualifiedName~TheDifferenceCovers_AreThePrototypes|FullyQualifiedName~TheCountedIo_AndTheScratchPeaks)&FullyQualifiedName~SaBuildDifferentialTests'
proj1 = ROOT + '/tests/SuffixTree/SuffixTree.Esa.Slim/SuffixTree.Esa.Slim.csproj'  # the External tests alone, see SuffixTree.Esa.Slim1b7d.csproj.template
proj2 = ROOT + '/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj'
mutants = json.load(open(mutants_file))
if sl:
    a, b = sl.split(':')
    mutants = mutants[int(a or 0):int(b) if b else None]
path = os.path.join(ROOT, mutants[0]['file'])
orig = open(path).read()
FLAGS = ['-p:RunAnalyzers=false', '-p:EnforceCodeStyleInBuild=false']


def sh(args, limit):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, start_new_session=True)
    try:
        txt, _ = p.communicate(timeout=limit)
        return txt, False
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        txt, _ = p.communicate()
        return txt, True


def build(proj):
    txt, hung = sh(['dotnet', 'build', proj, '-c', 'Release', '--nologo', '-v', 'q'] + FLAGS, BUILD_LIMIT)
    if hung:
        return 'BUILD-TIMEOUT'
    if re.search(r'error (CS|IDE|CA|S)\d+', txt) or 'Build FAILED' in txt:
        mm = re.search(r'error ((CS|IDE|CA|S)\d+)', txt)
        return 'BUILD-ERROR ' + (mm.group(1) if mm else '')
    return ''


def test(proj, filt):
    txt, hung = sh(['dotnet', 'test', proj, '-c', 'Release', '--nologo', '--no-build', '--filter', filt], TEST_LIMIT)
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
        r1, r2 = '', ''
        if '1' in STAGES:
            b = build(proj1)
            if b:
                results.append((m['id'], m['desc'], b, ''))
                open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
                continue
            r1 = test(proj1, m.get('filter', stage1_default))
        if '2' in STAGES and (r1 == 'survived' or STAGES == '2'):
            b2 = build(proj2)
            r2 = b2 if b2 else test(proj2, stage2)
        results.append((m['id'], m['desc'], r1 or '-', r2))
        open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
finally:
    open(path, 'w').write(orig)
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', m['file']], env=env)
