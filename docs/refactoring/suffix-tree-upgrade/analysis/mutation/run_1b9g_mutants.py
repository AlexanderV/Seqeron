#!/usr/bin/env python3
"""Single-site mutants of the code of work package 1b.9g (prefix 1b9g_), each built and run in the Release configuration with a WATCHDOG (a mutant that breaks the protocol may hang or crash the test host: the whole process
   group is killed on expiry, a hang or a host crash counts as killed and is re-run alone to classify). The runs do NOT use NUnit.StopOnError. Stage 1 runs the tests of the changed classes (the registry, the protocol of the gate,
   the tree-level tests, the lifetime and reader tests of the facade, the pin of the compilation of the mapped descent; then, if nothing failed, the child process of the Dispose race); a survivor of stage 1 runs stage 2, the whole default set of
   SuffixTree.Esa.Tests (the categories that `dotnet test` runs by default). Analyzers are off for the mutant builds. The source file is restored after every mutant (and at the end, also on an exception). The mutants run in a separate
   git worktree, never in the working checkout.
   usage: run_1b9g_mutants.py <worktree> <1b9g_<file>_listed.json> <result.tsv> [<start>:<end> slice]
   environment: MUTANT_BUILD_LIMIT (s, default 400), MUTANT_TEST_LIMIT (s, default 420), MUTANT_FULL_LIMIT (s, default 2700), STAGE2=0 to skip the whole set"""
import json, os, re, signal, subprocess, sys


def printable(s):
    """The text of a failure as the harness printed it may hold control characters (the text of a fuzz case): a TSV that carries them is a binary file for git, so they are written as \\xNN escapes."""
    return re.sub(r'[\x00-\x08\x0b-\x1f\x7f]', lambda m: '\\x%02x' % ord(m.group()), s)

ROOT = sys.argv[1]
BUILD_LIMIT = int(os.environ.get('MUTANT_BUILD_LIMIT', '400'))
TEST_LIMIT = int(os.environ.get('MUTANT_TEST_LIMIT', '420'))
FULL_LIMIT = int(os.environ.get('MUTANT_FULL_LIMIT', '2700'))
STAGE2 = os.environ.get('STAGE2', '1') != '0'
STAGE1 = ('FullyQualifiedName~ReaderCellTests|FullyQualifiedName~ReaderGateTests|FullyQualifiedName~FileBackedReaderProtocolTests|FullyQualifiedName~FacadeReaderCountTests|FullyQualifiedName~FacadeLifetimeTests'
          '|FullyQualifiedName~FacadeEsaDisposeTests|FullyQualifiedName~FacadeCallReachabilityTests|FullyQualifiedName~MappedSearchTiering1b9gTests')
STAGE1_RACE = 'FullyQualifiedName~FacadeDisposeRaceTests'  # the child process of the Dispose race (about 20 s): run only when the deterministic tests of stage 1 do not kill
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1')  # SEQERON_BENCH_DIR comes from the caller's environment
mutants_file, out = sys.argv[2], sys.argv[3]
sl = sys.argv[4] if len(sys.argv) > 4 else ''
PROJ = {'esa': ROOT + '/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj'}
mutants = json.load(open(mutants_file))
if sl:
    a, b = sl.split(':')
    mutants = mutants[int(a or 0):int(b) if b else None]
path = os.path.join(ROOT, mutants[0]['file'])
orig = open(path).read()
FLAGS = ['-p:RunAnalyzers=false', '-p:EnforceCodeStyleInBuild=false', '-p:TreatWarningsAsErrors=false']


def sh(args, limit, extra_env=None):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=dict(env, **(extra_env or {})), start_new_session=True)
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


def verdict(txt, hung, what):
    if hung:
        return 'killed(hang, %s)' % what
    if 'Test host process crashed' in txt or 'test run was aborted' in txt:
        mm = re.search(r'(System\.\w+Exception|Fatal error[^\n]*)', txt)
        return 'killed(host crash: %s, %s)' % (mm.group(1) if mm else '?', what)
    if 'Failed!' in txt:
        failed = re.findall(r'^\s+Failed (\S+)', txt, re.M)
        return 'killed(%s: %d failed, first %s)' % (what, len(failed), failed[0] if failed else '?')
    if 'Passed!' in txt:
        return 'survived'
    return 'NO-RESULT'


def unit(proj, filt, limit, what):
    cmd = ['dotnet', 'test', proj, '-c', 'Release', '--nologo', '--no-build'] + (['--filter', filt] if filt else [])
    txt, hung = sh(cmd, limit)
    return verdict(txt, hung, what)


results = []
try:
    for m in mutants:
        if os.path.join(ROOT, m['file']) != path:
            raise SystemExit('one source file per run')
        if orig.count(m['old']) != 1:
            results.append((m['id'], m['desc'], 'PATTERN-NOT-UNIQUE'))
            continue
        open(path, 'w').write(orig.replace(m['old'], m['new']))
        kind = m.get('project', 'esa')
        proj = PROJ[kind]
        b = build(proj)
        if b:
            results.append((m['id'], m['desc'], b))
        else:
            r = unit(proj, STAGE1, TEST_LIMIT, 'stage 1')
            if r == 'survived':
                r = unit(proj, STAGE1_RACE, TEST_LIMIT, 'stage 1, race')
            if r == 'survived' and STAGE2:
                r = unit(proj, '', FULL_LIMIT, 'stage 2')
                if r == 'survived':
                    r = 'survived (stage 1 and stage 2)'
            results.append((m['id'], m['desc'], r))
        open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
finally:
    open(path, 'w').write(orig)
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', mutants[0]['file']], env=env)
