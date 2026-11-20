#!/usr/bin/env python3
"""Single-site mutants of the production code of 1b.7d and 1b.7e (src/SuffixTree/Algorithms/SuffixTree.Esa/External), run by work package 1b.9c (the mutants that the audits of 1b.7d and 1b.7e did not run, the re-run of the kills
   that only a watchdog or a crash of the host had counted, and the mutants of the new tests). Each mutant is built and tested in the Release configuration with a WATCHDOG (a wrong cursor can make a loop non-terminating or crash the
   test host: the whole process group is killed on expiry); the runs do NOT use NUnit.StopOnError (a test host that must stop at the first failure can hang, so a time limit under it is not evidence). Analyzers are off for the
   mutant builds (a mutant that leaves a field unassigned or a local unused is still a valid mutant). The source file is restored after every mutant (and at the end, also on an exception, and at the start, in case an earlier
   run was killed in the middle). The mutants run in a separate git worktree, never in the working checkout.
   Stages (the project is always the whole SuffixTree.Esa.Tests, so that a test class added after the lists were written is built and can be selected). The test sets run one after the other and the first one that fails
   (an assertion, a hang or a crash of the host) ends the mutant: a kill is a kill whichever set found it, and a mutant that is killed by the first set costs one short run (no NUnit.StopOnError is used for that):
     1  the filter of the list as it is (the tests that the author of the list named for the file; narrow)
     1+ one run of the test classes that name a type declared in the mutated file and are not named by the filter of the list (found by reading the test sources, ``seeing_classes``); explicit tests are left out
     2  only for a survivor of those: the whole default set of SuffixTree.Esa.Tests (every test that is not [Explicit]; the tests of the stages before run again)
     3  only with MUTANT_STAGES containing 3 and only for a survivor of 2: the differential tests of SuffixTree.EmProto.Tests (they exercise PH1-PH3 and PH6 of the production stage)
     H  only with MUTANT_STAGES containing H and only for a survivor of the earlier stages: the [Explicit] heavy classes of the External stages (exhaustive and fuzz), a budget per mutant
   A mutant that survives every stage is "survived".
   The result file is APPENDED to (one line after every mutant, flushed and fsynced) and the mutants that it already holds are skipped, so a restart loses at most the mutant that was running.
   usage: run_1b9c_mutants.py <worktree> <listed.json> <result.tsv> [<start>:<end> slice]
   environment: MUTANT_BUILD_LIMIT (s, default 400), MUTANT_TEST_LIMIT (s, default 360), MUTANT_STAGES (default "1": stage 1 is the listed filter and the classes that see the file; "12" adds the whole default set, "2" runs the whole default set only (for the survivors of an earlier run of stage 1), "23", "2H" ...), MUTANT_FILTER_ONLY (replaces the filter of stage 1, used by the one-by-one re-runs),
   MUTANT_TMPDIR (a directory for the temporary files of the tests, e.g. on /dev/shm), MUTANT_REPEAT (default 1: how many times a mutant is run, the verdicts are all written)."""
import glob, json, os, re, shutil, signal, subprocess, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mutlib_1b9c import seeing_classes as _seeing

ROOT = sys.argv[1]
BUILD_LIMIT = int(os.environ.get('MUTANT_BUILD_LIMIT', '400'))
TEST_LIMIT = int(os.environ.get('MUTANT_TEST_LIMIT', '360'))
STAGES = os.environ.get('MUTANT_STAGES', '1')
FILTER_ONLY = os.environ.get('MUTANT_FILTER_ONLY', '')
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', SEQERON_BENCH_DIR='/tmp/stbench')
env.pop('DOTNET_GCHeapHardLimit', None)
TMP = os.environ.get('MUTANT_TMPDIR', '')  # a directory of a memory file system for the temporary files of the tests (the tests create and map many small files: a disk makes them 3x slower); emptied after every mutant
if TMP:
    os.makedirs(TMP, exist_ok=True)
    env['TMPDIR'] = TMP
mutants_file, out = sys.argv[2], sys.argv[3]
sl = sys.argv[4] if len(sys.argv) > 4 else ''
TESTDIR = os.path.join(ROOT, 'tests/SuffixTree/SuffixTree.Esa.Tests')
PROJ_ESA = TESTDIR + '/SuffixTree.Esa.Tests.csproj'
PROJ_EM = ROOT + '/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj'
FLAGS = ['-p:RunAnalyzers=false', '-p:EnforceCodeStyleInBuild=false', '-p:TreatWarningsAsErrors=false']
NOT_HEAVY = 'TestCategory!=Exhaustive&TestCategory!=Fuzz&TestCategory!=Cost'
STAGE3 = '(FullyQualifiedName~SaBuildDifferentialTests|FullyQualifiedName~ProductionPh6Tests)'
HEAVY = ('(FullyQualifiedName~SaExhaustiveTests|FullyQualifiedName~SesaExhaustiveTests|FullyQualifiedName~LcpChildExhaustiveTests|FullyQualifiedName~SaBuildFuzzTests|FullyQualifiedName~SesaBuildFuzzTests'
         '|FullyQualifiedName~LcpChildFuzzTests)')
HEAVY_ENV = {'ESA_FUZZ_BUDGET_S': os.environ.get('HEAVY_BUDGET_S', '60'), 'ESA_FUZZ_THREADS': '2', 'ESA_FUZZ_MAXFAIL': '1', 'ESA_FUZZ_HANG_S': '90'}


def sh(args, limit, extra_env=None):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace', env=dict(env, **(extra_env or {})), start_new_session=True)
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


def verdict(txt, hung):
    """killed(hang) / killed(host crash: <exception>) / killed(<n> failed: <first test names>) / survived; the first failing test names are kept (a kill by one flaky test must be visible)."""
    if hung:
        return 'killed(hang)'
    if 'Test host process crashed' in txt or 'test run was aborted' in txt:
        mm = re.search(r'(System\.\w+Exception|Fatal error[^\n]*|Stack overflow)', txt)
        return 'killed(host crash: %s)' % (mm.group(1)[:60] if mm else '?')
    if 'Failed!' in txt:
        failed = re.findall(r'^\s+Failed (\S+)', txt, re.M)
        names = ','.join(x.split('.')[-1][:60] for x in failed[:3])
        return 'killed(%d failed: %s)' % (len(failed), names)
    if 'Passed!' in txt:
        return 'survived'
    return 'NO-RESULT'


def test(proj, filt, extra_env=None):
    cmd = ['dotnet', 'test', proj, '-c', 'Release', '--nologo', '--no-build']
    if filt:
        cmd += ['--filter', filt]
    txt, hung = sh(cmd, TEST_LIMIT, extra_env)
    return verdict(txt, hung)


def main():
    mutants = json.load(open(mutants_file))
    if sl:
        a, b = sl.split(':')
        mutants = mutants[int(a or 0):int(b) if b else None]
    done = set()
    if os.path.exists(out):
        for line in open(out):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3:
                done.add(p[0])
    todo = [m for m in mutants if m['id'] not in done]
    if not todo:
        return
    rel = todo[0]['file']
    path = os.path.join(ROOT, rel)
    subprocess.run(['git', '-C', ROOT, 'checkout', '--', rel], env=env)  # a run that was killed in the middle may have left the file mutated
    orig = open(path).read()
    seeing = _seeing(path, TESTDIR)
    repeat = int(os.environ.get('MUTANT_REPEAT', '1'))
    try:
        for m in todo:
            if os.path.join(ROOT, m['file']) != path:
                raise SystemExit('one source file per run')
            if orig.count(m['old']) != 1:
                emit(m, 'PATTERN-NOT-UNIQUE(%d)' % orig.count(m['old']), '', '', 0)
                continue
            t0 = time.time()
            open(path, 'w').write(orig.replace(m['old'], m['new']))
            for _ in range(repeat):
                b = build(PROJ_ESA)
                if b:
                    emit(m, b, '', '', time.time() - t0)
                    continue
                listed = m.get('filter', '')
                own_cls = set(x for x in re.findall(r'FullyQualifiedName~([\w.]+)', listed) if '.' not in x)  # the classes that the listed filter names whole (a token with a dot names a method or a part)
                groups = []
                if FILTER_ONLY:
                    groups = [('only', FILTER_ONLY)]
                elif '1' in STAGES:
                    if listed and not os.environ.get('MUTANT_SKIP_LISTED'):  # MUTANT_SKIP_LISTED: the listed filter was run before (stage 1), stage 2 starts with the classes that name the file
                        groups.append(('listed', '%s&%s' % (listed, NOT_HEAVY)))
                    others = sorted(seeing - own_cls) if not os.path.exists('/tmp/1b9c_skip_seen') else []  # the flag file makes stage 1 the listed filter alone (stage 2, the whole default set, follows for every survivor)
                    if others:
                        groups.append(('seen', '(%s)&%s' % ('|'.join('FullyQualifiedName~' + c for c in others), NOT_HEAVY)))
                r1, r2, r3 = 'survived', '', ''
                for name, filt in groups:
                    r = test(PROJ_ESA, filt)
                    if r.startswith('killed'):
                        r1 = '%s[%s]' % (r, name)
                        break
                    if r != 'survived':  # NO-RESULT: the filter matched no test or the run did not report; the next set decides
                        r1 = 'NO-RESULT[%s]' % name
                if '1' not in STAGES and not FILTER_ONLY:
                    r1 = '-'  # stage 1 was run before
                open_ = r1 in ('survived', '-') or r1.startswith('NO-RESULT')
                if open_ and '2' in STAGES and not FILTER_ONLY:
                    r2 = test(PROJ_ESA, '')  # the whole default set of SuffixTree.Esa.Tests (every test that is not [Explicit]), the sets of the stages before included
                    open_ = r2 == 'survived'
                if open_ and '3' in STAGES:
                    b3 = build(PROJ_EM)
                    r3 = b3 if b3 else test(PROJ_EM, STAGE3)
                    open_ = r3 == 'survived'
                if open_ and 'H' in STAGES:
                    r3 = (r3 + ' ' if r3 else '') + 'heavy:' + test(PROJ_ESA, HEAVY, HEAVY_ENV)
                emit(m, r1, r2, r3, time.time() - t0)
    finally:
        open(path, 'w').write(orig)
        subprocess.run(['git', '-C', ROOT, 'checkout', '--', rel], env=env)


def clean_tmp():
    if TMP:
        for name in os.listdir(TMP):
            p = os.path.join(TMP, name)
            shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) and not os.path.islink(p) else (os.path.exists(p) and os.remove(p))


def emit(m, r1, r2, r3, secs):
    clean_tmp()
    with open(out, 'a') as f:
        f.write('\t'.join([m['id'], m['desc'].replace('\t', ' '), r1, r2, r3, '%d' % secs]) + '\n')
        f.flush()
        os.fsync(f.fileno())


main()
