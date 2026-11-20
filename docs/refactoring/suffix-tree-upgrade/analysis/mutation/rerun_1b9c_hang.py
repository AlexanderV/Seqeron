#!/usr/bin/env python3
"""Work package 1b.9c: the one-by-one re-run of the mutants that an earlier audit counted as killed only because a watchdog expired or the test host crashed (14 of 1b.7d, 3 of 1b.7e). Nothing else runs on the machine
   (the runs are made one after the other, not under load). One mutant: build, then ONE run of every test that can see it (the filter of its list, the test classes that name the types of its file and, with stage "A", the whole
   default set), with --blame (the sequence of the tests: the one that was running when the host died or hung), a time limit per test (--blame-hang-timeout) and no limit on the run itself but a large one, and no NUnit.StopOnError.
   The line written is the classification: assertion (a test fails on its own, the run ends), crash (the host dies: the exception and the last test that was running), hang (a test does not end: the last test running), none
   (nothing failed: the mutant survived this time, so the earlier kill was the environment), with the names of all failing tests.
   usage: rerun_1b9c_hang.py <worktree> <list.json> <result.tsv> [A]   (the list holds the mutants: id, desc, file, old, new, filter)
   environment: HANG_TEST_S (seconds per test before the host is killed, default 150), RUN_LIMIT_S (default 1800), MUTANT_TMPDIR, CHECKED=1 (the checked build and ESA_EXPECT_BOUNDS=true: the second pass over the crashes and hangs of the first)"""
import glob, json, os, re, signal, subprocess, sys, time, xml.etree.ElementTree as ET
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mutlib_1b9c import seeing_classes

ROOT, listfile, out = sys.argv[1], sys.argv[2], sys.argv[3]
ALL = len(sys.argv) > 4 and sys.argv[4] == 'A'
HANG_S = os.environ.get('HANG_TEST_S', '150')
RUN_LIMIT = int(os.environ.get('RUN_LIMIT_S', '1800'))
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', SEQERON_BENCH_DIR='/tmp/stbench')
env.pop('DOTNET_GCHeapHardLimit', None)
TMP = os.environ.get('MUTANT_TMPDIR', '')
if TMP:
    os.makedirs(TMP, exist_ok=True)
    env['TMPDIR'] = TMP
PROJ = ROOT + '/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj'
RESULTS = '/dev/shm/1b9c-blame'
NOT_HEAVY = 'TestCategory!=Exhaustive&TestCategory!=Fuzz&TestCategory!=Cost'
FLAGS = ['-p:RunAnalyzers=false', '-p:EnforceCodeStyleInBuild=false', '-p:TreatWarningsAsErrors=false']
if os.environ.get('CHECKED') == '1':  # the checked build (-p:EsaBounds=true, ESA_EXPECT_BOUNDS=true): the reads of the text and of the tables are bounds-checked, so a wrong cursor is an exception where the production build reads outside its mapping
    FLAGS += ['-p:EsaBounds=true']
    env['ESA_EXPECT_BOUNDS'] = 'true'


def sh(args, limit):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace', env=env, start_new_session=True)
    try:
        txt, _ = p.communicate(timeout=limit)
        return txt, False
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        txt, _ = p.communicate()
        return txt, True


def last_running(dirpath):
    """The test that was running when the host died or was killed: the last entry of the sequence file that has no completed marker."""
    for f in sorted(glob.glob(os.path.join(dirpath, '**', 'Sequence_*.xml'), recursive=True)):
        try:
            tests = ET.parse(f).getroot().findall('.//Test')
        except ET.ParseError:
            continue
        pending = [t.get('Name') for t in tests if t.get('Completed') in (None, 'False')]
        if pending:
            return pending[-1]
    return ''


def classify(txt, hung):
    failed = re.findall(r'^\s+Failed (\S+)', txt, re.M)
    crashed = 'Test host process crashed' in txt or 'test run was aborted' in txt.lower() or 'aborted' in txt.lower() and 'Failed!' not in txt and 'Passed!' not in txt
    mm = re.search(r'(Stack overflow|Fatal error\.?[^\n]*\n[^\n]{0,160}|Fatal error[^\n]*|System\.\w+Exception[^\n]{0,80}|OutOfMemory[^\n]{0,60}|Segmentation fault[^\n]*|Aborted[^\n]*)', txt)
    if hung:
        return 'run-limit', failed, mm.group(1).replace(chr(10), ' ').replace(chr(9), ' ') if mm else ''
    if 'hang' in txt.lower() and 'blame' in txt.lower() and 'testhost' in txt.lower() and ('was killed' in txt.lower() or 'timeout' in txt.lower() and 'hang' in txt.lower()):
        return 'hang', failed, ''
    if crashed:
        return 'crash', failed, mm.group(1).replace(chr(10), ' ').replace(chr(9), ' ') if mm else ''
    if 'Failed!' in txt:
        return 'assertion', failed, ''
    if 'Passed!' in txt:
        return 'none', [], ''
    return 'no-result', failed, ''


mutants = json.load(open(listfile))
done = set(l.split('\t')[0] for l in open(out)) if os.path.exists(out) else set()
for m in mutants:
    kid = (m['package'] + ':' if 'package' in m else '') + m['id']  # the ids of 1b.7d and 1b.7e overlap (S09)
    if kid in done:
        continue
    path = os.path.join(ROOT, m['file'])
    subprocess.run(['git', '-C', ROOT, 'checkout', '--', m['file']], env=env)
    orig = open(path).read()
    t0 = time.time()
    try:
        if orig.count(m['old']) != 1:
            line = [kid, m['desc'], 'pattern', '', '', '', '0']
        else:
            open(path, 'w').write(orig.replace(m['old'], m['new']))
            txt, _ = sh(['dotnet', 'build', PROJ, '-c', 'Release', '--nologo', '-v', 'q'] + FLAGS, 600)
            if re.search(r'error (CS|IDE|CA|S)\d+', txt) or 'Build FAILED' in txt:
                line = [kid, m['desc'], 'build-error', '', '', '', '0']
            else:
                classes = sorted(set(re.findall(r'FullyQualifiedName~(\w+)', m.get('filter', ''))) | seeing_classes(path, os.path.join(ROOT, 'tests/SuffixTree/SuffixTree.Esa.Tests')))
                filt = '' if ALL else '(%s)&%s' % ('|'.join('FullyQualifiedName~' + c for c in classes), NOT_HEAVY)
                shutil_dir = os.path.join(RESULTS, kid.replace(':', '-'))
                subprocess.run(['rm', '-rf', shutil_dir])
                args = ['dotnet', 'test', PROJ, '-c', 'Release', '--nologo', '--no-build', '--blame', '--blame-hang', '--blame-hang-timeout', HANG_S + 's', '--blame-hang-dump-type', 'none', '--results-directory', shutil_dir]
                if filt:
                    args += ['--filter', filt]
                txt, hung = sh(args, RUN_LIMIT)
                kind, failed, why = classify(txt, hung)
                where = last_running(shutil_dir) if kind in ('crash', 'hang', 'run-limit') else ''
                line = [kid, m['desc'], kind, ','.join(x.split('.')[-1][:70] for x in failed[:6]), why[:220], where[:110], '%d' % (time.time() - t0)]
    finally:
        open(path, 'w').write(orig)
        subprocess.run(['git', '-C', ROOT, 'checkout', '--', m['file']], env=env)
    with open(out, 'a') as f:
        f.write('\t'.join(line) + '\n')
        f.flush()
        os.fsync(f.fileno())
    if TMP:
        subprocess.run(['rm', '-rf', TMP], env=env)
        os.makedirs(TMP, exist_ok=True)
