#!/usr/bin/env python3
"""Single-site mutants of the code of work package 1b.9a (prefix 1b9a_), each built and run in the Release configuration with a WATCHDOG (a wrong cursor or a table that the code no longer bounds can make a loop endless or crash
   the test host: the whole process group is killed on expiry, a hang or a host crash counts as killed). The runs do NOT use NUnit.StopOnError. Round 1 runs the unit tests of the new code (MappedSearchTests, SesaChildCodecTests, the
   reader, builder, verifier and hostile tests); a survivor then runs the BUDGETED HARNESS of 1b.8 (the explicit class FullDifferentialFuzzTests, a short budget per mutant, as run_1b8_mutants.py does). The mutants of
   apps/SuffixTree.EmProto run ContainerVerifyTests of SuffixTree.EmProto.Tests. Analyzers are off for the mutant builds (a mutant that leaves a field unassigned or a local unused is still a mutant). The source file is
   restored after every mutant (and at the end, also on an exception). The mutants run in a separate git worktree, never in the working checkout.
   usage: run_1b9a_mutants.py <worktree> <1b9a_<file>_listed.json> <result.tsv> [<start>:<end> slice]
   environment: MUTANT_BUILD_LIMIT (s, default 400), MUTANT_TEST_LIMIT (s, default 420), FUZZ_BUDGET_S (default 60), FUZZ_THREADS (default 2), FUZZ_HEAP_LIMIT (default 0xC0000000), HARNESS=0 to skip the harness for survivors."""
import json, os, re, signal, subprocess, sys


def printable(s):
    """The text of a failure as the harness printed it may hold control characters (the text of a fuzz case): a TSV that carries them is a binary file for git, so they are written as \\xNN escapes."""
    return re.sub(r'[\x00-\x08\x0b-\x1f\x7f]', lambda m: '\\x%02x' % ord(m.group()), s)

ROOT = sys.argv[1]
BUILD_LIMIT = int(os.environ.get('MUTANT_BUILD_LIMIT', '400'))
TEST_LIMIT = int(os.environ.get('MUTANT_TEST_LIMIT', '420'))
BUDGET = os.environ.get('FUZZ_BUDGET_S', '60')
THREADS = os.environ.get('FUZZ_THREADS', '2')
HARNESS = os.environ.get('HARNESS', '1') != '0'
# the unit tests that can see a mutant of the read path (the descent, the mapped child table, the reader's views) and those that can see one of the format, the builder and the verifier
READ_PATH = ('FullyQualifiedName~MappedSearchTests|FullyQualifiedName~SesaChildCodecTests|FullyQualifiedName~SesaReaderTests|FullyQualifiedName~SesaReaderHostileTests|FullyQualifiedName~SesaReaderIsolationTests')
FORMAT = ('FullyQualifiedName~SesaChildCodecTests|FullyQualifiedName~SesaBuilderTests|FullyQualifiedName~SesaAuditTests|FullyQualifiedName~SesaVerifierTests|FullyQualifiedName~SesaReaderCorruptionTests|FullyQualifiedName~SesaSectionTests')
READ_FILES = ('MappedSearch.cs', 'ChildLayouts.cs', 'ChildTable.cs', 'SesaFile.cs')
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', SEQERON_BENCH_DIR='/tmp/stbench')
mutants_file, out = sys.argv[2], sys.argv[3]
sl = sys.argv[4] if len(sys.argv) > 4 else ''
PROJ = {'esa': ROOT + '/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj', 'emproto': ROOT + '/tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj'}
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


def unit(proj, filt):
    txt, hung = sh(['dotnet', 'test', proj, '-c', 'Release', '--nologo', '--no-build', '--filter', filt], TEST_LIMIT)
    return verdict(txt, hung, 'unit tests')


def harness(mid):
    rep = '%s.%s.rep' % (out, mid)
    for f in (rep, rep + '.log'):
        if os.path.exists(f):
            os.remove(f)
    xenv = {'ESA_FUZZ_BUDGET_S': BUDGET, 'ESA_FUZZ_THREADS': THREADS, 'ESA_FUZZ_MAXFAIL': '1', 'ESA_FUZZ_MINIMIZE_S': '4', 'ESA_FUZZ_HANG_S': '90',
            'DOTNET_GCHeapHardLimit': os.environ.get('FUZZ_HEAP_LIMIT', '0xC0000000'), 'ESA_FUZZ_REPORT': rep}
    txt, hung = sh(['dotnet', 'test', PROJ['esa'], '-c', 'Release', '--nologo', '--no-build', '--filter', 'FullyQualifiedName~FullDifferentialFuzzTests'], TEST_LIMIT, xenv)
    first = ''
    if os.path.exists(rep + '.log'):
        for line in open(rep + '.log'):
            if line.startswith('FAIL '):
                first = printable(re.sub(r'\s+', ' ', line.strip()))[:200]
                break
    if hung:
        return 'killed(hang, harness)'
    if 'HANG:' in txt:
        return 'killed(watchdog of the harness)'
    if 'Test host process crashed' in txt or 'test run was aborted' in txt:
        return 'killed(host crash, harness)'
    if 'Failed!' in txt:
        return 'killed(harness: %s)' % (first or 'test failed')
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
        kind = m.get('project', 'esa')
        proj = PROJ[kind]
        b = build(proj)
        if b:
            results.append((m['id'], m['desc'], b))
        else:
            r = unit(proj, (READ_PATH if m['file'].endswith(READ_FILES) else FORMAT) if kind == 'esa' else 'FullyQualifiedName~ContainerVerifyTests')
            if r == 'survived' and kind == 'esa' and HARNESS:
                r = harness(m['id'])
                if r == 'survived':
                    r = 'survived (unit tests and harness)'
            results.append((m['id'], m['desc'], r))
        open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
finally:
    open(path, 'w').write(orig)
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', mutants[0]['file']], env=env)
