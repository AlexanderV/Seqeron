#!/usr/bin/env python3
"""Single-site mutants of the test-support code of work package S2.1 (prefix s21_), each built and run in the Release configuration with a WATCHDOG (a mutant may hang the test host: the whole process group is killed on
   expiry, a hang or a host crash counts as killed). Stage 1 runs the test class that owns the mutated file (the filter of the mutant); a survivor of stage 1 runs stage 2, the four classes of the package
   (ClosedFormSelfTests, Crc32CCombineTests, SparseContainerTests, Sesa5EncoderTests); a survivor of both is reported as `survived` and is a gap to close or an equivalent mutant to argue.
   Analyzers are off for the mutant builds. The source file is restored after every mutant (and at the end, also on an exception). Run it in a separate git worktree, never in the working checkout.
   usage: run_s21_mutants.py <worktree> <s21_listed.json> <result.tsv> [<start>:<end> slice]
   environment: MUTANT_BUILD_LIMIT (s, default 400), MUTANT_TEST_LIMIT (s, default 300); paths are repo-relative in the list"""
import json, os, re, signal, subprocess, sys


def printable(s):
    """A failure text may hold control characters: a TSV that carries them is a binary file for git, so they are written as \\xNN escapes."""
    return re.sub(r'[\x00-\x08\x0b-\x1f\x7f]', lambda m: '\\x%02x' % ord(m.group()), s)


ROOT = sys.argv[1]
BUILD_LIMIT = int(os.environ.get('MUTANT_BUILD_LIMIT', '400'))
TEST_LIMIT = int(os.environ.get('MUTANT_TEST_LIMIT', '300'))
STAGE2 = ('FullyQualifiedName~Large.ClosedFormSelfTests|FullyQualifiedName~Large.Crc32CCombineTests|FullyQualifiedName~Large.SparseContainerTests|FullyQualifiedName~Large.Sesa5EncoderTests')
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1')
mutants_file, out = sys.argv[2], sys.argv[3]
sl = sys.argv[4] if len(sys.argv) > 4 else ''
PROJ = ROOT + '/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj'
mutants = json.load(open(mutants_file))
if sl:
    a, b = sl.split(':')
    mutants = mutants[int(a or 0):int(b) if b else None]
FLAGS = ['-p:RunAnalyzers=false', '-p:EnforceCodeStyleInBuild=false', '-p:TreatWarningsAsErrors=false']
originals = {}


def sh(args, limit):
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, start_new_session=True)
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


def verdict(txt, hung, what):
    if hung:
        return 'killed(hang, %s)' % what
    if 'Test host process crashed' in txt or 'test run was aborted' in txt:
        mm = re.search(r'(System\.\w+Exception|Fatal error[^\n]*)', txt)
        return 'killed(host crash: %s, %s)' % (mm.group(1) if mm else '?', what)
    if 'Failed!' in txt:
        failed = re.findall(r'^\s+Failed (\S+)', txt, re.M)
        return 'killed(%s: %d failed, first %s)' % (what, len(failed), printable(failed[0]) if failed else '?')
    if 'Passed!' in txt:
        return 'survived'
    return 'NO-RESULT'


def unit(filt, limit, what):
    txt, hung = sh(['dotnet', 'test', PROJ, '-c', 'Release', '--nologo', '--no-build'] + (['--filter', filt] if filt else []), limit)
    return verdict(txt, hung, what)


results = []
try:
    for m in mutants:
        path = os.path.join(ROOT, m['file'])
        if path not in originals:
            originals[path] = open(path).read()
        orig = originals[path]
        if orig.count(m['old']) != 1:
            results.append((m['id'], m['desc'], 'PATTERN-NOT-UNIQUE'))
            continue
        open(path, 'w').write(orig.replace(m['old'], m['new']))
        b = build()
        if b:
            results.append((m['id'], m['desc'], b))
        else:
            r = unit(m['filter'], TEST_LIMIT, 'stage 1, the owning class')
            if r == 'survived':
                r = unit(STAGE2, TEST_LIMIT, 'stage 2, the four classes')
                if r == 'survived':
                    r = 'survived (stages 1 and 2)'
            results.append((m['id'], m['desc'], r))
        open(path, 'w').write(orig)
        open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
finally:
    for path, orig in originals.items():
        open(path, 'w').write(orig)
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', 'tests'], env=env)
