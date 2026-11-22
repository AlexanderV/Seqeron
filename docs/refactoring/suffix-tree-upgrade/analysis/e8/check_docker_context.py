#!/usr/bin/env python3
"""Static check of an image build WITHOUT a Docker daemon (E8; extended for the package s1b-handoff, work package 1b.9f): applies a Dockerfile.dockerignore to the tracked files of the checkout with the
matching rules of the Docker client (patterns in order, a later pattern wins, matches through its parent directories, '!' re-includes, '**' crosses directories), copies exactly the surviving files into a clean
directory, puts the substitute Directory.Build.props of the Dockerfile in place, and runs `dotnet publish` as the build stage with an EMPTY NuGet cache and an EMPTY package source: the build then fails if it
needs a package, or a file that the allow-list does not let in. It then starts the result.
It does NOT test Docker itself (BuildKit's reading of the ignore file, the base images, the apt step, the scripts inside the image).
usage: python3 -I check_docker_context.py [<repo root>] [<work dir>] [options]        needs git and dotnet on the PATH
  no option: the E8 context (Dockerfile.dockerignore next to this script, publishes apps/SuffixTree.EmProto, rejects test files and .tsv files in the context)
  --package-dir DIR   the directory (relative to the repository root) whose Dockerfile.dockerignore is applied (default: docs/refactoring/suffix-tree-upgrade/analysis/e8)
  --publish PROJECT   a project (relative .csproj path) to publish in the chain build; repeatable (default: apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj)
  --allow-tests       the context may hold test projects (the s1b-handoff image runs the test sets); --allow-tsv  data files of the tests may be in the context
  --tests-source DIR  also build the test projects from a PRISTINE copy of the context (the repository's own props) restoring from the local package source DIR (the global packages folder of a machine
                      that ran the tests: a hierarchical feed), then run the new Handoff tests of SuffixTree.Esa.Tests: proves that the context holds everything the tests need
"""
import os, re, shutil, subprocess, sys, tempfile

_argv = sys.argv[1:]
OPT = {'publish': [], 'package-dir': 'docs/refactoring/suffix-tree-upgrade/analysis/e8', 'allow-tests': False, 'allow-tsv': False, 'tests-source': None}
_pos = []
_i = 0
while _i < len(_argv):
    a = _argv[_i]
    if a in ('--allow-tests', '--allow-tsv'):
        OPT[a[2:]] = True
    elif a in ('--publish', '--package-dir', '--tests-source'):
        _i += 1
        if a == '--publish':
            OPT['publish'].append(_argv[_i])
        else:
            OPT[a[2:]] = _argv[_i]
    else:
        _pos.append(a)
    _i += 1
if not OPT['publish']:
    OPT['publish'] = ['apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj']

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(_pos[0]) if len(_pos) > 0 else os.path.abspath(os.path.join(HERE, '../../../../..'))
WORK = _pos[1] if len(_pos) > 1 else tempfile.mkdtemp(prefix='e8ctx')


def glob_re(pat):
    out, i = '', 0
    while i < len(pat):
        c = pat[i]
        if pat.startswith('**/', i):
            out += '(?:.*/)?'; i += 3; continue
        if pat.startswith('**', i):
            out += '.*'; i += 2; continue
        out += {'*': '[^/]*', '?': '[^/]'}.get(c, re.escape(c)); i += 1
    return re.compile('^' + out + '$')


def load_patterns(path):
    pats = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        neg = line.startswith('!')
        pats.append((neg, glob_re(os.path.normpath(line[1:] if neg else line).replace(os.sep, '/').lstrip('/')), line))
    return pats


def excluded(path, pats):
    parents = ['/'.join(path.split('/')[:k]) for k in range(1, path.count('/') + 1)]
    state = False
    for neg, rx, _ in pats:
        if rx.match(path) or any(rx.match(p) for p in parents):
            state = not neg
    return state


PKG = os.path.join(ROOT, OPT['package-dir'])
pats = load_patterns(os.path.join(PKG, 'Dockerfile.dockerignore'))
files = subprocess.run(['git', '-C', ROOT, 'ls-files', '-z'], capture_output=True, check=True).stdout.decode().split('\0')
files = [f for f in files if f]
kept = [f for f in files if not excluded(f, pats)]
print('tracked files: %d, sent to the build context: %d (%.1f MB)' % (len(files), len(kept), sum(os.path.getsize(os.path.join(ROOT, f)) for f in kept) / 1e6))
tops = {}
for f in kept:
    tops[f.split('/')[0] + '/' + (f.split('/')[1] if '/' in f else '')] = tops.get(f.split('/')[0] + '/' + (f.split('/')[1] if '/' in f else ''), 0) + 1
for k in sorted(tops):
    print('   %-45s %d files' % (k, tops[k]))
bad = [f for f in kept if f.startswith(('.git/', 'wiki/') + (() if OPT['allow-tests'] else ('tests/',))) or '/bin/' in f or '/obj/' in f or f.endswith(('.pdb',) + (() if OPT['allow-tsv'] else ('.tsv',)))]
print('files that must not be in the context: %s' % (bad or 'none'))
ctx = os.path.join(WORK, 'context')
shutil.rmtree(ctx, ignore_errors=True)
for f in kept:
    dst = os.path.join(ctx, f)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(os.path.join(ROOT, f), dst)
shutil.copy2(os.path.join(ROOT, 'docs/refactoring/suffix-tree-upgrade/analysis/e8/docker/Directory.Build.props'), os.path.join(ctx, 'Directory.Build.props'))
if os.path.exists(os.path.join(ctx, 'Directory.Packages.props')):
    os.remove(os.path.join(ctx, 'Directory.Packages.props'))      # the chain stage of the Dockerfile removes it, too: nothing the chain builds uses a package
empty_src = os.path.join(WORK, 'empty-source'); os.makedirs(empty_src, exist_ok=True)
pkgs = os.path.join(WORK, 'empty-packages'); shutil.rmtree(pkgs, ignore_errors=True); os.makedirs(pkgs)
env = dict(os.environ, NUGET_PACKAGES=pkgs, DOTNET_CLI_HOME=os.path.join(WORK, 'home'), DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', NUGET_SCRATCH=os.path.join(WORK, 'nugetscratch'),
           NUGET_HTTP_CACHE_PATH=os.path.join(WORK, 'httpcache'))
env.pop('SEQERON_BENCH_DIR', None)
ok = True
outs = []
for proj in OPT['publish']:
    out = os.path.join(WORK, 'out-' + os.path.basename(os.path.dirname(proj)))
    shutil.rmtree(out, ignore_errors=True)
    cmd = ['dotnet', 'publish', proj, '-c', 'Release', '-o', out, '--nologo', '-p:UseAppHost=false', '-p:IsTrimmable=false', '-p:IsAotCompatible=false', '-p:EnableTrimAnalyzer=false',
           '-p:EnableAotAnalyzer=false', '-p:GenerateDocumentationFile=false', '-p:TreatWarningsAsErrors=false', '-p:EnforceCodeStyleInBuild=false', '-p:RestoreSources=' + empty_src, '-v:q']
    print('running: ' + ' '.join(cmd[:4]) + ' ... (empty NuGet cache and source)')
    p = subprocess.run(cmd, cwd=ctx, env=env, capture_output=True, text=True)
    print(p.stdout[-1500:] + p.stderr[-1500:])
    dll = os.path.join(out, os.path.basename(proj)[:-len('.csproj')] + '.dll')
    good = p.returncode == 0 and os.path.exists(dll)
    ok = ok and good
    print('publish %s: %s; packages downloaded: %s' % (proj, 'OK' if good else 'FAILED', os.listdir(pkgs) or 'none'))
    if good:
        outs.append(dll)
        print('output: %s' % ', '.join(sorted(f for f in os.listdir(out) if f.endswith(('.dll', '.json')))))
for dll in outs:
    if dll.endswith('SuffixTree.EmProto.dll'):
        r = subprocess.run(['dotnet', dll, 'e8-env', '--case', 'ctx'], capture_output=True, text=True, env=env)
        rows = [x.split('\t') for x in r.stdout.splitlines() if x.startswith('R\t')]
        print('e8-env exit code %d: %s' % (r.returncode, ' | '.join('%s=%s' % (c[3], c[4]) for c in rows if c[3] in ('env_framework', 'env_process_arch', 'env_lce_path'))))
    elif dll.endswith('SuffixTree.Esa.Bench.dll'):
        r = subprocess.run(['dotnet', dll, 'sesa', 'handoff', 'env'], capture_output=True, text=True, env=env)
        rows = [x.split('\t') for x in r.stdout.splitlines() if x.startswith('R\t')]
        print('sesa handoff env exit code %d: %s' % (r.returncode, ' | '.join('%s=%s' % (c[2], c[3]) for c in rows if c[2] in ('env_process_arch', 'env_crc32c_hardware', 'env_crc32c_check_value'))))
        ok = ok and r.returncode == 0
if OPT['tests-source']:
    # the test projects from a pristine copy of the context: the repository's own props (analyzer package included), the restore from a local feed only
    tctx = os.path.join(WORK, 'tests-context')
    shutil.rmtree(tctx, ignore_errors=True)
    for f in kept:
        dst = os.path.join(tctx, f)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(os.path.join(ROOT, f), dst)
    tenv = dict(env, NUGET_PACKAGES=os.path.join(WORK, 'tests-packages'))
    shutil.rmtree(tenv['NUGET_PACKAGES'], ignore_errors=True)
    flags = ['-p:RestoreSources=' + OPT['tests-source'], '-p:NuGetAudit=false', '-p:TreatWarningsAsErrors=false', '-p:EnforceCodeStyleInBuild=false', '-v:q', '--nologo']
    for proj in ('tests/SuffixTree/SuffixTree.Esa.Tests', 'tests/SuffixTree/SuffixTree.EmProto.Tests', 'tests/SuffixTree/SuffixTree.Persistent.Tests'):
        if not os.path.isdir(os.path.join(tctx, proj)):
            print('tests: %s is not in the context' % proj)
            ok = False
            continue
        p = subprocess.run(['dotnet', 'build', proj, '-c', 'Debug', *flags], cwd=tctx, env=tenv, capture_output=True, text=True)
        print('build %s: %s' % (proj, 'OK' if p.returncode == 0 else 'FAILED\n' + (p.stdout + p.stderr)[-1500:]))
        ok = ok and p.returncode == 0
    p = subprocess.run(['dotnet', 'test', 'tests/SuffixTree/SuffixTree.Esa.Tests', '--no-build', '--filter', 'FullyQualifiedName~Handoff', '--nologo'], cwd=tctx, env=tenv, capture_output=True, text=True)
    last = [x for x in p.stdout.splitlines() if 'Passed!' in x or 'Failed!' in x]
    print('Handoff tests from the context: %s' % (last[-1].strip() if last else (p.stdout + p.stderr)[-600:]))
    ok = ok and p.returncode == 0
print('work directory: ' + WORK)
sys.exit(0 if ok and not bad else 1)
