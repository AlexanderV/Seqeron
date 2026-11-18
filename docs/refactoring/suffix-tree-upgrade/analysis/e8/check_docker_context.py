#!/usr/bin/env python3
"""Static check of the image build WITHOUT a Docker daemon (E8): applies Dockerfile.dockerignore to the tracked files of the checkout with the matching rules of the Docker client (patterns in order, a later pattern wins, a file also
matches through its parent directories, '!' re-includes, '**' crosses directories), copies exactly the surviving files into a clean directory, puts the substitute Directory.Build.props of the Dockerfile in place, and runs the same
`dotnet publish` as the build stage with an EMPTY NuGet cache and an EMPTY package source: the build then fails if it needs a package, or a file that the allow-list does not let in. It then starts the result.
It does NOT test Docker itself (BuildKit's reading of the ignore file, the base images, the apt step, the scripts inside the image).
usage: python3 -I check_docker_context.py [<repo root>] [<work dir>]        needs git and dotnet on the PATH"""
import os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.abspath(os.path.join(HERE, '../../../../..'))
WORK = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix='e8ctx')


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


pats = load_patterns(os.path.join(HERE, 'Dockerfile.dockerignore'))
files = subprocess.run(['git', '-C', ROOT, 'ls-files', '-z'], capture_output=True, check=True).stdout.decode().split('\0')
files = [f for f in files if f]
kept = [f for f in files if not excluded(f, pats)]
print('tracked files: %d, sent to the build context: %d (%.1f MB)' % (len(files), len(kept), sum(os.path.getsize(os.path.join(ROOT, f)) for f in kept) / 1e6))
tops = {}
for f in kept:
    tops[f.split('/')[0] + '/' + (f.split('/')[1] if '/' in f else '')] = tops.get(f.split('/')[0] + '/' + (f.split('/')[1] if '/' in f else ''), 0) + 1
for k in sorted(tops):
    print('   %-45s %d files' % (k, tops[k]))
bad = [f for f in kept if f.startswith(('.git/', 'wiki/', 'tests/')) or '/bin/' in f or '/obj/' in f or f.endswith(('.pdb', '.tsv'))]
print('files that must not be in the context: %s' % (bad or 'none'))
ctx = os.path.join(WORK, 'context')
shutil.rmtree(ctx, ignore_errors=True)
for f in kept:
    dst = os.path.join(ctx, f)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(os.path.join(ROOT, f), dst)
shutil.copy2(os.path.join(ROOT, 'docs/refactoring/suffix-tree-upgrade/analysis/e8/docker/Directory.Build.props'), os.path.join(ctx, 'Directory.Build.props'))
empty_src = os.path.join(WORK, 'empty-source'); os.makedirs(empty_src, exist_ok=True)
pkgs = os.path.join(WORK, 'empty-packages'); shutil.rmtree(pkgs, ignore_errors=True); os.makedirs(pkgs)
out = os.path.join(WORK, 'out'); shutil.rmtree(out, ignore_errors=True)
env = dict(os.environ, NUGET_PACKAGES=pkgs, DOTNET_CLI_HOME=os.path.join(WORK, 'home'), DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', NUGET_SCRATCH=os.path.join(WORK, 'nugetscratch'),
           NUGET_HTTP_CACHE_PATH=os.path.join(WORK, 'httpcache'))
cmd = ['dotnet', 'publish', 'apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj', '-c', 'Release', '-o', out, '--nologo', '-p:UseAppHost=false', '-p:IsTrimmable=false', '-p:IsAotCompatible=false', '-p:EnableTrimAnalyzer=false',
       '-p:EnableAotAnalyzer=false', '-p:GenerateDocumentationFile=false', '-p:TreatWarningsAsErrors=false', '-p:EnforceCodeStyleInBuild=false', '-p:RestoreSources=' + empty_src, '--no-cache' if False else '-v:q']
print('running: ' + ' '.join(cmd[:4]) + ' ... (empty NuGet cache and source)')
p = subprocess.run(cmd, cwd=ctx, env=env, capture_output=True, text=True)
print(p.stdout[-1500:] + p.stderr[-1500:])
ok = p.returncode == 0 and os.path.exists(os.path.join(out, 'SuffixTree.EmProto.dll'))
print('publish: %s; packages downloaded: %s' % ('OK' if ok else 'FAILED', os.listdir(pkgs) or 'none'))
if ok:
    print('output: %s' % ', '.join(sorted(f for f in os.listdir(out) if f.endswith(('.dll', '.json')))))
    r = subprocess.run(['dotnet', os.path.join(out, 'SuffixTree.EmProto.dll'), 'e8-env', '--case', 'ctx'], capture_output=True, text=True, env=env)
    rows = [x.split('\t') for x in r.stdout.splitlines() if x.startswith('R\t')]
    print('e8-env exit code %d: %s' % (r.returncode, ' | '.join('%s=%s' % (c[3], c[4]) for c in rows if c[3] in ('env_framework', 'env_process_arch', 'env_lce_path'))))
print('work directory: ' + WORK)
sys.exit(0 if ok and not bad else 1)
