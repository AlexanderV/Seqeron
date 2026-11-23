#!/usr/bin/env python3
"""Static check of an image build WITHOUT a Docker daemon (E8; extended for the package s1b-handoff, work package 1b.9f): applies a Dockerfile.dockerignore to the tracked files of the checkout with the
matching rules of the Docker client (patterns in order, a later pattern wins, matches through its parent directories, '!' re-includes, '**' crosses directories), copies exactly the surviving files into a clean
directory, puts the substitute Directory.Build.props of the Dockerfile in place, and runs `dotnet publish` as the build stage with an EMPTY NuGet cache and an EMPTY package source: the build then fails if it
needs a package, or a file that the allow-list does not let in. It then starts the result.
Before the publish (work package 1b.9h) the CLOSURE is checked statically: every project reference of the projects that the image builds (and of the test projects with --allow-tests), the sources of those projects, the linked
sources outside their directories (`Compile Include` / `Link`) and every source of a `COPY` line of the Dockerfile (`COPY --from` and `COPY . .` aside) must be in the context that the allow-list leaves; what is not is printed
as MISSING with the file that needs it, and the check fails. A line missing from the allow-list used to be found only by the publish (CS0234), a COPY of a file that the list drops only by the image build itself.
It does NOT test Docker itself (BuildKit's reading of the ignore file, the base images, the apt step, the scripts inside the image).
usage: python3 -I check_docker_context.py [<repo root>] [<work dir>] [options]        needs git and dotnet on the PATH
  no option: the E8 context (Dockerfile.dockerignore next to this script, publishes apps/SuffixTree.EmProto, rejects test files and .tsv files in the context)
  --package-dir DIR   the directory (relative to the repository root) whose Dockerfile.dockerignore is applied (default: docs/refactoring/suffix-tree-upgrade/analysis/e8)
  --publish PROJECT   a project (relative .csproj path) to publish in the chain build; repeatable (default: apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj)
  --allow-tests       the context may hold test projects (the s1b-handoff image runs the test sets); --allow-tsv  data files of the tests may be in the context
  --static-only       only the static closure (below), no copy and no publish: needs git, not dotnet; exit code 1 when a needed path is not in the context
  --ignore-file FILE  apply this allow-list instead of the Dockerfile.dockerignore of the package (the tests of the closure use it); --dockerfile FILE  the Dockerfile whose COPY lines are checked (default: the one next to the allow-list)
  --tests-source DIR  also build the test projects from a PRISTINE copy of the context (the repository's own props) restoring from the local package source DIR (the global packages folder of a machine
                      that ran the tests: a hierarchical feed), then run the new Handoff tests of SuffixTree.Esa.Tests: proves that the context holds everything the tests need
"""
import os, re, shutil, subprocess, sys, tempfile

_argv = sys.argv[1:]
OPT = {'publish': [], 'package-dir': 'docs/refactoring/suffix-tree-upgrade/analysis/e8', 'allow-tests': False, 'allow-tsv': False, 'tests-source': None, 'static-only': False, 'ignore-file': None, 'dockerfile': None}
_pos = []
_i = 0
while _i < len(_argv):
    a = _argv[_i]
    if a in ('--allow-tests', '--allow-tsv', '--static-only'):
        OPT[a[2:]] = True
    elif a in ('--publish', '--package-dir', '--tests-source', '--ignore-file', '--dockerfile'):
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


def _norm(base_dir, rel):
    """A path of a project file (backslashes allowed), relative to the directory `base_dir`, as a normalized repository path with forward slashes."""
    return os.path.normpath(os.path.join(base_dir, rel.replace('\\', '/'))).replace(os.sep, '/')


def _wild(path):
    """The regular expression of a path with `*`, `**` and `?` (the rules of glob_re) for the tracked files."""
    return glob_re(path)


def project_needs(tracked, proj, seen, needs):
    """Walks the project `proj` (a repository path of a .csproj) and what it pulls in: the project file, its own sources (every tracked .cs / .resx below its directory that is not build output), the project references (recursively),
    the `Compile` / `EmbeddedResource` / `None` / `Content` items that name a file outside its directory (links) or a wildcard. `needs` maps a needed repository path to the text 'referenced by ...' / 'linked by ...' / 'a source of ...'."""
    import xml.etree.ElementTree as ET
    if proj in seen:
        return
    seen.add(proj)
    d = os.path.dirname(proj)
    needs.setdefault(proj, 'the project to build')
    try:
        tree = ET.parse(os.path.join(ROOT, proj)).getroot()
    except (OSError, ET.ParseError):
        return
    for f in tracked:
        if f.startswith(d + '/') and f.endswith(('.cs', '.resx')) and '/bin/' not in f and '/obj/' not in f and not f.startswith(d + '/bin/') and not f.startswith(d + '/obj/'):
            needs.setdefault(f, 'a source of ' + proj)
    for el in tree.iter():
        tag = el.tag.split('}')[-1]
        inc = el.get('Include')
        if not inc or '$(' in inc:
            continue
        if tag == 'ProjectReference':
            target = _norm(d, inc)
            needs.setdefault(target, 'referenced by ' + proj)
            project_needs(tracked, target, seen, needs)
        elif tag in ('Compile', 'EmbeddedResource', 'None', 'Content'):
            pat = _norm(d, inc)
            outside = not pat.startswith(d + '/')
            if outside or el.get('Link'):
                if '*' in pat or '?' in pat:
                    rx = _wild(pat)
                    hits = [f for f in tracked if rx.match(f)]
                    for f in hits:
                        needs.setdefault(f, 'linked by ' + proj)
                    if not hits:
                        needs.setdefault(pat, 'linked by ' + proj + ' (a wildcard that matches no tracked file)')
                else:
                    needs.setdefault(pat, 'linked by ' + proj)


def copy_sources(dockerfile):
    """The (line, source) pairs of the COPY instructions of a Dockerfile that copy from the build context: not `COPY --from=...`, and not the `.` that sends the whole context."""
    text = open(dockerfile).read().replace('\\\n', ' ')
    out = []
    for line in text.splitlines():
        parts = line.strip().split()
        if len(parts) < 3 or parts[0].upper() != 'COPY':
            continue
        flags = [x for x in parts[1:] if x.startswith('--')]
        args = [x for x in parts[1:] if not x.startswith('--')]
        if any(x.startswith('--from') for x in flags) or len(args) < 2:
            continue
        for src in args[:-1]:
            if src not in ('.', './'):
                out.append((line.strip(), src[2:] if src.startswith('./') else src))
    return out


def closure_missing(tracked, kept, roots, dockerfile):
    """What the image build needs and the context does not have: a list of (path or pattern, why)."""
    keptset = set(kept)
    needs = {}
    seen = set()
    tracked_set = set(tracked)
    for r in roots:
        if r not in tracked_set:
            needs.setdefault(r, 'a project to build that is not tracked')
        project_needs(tracked, r, seen, needs)
    missing = [(path, why) for path, why in sorted(needs.items()) if path in tracked_set and path not in keptset]
    missing += [(path, why) for path, why in sorted(needs.items()) if path not in tracked_set]
    if dockerfile and os.path.isfile(dockerfile):
        for line, src in copy_sources(dockerfile):
            if '*' in src or '?' in src:
                rx = _wild(src)
                if not any(rx.match(f) for f in kept):
                    missing.append((src, 'COPY source of "%s": no file of the context matches it' % line[:90]))
            elif src in keptset or any(f.startswith(src.rstrip('/') + '/') for f in kept):
                continue
            else:
                missing.append((src, 'COPY source of "%s": not in the context' % line[:90]))
    return missing


PKG = os.path.join(ROOT, OPT['package-dir'])
pats = load_patterns(OPT['ignore-file'] or os.path.join(PKG, 'Dockerfile.dockerignore'))
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
roots = list(OPT['publish']) + (['tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj', 'tests/SuffixTree/SuffixTree.EmProto.Tests/SuffixTree.EmProto.Tests.csproj',
                                 'tests/SuffixTree/SuffixTree.Persistent.Tests/SuffixTree.Persistent.Tests.csproj'] if OPT['allow-tests'] else [])
missing = closure_missing(files, kept, roots, OPT['dockerfile'] or os.path.join(PKG, 'Dockerfile'))
if missing:
    print('context closure: MISSING from the context (the allow-list does not let in what the image build needs): %d paths' % len(missing))
    sources = {}
    for path, why in missing:
        if why.startswith('a source of '):
            sources.setdefault(why[len('a source of '):], []).append(path)
        else:
            print('   MISSING %s   (%s)' % (path, why))
    for proj, paths in sorted(sources.items()):
        for path in paths[:3]:
            print('   MISSING %s   (a source of %s)' % (path, proj))
        if len(paths) > 3:
            print('   ... and %d more sources of %s (%d in all)' % (len(paths) - 3, proj, len(paths)))
else:
    print('context closure: OK (the projects that the image builds, their references and links, and every COPY source of the Dockerfile are in the context)')
if OPT['static-only']:
    sys.exit(0 if not missing else 1)
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
sys.exit(0 if ok and not bad and not missing else 1)
