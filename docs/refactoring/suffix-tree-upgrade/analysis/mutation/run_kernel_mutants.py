#!/usr/bin/env python3
"""Single-line mutants of an optimized kernel file (acceptance artefact of every kernel iteration, REVIEW-S1a section 6): each mutant is built and tested in the Release (raw pointer) configuration and, if it survives, in the checked configuration (-p:EsaBounds=true).
   usage: run_kernel_mutants.py <kernel .cs path relative to the repo root> <mutants.json> <result.tsv> [nunit filter, default FullyQualifiedName~Kernels]
   The source file is restored afterwards (a try/finally; check `git diff` anyway)."""
import json, os, subprocess, sys, shutil
ROOT = os.environ.get('ESA_ROOT', '/home/user/Seqeron')  # a git worktree can be mutated while the main tree stays editable
BOTH = os.environ.get('MUTANT_BOTH') == '1'  # run the checked configuration for every mutant, not only for the Release survivors
RELEASE_ONLY = os.environ.get('MUTANT_RELEASE_ONLY') == '1'  # never run the checked configuration (long exhaustive filters: the equivalence evidence needs one configuration)
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1')
kernel, mutants_file, out = sys.argv[1], sys.argv[2], sys.argv[3]
filt = sys.argv[4] if len(sys.argv) > 4 else 'FullyQualifiedName~Kernels'
path = os.path.join(ROOT, kernel)
orig = open(path).read()
mutants = json.load(open(mutants_file))
proj = ROOT + '/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj'
def run(bounds):
    e = dict(env)
    args = ['dotnet', 'test', proj, '-c', 'Release', '--nologo', '-p:EsaBounds=' + ('true' if bounds else 'false'), '--filter', filt]
    if bounds: e['ESA_EXPECT_BOUNDS'] = 'true'
    r = subprocess.run(args, capture_output=True, text=True, env=e)
    txt = r.stdout + r.stderr
    if 'error CS' in txt or 'error IDE' in txt or 'error CA' in txt or 'error S' in txt:
        return 'BUILD-ERROR'
    if 'Test host process crashed' in txt or 'test run was aborted' in txt: return 'killed(host crash)'  # a mutant that corrupts memory may take the test host down
    if 'Failed!' in txt: pass
    elif 'Passed!' in txt: return 'survived'
    else: return 'NO-RESULT'
    import re
    m = re.search(r'Failed!.*?Failed:\s+(\d+)', txt)
    return 'killed(%s)' % (m.group(1) if m else '?')
results = []
try:
    for m in mutants:
        src = orig
        pairs = [(m['old'], m['new'])] + [tuple(x) for x in m.get('also', [])]  # 'also': further replacements of the same mutant
        if any(src.count(o) != 1 for o, _ in pairs):
            results.append((m['id'], m['desc'], 'PATTERN-NOT-UNIQUE', '')); continue
        for o, nw in pairs: src = src.replace(o, nw)
        open(path, 'w').write(src)
        r1 = run(False)
        r2 = run(True) if ((r1 == 'survived' or BOTH) and not RELEASE_ONLY) else ''
        results.append((m['id'], m['desc'], r1, r2))
        open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
finally:
    open(path, 'w').write(orig)
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', kernel], env=env)
