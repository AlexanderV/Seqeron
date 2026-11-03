#!/usr/bin/env python3
"""Single-line mutants of an optimized kernel file (acceptance artefact of every kernel iteration, REVIEW-S1a section 6): each mutant is built and tested in the Release (raw pointer) configuration and, if it survives, in the checked configuration (-p:EsaBounds=true).
   usage: run_kernel_mutants.py <kernel .cs path relative to the repo root> <mutants.json> <result.tsv> [nunit filter, default FullyQualifiedName~Kernels]
   The source file is restored afterwards (a try/finally; check `git diff` anyway)."""
import json, os, subprocess, sys, shutil
ROOT = '/home/user/Seqeron'
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
    if 'Passed!' in txt: return 'survived'
    import re
    m = re.search(r'Failed!.*?Failed:\s+(\d+)', txt)
    return 'killed(%s)' % (m.group(1) if m else '?')
results = []
try:
    for m in mutants:
        src = orig
        if src.count(m['old']) != 1:
            results.append((m['id'], m['desc'], 'PATTERN-NOT-UNIQUE', '')); continue
        open(path, 'w').write(src.replace(m['old'], m['new']))
        r1 = run(False)
        r2 = run(True) if r1 == 'survived' else ''
        results.append((m['id'], m['desc'], r1, r2))
        open(out, 'w').write('\n'.join('\t'.join(map(str, r)) for r in results) + '\n')
finally:
    open(path, 'w').write(orig)
    subprocess.run(['git', '-C', ROOT, 'diff', '--stat', '--', kernel], env=env)
