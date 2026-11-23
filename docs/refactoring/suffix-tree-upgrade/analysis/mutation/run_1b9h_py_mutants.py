#!/usr/bin/env python3
"""Mutants of the static closure check of the image build context (work package 1b.9h, item 5): single-site mutants of
docs/refactoring/suffix-tree-upgrade/analysis/e8/check_docker_context.py, each run against its own tests (test_docker_context.py, 9 tests; the checker and the tests run in a separate git worktree, never in the working checkout).
The list is inside this file (Python needs no build): a mutant is a pair (old text, new text) with `old` occurring exactly once in the checker at the head. A mutant is KILLED when the tests fail (a failed assertion) or the checker no
longer runs (a crash, which the tests see as a failed assertion), SURVIVED when all 9 tests pass; a run that does not end in 300 s is a HANG (counted as killed, said so). The checker is restored after every mutant.
usage: run_1b9h_py_mutants.py <worktree> <result.tsv>      (the result is appended to; mutants it already holds are skipped)"""
import os
import re
import subprocess
import sys

ROOT, OUT = sys.argv[1], sys.argv[2]
E8 = os.path.join(ROOT, 'docs/refactoring/suffix-tree-upgrade/analysis/e8')
CHECK = os.path.join(E8, 'check_docker_context.py')
TESTS = os.path.join(E8, 'test_docker_context.py')

M = []


def m(i, desc, old, new):
    M.append((i, desc, old, new))


m('P01', "glob: '**/' no longer matches zero directories", "out += '(?:.*/)?'; i += 3; continue", "out += '(?:.*/)'; i += 3; continue")
m('P02', "glob: '**' stops at a slash", "out += '.*'; i += 2; continue", "out += '[^/]*'; i += 2; continue")
m('P03', "glob: '*' crosses a slash", "{'*': '[^/]*', '?': '[^/]'}", "{'*': '.*', '?': '[^/]'}")
m('P04', "glob: '?' crosses a slash", "'?': '[^/]'}", "'?': '.'}")
m('P05', "allow-list: '!' is not a negation", "neg = line.startswith('!')", "neg = False")
m('P06', "allow-list: a comment line is a pattern", "if not line or line.startswith('#'):", "if not line:")
m('P07', "allow-list: the patterns are read with a leading slash", ".replace(os.sep, '/').lstrip('/')", ".replace(os.sep, '/')")
m('P08', "excluded: a file below an excluded directory is not excluded", "if rx.match(path) or any(rx.match(p) for p in parents):", "if rx.match(path):")
m('P09', "excluded: a '!' line excludes", "state = not neg", "state = True")
m('P10', "excluded: a '!' line is not looked at", "state = not neg", "state = state or not neg")
m('P11', "excluded: the first matching pattern wins", "            state = not neg\n    return state", "            state = not neg\n            break\n    return state")
m('P12', "project paths: a backslash is not a separator", "os.path.join(base_dir, rel.replace('\\\\', '/'))", "os.path.join(base_dir, rel)")
m('P13', "closure: a project reference is not followed", "            project_needs(tracked, target, seen, needs)\n        elif tag in", "            pass\n        elif tag in")
m('P14', "closure: the project referenced is not needed itself", "            needs.setdefault(target, 'referenced by ' + proj)\n", "")
m('P15', "closure: the sources of a project are not needed", "if f.startswith(d + '/') and f.endswith(('.cs', '.resx'))", "if False and f.endswith(('.cs', '.resx'))")
m('P16', "closure: the project file is not needed", "    needs.setdefault(proj, 'the project to build')\n", "")
m('P17', "closure: build output is a needed source (bin)", "and '/bin/' not in f and '/obj/' not in f and not f.startswith(d + '/bin/') and not f.startswith(d + '/obj/'):", "and '/obj/' not in f and not f.startswith(d + '/obj/'):")
m('P18', "closure: build output is a needed source (obj)", "and '/bin/' not in f and '/obj/' not in f and not f.startswith(d + '/bin/') and not f.startswith(d + '/obj/'):", "and '/bin/' not in f and not f.startswith(d + '/bin/'):")
m('P19', "closure: a link inside the project directory is a needed path", "outside = not pat.startswith(d + '/')", "outside = True")
m('P20', "closure: a link outside the project directory is not looked at", "outside = not pat.startswith(d + '/')", "outside = False")
m('P21', "closure: a linked wildcard is taken literally", "if '*' in pat or '?' in pat:\n                    rx = _wild(pat)", "if False:\n                    rx = _wild(pat)")
m('P22', "closure: a wildcard that matches no tracked file is not named", "                    if not hits:\n                        needs.setdefault(pat, 'linked by ' + proj + ' (a wildcard that matches no tracked file)')\n", "")
m('P23', "closure: an item with a variable in its path is followed", "if not inc or '$(' in inc:", "if not inc:")
m('P24', "closure: a project already walked is walked again (no memo)", "    if proj in seen:\n        return\n    seen.add(proj)\n", "    seen.add(proj)\n")
m('P25', "closure: only Compile items are links", "elif tag in ('Compile', 'EmbeddedResource', 'None', 'Content'):", "elif tag in ('Compile',):")
m('P26', "closure: a project that is not tracked is not named", "        if r not in tracked_set:\n            needs.setdefault(r, 'a project to build that is not tracked')\n", "")
m('P27', "closure: a needed path that the list lets in is missing", "if path in tracked_set and path not in keptset]", "if path in tracked_set and path in keptset]")
m('P28', "closure: a needed path that the list drops is not missing", "if path in tracked_set and path not in keptset]", "if False]")
m('P29', "closure: a needed path that is not tracked is not named", "    missing += [(path, why) for path, why in sorted(needs.items()) if path not in tracked_set]\n", "")
m('P30', "COPY: a source taken from another stage is checked", "if any(x.startswith('--from') for x in flags) or len(args) < 2:", "if len(args) < 2:")
m('P31', "COPY: the whole context is a source", "if src not in ('.', './'):", "if True:")
m('P32', "COPY: a source with a leading ./ keeps it", "src[2:] if src.startswith('./') else src", "src")
m('P33', "COPY: only the first source of a line is checked", "for src in args[:-1]:", "for src in args[:1]:")
m('P34', "COPY: the last argument (the destination) is a source", "for src in args[:-1]:", "for src in args:")
m('P35', "COPY: a continued line is read as two", ".read().replace('\\\\\\n', ' ')", ".read()")
m('P36', "COPY: a directory source that has files in the context is missing", "elif src in keptset or any(f.startswith(src.rstrip('/') + '/') for f in kept):", "elif src in keptset:")
m('P37', "COPY: a source that is not in the context is not named", "                missing.append((src, 'COPY source of \"%s\": not in the context' % line[:90]))\n", "                pass\n")
m('P38', "COPY: a wildcard source is looked for among the tracked files, not the context", "if not any(rx.match(f) for f in kept):", "if not any(rx.match(f) for f in tracked):")
m('P39', "COPY: a wildcard that matches nothing is not named", "                    missing.append((src, 'COPY source of \"%s\": no file of the context matches it' % line[:90]))\n", "                    pass\n")
m('P40', "COPY: the Dockerfile is not read", "if dockerfile and os.path.isfile(dockerfile):", "if False:")
m('P41', "exit code: a missing path does not fail the static check", "    sys.exit(0 if not missing else 1)", "    sys.exit(0)")
m('P42', "exit code: a clean closure fails the static check", "    sys.exit(0 if not missing else 1)", "    sys.exit(1)")
m('P43', "options: --static-only is not an option", "if a in ('--allow-tests', '--allow-tsv', '--static-only'):", "if a in ('--allow-tests', '--allow-tsv'):")
m('P44', "options: --ignore-file is ignored", "pats = load_patterns(OPT['ignore-file'] or os.path.join(PKG, 'Dockerfile.dockerignore'))", "pats = load_patterns(os.path.join(PKG, 'Dockerfile.dockerignore'))")
m('P45', "options: --dockerfile is ignored", "missing = closure_missing(files, kept, roots, OPT['dockerfile'] or os.path.join(PKG, 'Dockerfile'))", "missing = closure_missing(files, kept, roots, os.path.join(PKG, 'Dockerfile'))")
m('P46', "roots: the test projects are not roots with --allow-tests", "if OPT['allow-tests'] else [])\nmissing", "if False else [])\nmissing")
m('P47', "roots: the projects to publish are not roots", "roots = list(OPT['publish']) + (", "roots = [] + (")

done = set()
if os.path.exists(OUT):
    done = {line.split('\t')[0] for line in open(OUT) if line.strip()}
orig = open(CHECK).read()
for i, desc, old, new in M:
    n = orig.count(old)
    if n != 1:
        raise SystemExit('%s: the pattern occurs %d times: %r' % (i, n, old[:80]))
    if old == new:
        raise SystemExit('%s: not a change' % i)
base = subprocess.run([sys.executable, '-I', TESTS], capture_output=True, text=True, timeout=600)
if base.returncode != 0:
    raise SystemExit('the tests fail at the head:\n' + base.stderr[-1500:])
try:
    for i, desc, old, new in M:
        if i in done:
            continue
        open(CHECK, 'w').write(orig.replace(old, new))
        try:
            p = subprocess.run([sys.executable, '-I', TESTS], capture_output=True, text=True, timeout=300)
            txt = p.stderr + p.stdout
            if p.returncode == 0:
                verdict = 'survived'
            else:
                names = re.findall(r'^(?:FAIL|ERROR): (\S+)', txt, re.M)
                verdict = 'killed(%d: %s)' % (len(names), ','.join(n[:60] for n in names[:4]))
        except subprocess.TimeoutExpired:
            verdict = 'killed(hang)'
        with open(OUT, 'a') as f:
            f.write('%s\t%s\t%s\n' % (i, desc, verdict))
            f.flush()
            os.fsync(f.fileno())
        print(i, verdict, flush=True)
finally:
    open(CHECK, 'w').write(orig)
