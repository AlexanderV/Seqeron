#!/usr/bin/env python3
"""The part of the runners of work package 1b.9c that finds the test classes that can see a source file by name (shared by run_1b9c_mutants.py and rerun_1b9c_hang.py)."""
import glob, os, re

GENERIC = {'Options', 'Result', 'State', 'Tuning', 'Kind', 'Mode', 'Context', 'Stats', 'Task'}
POPULAR = 20  # a helper type that more than this many test files mention (TempDir, ...) says nothing about the code under test


def project_files(testdir):
    """Every source of the Esa test project: (text, declared types, public classes, is a default test file). A file that holds an [Explicit] class is left out of the test classes (the heavy classes run in stage H)."""
    out = []
    for f in sorted(glob.glob(os.path.join(testdir, '**', '*.cs'), recursive=True)):
        if '/obj/' in f or '/bin/' in f:
            continue
        t = open(f).read()
        types = set(re.findall(r'\b(?:class|struct|enum|interface|record)\s+([A-Z]\w+)', t))
        classes = re.findall(r'\bpublic (?:sealed |static |partial )*class (\w+)', t)
        out.append((t, types, classes, '[Test' in t and '[Explicit' not in t))
    return out


def seeing_classes(src_path, testdir):
    """The test classes (not [Explicit]) that can see the mutated file by name: those whose source names a type that the file declares, and those that name a type of a helper file (a file without tests: SesaSupport, SaSupport ...) that
    names such a type, unless that helper type is used everywhere (TempDir). Stage 2 (the whole default set) covers what no name shows."""
    types = set(re.findall(r'\b(?:class|struct|enum|interface|record)\s+([A-Z]\w+)', open(src_path).read())) - GENERIC
    files = project_files(testdir)

    def mentions(text, ts):
        return any(re.search(r'\b%s\b' % re.escape(x), text) for x in ts)

    helper_types = set()
    for t, declared, _, is_test in files:
        if not is_test and '[Explicit' not in t and mentions(t, types):
            helper_types |= declared
    test_files = [x for x in files if x[3]]
    helper_types = set(x for x in helper_types if sum(1 for tf in test_files if re.search(r'\b%s\b' % re.escape(x), tf[0])) <= POPULAR) - types
    found = set()
    for t, _, classes, is_test in files:
        if is_test and mentions(t, types | helper_types):
            found.update(classes)
    return found
