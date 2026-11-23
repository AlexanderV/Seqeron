#!/usr/bin/env python3
"""Tests of the static closure check of check_docker_context.py (work package 1b.9h, item 5; the open point of 1b.9f): an image build context is an ALLOW-LIST, and a path that the build needs and the list does not let
in fails the image build at its first use (CS0234 for a project that was not sent) and, in the checker, used to fail only after a whole `dotnet publish`. The static check needs no Docker and no dotnet: it applies the allow-list to
the tracked files and requires that every project reference, every linked source and every COPY source of the Dockerfile that the image build uses (`COPY --from` and `COPY . .` aside) is in the context that results.
What is pinned here: the two real allow-lists (E8, hand-off package) pass; the E8 list without its `SuffixTree.Persistent` line, the hand-off list without its `src/SuffixTree/Algorithms` line, and a list without the
directory of the scripts that a COPY names, each FAIL and name what is missing; a COPY of a file that is not tracked and a wildcard that matches nothing are named, too.
usage: python3 -I test_docker_context.py   (needs the git checkout; skips itself outside one, e.g. inside an image)"""
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
CHECK = os.path.join(HERE, 'check_docker_context.py')
ROOT = os.path.abspath(os.path.join(HERE, '../../../../..'))
E8 = 'docs/refactoring/suffix-tree-upgrade/analysis/e8'
S1B = 'docs/refactoring/suffix-tree-upgrade/analysis/s1b-handoff'
PERSISTENT = '!src/SuffixTree/Algorithms/SuffixTree.Persistent'
S1B_ARGS = ['--package-dir', S1B, '--allow-tests', '--allow-tsv', '--publish', 'apps/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.csproj', '--publish', 'apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj']


def in_checkout():
    r = subprocess.run(['git', '-C', ROOT, 'rev-parse', '--is-inside-work-tree'], capture_output=True, text=True)
    return r.returncode == 0 and r.stdout.strip() == 'true' and os.path.isfile(os.path.join(ROOT, E8, 'Dockerfile.dockerignore'))


def static_check(*args, ignore_lines=None):
    """Runs the checker in its static mode; with ignore_lines the allow-list is that text instead of the package's own file."""
    cmd = [sys.executable, '-I', CHECK, ROOT, '--static-only', *args]
    tmp = None
    try:
        if ignore_lines is not None:
            tmp = tempfile.NamedTemporaryFile('w', suffix='.dockerignore', delete=False)
            tmp.write('\n'.join(ignore_lines) + '\n')
            tmp.close()
            cmd += ['--ignore-file', tmp.name]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return p.returncode, p.stdout + p.stderr
    finally:
        if tmp is not None:
            os.unlink(tmp.name)


def allow_list(package_dir):
    with open(os.path.join(ROOT, package_dir, 'Dockerfile.dockerignore')) as f:
        return f.read().splitlines()


@unittest.skipUnless(in_checkout(), 'needs the git checkout of the repository')
class StaticClosure(unittest.TestCase):
    def test_the_e8_allow_list_holds_everything_the_image_build_needs(self):
        rc, out = static_check()
        self.assertEqual(rc, 0, out)
        self.assertIn('context closure: OK', out)
        self.assertNotIn('MISSING', out)

    def test_the_e8_allow_list_has_the_persistent_line(self):
        self.assertIn(PERSISTENT, [x.strip() for x in allow_list(E8)])

    def test_the_e8_allow_list_without_the_persistent_line_fails_and_names_the_project(self):
        lines = [x for x in allow_list(E8) if x.strip() != PERSISTENT]
        self.assertEqual(len(lines), len(allow_list(E8)) - 1, 'the line occurs once')
        rc, out = static_check(ignore_lines=lines)
        self.assertEqual(rc, 1, out)
        self.assertIn('MISSING', out)
        self.assertIn('src/SuffixTree/Algorithms/SuffixTree.Persistent/SuffixTree.Persistent.csproj', out)
        self.assertIn('referenced by src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj', out)
        self.assertRegex(out, r'src/SuffixTree/Algorithms/SuffixTree\.Persistent/\S+\.cs')   # its sources are named too, not only its project file

    def test_the_missing_list_of_the_e8_context_without_the_line_is_the_persistent_project_alone(self):
        lines = [x for x in allow_list(E8) if x.strip() != PERSISTENT]
        _, out = static_check(ignore_lines=lines)
        missing = [m.group(1) for m in re.finditer(r'^\s+MISSING (\S+)', out, re.M)]
        self.assertTrue(missing, out)
        self.assertTrue(all(m.startswith('src/SuffixTree/Algorithms/SuffixTree.Persistent/') for m in missing), missing)

    def test_the_handoff_allow_list_holds_everything_the_image_build_needs(self):
        rc, out = static_check(*S1B_ARGS)
        self.assertEqual(rc, 0, out)
        self.assertIn('context closure: OK', out)

    def test_the_handoff_allow_list_without_the_algorithms_line_fails(self):
        lines = [x for x in allow_list(S1B) if x.strip() != '!src/SuffixTree/Algorithms']
        self.assertEqual(len(lines), len(allow_list(S1B)) - 1)
        rc, out = static_check(*S1B_ARGS, ignore_lines=lines)
        self.assertEqual(rc, 1, out)
        self.assertIn('src/SuffixTree/Algorithms/SuffixTree.Esa/SuffixTree.Esa.csproj', out)
        self.assertIn('src/SuffixTree/Algorithms/SuffixTree.Persistent/SuffixTree.Persistent.csproj', out)

    def test_a_copy_source_that_the_list_does_not_let_in_fails(self):
        # the E8 Dockerfile copies the scripts of its own directory and the substitute props: without the directory line they are not in the context
        lines = [x for x in allow_list(E8) if x.strip() != '!' + E8]
        self.assertEqual(len(lines), len(allow_list(E8)) - 1)
        rc, out = static_check(ignore_lines=lines)
        self.assertEqual(rc, 1, out)
        self.assertIn('COPY', out)
        self.assertIn(E8 + '/docker/Directory.Build.props', out)
        self.assertIn(E8 + '/*.sh', out)

    def test_a_linked_source_outside_the_project_directory_is_a_needed_path(self):
        # EmProto links four files of the Bench project: without the lines the publish would fail at the first of them
        lines = [x for x in allow_list(E8) if x.strip() != '!apps/SuffixTree.Esa.Bench/Links']
        self.assertEqual(len(lines), len(allow_list(E8)) - 1)
        rc, out = static_check(ignore_lines=lines)
        self.assertEqual(rc, 1, out)
        self.assertIn('apps/SuffixTree.Esa.Bench/Links/LinkTables.cs', out)
        self.assertIn('linked by apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj', out)

    def test_a_build_output_directory_is_not_a_needed_path(self):
        # the closure is over tracked files: the bin and obj rules of the list must not make the check fail
        rc, out = static_check()
        self.assertEqual(rc, 0, out)
        self.assertNotIn('/bin/', out)
        self.assertNotIn('/obj/', out)


if __name__ == '__main__':
    unittest.main(verbosity=1)
