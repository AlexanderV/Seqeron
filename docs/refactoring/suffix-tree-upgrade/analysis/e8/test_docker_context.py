#!/usr/bin/env python3
"""Tests of the static closure check of check_docker_context.py (work package 1b.9h, item 5; the open point of 1b.9f): an image build context is an ALLOW-LIST, and a path that the build needs and the list does not let
in fails the image build at its first use (CS0234 for a project that was not sent) and, in the checker, used to fail only after a whole `dotnet publish`. The static check needs no Docker and no dotnet: it applies the allow-list to
the tracked files and requires that every project reference, every linked source and every COPY source of the Dockerfile that the image build uses (`COPY --from` and `COPY . .` aside) is in the context that results.
What is pinned here: the two real allow-lists (E8, hand-off package) pass; the E8 list without its `SuffixTree.Persistent` line, the hand-off list without its `src/SuffixTree/Algorithms` line, and a list without the
directory of the scripts that a COPY names, each FAIL and name what is missing; a COPY of a file that is not tracked and a wildcard that matches nothing are named, too.
The second half (`Synthetic*`, 26 tests) builds small git repositories and pins each thing that the walk follows and leaves alone: the project to build, its sources, project references (two steps, an untracked one, a cycle), build output, items with and without a `Link`
and outside the project directory (`Compile`, `EmbeddedResource`, `None`, `Content`), wildcards, properties in a path, the matching rules of the allow-list, every form of a `COPY` source, and the test projects of `--allow-tests`. These need git only.
usage: python3 -I test_docker_context.py   (needs the git checkout; skips itself outside one, e.g. inside an image)"""
import os
import re
import shutil
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
        self.assertRegex(out, r'src/SuffixTree/Algorithms/SuffixTree\.Persistent/\S+\.cs\s')   # its sources are named too, not only its project file

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


def make_repo(files):
    """A small git repository (tracked files only: the closure is over `git ls-files`) in a temporary directory: `files` maps a repository path to its text."""
    root = tempfile.mkdtemp(prefix='ctx1b9h')
    for rel, text in files.items():
        path = os.path.join(root, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            f.write(text)
    subprocess.run(['git', '-C', root, 'init', '-q'], check=True, capture_output=True)
    subprocess.run(['git', '-C', root, 'add', '-f', '-A'], check=True, capture_output=True)
    return root


def proj(*items):
    return '<Project Sdk="Microsoft.NET.Sdk">\n  <ItemGroup>\n' + ''.join('    %s\n' % i for i in items) + '  </ItemGroup>\n</Project>\n'


def closure_of(root, allow, *args, dockerfile=None, publish='apps/App/App.csproj'):
    """The static closure of the synthetic repository under the allow-list `allow` (lines): (exit code, output)."""
    ign = os.path.join(root, '..', os.path.basename(root) + '.dockerignore')
    with open(ign, 'w') as f:
        f.write('\n'.join(allow) + '\n')
    cmd = [sys.executable, '-I', CHECK, root, '--static-only', '--ignore-file', ign, '--publish', publish, *args]
    if dockerfile is not None:
        df = os.path.join(root, '..', os.path.basename(root) + '.Dockerfile')
        with open(df, 'w') as f:
            f.write(dockerfile)
        cmd += ['--dockerfile', df]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return p.returncode, p.stdout + p.stderr
    finally:
        os.unlink(ign)
        if dockerfile is not None:
            os.unlink(df)


def missing_of(out):
    return [m.group(1) for m in re.finditer(r'^\s+MISSING (\S+)', out, re.M)]


BASE = {
    'apps/App/App.csproj': proj(r'<ProjectReference Include="..\..\src\Lib\Lib.csproj" />'),
    'apps/App/Program.cs': 'class P {}\n',
    'src/Lib/Lib.csproj': proj(r'<ProjectReference Include="..\Other\Other.csproj" />'),
    'src/Lib/Lib.cs': 'class L {}\n',
    'src/Other/Other.csproj': proj(),
    'src/Other/Other.cs': 'class O {}\n',
}
ALL = ['*', '!apps/App', '!src/Lib', '!src/Other']


class SyntheticClosure(unittest.TestCase):
    """The closure on small repositories that these tests build: every kind of item that the walk follows, and what it must leave alone. Needs git only."""

    def check(self, files, allow, *args, **kw):
        root = make_repo({**BASE, **files})
        try:
            return closure_of(root, allow, *args, **kw)
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_a_list_that_lets_everything_in_passes(self):
        rc, out = self.check({}, ALL)
        self.assertEqual(rc, 0, out)
        self.assertIn('context closure: OK', out)

    def test_the_project_to_build_is_needed_itself(self):
        rc, out = self.check({}, ['*', '!apps/App/Program.cs', '!src/Lib', '!src/Other'])
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['apps/App/App.csproj'])
        self.assertIn('the project to build', out)

    def test_the_sources_of_a_project_are_needed(self):
        rc, out = self.check({}, ['*', '!apps/App', '!src/Lib/Lib.csproj', '!src/Other'])
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['src/Lib/Lib.cs'])
        self.assertIn('a source of src/Lib/Lib.csproj', out)

    def test_a_referenced_project_is_needed_and_so_are_the_references_of_that_project(self):
        rc, out = self.check({}, ['*', '!apps/App', '!src/Lib'])          # Other is referenced by Lib, which App references: two steps
        self.assertEqual(rc, 1, out)
        self.assertEqual(sorted(missing_of(out)), ['src/Other/Other.cs', 'src/Other/Other.csproj'])
        self.assertIn('referenced by src/Lib/Lib.csproj', out)

    def test_a_reference_to_a_project_that_is_not_tracked_is_named(self):
        rc, out = self.check({'src/Other/Other.csproj': proj(r'<ProjectReference Include="..\Ghost\Ghost.csproj" />')}, ALL)
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['src/Ghost/Ghost.csproj'])
        self.assertIn('referenced by src/Other/Other.csproj', out)

    def test_a_project_to_build_that_is_not_tracked_is_named(self):
        rc, out = self.check({}, ALL, publish='apps/Nope/Nope.csproj')
        self.assertEqual(rc, 1, out)
        self.assertIn('a project to build that is not tracked', out)
        self.assertEqual(missing_of(out), ['apps/Nope/Nope.csproj'])

    def test_projects_that_reference_each_other_end_the_walk(self):
        rc, out = self.check({'src/Other/Other.csproj': proj(r'<ProjectReference Include="..\..\apps\App\App.csproj" />')}, ALL)
        self.assertEqual(rc, 0, out)
        self.assertNotIn('Traceback', out)

    def test_build_output_below_a_project_is_not_a_needed_source(self):
        files = {'src/Lib/bin/Debug/Gen.cs': 'x', 'src/Lib/obj/Gen2.cs': 'x', 'apps/App/obj/Gen3.cs': 'x', 'apps/App/bin/Gen4.cs': 'x'}
        rc, out = self.check(files, ALL + ['**/bin', '**/obj'])
        self.assertEqual(rc, 0, out)
        self.assertNotIn('MISSING', out)

    def test_a_resource_file_of_a_project_is_a_needed_source(self):
        rc, out = self.check({'src/Lib/Strings.resx': '<root/>'}, ALL[:2] + ['!src/Lib/Lib.csproj', '!src/Lib/Lib.cs', '!src/Other'])
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['src/Lib/Strings.resx'])

    def test_build_output_in_a_directory_below_the_project_directory_is_not_a_needed_source(self):
        files = {'src/Lib/sub/bin/Gen.cs': 'x', 'src/Lib/sub/obj/Gen2.cs': 'x'}
        rc, out = self.check(files, ALL + ['**/bin', '**/obj'])
        self.assertEqual(rc, 0, out)
        self.assertNotIn('MISSING', out)

    def test_a_question_mark_in_a_linked_path_matches_one_character(self):
        files = {'shared/A.cs': 'x', 'shared/B.cs': 'x', 'shared/CC.cs': 'x', 'apps/App/App.csproj': proj(r'<ProjectReference Include="..\..\src\Lib\Lib.csproj" />', r'<Compile Include="..\..\shared\?.cs" />')}
        rc, out = self.check(files, ALL + ['!shared/A.cs'])
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['shared/B.cs'], 'B.cs matches ?.cs and is not in the list; CC.cs does not match; a ? that was read as a character would name shared/?.cs')

    def test_an_item_inside_the_project_directory_without_a_link_is_not_a_needed_path(self):
        rc, out = self.check({'apps/App/notes.txt': 'x', 'apps/App/App.csproj': proj(r'<ProjectReference Include="..\..\src\Lib\Lib.csproj" />', '<None Include="notes.txt" />')},
                             ['*', '!apps/App/App.csproj', '!apps/App/Program.cs', '!src/Lib', '!src/Other'])
        self.assertEqual(rc, 0, out)

    def test_an_item_inside_the_project_directory_with_a_link_is_a_needed_path(self):
        rc, out = self.check({'apps/App/notes.txt': 'x', 'apps/App/App.csproj': proj(r'<ProjectReference Include="..\..\src\Lib\Lib.csproj" />', '<None Include="notes.txt" Link="n.txt" />')},
                             ['*', '!apps/App/App.csproj', '!apps/App/Program.cs', '!src/Lib', '!src/Other'])
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['apps/App/notes.txt'])
        self.assertIn('linked by apps/App/App.csproj', out)

    def test_an_item_outside_the_project_directory_is_a_needed_path_with_or_without_a_link_attribute(self):
        for tag in ('Compile', 'EmbeddedResource', 'None', 'Content'):
            with self.subTest(tag=tag):
                files = {'shared/Shared.txt': 'x', 'apps/App/App.csproj': proj(r'<ProjectReference Include="..\..\src\Lib\Lib.csproj" />', r'<%s Include="..\..\shared\Shared.txt" />' % tag)}
                rc, out = self.check(files, ALL)
                self.assertEqual(rc, 1, out)
                self.assertEqual(missing_of(out), ['shared/Shared.txt'])
                rc, out = self.check(files, ALL + ['!shared'])
                self.assertEqual(rc, 0, out)

    def test_a_wildcard_item_outside_the_project_directory_needs_every_file_it_matches(self):
        files = {'shared/A.cs': 'x', 'shared/B.cs': 'x', 'shared/C.txt': 'x', 'apps/App/App.csproj': proj(r'<ProjectReference Include="..\..\src\Lib\Lib.csproj" />', r'<Compile Include="..\..\shared\*.cs" />')}
        rc, out = self.check(files, ALL + ['!shared/A.cs'])
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['shared/B.cs'])

    def test_a_wildcard_item_that_matches_no_tracked_file_is_named(self):
        files = {'apps/App/App.csproj': proj(r'<ProjectReference Include="..\..\src\Lib\Lib.csproj" />', r'<Compile Include="..\..\shared\*.cs" />')}
        rc, out = self.check(files, ALL)
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['shared/*.cs'])
        self.assertIn('a wildcard that matches no tracked file', out)

    def test_an_item_with_a_property_in_its_path_is_not_followed(self):
        files = {'apps/App/App.csproj': proj(r'<ProjectReference Include="..\..\src\Lib\Lib.csproj" />', r'<Compile Include="$(Shared)\Nope.cs" />', r'<ProjectReference Include="$(Libs)\Nope.csproj" />')}
        rc, out = self.check(files, ALL)
        self.assertEqual(rc, 0, out)

    def test_an_allow_list_with_a_leading_slash_and_globs_means_what_the_docker_client_means(self):
        rc, out = self.check({}, ['*', '!/apps/App', '!/src/Lib', '!/src/Other'])
        self.assertEqual(rc, 0, out)
        # '*' does not cross a slash: 'apps/*.cs' is not 'apps/App/Program.cs' ('?' neither)
        rc, out = self.check({}, ['*', '!apps/App/App.csproj', '!apps/*.cs', '!src/Lib', '!src/Other'])
        self.assertEqual(missing_of(out), ['apps/App/Program.cs'], out)
        rc, out = self.check({}, ['*', '!apps/App/App.csproj', '!apps/App?Program.cs', '!src/Lib', '!src/Other'])
        self.assertEqual(missing_of(out), ['apps/App/Program.cs'], out)
        # '**' crosses slashes: at the end of a pattern and before a slash (also for a top-level name)
        rc, out = self.check({'obj/Top.cs': 'x', 'src/Lib/obj/Deep.cs': 'x'}, ['*', '!apps/App/**', '!src/Lib/**', '!src/Other', '!obj/Top.cs', '**/obj'])
        self.assertEqual(rc, 0, out)
        rc, out = self.check({}, ['*', '!apps/App/**/Program.cs', '!apps/App/App.csproj', '!src/Lib', '!src/Other'])
        self.assertEqual(rc, 0, out)

    def test_a_later_line_of_the_list_wins_and_a_bang_line_lets_a_file_in(self):
        rc, out = self.check({}, ALL + ['apps/App/Program.cs'])
        self.assertEqual(missing_of(out), ['apps/App/Program.cs'], out)
        rc, out = self.check({}, ALL + ['apps/App/Program.cs', '!apps/App/Program.cs'])
        self.assertEqual(rc, 0, out)
        rc, out = self.check({}, ['*', '!apps/App', '!src/Lib', '!src/Other', '# apps', '', '   '])         # comments and blank lines are no patterns
        self.assertEqual(rc, 0, out)


DOCKERFILE_COPY_OK = 'FROM x\nCOPY scripts/run.sh /app/\nCOPY ./scripts/two.sh /app/\nCOPY scripts/a.sh scripts/b.sh /app/\nCOPY --from=build /out /app\nCOPY . .\n'


class SyntheticCopySources(unittest.TestCase):
    """The COPY lines of the Dockerfile that the image build uses: the sources that come from the context must be in it."""
    FILES = {'scripts/run.sh': 'x', 'scripts/two.sh': 'x', 'scripts/a.sh': 'x', 'scripts/b.sh': 'x'}

    def check(self, allow, dockerfile, extra=None):
        root = make_repo({**BASE, **self.FILES, **(extra or {})})
        try:
            return closure_of(root, allow, dockerfile=dockerfile)
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_every_source_that_is_in_the_context_passes_and_the_others_forms_are_not_sources(self):
        rc, out = self.check(ALL + ['!scripts'], DOCKERFILE_COPY_OK)           # './x', several sources, COPY --from and COPY . . are not looked at
        self.assertEqual(rc, 0, out)

    def test_a_source_with_a_dot_slash_is_looked_up_without_it(self):
        rc, out = self.check(ALL + ['!scripts/run.sh', '!scripts/a.sh', '!scripts/b.sh'], DOCKERFILE_COPY_OK)
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['scripts/two.sh'])

    def test_each_source_of_a_line_is_looked_at_and_the_destination_is_not(self):
        rc, out = self.check(ALL + ['!scripts/run.sh', '!scripts/two.sh', '!scripts/a.sh'], DOCKERFILE_COPY_OK)
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['scripts/b.sh'])
        self.assertIn('COPY source of "COPY scripts/a.sh scripts/b.sh /app/"', out)

    def test_a_line_that_is_continued_is_one_instruction(self):
        rc, out = self.check(ALL + ['!scripts/run.sh'], 'FROM x\nCOPY scripts/run.sh \\\n     scripts/a.sh \\\n     /app/\n')
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['scripts/a.sh'])

    def test_a_directory_source_is_there_when_a_file_below_it_is(self):
        rc, out = self.check(ALL + ['!scripts/a.sh'], 'FROM x\nCOPY scripts /app/\n')
        self.assertEqual(rc, 0, out)
        rc, out = self.check(ALL, 'FROM x\nCOPY scripts /app/\n')
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['scripts'])

    def test_a_wildcard_source_is_matched_against_the_context_not_against_the_tracked_files(self):
        rc, out = self.check(ALL + ['!scripts/a.sh'], 'FROM x\nCOPY scripts/*.sh /app/\n')
        self.assertEqual(rc, 0, out)
        rc, out = self.check(ALL, 'FROM x\nCOPY scripts/*.sh /app/\n')            # the files exist in the repository, the list does not let them in
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['scripts/*.sh'])
        self.assertIn('no file of the context matches it', out)

    def test_a_source_that_is_not_tracked_is_not_in_the_context(self):
        rc, out = self.check(ALL + ['!scripts'], 'FROM x\nCOPY scripts/ghost.sh /app/\n')
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['scripts/ghost.sh'])
        self.assertIn('not in the context', out)

    def test_the_dockerfile_that_is_given_is_the_one_that_is_read(self):
        # the synthetic repository has no package directory with a Dockerfile of its own: only the one given by --dockerfile can make this fail
        rc, out = self.check(ALL, 'FROM x\nCOPY scripts/run.sh /app/\n')
        self.assertEqual(rc, 1, out)


class SyntheticTestProjects(unittest.TestCase):
    """--allow-tests: the three test projects are roots of the closure, as the projects to publish are."""
    TESTS = {'tests/SuffixTree/%s/%s.csproj' % (n, n): proj(r'<ProjectReference Include="..\..\..\src\Lib\Lib.csproj" />') for n in ('SuffixTree.Esa.Tests', 'SuffixTree.EmProto.Tests', 'SuffixTree.Persistent.Tests')}

    def check(self, allow, *args):
        root = make_repo({**BASE, **self.TESTS})
        try:
            return closure_of(root, allow, *args)
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_a_test_project_that_the_list_leaves_out_is_missing_with_the_option_and_not_without(self):
        allow = ALL + ['!tests/SuffixTree/SuffixTree.EmProto.Tests', '!tests/SuffixTree/SuffixTree.Persistent.Tests']
        rc, out = self.check(allow, '--allow-tests')
        self.assertEqual(rc, 1, out)
        self.assertEqual(missing_of(out), ['tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj'])
        rc, out = self.check(allow)
        self.assertEqual(rc, 0, out)

    def test_all_three_test_projects_in_the_list_pass(self):
        rc, out = self.check(ALL + ['!tests/SuffixTree'], '--allow-tests')
        self.assertEqual(rc, 0, out)


@unittest.skipUnless(in_checkout() and shutil.which('dotnet'), 'needs the git checkout of the repository and the dotnet SDK (the full check publishes the program with an empty package cache; about 10 s)')
class FullCheck(unittest.TestCase):
    """The checker without --static-only: its exit code carries the closure as well, not only the copy and the publish. A COPY source that the context does not hold passes the publish (it is not read by it) and
    must still fail the check; the unit tests above stop at the static mode, which has an exit code of its own."""

    def test_a_copy_source_that_is_not_in_the_context_fails_the_full_check_although_the_publish_passes(self):
        dockerfile = os.path.join(ROOT, E8, 'Dockerfile')
        with open(dockerfile) as f:
            text = f.read()
        tmp = tempfile.NamedTemporaryFile('w', suffix='.Dockerfile', delete=False)
        try:
            tmp.write(text + '\nCOPY a-file-that-the-context-does-not-hold.txt /x\n')
            tmp.close()
            p = subprocess.run([sys.executable, '-I', CHECK, ROOT, '--dockerfile', tmp.name], capture_output=True, text=True, timeout=600)
        finally:
            os.unlink(tmp.name)
        out = p.stdout + p.stderr
        self.assertIn('publish apps/SuffixTree.EmProto/SuffixTree.EmProto.csproj: OK', out)
        self.assertEqual(missing_of(out), ['a-file-that-the-context-does-not-hold.txt'], out[-800:])
        self.assertEqual(p.returncode, 1, out[-800:])

    def test_the_full_check_of_the_real_list_passes(self):
        p = subprocess.run([sys.executable, '-I', CHECK, ROOT], capture_output=True, text=True, timeout=600)
        out = p.stdout + p.stderr
        self.assertEqual(p.returncode, 0, out[-800:])
        self.assertIn('context closure: OK', out)


if __name__ == '__main__':
    unittest.main(verbosity=1)
