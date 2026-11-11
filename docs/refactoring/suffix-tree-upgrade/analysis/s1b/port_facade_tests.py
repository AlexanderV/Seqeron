#!/usr/bin/env python3
"""
Work package 1b.6: derives the backend-parameterised copies of the existing ISuffixTree-level tests.

Reads the original test files (they are NOT edited), rewrites only the construction of the tree
(`SuffixTree.Build(x)` -> `Tree(x)`, `PersistentSuffixTreeFactory.Create(new StringTextSource(x))` -> `Tree(x)`),
the namespace and the class header (fixture runs once per backend of `SuffixTreeBackends`), removes the methods that are not
backend-agnostic (config in port_facade_config.py, with the reason of every removal) and writes the derived files below
tests/SuffixTree/SuffixTree.Esa.Tests/Facade/Ported. Assertions are copied verbatim.

Also writes docs/.../analysis/s1b/facade_traceability.tsv: one row per test class (out of scope classes) or per test method.

Usage: python3 port_facade_tests.py            (from anywhere; paths are derived from this file)
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..', '..'))
sys.path.insert(0, HERE)
import port_facade_config as cfg  # noqa: E402

OUT_DIR = os.path.join(ROOT, 'tests', 'SuffixTree', 'SuffixTree.Esa.Tests', 'Facade', 'Ported')
TSV = os.path.join(HERE, 'facade_traceability.tsv')

ORIGINS = {
    'Ukkonen': ('tests/SuffixTree/SuffixTree.Tests', 'SuffixTree.Tests'),
    'Persistent': ('tests/SuffixTree/SuffixTree.Persistent.Tests', 'SuffixTree.Persistent.Tests'),
    'McpCore': ('tests/SuffixTree/SuffixTree.Mcp.Core.Tests', 'SuffixTree.Mcp.Core.Tests'),
}

TEST_ATTR = re.compile(r'^\s*\[(Test|TestCase|TestCaseSource|Theory|Property|Repeat|Combinatorial|Sequential)\b')
METHOD = re.compile(r'^\s*(public|internal)\s+(?:async\s+)?(?:static\s+)?[\w<>\[\],\.\?\(\) ]+?\s+(\w+)\s*\(')
CLASS = re.compile(r'^\s*(?:public|internal)\s+(?:sealed\s+|static\s+|abstract\s+)*class\s+(\w+)')


def read(path):
    with open(path, encoding='utf-8') as f:
        return f.read()


def skip_literal(s, i):
    """Returns the index after the string/char literal that starts at s[i] (a quote)."""
    q = s[i]
    verbatim = i > 0 and s[i - 1] == '@' or (i > 1 and s[i - 2:i] in ('$@', '@$'))
    i += 1
    while i < len(s):
        c = s[i]
        if verbatim:
            if c == '"':
                if i + 1 < len(s) and s[i + 1] == '"':
                    i += 2
                    continue
                return i + 1
        else:
            if c == '\\':
                i += 2
                continue
            if c == q:
                return i + 1
        i += 1
    return i


def code_chars(line):
    """Yield (index, char) of the characters of a line outside string literals and comments."""
    i = 0
    while i < len(line):
        c = line[i]
        if c == '/' and i + 1 < len(line) and line[i + 1] == '/':
            return
        if c in '"\'':
            i = skip_literal(line, i)
            continue
        yield i, c
        i += 1


def brace_delta(line):
    d = 0
    for _, c in code_chars(line):
        if c == '{':
            d += 1
        elif c == '}':
            d -= 1
    return d


def method_span(lines, idx):
    """Span [start, end) of the method whose signature is on lines[idx] (attributes and doc comments before it included)."""
    start = idx
    while start > 0 and (lines[start - 1].strip().startswith('[') or lines[start - 1].strip().startswith('///')):
        start -= 1
    # an expression-bodied member ends at the first line that ends with ';' (outside literals)
    depth = 0
    seen_open = False
    j = idx
    while j < len(lines):
        d = brace_delta(lines[j])
        if d != 0 or '{' in lines[j]:
            seen_open = True
        depth += d
        if not seen_open and lines[j].rstrip().endswith(';'):
            return start, j + 1
        if seen_open and depth <= 0:
            return start, j + 1
        j += 1
    raise RuntimeError('unterminated method at line %d' % (idx + 1))


def find_matching_paren(s, open_idx):
    depth = 0
    i = open_idx
    while i < len(s):
        c = s[i]
        if c in '"\'':
            i = skip_literal(s, i)
            continue
        if c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise RuntimeError('unbalanced parentheses')


def rewrite_calls(s, callee, replacement, arg_unwrap=None):
    """Replaces `callee(ARGS)` by `replacement(ARGS')`; with arg_unwrap, `new StringTextSource(X)` inside ARGS becomes X."""
    out = []
    i = 0
    while True:
        k = s.find(callee, i)
        if k < 0:
            out.append(s[i:])
            break
        p = k + len(callee)
        if k > 0 and (s[k - 1] in '.:_' or s[k - 1].isalnum()):
            # part of a longer name (global::SuffixTree.SuffixTree.Build is the oracle and stays)
            out.append(s[i:p])
            i = p
            continue
        if p >= len(s) or s[p] != '(':
            out.append(s[i:p])
            i = p
            continue
        close = find_matching_paren(s, p)
        args = s[p + 1:close]
        if arg_unwrap:
            m = re.fullmatch(r'\s*new StringTextSource\((.*)\)\s*(,.*)?', args, re.S)
            if m:
                inner = m.group(1)
                rest = m.group(2) or ''
                if rest.strip():
                    out.append(s[i:k] + replacement + '(' + inner + rest + ')')
                else:
                    out.append(s[i:k] + replacement + '(' + inner + ')')
                i = close + 1
                continue
        out.append(s[i:k] + replacement + '(' + args + ')')
        i = close + 1
    return ''.join(out)


def list_tests(text):
    """[(class, method, line)] of the test methods of a file (top-level and nested classes)."""
    lines = text.split('\n')
    cls = None
    tests = []
    pending = False
    for n, line in enumerate(lines):
        m = CLASS.match(line)
        if m:
            cls = m.group(1)
        if TEST_ATTR.match(line):
            pending = True
            continue
        if pending:
            if line.strip().startswith('[') or line.strip().startswith('///') or not line.strip():
                continue
            mm = METHOD.match(line)
            if mm:
                tests.append((cls, mm.group(2), n + 1))
            pending = False
    return tests


def derive(origin, rel, spec):
    folder, ns_root = ORIGINS[origin]
    path = os.path.join(ROOT, folder, rel)
    text = read(path)
    original_tests = list_tests(text)
    lines = text.split('\n')
    removed = []
    for method, reason in spec.get('exclude', {}).items():
        idxs = [n for n, l in enumerate(lines) if METHOD.match(l) and METHOD.match(l).group(2) == method]
        if not idxs:
            raise RuntimeError('%s: excluded method %s not found' % (rel, method))
        # remove from the last occurrence upwards (overloads of one name are removed together)
        for idx in reversed(idxs):
            a, b = method_span(lines, idx)
            lines[a:b] = ['        // not ported: %s - %s' % (method, reason)] if idx == idxs[0] else []
        removed.append(method)
    # per-backend skips: the assertion is the original one, but a backend that differs on purpose skips it with the reason
    for method, skips in spec.get('backend_skip', {}).items():
        idxs = [n for n, l in enumerate(lines) if METHOD.match(l) and METHOD.match(l).group(2) == method]
        if not idxs:
            raise RuntimeError('%s: backend_skip method %s not found' % (rel, method))
        for idx in reversed(idxs):
            j = idx
            while '{' not in lines[j]:
                j += 1
            inserted = []
            for backend, reason in skips.items():
                inserted.append('            if (Backend.Name == "%s") Assert.Ignore("%s");' % (backend, reason.replace('"', "'")))
            lines[j + 1:j + 1] = inserted
    text = '\n'.join(lines)

    # construction
    if origin == 'Ukkonen':
        text = rewrite_calls(text, 'SuffixTree.Build', 'Tree')
        text = rewrite_calls(text, 'SuffixTree.TryBuild', 'TryBuild')
    text = rewrite_calls(text, 'PersistentSuffixTreeFactory.Create', 'Tree', arg_unwrap=True)
    text = rewrite_calls(text, 'PersistentSuffixTreeFactory.CreatePersistent', 'Tree', arg_unwrap=True)
    if origin == 'Persistent':
        # the persistent suites construct through the factory and dispose through using; the backend factory trees are disposed by the base class
        text = re.sub(r'\(PersistentSuffixTree\)\s*Tree\(', 'Tree(', text)
        text = re.sub(r'\(PersistentSuffixTree\)\s*(Build|Create)Tree\(', r'\1Tree(', text)
        text = re.sub(r'using\s+\(var\s+(\w+)\s*=\s*(Tree|CreateTree|BuildTree)\((.*?)\)\s+as\s+IDisposable\)(\s*)\{(\s*)var\s+(\w+)\s*=\s*\(ISuffixTree\)\1!;',
                      lambda m: '{' + m.group(5) + 'var ' + m.group(6) + ' = ' + m.group(2) + '(' + m.group(3) + ');', text, flags=re.S)
        text = re.sub(r'using\s+var\s+(\w+)\s*=\s*(Tree|CreateTree|BuildTree)\((.*?)\)\s+as\s+IDisposable;(\s*)var\s+(\w+)\s*=\s*\(ISuffixTree\)\1!;',
                      lambda m: 'var ' + m.group(5) + ' = ' + m.group(2) + '(' + m.group(3) + ');', text, flags=re.S)
        text = re.sub(r'using\s+var\s+(\w+)\s*=\s*(?:\(IDisposable\))?\s*(Tree|BuildTree|CreateTree)\(', r'var \1 = \2(', text)
        text = re.sub(r'^([ \t]*)using\s+\(var\s+(\w+)\s*=\s*(Tree|BuildTree|CreateTree)\(([^\n]*)\)\)[ \t]*$', r'\1var \2 = \3(\4);', text, flags=re.M)
    for pat, repl in spec.get('subs', []):
        text, count = re.subn(pat, repl, text, flags=re.S | re.M)
        if count == 0:
            raise RuntimeError('%s: substitution %r matched nothing' % (rel, pat))
    if origin == 'Ukkonen':
        text = text.replace('SuffixTree.FindLongestCommonSubstrings(', 'FindLongestCommonSubstrings(')
        text = text.replace('SuffixTree.LongestCommonSubstringLengthsBySupport(', 'LongestCommonSubstringLengthsBySupport(')

    # namespace
    sub_ns = '.'.join(ns_root.split('.')[1:])
    def ns_sub(m):
        rest = m.group(1)
        return 'namespace SuffixTree.Esa.Tests.Facade.Ported.%s%s' % (origin, rest)
    text = re.sub(r'namespace\s+' + re.escape(ns_root) + r'((?:\.\w+)*)', ns_sub, text, count=1)

    # class headers: every [TestFixture] class
    def cls_sub(m):
        attrs = m.group(1)
        name = m.group(3)
        base = m.group(4) or ''
        if base.strip():
            raise RuntimeError('%s: class %s has a base list' % (rel, name))
        return (attrs + '[TestFixtureSource(typeof(SuffixTreeBackends), nameof(SuffixTreeBackends.PortedFixtures))]\n'
                + m.group(2) + 'class ' + name + ' : PortedBackendTestBase\n')
    pat = re.compile(r'((?:^[ \t]*\[(?:TestFixture|Category|NonParallelizable|Parallelizable)[^\]\n]*\][ \t]*\n)*?^[ \t]*\[TestFixture\][ \t]*\n(?:^[ \t]*\[[^\]\n]*\][ \t]*\n)*)([ \t]*public (?:sealed )?)class (\w+)(\s*:\s*[^\n{]+)?\n', re.M)
    text, nfix = pat.subn(cls_sub, text)
    if nfix == 0 and not spec.get('no_fixture'):
        raise RuntimeError('%s: no TestFixture class found' % rel)
    # constructor after the opening brace of each fixture class
    def ctor_sub(m):
        indent = m.group(2)
        return m.group(1) + m.group(2) + '{\n' + indent + '    public ' + m.group(3) + '(SuffixTreeBackend backend) : base(backend) { }\n'
    text = re.sub(r'(class (\w+) : PortedBackendTestBase\n)([ \t]*)\{\n', lambda m: m.group(1) + m.group(3) + '{\n' + m.group(3) + '    public ' + m.group(2) + '(SuffixTreeBackend backend) : base(backend) { }\n\n', text)

    header = ('// <auto-generated>\n'
              '// Derived from %s/%s by docs/refactoring/suffix-tree-upgrade/analysis/s1b/port_facade_tests.py (work package 1b.6).\n'
              '// Assertions are the original ones; only the construction of the tree goes through the backend factory. Do not edit: regenerate.\n'
              '// </auto-generated>\n#nullable enable\n') % (folder, rel.replace('\\', '/'))
    text = header + text
    dest_rel = os.path.join(origin, rel)
    dest = os.path.join(OUT_DIR, dest_rel)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, 'w', encoding='utf-8', newline='\n') as f:
        f.write(text)
    return original_tests, removed, nfix


def main():
    rows = []  # (origin, file, class, method, disposition, target, reason)
    for origin, files in cfg.PORTED.items():
        for rel, spec in files.items():
            tests, removed, nfix = derive(origin, rel, spec)
            exclude = spec.get('exclude', {})
            cls_name = os.path.splitext(os.path.basename(rel))[0]
            for cls, method, line in tests:
                if method in exclude:
                    rows.append((origin, rel, cls, method, 'not backend-agnostic', '', exclude[method]))
                else:
                    target = 'Facade/Ported/%s/%s :: %s.%s' % (origin, rel.replace('\\', '/'), cls, method)
                    note = spec.get('note', 'verbatim assertions through the backend factory')
                    for backend, reason in spec.get('backend_skip', {}).get(method, {}).items():
                        note += '; skipped on %s: %s' % (backend, reason)
                    rows.append((origin, rel, cls, method, 'ported', target, note))
    # classes and files that are not ported
    for origin, files in cfg.NOT_PORTED.items():
        folder, _ = ORIGINS[origin]
        for rel, spec in files.items():
            path = os.path.join(ROOT, folder, rel)
            text = read(path)
            tests_here = list_tests(text)
            if not tests_here:
                m = CLASS.search(text)
                tests_here = [(m.group(1) if m else os.path.splitext(os.path.basename(rel))[0], '(no test methods of its own)', 0)]
            for cls, method, line in tests_here:
                own = spec.get('methods', {}).get(method)
                if own is None:
                    disp, target, reason = spec.get('disposition', 'not backend-agnostic'), spec.get('target', ''), spec['reason']
                elif isinstance(own, str):
                    disp, target, reason = 'not backend-agnostic', '', own
                else:
                    disp, target, reason = own
                if disp != 'not backend-agnostic' and not target:
                    raise RuntimeError('%s %s: disposition %s needs a target' % (rel, method, disp))
                rows.append((origin, rel, cls, method, disp, target, reason))
    # every file with tests in the three projects must be accounted for
    for origin, (folder, _) in ORIGINS.items():
        base = os.path.join(ROOT, folder)
        for dirpath, dirs, files in os.walk(base):
            dirs[:] = [d for d in dirs if d not in ('obj', 'bin')]
            for name in files:
                if not name.endswith('.cs'):
                    continue
                full = os.path.join(dirpath, name)
                rel = os.path.relpath(full, base).replace(os.sep, '/')
                if not list_tests(read(full)):
                    continue
                if rel not in cfg.PORTED.get(origin, {}) and rel not in cfg.NOT_PORTED.get(origin, {}):
                    raise RuntimeError('%s/%s has tests but is in neither PORTED nor NOT_PORTED' % (origin, rel))
    rows.sort(key=lambda r: (r[0], r[1], r[2] or '', r[3]))
    with open(TSV, 'w', encoding='utf-8', newline='\n') as f:
        f.write('origin\tfile\tclass\ttest\tdisposition\tported_to\treason_or_note\n')
        for r in rows:
            for field in r:
                if '\t' in field:
                    raise RuntimeError('tab in field: %r' % (r,))
            f.write('\t'.join(r) + '\n')
    n_port = sum(1 for r in rows if r[4] == 'ported')
    print('rows %d, ported %d, other %d' % (len(rows), n_port, len(rows) - n_port))


if __name__ == '__main__':
    main()
