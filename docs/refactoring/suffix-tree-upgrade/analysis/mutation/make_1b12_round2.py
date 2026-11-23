#!/usr/bin/env python3
"""Second round of the audit of work package 1b.12 (prefix 1b12_r2_): mutants that the first round (make_1b12_mutants.py) did not contain, written after its result (37 of 37 killed) to look for a gap on purpose:
   a shortcut for the empty pattern that reads the locus after the write, a write to the result in front of the call in the two public overloads of TryExtendRight, a partial result on each failure path of the extension
   loop, the positions of Locate cleared before the locus is checked (the destination of Locate may overlap the locus), and two mutants that should SURVIVE because they are equivalent (a copy of the locus made at the
   facade, or in the symbol overload: the answer is the same by construction, nothing can tell them from the original). Every `old` occurs exactly once.
   usage: make_1b12_round2.py <repo root>      (writes 1b12_r2_<file>_listed.json next to this script)"""
import json, os, sys

ROOT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'
NEW = 'FullyQualifiedName~IndexAliasingTests|FullyQualifiedName~IndexAliasingPropertyTests|FullyQualifiedName~FacadeAliasingTests'
FILES = {'esasuffixindex': SRC + 'EsaSuffixIndex.cs', 'esasuffixtree.navigation': SRC + 'EsaSuffixTree.Navigation.cs'}
G = {k: [] for k in FILES}


def m(key, i, desc, old, new):
    G[key].append({'id': '1b12_r2_' + i, 'desc': desc, 'old': old, 'new': new})


A = lambda *a: m('esasuffixindex', *a)  # noqa: E731
A('D01', 'Extend: an empty pattern returns a locus built from the locus after the write (the same for two variables)', '        result = default;\n        EsaText text = _index.TextView;\n        foreach (char symbol in pattern)',
  '        result = default;\n        if (pattern.IsEmpty) { result = new IndexLocus(node.Lb, node.Rb, locus.Length, node.Depth); return true; }\n        EsaText text = _index.TextView;\n        foreach (char symbol in pattern)')
A('D02', 'TryExtendRight(symbol): the result is cleared before the call of Extend', '        ReadOnlySpan<char> one = [(char)symbol];\n        return Extend(locus, one, out result);', '        ReadOnlySpan<char> one = [(char)symbol];\n        result = default;\n        return Extend(locus, one, out result);')
A('D03', 'TryExtendRight(span): the result is cleared before the call of Extend', 'public bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result) => Extend(locus, pattern, out result);',
  'public bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result)\n    {\n        result = default;\n        return Extend(locus, pattern, out result);\n    }')
A('E01', 'Extend: a mismatch inside an edge leaves the partial locus in the result', 'if (at >= (uint)_n || text[(int)at] != symbol) return false;', 'if (at >= (uint)_n || text[(int)at] != symbol) { result = new IndexLocus(node.Lb, node.Rb, length, node.Depth); return false; }')
A('E02', 'Extend: a missing child leaves the partial locus in the result', 'if (!FindChild(node, symbol, out EsaNode next)) return false;', 'if (!FindChild(node, symbol, out EsaNode next)) { result = new IndexLocus(node.Lb, node.Rb, length, node.Depth); return false; }')
A('L01', 'Locate: the destination is cleared before the locus is checked (the same unless the destination overlaps the locus)', '        ArgumentOutOfRangeException.ThrowIfNegative(skip);\n        CheckLocus(locus, nameof(locus));\n        int from', '        ArgumentOutOfRangeException.ThrowIfNegative(skip);\n        positions.Clear();\n        CheckLocus(locus, nameof(locus));\n        int from')
A('Q02', 'EQUIVALENT: TryExtendRight(symbol) hands Extend a copy of the locus', '        ReadOnlySpan<char> one = [(char)symbol];\n        return Extend(locus, one, out result);', '        ReadOnlySpan<char> one = [(char)symbol];\n        IndexLocus copy = locus;\n        return Extend(copy, one, out result);')
C = lambda *a: m('esasuffixtree.navigation', *a)  # noqa: E731
C('L02', 'facade Locate: the destination is cleared before the call', '        using ReadScope scope = Enter();\n        return _contract.Locate(in locus, skip, positions);', '        using ReadScope scope = Enter();\n        positions.Clear();\n        return _contract.Locate(in locus, skip, positions);')
C('Q01', 'EQUIVALENT: facade TryParent hands the index a copy of the locus', '        using ReadScope scope = Enter();\n        return _contract.TryParent(in locus, out parent);', '        using ReadScope scope = Enter();\n        IndexLocus copy = locus;\n        return _contract.TryParent(in copy, out parent);')

out = []
for key, rel in FILES.items():
    src = open(os.path.join(ROOT, rel), encoding='utf-8').read()
    lst = []
    for x in G[key]:
        n = src.count(x['old'])
        if n != 1:
            raise SystemExit('%s: pattern occurs %d times: %r' % (x['id'], n, x['old'][:90]))
        if x['old'] == x['new']:
            raise SystemExit('%s: not a change' % x['id'])
        d = dict(x)
        d['file'] = rel
        d['filter'] = '(' + NEW + ')'
        lst.append(d)
    json.dump(lst, open(os.path.join(HERE, '1b12_r2_%s_listed.json' % key), 'w'), indent=1)
    out.append((key, len(lst)))
print('listed', out, sum(n for _, n in out), 'mutants')
