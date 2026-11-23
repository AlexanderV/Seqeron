#!/usr/bin/env python3
"""Mutants of the audit of work package 1b.12 (prefix 1b12_): the members of the index contract that take a locus by `in` and hand a locus (or children) back, the fix of finding R01 of the review of 1b.10.
   Three files, all single-site, every `old` occurs exactly once:
     EsaSuffixIndex.cs (A)            Extend, TryExtendRight (symbol), TryExtendLeft: the three members that the fix changed. Mutants restore the read after the write in each place, leave a local uncopied, read the wrong local,
                                      leave a non-default result on a failure path, and move the extension depth by one.
     EsaSuffixIndex.Tree.cs (B)       TryParent, TryChild, TrySuffixLink, GetChildren: not changed by the fix (they resolve the locus into a handle before they write). Mutants put the write before the read, or read
                                      the locus again after the write in a way that changes the answer only when result and locus are one variable.
     EsaSuffixTree.Navigation.cs (C)  the facade: each pass-through gets a write to the result (or to the first slot of the children) in front of the call.
   usage: make_1b12_mutants.py <repo root>      (writes 1b12_<file>_listed.json next to this script)"""
import json, os, sys

ROOT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'
NEW = 'FullyQualifiedName~IndexAliasingTests|FullyQualifiedName~IndexAliasingPropertyTests|FullyQualifiedName~FacadeAliasingTests'
FILES = {
    'esasuffixindex': SRC + 'EsaSuffixIndex.cs',
    'esasuffixindex.tree': SRC + 'EsaSuffixIndex.Tree.cs',
    'esasuffixtree.navigation': SRC + 'EsaSuffixTree.Navigation.cs',
}
G = {k: [] for k in FILES}


def m(key, i, desc, old, new):
    G[key].append({'id': '1b12_' + i, 'desc': desc, 'old': old, 'new': new})


A = lambda *a: m('esasuffixindex', *a)  # noqa: E731
EXT_HEAD = '        long length = locus.Length;\n        result = default;\n        EsaText text = _index.TextView;'
A('A01', 'Extend: the depth is read after the write to the result (the original defect)', EXT_HEAD, '        result = default;\n        long length = locus.Length;\n        EsaText text = _index.TextView;')
A('A02', 'Extend: the locus is resolved after the write to the result', '        EsaNode node = Resolve(locus, nameof(locus));\n        long length = locus.Length;\n        result = default;',
  '        result = default;\n        EsaNode node = Resolve(locus, nameof(locus));\n        long length = locus.Length;')
A('A03', 'Extend: the depth is taken from the (cleared) result, the wrong local', EXT_HEAD, '        result = default;\n        long length = result.Length;\n        EsaText text = _index.TextView;')
A('A04', 'Extend: the depth local is not copied from the locus (starts at zero)', EXT_HEAD, '        result = default;\n        long length = 0;\n        EsaText text = _index.TextView;')
A('A05', 'Extend: the depth moves by two per symbol', '            length++;\n        }\n        result = new IndexLocus(node.Lb, node.Rb, length, node.Depth);', '            length += 2;\n        }\n        result = new IndexLocus(node.Lb, node.Rb, length, node.Depth);')
A('A06', 'Extend: the depth of the result is one short', 'result = new IndexLocus(node.Lb, node.Rb, length, node.Depth);', 'result = new IndexLocus(node.Lb, node.Rb, length - 1, node.Depth);')
A('A07', 'Extend: the aux of the result is the matched length, not the depth of the node below', 'result = new IndexLocus(node.Lb, node.Rb, length, node.Depth);', 'result = new IndexLocus(node.Lb, node.Rb, length, length);')
A('A08', 'Extend: a failed extension leaves the locus in the result, not the default', EXT_HEAD, '        result = locus;\n        long length = locus.Length;\n        EsaText text = _index.TextView;')
A('A09', 'TryExtendRight(symbol): the locus is checked after the write to the result (a symbol that is no code unit)', '            Resolve(locus, nameof(locus));\n            result = default;\n            return false;\n        }\n        ReadOnlySpan<char> one',
  '            result = default;\n            Resolve(locus, nameof(locus));\n            return false;\n        }\n        ReadOnlySpan<char> one')
A('A10', 'TryExtendRight(symbol): a symbol that is no code unit leaves the locus in the result', '            Resolve(locus, nameof(locus));\n            result = default;\n            return false;\n        }\n        ReadOnlySpan<char> one',
  '            Resolve(locus, nameof(locus));\n            result = locus;\n            return false;\n        }\n        ReadOnlySpan<char> one')
A('A11', 'TryExtendRight(symbol): the code unit U+FFFF is refused', 'if ((uint)symbol > char.MaxValue)\n        {\n            Resolve(locus, nameof(locus));', 'if ((uint)symbol >= char.MaxValue)\n        {\n            Resolve(locus, nameof(locus));')
LEFT_HEAD = '        IndexLocus at = locus; // a copy: result may be the same variable as locus (the contract), and result is written before the locus is read again\n        Resolve(at, nameof(locus));\n        result = default;\n'
A('A12', 'TryExtendLeft: the copy of the locus is taken after the write to the result', LEFT_HEAD, '        Resolve(locus, nameof(locus));\n        result = default;\n        IndexLocus at = locus;\n')
A('A13', 'TryExtendLeft: the empty-string branch reads the locus, not the copy', 'if (at.Length == 0) return TryExtendRight(at, symbol, out result);', 'if (locus.Length == 0) return TryExtendRight(locus, symbol, out result);')
A('A14', 'TryExtendLeft: the first rank is read from the locus, not the copy', 'for (int rank = (int)at.Lo, last = (int)at.Hi;', 'for (int rank = (int)locus.Lo, last = (int)at.Hi;')
A('A15', 'TryExtendLeft: the last rank is read from the locus, not the copy', 'for (int rank = (int)at.Lo, last = (int)at.Hi;', 'for (int rank = (int)at.Lo, last = (int)locus.Hi;')
A('A16', 'TryExtendLeft: the depth of the result is read from the locus, not the copy', 'result = new IndexLocus(low, high, at.Length + 1, aux);', 'result = new IndexLocus(low, high, locus.Length + 1, aux);')
A('A17', 'TryExtendLeft: the depth of the result is one too long', 'result = new IndexLocus(low, high, at.Length + 1, aux);', 'result = new IndexLocus(low, high, at.Length + 2, aux);')
A('A18', 'TryExtendLeft: no preceding symbol leaves the locus in the result, not the default', '        if (high < 0) return false;', '        if (high < 0) { result = at; return false; }')
A('A19', 'TryExtendLeft: a symbol that is no code unit leaves the locus in the result (the result is set to the copy at the start)', '        Resolve(at, nameof(locus));\n        result = default;\n        if ((uint)symbol > char.MaxValue) return false;', '        Resolve(at, nameof(locus));\n        result = at;\n        if ((uint)symbol > char.MaxValue) return false;')
A('A20', 'TryExtendLeft: the locus is not checked', '        Resolve(at, nameof(locus));\n        result = default;\n        if ((uint)symbol > char.MaxValue) return false;', '        result = default;\n        if ((uint)symbol > char.MaxValue) return false;')

B = lambda *a: m('esasuffixindex.tree', *a)  # noqa: E731
B('B01', 'TryParent: the write to the result is before the read of the locus', '        EsaNode node = Resolve(locus, nameof(locus)); // the locus is copied into the handle before parent is written: parent may be the same variable as locus\n        parent = default;',
  '        parent = default;\n        EsaNode node = Resolve(locus, nameof(locus));')
B('B02', 'TryParent: the parent depth is bounded by the matched length of the locus, read after the write (the same for two variables)', '        int parentDepth = node.ParentDepth;', '        int parentDepth = (int)Math.Min(node.ParentDepth, locus.Length);')
B('B03', 'TryParent: the root test reads the locus after the write (the same for two variables)', '        if (node.Depth == 0) return false;', '        if (locus.Aux == 0) return false;')
B('B04', 'TryChild: the write to the result is before the read of the node', '        EsaNode resolved = ResolveNode(node, nameof(node)); // copied before child is written: child may be the same variable as node\n        child = default;',
  '        child = default;\n        EsaNode resolved = ResolveNode(node, nameof(node));')
B('B05', 'TryChild: the node is compared with its handle after the write (the same for two variables)', 'if (!FindChild(resolved, firstSymbol, out EsaNode found)) return false;', 'if (node.Aux != resolved.Depth || !FindChild(resolved, firstSymbol, out EsaNode found)) return false;')
B('B06', 'TrySuffixLink: the write to the result is before the read of the node', '        EsaNode resolved = ResolveNode(node, nameof(node)); // copied before link is written: link may be the same variable as node\n        link = default;',
  '        link = default;\n        EsaNode resolved = ResolveNode(node, nameof(node));')
B('B07', 'TrySuffixLink: the node is compared with its handle after the write (the same for two variables)', 'if (resolved.Depth == 0 || resolved.IsLeaf) return false;', 'if (resolved.Depth == 0 || resolved.IsLeaf || node.Aux != resolved.Depth) return false;')
B('B08', 'GetChildren: the first slot of the destination is cleared before the node is read', '        EsaNode resolved = ResolveNode(node, nameof(node)); // copied before the first child is written: children may hold the variable that holds node',
  '        if (children.Length > 0) children[0] = default;\n        EsaNode resolved = ResolveNode(node, nameof(node));')
B('B09', 'GetChildren: the number is returned only if the node still has its depth after the writes (the same for a separate array)', '        };\n        return count;\n    }', '        };\n        return node.Length == resolved.Depth ? count : 0;\n    }')
B('B10', 'GetChildren: an exact fit is taken for too small', 'if (count > children.Length) return count;', 'if (count >= children.Length) return count;')

C = lambda *a: m('esasuffixtree.navigation', *a)  # noqa: E731
C('C01', 'facade TryExtendRight(symbol): the result is cleared before the call', '        using ReadScope scope = Enter();\n        return _contract.TryExtendRight(in locus, symbol, out result);', '        using ReadScope scope = Enter();\n        result = default;\n        return _contract.TryExtendRight(in locus, symbol, out result);')
C('C02', 'facade TryExtendRight(span): the result is cleared before the call', '        using ReadScope scope = Enter();\n        return _contract.TryExtendRight(in locus, pattern, out result);', '        using ReadScope scope = Enter();\n        result = default;\n        return _contract.TryExtendRight(in locus, pattern, out result);')
C('C03', 'facade TryExtendLeft: the result is cleared before the call', '        using ReadScope scope = Enter();\n        return _contract.TryExtendLeft(in locus, symbol, out result);', '        using ReadScope scope = Enter();\n        result = default;\n        return _contract.TryExtendLeft(in locus, symbol, out result);')
C('C04', 'facade TryParent: the result is cleared before the call', '        using ReadScope scope = Enter();\n        return _contract.TryParent(in locus, out parent);', '        using ReadScope scope = Enter();\n        parent = default;\n        return _contract.TryParent(in locus, out parent);')
C('C05', 'facade TryChild: the result is cleared before the call', '        using ReadScope scope = Enter();\n        return _contract.TryChild(in node, firstSymbol, out child);', '        using ReadScope scope = Enter();\n        child = default;\n        return _contract.TryChild(in node, firstSymbol, out child);')
C('C06', 'facade TrySuffixLink: the result is cleared before the call', '        using ReadScope scope = Enter();\n        return _contract.TrySuffixLink(in node, out link);', '        using ReadScope scope = Enter();\n        link = default;\n        return _contract.TrySuffixLink(in node, out link);')
C('C07', 'facade GetChildren: the first slot of the destination is cleared before the call', '        using ReadScope scope = Enter();\n        return _contract.GetChildren(in node, children);', '        using ReadScope scope = Enter();\n        if (children.Length > 0) children[0] = default;\n        return _contract.GetChildren(in node, children);')

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
    json.dump(lst, open(os.path.join(HERE, '1b12_%s_listed.json' % key), 'w'), indent=1)
    out.append((key, len(lst)))
print('listed', out, sum(n for _, n in out), 'mutants')
