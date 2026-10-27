#!/usr/bin/env python3
"""Executes the GENUINE CRISPOR scoring functions (crispor.py of maximilianh/crisporWebsite) as an oracle (B08 H4/H5).

Every top-level function/class definition, literal assignment and stdlib import of crispor.py is executed in a fresh
namespace by AST; nothing else (no main, no side effects). `annotateOfftargets`, `calcHitScore`, `calcCfdScore`, ...
then run unchanged. Usage: crispor_exec.py PATH_TO_crispor.py   (smoke test)
"""
import ast, sys, os
def load(path):
    src = open(path, encoding='utf-8').read()
    tree = ast.parse(src)
    ns = {'__name__': 'crispor_ns', '__file__': path}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            exec(compile(ast.Module([node], []), path, 'exec'), ns)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            try: exec(compile(ast.Module([node], []), path, 'exec'), ns)
            except Exception: pass
        elif isinstance(node, ast.Assign):
            try:
                ast.literal_eval(node.value)
                exec(compile(ast.Module([node], []), path, 'exec'), ns)
            except Exception:
                pass
    return ns
if __name__ == '__main__':
    ns = load(sys.argv[1])
    print(len(ns), 'names;', 'annotateOfftargets' in ns, 'calcHitScore' in ns, ns.get('MAXOCC'))
