#!/usr/bin/env python3
"""ADR-005 section 6: every place where a 32-bit type limits a size or a position in the suffix-tree projects (file:line, category, the line), by pattern.
usage: int_inventory.py [repo root] > results/int_inventory.tsv        (read only; the categories are grep patterns, the curated reading is in ADR-005 section 6)
Categories: LIMIT (a named bound), VIEW (a 32-bit array, pointer or span of positions or ranks), CAST (a narrowing cast of a position), SIG (a public or internal member that takes or returns int for a position / count),
B32 (the builder's 32-bit positions: uint arrays, DistTarget<uint>, ...)."""
import os, re, sys
root = sys.argv[1] if len(sys.argv) > 1 else os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../../..'))
dirs = ['src/SuffixTree/Algorithms/SuffixTree.Esa', 'src/SuffixTree/Algorithms/SuffixTree.Core']
pats = [
    ('LIMIT', re.compile(r'int\.MaxValue|Array\.MaxLength|uint\.MaxValue|MaxTextLength|TaggedChild\.MaxLength|MaxSymbols|MaxLength\b')),
    ('VIEW',  re.compile(r'\bint\[\]|\bint\*|Buf<int>|PinnedArray<int>|Span<int>|ReadOnlySpan<int>|IReadOnlyList<int>|List<int>|uint\[\]|uint\*|Buf<uint>|DistTarget<uint>')),
    ('CAST',  re.compile(r'\(int\)')),
    ('SIG',   re.compile(r'^\s*(?:public|internal|private|protected)[^=;{(]*\b(?:int|IReadOnlyList<int>)\s+\w+(?:\s*\(|\s*\{|\s*=>)')),
]
print('file\tline\tcategory\ttext')
seen = set()
for d in dirs:
    for dp, dn, fn in os.walk(os.path.join(root, d)):
        dn[:] = [x for x in dn if x not in ('bin', 'obj')]
        for f in sorted(fn):
            if not f.endswith('.cs'): continue
            p = os.path.join(dp, f)
            rel = os.path.relpath(p, root)
            for i, line in enumerate(open(p, encoding='utf-8'), 1):
                s = line.strip()
                if s.startswith('///') or s.startswith('//') or s.startswith('*'): continue
                for cat, rx in pats:
                    if rx.search(line):
                        key = (rel, i, cat)
                        if key in seen: continue
                        seen.add(key)
                        print(f"{rel}\t{i}\t{cat}\t{s[:170]}")
