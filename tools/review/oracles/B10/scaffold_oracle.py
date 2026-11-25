"""Independent oracle for ASSEMBLY-SCAFFOLD-001.
Greedy (Huson, Reinert & Myers 2002, Alg. 5.1, end-to-end merges; Bambus 'most links first'):
mate-edges in descending weight (stable on input order); an edge joins two contig END nodes;
accepted iff both end nodes still have no selected mate-edge and they lie in different selected
paths. Paths are reverse-complemented with Biopython. Negative distance d: ABySS MergeContigs
(the -d bases of the tail must agree, N matching anything) -> collapse; otherwise AGP 'U' 100.
"""
from Bio.Seq import Seq

UNKNOWN = 100

def rc(s): return str(Seq(s).reverse_complement())

def scaffold(contigs, links, gap='N'):
    n = len(contigs)
    # link = (c1, c2, gap, rev1, rev2, weight)
    L = [l if len(l) == 6 else (l[0], l[1], l[2], False, False, 1) for l in links]
    order = sorted(range(len(L)), key=lambda i: -L[i][5])
    partner = {}  # (contig,'L'/'R') -> ((contig,end), gap)
    comp = list(range(n))
    def find(x):
        while comp[x] != x: x = comp[x]
        return x
    for i in order:
        c1, c2, g, r1, r2, w = L[i]
        if not (0 <= c1 < n and 0 <= c2 < n) or c1 == c2: continue
        a = (c1, 'L' if r1 else 'R'); b = (c2, 'R' if r2 else 'L')
        if a in partner or b in partner: continue
        if find(c1) == find(c2): continue
        partner[a] = (b, g); partner[b] = (a, g); comp[find(c1)] = find(c2)
    groups = {}
    for c in range(n): groups.setdefault(find(c), []).append(c)
    out = []
    for members in groups.values():
        # terminal: contig with a free end; walk from it with free end on the left
        def walk(start):
            path = []; c = start
            o = '+' if (c, 'L') not in partner else '-'
            g = None
            while True:
                path.append((c, o, g))
                ex = (c, 'R' if o == '+' else 'L')
                if ex not in partner: break
                (nc, ne), g = partner[ex]
                c, o = nc, ('+' if ne == 'L' else '-')
            return path
        terms = [c for c in members if (c, 'L') not in partner or (c, 'R') not in partner]
        cands = [walk(t) for t in terms]
        lo = min(members)
        def key(p):
            fwd = sum(1 for (_, o, _) in p if o == '+')
            lo_fwd = any(c == lo and o == '+' for (c, o, _) in p)
            return (-fwd, 0 if lo_fwd else 1)
        best = min(cands, key=key)
        s = ''; since = 0
        for (c, o, g) in best:
            seq = contigs[c] if o == '+' else rc(contigs[c])
            if g is None:
                s = seq; since = len(seq); continue
            if g > 0:
                s += gap * g + seq; since = len(seq); continue
            ov = -g if g < 0 else 0
            if g < 0 and ov <= since and ov <= len(seq) and all(
                    x.upper() == y.upper() or x.upper() == 'N' or y.upper() == 'N'
                    for x, y in zip(s[len(s)-ov:], seq[:ov])):
                tail = s[len(s)-ov:]
                merged = ''.join(y if x.upper() == 'N' else x for x, y in zip(tail, seq[:ov]))
                s = s[:len(s)-ov] + merged + seq[ov:]; since = since - ov + len(seq)
            else:
                s += gap * UNKNOWN + seq; since = len(seq)
        out.append((lo, s))
    return [s for _, s in sorted(out)]

if __name__ == '__main__':
    cases = [
      ("back-link", ["AA","TT"], [(1,0,2)]),
      ("pred-first", ["AA","CC","GG"], [(2,0,1),(0,1,1)]),
      ("branch-out", ["AAAA","CCCC","GGGG"], [(0,1,2),(0,2,2)]),
      ("heavier-wins", ["AAAA","CCCC","GGGG"], [(0,1,2,False,False,1),(0,2,3,False,False,5)]),
      ("convergent", ["AAAA","CCCC","GGGG"], [(0,2,2),(1,2,2)]),
      ("cycle3", ["AA","CC","GG"], [(0,1,1),(1,2,1),(2,0,1)]),
      ("skip-used-old", ["AA","CC","GG"], [(1,0,1),(1,2,1)]),
      ("rc-second", ["AACG","TTGC"], [(0,1,3,False,True,1)]),
      ("rc-first", ["AACG","TTGC"], [(0,1,3,True,False,1)]),
      ("both-rc", ["AACG","TTGC"], [(0,1,3,True,True,1)]),
      ("rc-chain", ["ACCA","GGTT","CATG"], [(0,1,2,False,True,1),(1,2,4,True,False,1)]),
      ("rc-inconsistent", ["ACCA","GGTT","CATG"], [(0,1,2,False,True,1),(1,2,4,False,False,1)]),
      ("overlap-merge", ["ACGTAC","TACGGA"], [(0,1,-3)]),
      ("overlap-mismatch", ["ACGTAC","TTCGGA"], [(0,1,-3)]),
      ("overlap-N", ["ACGTNC","TACGGA"], [(0,1,-3)]),
      ("overlap-too-long", ["ACG","ACGT"], [(0,1,-5)]),
      ("overlap-rc", ["GGATCC","AAACCG"], [(0,1,-2,False,True,1)]),
      ("zero", ["AAAA","TTTT"], [(0,1,0)]),
      ("singletons", ["AAA","CCC","GGG"], [(2,1,4)]),
      ("chain-overlap-then-gap", ["AAACG","CGTTT","GGGG"], [(0,1,-2),(1,2,3)]),
    ]
    for name, c, l in cases:
        print(name, scaffold(c, l))
