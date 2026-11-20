#!/usr/bin/env python3
"""Work package 1b.9c: a second implementation of DifferenceCover.Pruned (SplitMix64 stream, greedy removal), written outside the production code; it gave the figures of SaCoverGap1b9cTests (|D| and the sum of the residues for the periods that the test pins)
   and showed which periods the mutant D14 (a 13th round for 256 < X <= 1024) changes: 323, 341, 380, 443, 672, 754. usage: pruned_cover_1b9c.py <lo> <hi> (the sizes of the covers of the periods lo..hi)."""
import sys
M=(1<<64)-1
class SM:
    def __init__(s,seed): s.st=seed&M
    def next(s):
        s.st=(s.st+0x9E3779B97F4A7C15)&M
        z=s.st
        z=((z^(z>>30))*0xBF58476D1CE4E5B9)&M
        z=((z^(z>>27))*0x94D049BB133111EB)&M
        return z^(z>>31)
def pruned(x, rounds_fn=None, largest=False, k0=False):
    m=1
    while m*m<x: m+=1
    start=set()
    for j in range(m): start.add(j%x)
    for k in range(0 if k0 else 1, m+1): start.add(m*k%x)
    start=sorted(start)
    best=list(start)
    rng=SM(0x9E3779B97F4A7C15 ^ x)
    rounds=6
    if x<=256: rounds=40
    elif x<=1024: rounds=12
    if rounds_fn: rounds=rounds_fn(x,rounds)
    for rd in range(rounds):
        d=list(start)
        keys=[rng.next() for _ in d]
        order=[d[i] for i in sorted(range(len(d)), key=lambda i: keys[i])]
        count=[0]*x
        for a in d:
            for b in d: count[(b-a+x)%x]+=1
        for r in order:
            delta=[]
            for b in d:
                delta.append((b-r+x)%x)
                if b!=r: delta.append((r-b+x)%x)
            ok=True
            for dd in delta: count[dd]-=1
            for dd in delta:
                if count[dd]<=0:
                    ok=False;break
            if ok: d.remove(r)
            else:
                for dd in delta: count[dd]+=1
        if (len(d)>len(best)) if largest else (len(d)<len(best)): best=d
    return sorted(best)
if __name__=='__main__':
    import json
    lo,hi=int(sys.argv[1]),int(sys.argv[2])
    out={}
    for x in range(lo,hi+1):
        # skip singer periods
        out[x]=len(pruned(x))
    print(out)
