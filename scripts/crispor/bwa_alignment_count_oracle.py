#!/usr/bin/env python3
"""BWA oracle for CRISPOR's `totalAlnCount` (B08 H4): BWA's X0+X1 for a guide equals the number of guide-length windows,
on either strand, within k mismatches (PAM-agnostic).  Runs the reference's own command line
(`bwa aln -o 0 -m 2000000 -n k -k k -N -l 20`, `bwa samse -n 60000`; BWA 0.7.17) on random genomes with planted
0..k+1-mismatch copies on both strands (every fifth guide a palindrome) and compares.  Needs `bwa` on PATH and writes
g.fa / q.fa / q.sai (+ BWA index files) into the CURRENT directory - run it in an empty scratch directory.
Usage: bwa_alignment_count_oracle.py NUMBER_OF_SEEDS     (the B08 run: 40 seeds, 0 differences)
"""
import random, subprocess, sys
COMP={'A':'T','C':'G','G':'C','T':'A'}
def rc(s): return ''.join(COMP[c] for c in reversed(s))
def count_windows(guide, genome, k):
    L=len(guide); n=0
    for strand in (genome, rc(genome)):
        for i in range(len(strand)-L+1):
            mm=0
            for a,b in zip(guide, strand[i:i+L]):
                if a!=b:
                    mm+=1
                    if mm>k: break
            if mm<=k: n+=1
    return n
def mutate(s,k,rng):
    s=list(s)
    for p in rng.sample(range(len(s)),k): s[p]=rng.choice([c for c in 'ACGT' if c!=s[p]])
    return ''.join(s)
def run(seed, k=4, pal=False):
    rng=random.Random(seed)
    if pal:
        h=''.join(rng.choice('ACGT') for _ in range(10)); guide=h+rc(h)
    else:
        guide=''.join(rng.choice('ACGT') for _ in range(20))
    parts=[]
    for _ in range(rng.randint(8,25)):
        parts.append(''.join(rng.choice('ACGT') for _ in range(rng.randint(30,300))))
        c=mutate(guide, rng.randint(0,k+1), rng)
        if rng.random()<0.5: c=rc(c)
        parts.append(c)
    parts.append(''.join(rng.choice('ACGT') for _ in range(100)))
    genome=''.join(parts)
    open('g.fa','w').write('>g\n'+genome+'\n'); open('q.fa','w').write('>q\n'+guide+'\n')
    subprocess.run(['bwa','index','g.fa'],capture_output=True,check=True)
    with open('q.sai','wb') as f:
        subprocess.run(['bwa','aln','-o','0','-m','2000000','-n',str(k),'-k',str(k),'-N','-l','20','g.fa','q.fa'],stdout=f,stderr=subprocess.DEVNULL,check=True)
    sam=subprocess.run(['bwa','samse','-n','60000','g.fa','q.sai','q.fa'],capture_output=True,text=True,stderr=subprocess.DEVNULL if False else None).stdout
    rec=[l for l in sam.splitlines() if not l.startswith('@')][0].split('\t')
    tags=dict((t[:2],t[5:]) for t in rec[11:])
    x0=int(tags.get('X0',0)); x1=int(tags.get('X1',0)); xa=len(tags['XA'].rstrip(';').split(';')) if 'XA' in tags else 0
    return guide, count_windows(guide,genome,k), x0, x1, xa, rec[1]
if __name__=='__main__':
    bad=0
    for seed in range(int(sys.argv[1])):
        g,mine,x0,x1,xa,flag=run(seed, 4, pal=(seed%5==0))
        ok = (mine==x0+x1)
        bad+= (not ok)
        if not ok or seed<4: print(seed, 'palindrome' if seed%5==0 else '', 'windows=%d X0=%d X1=%d XA=%d  X0+X1=%d  %s'%(mine,x0,x1,xa,x0+x1,'OK' if ok else 'DIFF'))
    print('seeds with DIFF:',bad,'of',sys.argv[1])
