#!/usr/bin/env python3
"""Two-sided differential of CrisprDesigner.FindPamSites(string, system) against CRISPOR's genuine PAM search (B08 H12).

Executes crispor.py's own `findAllPams` + `flankSeqIter(doFilterNs=True)` (functions loaded by AST, see crispor_exec.py)
on random strings that contain `N`, for seven systems, and compares the sites: (PAM start, PAM sequence, protospacer,
protospacer start, strand).  B08 run: 10 seeds x 500 cases, 0 sites on either side only (the pre-fix scan had 774-824
Seqeron-only and 43-48 CRISPOR-only sites per 400 cases).

Usage: pam_sites_differential.py CRISPOR_PY DRIVER_DLL SEED CASES N_PROBABILITY
DRIVER_DLL is a 12-line console project referencing Seqeron.Genomics.MolTools that reads "sequence system" lines from
stdin and prints "#case n" followed by "Position<TAB>PamSequence<TAB>TargetSequence<TAB>TargetStart<TAB>fwd(1/0)" per
site of CrisprDesigner.FindPamSites(sequence, Enum.Parse<CrisprSystemType>(system)).
"""
import os, sys, random, subprocess, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import crispor_exec
ns = crispor_exec.load(sys.argv[1])
SYS = {'SpCas9': ('NGG', 20, False), 'SaCas9': ('NNGRRT', 21, False), 'Cas12a': ('TTTV', 23, True), 'CasX': ('TTCN', 20, True),
       'SpCas9_NAG': ('NAG', 20, False), 'SpCas9_NGA': ('NGA', 20, False), 'SaCas9_NNGRRN': ('NNGRRN', 21, False)}
def crispor_sites(seq, pam, guidelen, first):
    ns['pamIsFirst'] = first; ns['GUIDELEN'] = guidelen
    startDict, _ = ns['findAllPams'](seq, pam)
    out = set()
    for name, pamStart, guideStart, strand, flank, pamSeq, plus in ns['flankSeqIter'](seq, startDict, len(pam), True):
        out.add((pamStart, pamSeq, flank, guideStart, 1 if strand == '+' else 0))
    return out
def rand_seq(rng, n, pam, pN):
    alphabet = 'ACGT'
    s = [rng.choice(alphabet) for _ in range(n)]
    for _ in range(rng.randint(0, 3)):           # plant PAM-ish motifs
        p = rng.randint(0, max(0, n - 8)); m = ''.join(rng.choice('ACGT') if c in 'NRV' else c for c in pam)
        s[p:p+len(m)] = list(m)
    for i in range(n):
        if rng.random() < pN: s[i] = 'N'
    return ''.join(s)
def main(seed, ncases, pN):
    rng = random.Random(seed)
    cases = []
    for _ in range(ncases):
        name = rng.choice(list(SYS)); pam, gl, first = SYS[name]
        n = rng.randint(gl + len(pam), 120)
        cases.append((rand_seq(rng, n, pam, pN), name))
    inp = '\n'.join(f'{s} {n}' for s, n in cases) + '\n'
    res = subprocess.run(['dotnet', sys.argv[2]], input=inp, capture_output=True, text=True)
    if res.returncode: print(res.stderr); sys.exit(1)
    got = collections.defaultdict(set); cur = None
    for ln in res.stdout.splitlines():
        if ln.startswith('#case'): cur = int(ln.split()[1]); got[cur]
        else:
            f = ln.split('\t'); got[cur].add((int(f[0]), f[1], f[2], int(f[3]), int(f[4])))
    extra = missing = total = 0; exs = []
    for i, (seq, name) in enumerate(cases):
        pam, gl, first = SYS[name]
        want = crispor_sites(seq, pam, gl, first); have = got[i]; total += len(want)
        extra += len(have - want); missing += len(want - have)
        if (have != want) and len(exs) < 3: exs.append((seq, name, sorted(have - want)[:2], sorted(want - have)[:2]))
    print(f'seed={seed} cases={ncases} pN={pN}: CRISPOR sites={total} Seqeron-only={extra} CRISPOR-only={missing}')
    for e in exs: print('  ', e)
if __name__ == '__main__':
    main(int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]))
