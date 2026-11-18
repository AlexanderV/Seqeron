#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-15b) oracle: ClustalW 2.1 FAST pairwise stage (`-QUICKTREE`, FastPairwiseAlign.cpp).

Seqeron API: SequenceAligner.ClustalWQuickPairwiseAlign(a, b, new ClustalWQuickOptions {...}, maxSequenceLength) and
MsaOptions { Distance = MsaDistanceMode.ClustalWQuick, QuickPairwise = ... }.

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7), `-TYPE=DNA -QUICKTREE -ALIGN [-KTUPLE -WINDOW -TOPDIAGS
-PAIRGAP -SCORE]`; for every pair it prints "Sequences (i:j) Aligned. Score: X" (X = %lg of the double percent score) or
"Sequences (i:j) Not Aligned" (score <= 0.1). The fragment list is bounded by 2 * maxAlnLength where maxAlnLength = 2 * (longest
sequence of the WHOLE input), so every row records the set's longest length (maxLen) - repetitive sets hit the
"(Partial alignment)" exit, which the oracle therefore covers.

ClustalW's command line only accepts -WINDOW/-TOPDIAGS/-PAIRGAP values > KTUPLE and sets PAIRGAP = KTUPLE + 4 with -KTUPLE
(CommandLineParser.cpp); the generator only produces such parameter sets and records the EFFECTIVE values.

Output TSV (tests/.../TestData/Alignment/clustalw_quick_oracle.tsv):
  kind ktuple window topdiags pairgap percent maxLen seq1 seq2 printed   (printed = %lg text or NA for "Not Aligned")
Usage: python3 scripts/alignment/generate_clustalw_quick_oracle.py OUTDIR      (needs `clustalw` on PATH)
"""
import os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261121)
ALPHA = "ACGT"


def gen(kind, L, alpha):
    if kind == 0: return ''.join(rnd.choice(alpha) for _ in range(L))
    if kind == 1: return ''.join(rnd.choice(alpha[:2]) for _ in range(L))
    if kind == 2: return ''.join(rnd.choice([alpha[0]] * 3 + [alpha[1]]) for _ in range(L))
    u = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 4)))
    return (u * L)[:L]


def mutate(s, rate, alpha):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: continue
        elif r < rate: out.append(ch); out.append(rnd.choice(alpha))
        else: out.append(ch)
    return ''.join(out) or alpha[0]


def make_set(maxlen, alpha, kmin=3, kmax=8):
    k = rnd.randint(kmin, kmax)
    base = gen(rnd.randint(0, 3), rnd.randint(5, maxlen), alpha)
    seqs = []
    for _ in range(k):
        if rnd.random() < 0.55: seqs.append(mutate(base, rnd.choice([0.05, 0.15, 0.3, 0.5]), alpha))
        else: seqs.append(gen(rnd.choice([0, 0, 1, 2, 3]), rnd.randint(3, maxlen), alpha))
    return seqs


def params():
    r = rnd.random()
    if r < 0.4: return dict(ktuple=None, window=None, topdiags=None, pairgap=None, percent=True)
    kt = rnd.choice([None, 1, 2, 3, 4])
    k = kt or 2
    return dict(ktuple=kt,
                window=rnd.choice([None, rnd.randint(k + 1, k + 9)]),
                topdiags=rnd.choice([None, rnd.randint(k + 1, k + 10)]),
                pairgap=rnd.choice([None, rnd.randint(k + 5, k + 12)]),   # CLI: ignored unless > ktup; effective default ktup + 4
                percent=rnd.random() < 0.8)


def effective(p):
    k = p['ktuple'] or 2
    return dict(ktuple=k,
                window=p['window'] if p['window'] and p['window'] > k else 4,
                topdiags=p['topdiags'] if p['topdiags'] and p['topdiags'] > k else 4,
                pairgap=p['pairgap'] if p['pairgap'] and p['pairgap'] > k else (k + 4 if p['ktuple'] else 5),
                percent=p['percent'])


def run_clustalw(seqs, p):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'x.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs): f.write(f'>s{i}\n{s}\n')
        args = ['clustalw', f'-INFILE={fa}', '-TYPE=DNA', '-QUICKTREE', '-ALIGN', f'-OUTFILE={d}/o.aln']
        if p['ktuple']: args.append(f"-KTUPLE={p['ktuple']}")
        if p['window']: args.append(f"-WINDOW={p['window']}")
        if p['topdiags']: args.append(f"-TOPDIAGS={p['topdiags']}")
        if p['pairgap']: args.append(f"-PAIRGAP={p['pairgap']}")
        if not p['percent']: args.append('-SCORE=ABSOLUTE')
        out = subprocess.run(args, capture_output=True, text=True, cwd=d).stdout
    res = {}
    for a, b, sc in re.findall(r'Sequences \((\d+):(\d+)\) Aligned\. Score:\s+(\S+)', out): res[(int(a) - 1, int(b) - 1)] = sc
    for a, b in re.findall(r'Sequences \((\d+):(\d+)\) Not Aligned', out): res[(int(a) - 1, int(b) - 1)] = 'NA'
    assert len(res) == len(seqs) * (len(seqs) - 1) // 2, out[-500:]
    return res


def main(outdir):
    rows = []
    sets = 0
    # kind: random/mutated sets of increasing length, IUB N-containing sets, long sets, repeat-rich sets (partial alignment)
    plan = [('acgt', 260, 60, ALPHA), ('acgt', 100, 160, ALPHA), ('long', 5, 700, ALPHA), ('rep', 120, 70, "AC"),
            ('iupac', 120, 70, "ACGTN")]
    for kind, n, maxlen, alpha in plan:
        for _ in range(n):
            p = params()
            seqs = make_set(maxlen, alpha)
            if kind == 'rep':
                seqs = [gen(rnd.choice([2, 3, 1]), rnd.randint(20, maxlen), alpha) for _ in range(rnd.randint(2, 5))] + seqs[:2]
            if kind == 'iupac':
                seqs = [''.join(c if rnd.random() > 0.12 else rnd.choice("NRYKMSWBDHV") for c in s) for s in seqs]
            maxlen_set = max(len(s) for s in seqs)
            res = run_clustalw(seqs, p)
            e = effective(p)
            sets += 1
            for (i, j), sc in sorted(res.items()):
                rows.append((kind, e['ktuple'], e['window'], e['topdiags'], e['pairgap'], int(e['percent']), maxlen_set,
                             seqs[i], seqs[j], sc))
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'clustalw_quick_oracle.tsv'), 'w') as f:
        f.write('kind\tktuple\twindow\ttopdiags\tpairgap\tpercent\tmaxLen\tseq1\tseq2\tprinted\n')
        for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
    print(f'{sets} sets, {len(rows)} pairs; not aligned: {sum(1 for r in rows if r[-1] == "NA")}')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
