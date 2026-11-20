#!/usr/bin/env python3
"""A4-4 (B09 F156) oracle: the option boundaries of the new ports, measured on the installed programs (minimap2 2.26, mgaps 3.23 and EMBOSS etandem 6.6.0).

The option objects of Minimap2Chaining (MinimapChainOptions), NucmerClusters (NucmerClusterOptions) and FindTandemRepeatsEtandem (-threshold) accepted values the programs reject or treat in a defined way;
this file records what each program does with them, so that the tests lock the program's behaviour and not an assumption:
  minimap2  : a fixed pair (1 800-base reference, 1 100-base query with three diverged copies, one reverse-complemented); for every option value the chains the program prints (`--print-chains` PAF s1 score and
              cm anchor count of every hit, `-N 1000 -p 0`), or rc != 0.   Rows: option property, value text, chains [[score, anchors], ...] sorted descending.
  mgaps     : five exact matches (1 1 20, 31 31 20, 61 70 20, 130 130 25, 400 410 30); the clusters `mgaps -l 1 -s 90 -d 5 -f 0.12` prints with one option changed.
  etandem   : -threshold int.MinValue / int.MinValue + 1 / -1 / 0 on short sequences: the report rows, or CRASH (the program dies with SIGFPE, exit status 136).
Output: tests/.../TestData/Alignment/option_boundary_oracle.json
Usage: python3 scripts/alignment/generate_option_boundary_oracle.py OUTFILE           (needs minimap2, mgaps and etandem on PATH)
"""
import json, os, random, subprocess, sys, tempfile

rnd = random.Random(2026)
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def rseq(n): return ''.join(rnd.choice('ACGT') for _ in range(n))


def mut(s, r):
    o = []
    for c in s:
        x = rnd.random()
        if x < r: o.append(rnd.choice('ACGT'))
        elif x < r * 1.3: pass
        else: o.append(c)
    return ''.join(o)


def rc(s): return ''.join(COMP[c] for c in reversed(s))


# `-r 500,500` in every run (bw = bw_long): the long-join bandwidth (-r b,B with B > b, default 20000) is the registered A4-7 and not part of the chaining port; the oracle generator of the chains does the same
MM_ROWS = [  # (property, value text, minimap2 arguments)
    ('Bandwidth', '500', []), ('Bandwidth', '-5', ['-r', '-5,-5']), ('Bandwidth', '0', ['-r', '0,0']), ('Bandwidth', '100', ['-r', '100,100']),
    ('MaxTargetGap', '-1', ['-g', '-1']), ('MaxTargetGap', '0', ['-g', '0']), ('MaxTargetGap', '100', ['-g', '100']),
    ('MinCount', '0', ['-n', '0']), ('MinCount', '-1', ['-n', '-1']), ('MinCount', '10', ['-n', '10']),
    ('MinScore', '-5', ['-m', '-5']), ('MinScore', '0', ['-m', '0']), ('MinScore', '100', ['-m', '100']),
    ('ChainGapScale', '-1', ['--chain-gap-scale', '-1']), ('ChainGapScale', '0', ['--chain-gap-scale', '0']), ('ChainGapScale', '3', ['--chain-gap-scale', '3']),
    ('ChainGapScale', 'NaN', ['--chain-gap-scale', 'nan']), ('ChainGapScale', 'Infinity', ['--chain-gap-scale', 'inf']), ('ChainGapScale', '-Infinity', ['--chain-gap-scale', '-inf']),
    ('ChainGapScale', '1e30', ['--chain-gap-scale', '1e30']), ('ChainGapScale', '1e9', ['--chain-gap-scale', '1e9']),
    ('ChainSkipScale', '1', ['--chain-skip-scale', '1']), ('ChainSkipScale', 'NaN', ['--chain-skip-scale', 'nan']), ('ChainSkipScale', 'Infinity', ['--chain-skip-scale', 'inf']),
    ('MaxChainSkip', '-1', ['--max-chain-skip', '-1']), ('MaxChainSkip', '0', ['--max-chain-skip', '0']), ('MaxChainSkip', '3', ['--max-chain-skip', '3']),
    ('MaxChainIterations', '0', ['--max-chain-iter', '0']), ('MaxChainIterations', '-1', ['--max-chain-iter', '-1']), ('MaxChainIterations', '3', ['--max-chain-iter', '3']),
]
MG_ROWS = [('MaxSeparation', '-1', ['-s', '-1']), ('MaxSeparation', '0', ['-s', '0']), ('MaxSeparation', '10', ['-s', '10']), ('MinClusterLength', '-1', ['-l', '-1']), ('MinClusterLength', '0', ['-l', '0']),
           ('MinClusterLength', '1000', ['-l', '1000']), ('DiagonalDifference', '-1', ['-d', '-1']), ('DiagonalDifference', '0', ['-d', '0']), ('DiagonalFraction', '-1', ['-f', '-1']),
           ('DiagonalFraction', '0', ['-f', '0']), ('DiagonalFraction', '2', ['-f', '2']), ('DiagonalFraction', 'NaN', ['-f', 'nan']), ('DiagonalFraction', 'Infinity', ['-f', 'inf']),
           ('DiagonalFraction', '1e10', ['-f', '1e10']), ('DiagonalFraction', '1e30', ['-f', '1e30']), ('DiagonalFraction', '-1e10', ['-f', '-1e10']), ('UseExtent', 'true', ['-e'])]
ET_ROWS = [('ACACACACACACGTGTGTGTGTGTAAA', 2, 4, t) for t in (-2147483648, -2147483647, -2147483646, -100, -1, 0)] + \
          [(s, 2, 4, -2147483648) for s in ('ACGT', 'AAAAAAAA', 'ACGTACGTAC', 'A', 'ACACAC', 'GGGGCCCCGGGGCCCC')]


def main(out):
    ref = rseq(1800)
    q = mut(ref[100:600], 0.03) + rseq(40) + rc(mut(ref[900:1300], 0.03)) + rseq(30) + mut(ref[1400:1700], 0.03)
    doc = {'minimap2': {'reference': ref, 'query': q, 'rows': []}, 'mgaps': {'matches': [[1, 1, 20], [31, 31, 20], [61, 70, 20], [130, 130, 25], [400, 410, 30]], 'rows': []}, 'etandem': []}
    with tempfile.TemporaryDirectory() as d:
        open(d + '/r.fa', 'w').write('>r0\n' + ref + '\n'); open(d + '/q.fa', 'w').write('>q\n' + q + '\n')
        for prop, val, args in MM_ROWS:
            p = subprocess.run(['minimap2', '-k', '15', '-w', '10', '-r', '500,500'] + args + ['-N', '1000', '-p', '0', d + '/r.fa', d + '/q.fa'], capture_output=True, text=True)
            chains = []
            for line in p.stdout.split('\n'):
                f = line.split('\t')
                if len(f) < 12: continue
                tags = {t.split(':')[0]: t.split(':')[2] for t in f[12:]}
                chains.append([int(tags['s1']), int(tags['cm'])])
            doc['minimap2']['rows'].append({'property': prop, 'value': val, 'rc': p.returncode, 'chains': sorted(chains, reverse=True)})
        inp = '> q\n' + ''.join(' %d %d %d\n' % tuple(m) for m in doc['mgaps']['matches'])
        for prop, val, args in MG_ROWS:
            p = subprocess.run(['mgaps', '-l', '1', '-s', '90', '-d', '5', '-f', '0.12'] + args, input=inp, capture_output=True, text=True)
            clusters = []; cur = []
            for line in p.stdout.split('\n'):
                t = line.split()
                if not t or t[0] == '>': continue
                if t[0] == '#':
                    if cur: clusters.append(cur); cur = []
                    continue
                cur.append([int(t[0]), int(t[1]), int(t[2])])
            if cur: clusters.append(cur)
            doc['mgaps']['rows'].append({'property': prop, 'value': val, 'rc': p.returncode, 'clusters': clusters})
        for seq, mn, mx, th in ET_ROWS:
            open(d + '/s.fa', 'w').write('>s\n' + seq + '\n')
            for f in ('o.rep', 'o.orig'):
                if os.path.exists(d + '/' + f): os.remove(d + '/' + f)
            p = subprocess.run(['etandem', '-sequence', d + '/s.fa', '-outfile', d + '/o.rep', '-origfile', d + '/o.orig', '-auto', '-threshold', str(th), '-minrepeat', str(mn), '-maxrepeat', str(mx)], capture_output=True, text=True)
            rows = None
            if p.returncode == 0 and os.path.exists(d + '/o.orig'):
                rows = ['|'.join(','.join(l.split()[:7]) for l in open(d + '/o.orig') if len(l.split()) >= 7)] if os.path.getsize(d + '/o.orig') else ['']
            doc['etandem'].append({'sequence': seq, 'minRepeat': mn, 'maxRepeat': mx, 'threshold': th, 'rc': p.returncode, 'rows': rows[0] if rows else None})
    with open(out, 'w') as fh:
        json.dump(doc, fh, indent=1)
        fh.write('\n')
    print({k: len(v['rows']) if isinstance(v, dict) else len(v) for k, v in doc.items()})


if __name__ == '__main__':
    main(sys.argv[1])
