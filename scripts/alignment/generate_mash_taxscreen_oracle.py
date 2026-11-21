#!/usr/bin/env python3
"""A5-10b2 (B09 F215) oracle: `mash taxscreen [-m MAP] -t DIR queries.msh pool...` of Mash 2.3 (the Kraken-style taxonomic report: clade / taxon hash counts per taxon, LCA of the references that share a hash).
The taxonomy is SYNTHETIC but in the format of the NCBI dump (nodes.dmp / names.dmp, random trees of 4 - 300 taxa with up to 40 children per node, ranks, scientific names with blanks / punctuation, synonyms, unique names), so no download is needed; the
program only parses the two files. References carry the taxonomy ID in the comment ('taxid N' among other words) or in a mapping file (the file wins), some have none (the program then counts their hashes under taxon 0) or an ID that is not in the dump.
JSON lines: {"id", "args", "files": [{"name","fasta"}] (the pool), "qmsh": base64 query sketch, "names": names.dmp, "nodes": nodes.dmp, "map": mapping file or null, "out": stdout}.
Usage: python3 scripts/alignment/generate_mash_taxscreen_oracle.py OUT.jsonl.gz [n=300] [seed]       (needs mash 2.3 on PATH)
"""
import base64, gzip, json, os, random, subprocess, sys, tempfile

AA = 'ACDEFGHIKLMNPQRSTVWY'
RANKS = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies', 'strain', 'clade', 'no rank', 'subgenus', 'tribe']
CODONS = {}
_code = 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF'
for i, x in enumerate('ACGT'):
    for j, y in enumerate('ACGT'):
        for k, z in enumerate('ACGT'):
            CODONS.setdefault(_code[16 * i + 4 * j + k], []).append(x + y + z)


def main(out, n, seed):
    rnd = random.Random(seed)

    def taxonomy(size, fan, shape='random'):
        ids = [1]
        parent = {1: 1}
        rank = {1: 'no rank'}
        pool = rnd.sample(range(2, 5000000), size - 1)
        kids = {1: 0}
        for t in pool:
            cands = [x for x in ids if kids[x] < fan]
            if shape == 'deep': p = ids[-rnd.randint(1, min(3, len(ids)))]
            elif shape == 'star': p = 1 if rnd.random() < 0.9 else rnd.choice(ids)
            else: p = rnd.choice(cands[-8:]) if rnd.random() < 0.6 and cands else rnd.choice(cands)
            parent[t] = p; kids[p] += 1; kids[t] = 0; ids.append(t)
            rank[t] = rnd.choice(RANKS)
        words = ['Bacteria', 'virus', 'phage', 'sp.', 'str.', 'ATCC', 'K-12', 'subsp.', 'x', 'alpha', 'beta', "O'Brien", 'a/b', 'MG1655']
        nodes = []; names = []
        order = ids[:]
        rnd.shuffle(order[1:])
        for t in order:
            nodes.append('%d\t|\t%d\t|\t%s\t|\t\t|\t%d\t|\t0\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t\t|\n' % (t, parent[t], rank[t], rnd.choice([0, 1, 2, 8])))
            nm = ' '.join(rnd.choice(words) for _ in range(rnd.randint(1, 3))) if t != 1 else 'root'
            if rnd.random() < 0.3: names.append('%d\t|\t%s\t|\t\t|\tsynonym\t|\n' % (t, rnd.choice(words)))
            uniq = nm + ' <x>' if rnd.random() < 0.2 else ''
            names.append('%d\t|\t%s\t|\t%s\t|\tscientific name\t|\n' % (t, nm, uniq))
            if rnd.random() < 0.15: names.append('%d\t|\t%s\t|\t\t|\tcommon name\t|\n' % (t, rnd.choice(words)))
        return ids, parent, ''.join(nodes), ''.join(names)

    def mutate(s, rate, alpha):
        return ''.join(rnd.choice(alpha) if rnd.random() < rate else c for c in s)

    def back(p):
        return ''.join(rnd.choice(CODONS.get(c, ['NNN'])) for c in p)

    def revcomp(s): return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]

    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            shape = rnd.choice(['random'] * 5 + ['deep'] * 3 + ['star'] * 2)
            big = rnd.random() < 0.04 or shape != 'random'
            size = rnd.randint(150, 300) if big else rnd.choice([4, 6, 10, 15, 25, 40, 60, 90])
            if shape == 'deep': size = rnd.choice([20, 40, 70, 140, 260])
            if shape == 'star': size = rnd.choice([25, 40, 70])
            fan = rnd.choice([40, 25, 8, 4, 3]) if big else rnd.choice([2, 3, 5, 40])
            if shape == 'star': fan = 100
            ids, parent, nodes, names = taxonomy(size, fan, shape)
            protein = rnd.random() < 0.3
            alpha = AA if protein else 'ACGT'
            nref = rnd.randint(25, 40) if big else rnd.choice([1, 2, 3, 5, 8, 12, 18])
            if shape == 'deep': nref = rnd.choice([3, 6, 9, 12, 15, 20])
            if shape == 'star': nref = rnd.choice([18, 24, 32, 45])
            root = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(120, 260)))
            # the references: a family, the taxa follow the family tree loosely (siblings are close)
            refs = []
            for i in range(nref):
                if refs and rnd.random() < 0.8:
                    base = rnd.choice(refs)
                    g = mutate(base['g'], rnd.choice([0.02, 0.1, 0.3]), alpha)
                    t = rnd.choice(ids) if rnd.random() < 0.5 else base['t']
                else:
                    g = mutate(root, 0.5, alpha) if rnd.random() < 0.5 else ''.join(rnd.choice(alpha) for _ in range(rnd.randint(150, 400)))
                    t = rnd.choice(ids)
                if shape == 'deep': t = rnd.choice(ids[len(ids) // 2:])
                if shape == 'star':
                    t = ids[1 + i % (len(ids) - 1)]
                    g = ''.join(rnd.choice(alpha) for _ in range(160))
                refs.append({'g': g, 't': t})
            files = []
            mapping = []
            individual = rnd.random() < 0.5
            qfiles = []
            for i, r in enumerate(refs):
                how = rnd.random()
                words = rnd.choice(['', 'some words', 'strain X']).split()
                tid = r['t'] if rnd.random() > 0.05 else rnd.randint(6000000, 7000000)       # an ID that is not in the dump
                comment = ''
                if how < 0.55: comment = ' '.join(words + ['taxid', str(tid)] + (['extra'] if rnd.random() < 0.3 else []))
                elif how < 0.75: mapping.append((tid, 'ref%d' % i if individual else 'r%d.fa' % i))
                elif how < 0.85:
                    comment = 'taxid %d' % tid; mapping.append((rnd.choice(ids), 'ref%d' % i if individual else 'r%d.fa' % i))
                # else: no taxonomy at all
                qfiles.append((i, comment, r['g']))
            if individual:
                fa = ''.join('>ref%d%s\n%s\n' % (i, (' ' + c) if c else '', g) for i, c, g in qfiles)
                qf = [{'name': 'refs.fa', 'fasta': fa}]
            else:
                qf = [{'name': 'r%d.fa' % i, 'fasta': '>ref%d%s\n%s\n' % (i, (' ' + c) if c else '', g)} for i, c, g in qfiles]
            # the pool: reads (or contigs) of some of the references
            pool = []
            for m in range(rnd.choice([1, 1, 2])):
                recs = []
                for r in range(rnd.randint(3, 60) if shape != 'star' else 0):
                    src = rnd.choice(refs)['g'] if rnd.random() < 0.85 else ''.join(rnd.choice(alpha) for _ in range(rnd.randint(80, 300)))
                    st = rnd.randint(0, max(0, len(src) - 60)); piece = src[st:st + rnd.choice([30, 60, 90, 150])]
                    if protein:
                        dn = back(piece)
                        if rnd.random() < 0.5: dn = revcomp(dn)
                        if rnd.random() < 0.3: dn = dn[rnd.randint(0, 2):]
                    else:
                        dn = revcomp(piece) if rnd.random() < 0.5 else piece
                        if rnd.random() < 0.1: dn = dn.lower()
                    recs.append(('p%d_%d' % (m, r), dn))
                if shape == 'star':               # whole genomes of most references: every child of the root holds the same number of hashes (ties for std::sort)
                    for r, ref in enumerate(refs):
                        if rnd.random() < 0.9:
                            dn = back(ref['g']) if protein else ref['g']
                            recs.append(('p%d_%d' % (m, r), revcomp(dn) if rnd.random() < 0.5 else dn))
                pool.append({'name': 'pool%d.fa' % (m + 1), 'fasta': ''.join('>%s\n%s\n' % x for x in recs)})
            for f in qf + pool:
                open(os.path.join(d, f['name']), 'w', newline='\n').write(f['fasta'])
            open(os.path.join(d, 'nodes.dmp'), 'w', newline='\n').write(nodes)
            open(os.path.join(d, 'names.dmp'), 'w', newline='\n').write(names)
            map_text = None
            if mapping and rnd.random() < 0.95:
                sep = rnd.choice(['\t', ' ', '\t', '  '])
                map_text = ''.join('%d%s%s%s\n' % (t, sep, nm, rnd.choice(['', '', ' '])) for t, nm in mapping)
                open(os.path.join(d, 'map.txt'), 'w', newline='\n').write(map_text)
            sk = ['sketch', '-o', 'q'] + (['-a'] if protein else []) + ([] if protein else ['-k', str(rnd.choice([8, 12, 16, 21]))]) + ['-s', str(rnd.choice([10, 30, 60, 120]))] + (['-i'] if individual else []) + (['-n'] if (not protein and rnd.random() < 0.2) else [])
            p = subprocess.run(['mash'] + sk + [f['name'] for f in qf], cwd=d, capture_output=True, text=True, timeout=120)
            if p.returncode != 0 or not os.path.exists(os.path.join(d, 'q.msh')): continue
            args = ['taxscreen'] + (['-m', 'map.txt'] if map_text is not None else []) + ['-t', '.', 'q.msh'] + [f['name'] for f in pool]
            try:
                p = subprocess.run(['mash'] + args, cwd=d, capture_output=True, text=True, timeout=60)
            except subprocess.TimeoutExpired:
                print('timeout'); continue
            if p.returncode != 0: continue
            rows.append({'id': len(rows), 'args': args, 'files': pool, 'qmsh': base64.b64encode(open(os.path.join(d, 'q.msh'), 'rb').read()).decode(), 'names': names, 'nodes': nodes, 'map': map_text, 'out': p.stdout})
            os.remove(os.path.join(d, 'q.msh'))
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    print(len(rows), 'rows', sum(1 for r in rows if r['out'].count('\n') > 3), 'with taxa')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 20261221)
