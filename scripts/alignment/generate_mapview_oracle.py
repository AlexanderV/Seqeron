#!/usr/bin/env python3
"""A3-9e4 oracle: `mapview` (perl script of MUMmer 3.23, Debian 3.23+dfsg-8).
Rows (JSON lines):
  {"kind":"mv","coords","utr","cds","args":[...],"figs":{name:text},"out","err","rc","side":{name:text},"commands":[...]}
        mapview <args> n.coords [utr.gff [cds.gff]] in a scratch directory: the FIG files it leaves (with -f pdf / ps the stub fig2dev that stands in for the real one records its command line and the FIG file it is given:
        "figs" are those, "commands" the command lines; the real fig2dev is not needed to compare the script), standard output / error, exit status, the side file of -g runs (<coords>coords).
  {"kind":"mv","ub":true,...}   hash-order dependent: the output of the script differs between two runs (perl randomises the order of keys: equal 5' UTR ends of two genes are listed in hash order) - recorded, not compared.
Inputs: coords files of `show-coords -rl` (nucmer and promer deltas of random related sequences; one or several references and queries), mgaps output (-g), GFF files of random genes on both strands
(UTR file = transcripts, CDS file = the coding parts); options -d -f -m -n -p -v -x1 -x2 -g -I -Ir.
Usage: python3 scripts/alignment/generate_mapview_oracle.py [seed=20261206] [n=330] | gzip > mapview_oracle.jsonl.gz"""
import json, os, random, shutil, stat, subprocess, sys

WORK = '/tmp/b09mv_%d' % os.getpid()
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
STUB = '''#!/usr/bin/env python3
import os, shutil, sys
log = os.environ['FIG2DEV_LOG']
args = sys.argv[1:]
open(log, 'a').write('fig2dev ' + ' '.join(args) + '\\n')
os.makedirs(os.path.dirname(log) + '/captured', exist_ok=True)
shutil.copy(args[-2], os.path.dirname(log) + '/captured/' + os.path.basename(args[-2]))
open(args[-1], 'w').write('stub\\n')
'''


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261206
    n_rows = int(sys.argv[2]) if len(sys.argv) > 2 else 330
    rnd = random.Random(seed)

    def rseq(k): return ''.join(rnd.choice('acgt') for _ in range(k))
    def rc(s): return s.translate(COMP)[::-1]
    def mutate(s, rate):
        o = []
        for c in s:
            r = rnd.random()
            if r < rate * .6: o.append(rnd.choice('acgt'))
            elif r < rate * .75: pass
            elif r < rate: o.append(c + rnd.choice('acgt'))
            else: o.append(c)
        return ''.join(o)
    def fasta(recs, width=60):
        return ''.join('>%s\n%s\n' % (nm, '\n'.join(sq[i:i + width] for i in range(0, len(sq), width))) for nm, sq in recs)
    def reset():
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
    def run(cmd, **kw): return subprocess.run(cmd, cwd=WORK, capture_output=True, **kw)
    def lat(b): return b.decode('latin-1')

    os.makedirs(WORK + '_bin', exist_ok=True)
    stub = WORK + '_bin/fig2dev'
    open(stub, 'w').write(STUB); os.chmod(stub, 0o755)

    def make_coords(nrefs, promer=False):
        refs = [('chr%s' % 'ABC'[i], rseq(rnd.randint(1500, 6000))) for i in range(nrefs)]
        qrys = []
        for k in range(rnd.randint(1, 5)):
            nm, base = rnd.choice(refs)
            a = rnd.randint(0, len(base) // 2); b = rnd.randint(a + 300, len(base))
            piece = mutate(base[a:b], rnd.choice([.005, .02, .05, .1]))
            if rnd.random() < .4: piece = rc(piece)
            qrys.append(('q%d' % k, piece))
        reset()
        open(WORK + '/ref.fa', 'w').write(fasta(refs)); open(WORK + '/qry.fa', 'w').write(fasta(qrys))
        tool = ['promer', '-p', 'p', '-l', '5', '-c', '10'] if promer else ['nucmer', '-p', 'p', '-l', '12', '-c', '30']
        if run(tool + ['ref.fa', 'qry.fa']).returncode != 0: return None
        p = run(['show-coords', '-rl', 'p.delta'])
        if p.returncode != 0: return None
        txt = lat(p.stdout)
        lines = [l for l in txt.split('\n')[5:] if l.strip()]
        if not lines: return None
        return txt, {nm: len(sq) for nm, sq in refs}

    def make_gff(lens, ref_ids):
        """UTR and CDS GFF files for random genes on the references"""
        utr, cds = [], []
        for rid in ref_ids:
            L = lens[rid]
            pos = rnd.randint(50, 400)
            for g in range(rnd.randint(1, 5)):
                nex = rnd.randint(1, 4)
                exons = []
                p = pos
                for _ in range(nex):
                    s = p + rnd.randint(5, 80); e = s + rnd.randint(30, 300); exons.append((s, e)); p = e + rnd.randint(30, 200)
                if exons[-1][1] > L - 20: break
                pos = exons[-1][1] + rnd.randint(-100, 400)      # genes may overlap a little
                strand = rnd.choice('+-')
                name = '%s_G%d' % (rid, g)
                # coding part: strip some of the first and last exon
                c = [list(x) for x in exons]
                c[0][0] += rnd.randint(0, min(40, c[0][1] - c[0][0] - 2))
                c[-1][1] -= rnd.randint(0, min(40, c[-1][1] - c[-1][0] - 2))
                def lines(ex, store):
                    seq = ex if strand == '+' else ex[::-1]
                    for i, (s, e) in enumerate(seq):
                        t = 'single-exon' if len(seq) == 1 else ('initial-exon' if i == 0 else ('last-exon' if i == len(seq) - 1 else 'internal-exon'))
                        store.append('%s\tsrc\t%s\t%d\t%d\t.\t%s\t.\t%s\n' % (rid, t, s, e, strand, name))
                lines(exons, utr); lines(c, cds)
        return ''.join(utr), ''.join(cds)

    rows = []

    def mv_row(coords, utr, cds, args, coordsname='n.coords', same=False):
        reset()
        open(WORK + '/' + coordsname, 'w').write(coords)
        files = [coordsname]
        if utr is not None: open(WORK + '/utr.gff', 'w').write(utr); files.append('utr.gff')
        if cds is not None and not same: open(WORK + '/cds.gff', 'w').write(cds); files.append('cds.gff')
        positional = [coordsname] + (['utr.gff'] if utr is not None else []) + (['utr.gff' if same else 'cds.gff'] if cds is not None else [])
        def go(env_extra):
            env = dict(os.environ, PATH=WORK + '_bin:' + os.environ['PATH'], FIG2DEV_LOG=WORK + '/fig2dev.log')
            env.update(env_extra)
            for f in os.listdir(WORK):
                if f not in files:
                    fp = os.path.join(WORK, f)
                    shutil.rmtree(fp) if os.path.isdir(fp) else os.remove(fp)
            p = run(['mapview'] + args + positional, env=env, timeout=60)
            figs, side = {}, {}
            cap = WORK + '/captured'
            for f in sorted(os.listdir(WORK)):
                if f in files or f in ('fig2dev.log',): continue
                if f == 'captured':
                    for g in sorted(os.listdir(cap)): figs[g] = open(os.path.join(cap, g), encoding='latin-1', newline='').read()
                    continue
                t = open(os.path.join(WORK, f), encoding='latin-1', newline='').read()
                if f.endswith('.fig'): figs[f] = t
                elif f.endswith('coords') and f != coordsname: side[f] = t
            cmds = open(WORK + '/fig2dev.log').read().split('\n')[:-1] if os.path.exists(WORK + '/fig2dev.log') else []
            return p, figs, side, cmds
        p, figs, side, cmds = go({})
        row = {'kind': 'mv', 'coords': coords, 'utr': utr, 'cds': cds, 'same': same, 'args': args, 'name': coordsname, 'figs': figs, 'out': lat(p.stdout), 'err': lat(p.stderr), 'rc': p.returncode, 'side': side, 'commands': cmds}
        base = (lat(p.stdout), figs, side, cmds)
        for env in ({}, {}, {}, {}, {'PERL_HASH_SEED': '5000', 'PERL_PERTURB_KEYS': '0'}, {'PERL_HASH_SEED': '987654321', 'PERL_PERTURB_KEYS': '0'}):
            q, f2, s2, c2 = go(env)
            if (lat(q.stdout), f2, s2, c2) != base:
                row['ub'] = True
                break
        return row

    ARGS = lambda: []
    for k in range(n_rows):
        multi = rnd.random() < .4
        promer = k % 11 == 5
        got = make_coords(rnd.choice([2, 3]) if multi else 1, promer)
        if got is None: continue
        coords, lens = got
        ref_ids = sorted({l.split()[-2] for l in coords.split('\n')[5:] if len(l.split()) >= 2})
        args = []
        if rnd.random() < .3: args += ['-n', str(rnd.choice([1, 2, 3, 5, 7, 10, 15]))]
        if rnd.random() < .25: args += ['-p', rnd.choice(['plot', 'x', 'sub_out'])]
        if rnd.random() < .3: args.append('-v')
        if rnd.random() < .15: args += ['-d', str(rnd.choice([100, 500, 2000, 50000]))]
        if rnd.random() < .12: args.append('-I')
        r = rnd.random()
        if r < .12: args += ['-f', 'pdf']
        elif r < .2: args += ['-f', 'ps']
        if rnd.random() < .1: args += ['-m', str(rnd.choice([0.5, 1.5, 2]))]
        utr = cds = None
        same = False
        if rnd.random() < .5:
            utr, cds = make_gff({**lens}, ref_ids)
            if not utr: utr = cds = None
            else:
                if rnd.random() < .2: same = True
                if rnd.random() < .3: args.append('-Ir')
        if utr is not None and not multi and rnd.random() < .3:
            args += ['-x1', str(rnd.randint(1, 1500)), '-x2', str(rnd.randint(1500, 5000))]
        elif rnd.random() < .06:
            args += ['-x1', str(rnd.randint(1, 1500)), '-x2', str(rnd.randint(1, 3000))]
        rows.append(mv_row(coords, utr, cds, args, same=same))
    # mgaps input (-g): clusters of the mummer | mgaps pipeline of random pairs
    for k in range(max(30, n_rows // 8)):
        reset()
        ref = rseq(rnd.randint(2000, 6000))
        qrys = []
        for i in range(rnd.randint(1, 4)):
            a = rnd.randint(0, len(ref) // 3); piece = mutate(ref[a:rnd.randint(a + 400, len(ref))], rnd.choice([.01, .04]))
            if rnd.random() < .5: piece = rc(piece)
            qrys.append(('q%d' % i, piece))
        open(WORK + '/ref.fa', 'w').write(fasta([('REFX', ref)])); open(WORK + '/qry.fa', 'w').write(fasta(qrys))
        mm = subprocess.run('mummer -mumreference -b -l 15 ref.fa qry.fa 2>/dev/null | mgaps -l %d -f .12 -s %d' % (rnd.choice([30, 60, 100]), rnd.choice([200, 600, 1000])), shell=True, cwd=WORK, capture_output=True)
        txt = lat(mm.stdout)
        if not any(l.strip() and not l.startswith(('>', '#')) for l in txt.split('\n')): continue
        args = ['-g', 'REFX']
        if rnd.random() < .5: args += ['-n', str(rnd.choice([1, 2, 4]))]
        if rnd.random() < .3: args.append('-v')
        if rnd.random() < .3: args += ['-p', 'mg']
        utr = cds = None
        if rnd.random() < .4:
            utr, cds = make_gff({'REFX': len(ref)}, ['REFX'])
            utr = utr or None; cds = cds or None
        rows.append(mv_row(txt, utr, cds, args))
    # special inputs
    SIMPLE = '/dev/null'
    for coords, args in [('', []), ('x\nNUCMER\n\n[S1]\n=\n', []), ('\n\n\n\n\n', ['-v'])]:
        rows.append(mv_row(coords, None, None, args))
    shutil.rmtree(WORK, ignore_errors=True); shutil.rmtree(WORK + '_bin', ignore_errors=True)
    for r in rows: print(json.dumps(r))
    kinds = {}
    for r in rows:
        key = r['kind'] + (':rc%d' % r['rc'] if r['rc'] else '') + (':ub' if r.get('ub') else '')
        kinds[key] = kinds.get(key, 0) + 1
    sys.stderr.write('rows=%d %s\n' % (len(rows), kinds))


main()
