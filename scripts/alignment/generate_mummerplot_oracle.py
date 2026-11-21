#!/usr/bin/env python3
"""A3-9e4 oracle: `mummerplot` (perl script 3.5 of MUMmer 3.23, Debian 3.23+dfsg-8).
Rows (JSON lines):
  {"kind":"mp","match","name","files":{name:text},"args":[...],"out","err","rc","outfiles":{name:text},"cmd":[argv of the gnuplot rendering call] | null}
        mummerplot <args> <name> in a scratch directory that holds the match file and `files` (the sequence files named in the first line of a delta file, the -R / -Q files); a stub `gnuplot` that prints
        "gnuplot 6.0 patchlevel 0" for --version and records the arguments of any other call (and a stub xclip) stands in for gnuplot, so the script text is compared, not the rendering (the real gnuplot 6.0 and
        the real fig2dev of mapview are installed on the machine that generated this file and were run by hand on samples: .ps / .png / .pdf come out); "outfiles" = every file the script wrote
        (.fplot .rplot .hplot .gp .filter); "cmd" = the arguments of the gnuplot call.
  {"kind":"mp","ub":true,...}   hash-order dependent: the output differs between runs (perl randomises the key order: layouts of sequences with equal lengths / offsets) - recorded, not compared.
Inputs: deltas of nucmer and promer (one or several references and queries, forward and reverse alignments), clusters (--nodelta), mummer -b -c -L output, show-tiling output; options -b --[no]color
-c --[no]coverage -f -l --fat -p -rv -r -q -R -Q -rport -qport -s --small --medium --large -S -t --x11 --postscript --png --title -x -y.
Usage: python3 scripts/alignment/generate_mummerplot_oracle.py [seed=20261207] [n=330] | gzip > mummerplot_oracle.jsonl.gz"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09mp_%d' % os.getpid()
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
GNUPLOT = '''#!/bin/sh
if [ "$1" = "--version" ]; then echo "gnuplot 6.0 patchlevel 0"; exit 0; fi
for a in "$@"; do printf '%s\\n' "$a" >> "$GP_LOG"; done
printf '\\n' >> "$GP_LOG"
exit 0
'''
XCLIP = '#!/bin/sh\nexit 0\n'


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261207
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

    binp = WORK + '_bin'
    os.makedirs(binp, exist_ok=True)
    for nm, txt in (('gnuplot', GNUPLOT), ('xclip', XCLIP)):
        open(binp + '/' + nm, 'w').write(txt); os.chmod(binp + '/' + nm, 0o755)

    rows = []

    def make_inputs(kind):
        """returns (match text, files, ref ids, qry ids)"""
        nref = rnd.choice([1, 1, 2, 3])
        refs = [('chr%s' % 'ABC'[i], rseq(rnd.randint(1500, 5000))) for i in range(nref)]
        qrys = []
        for k in range(rnd.randint(1, 5)):
            nm, base = rnd.choice(refs)
            a = rnd.randint(0, len(base) // 2); b = rnd.randint(a + 400, len(base))
            piece = mutate(base[a:b], rnd.choice([.005, .02, .05, .1]))
            if rnd.random() < .45: piece = rc(piece)
            qrys.append(('q%d' % k, piece))
        reset()
        open(WORK + '/ref.fa', 'w').write(fasta(refs)); open(WORK + '/qry.fa', 'w').write(fasta(qrys))
        files = {'ref.fa': fasta(refs), 'qry.fa': fasta(qrys)}
        if kind in ('delta', 'cluster', 'tiling'):
            tool = ['promer', '-p', 'p', '-l', '5', '-c', '10'] if kind == 'delta' and rnd.random() < .15 else ['nucmer', '-p', 'p', '-l', '12', '-c', '30']
            if kind == 'cluster': tool.append('--nodelta')
            if run(tool + ['ref.fa', 'qry.fa']).returncode != 0: return None
            if kind == 'cluster':
                txt = open(WORK + '/p.cluster', encoding='latin-1', newline='').read()
            else:
                txt = open(WORK + '/p.delta', encoding='latin-1', newline='').read()
            if kind != 'tiling':
                ls = txt.split('\n'); ls[0] = 'ref.fa qry.fa'; txt = '\n'.join(ls)       # the delta header names the files by their absolute path
            if kind == 'tiling':
                p = run(['show-tiling', 'p.delta'])
                txt = lat(p.stdout)
                if '>' not in txt or len(txt.split('\n')) < 3: return None
        else:    # mummer
            p = run(['mummer', '-mum', '-b', '-c', '-L', '-l', '14', 'ref.fa', 'qry.fa'])
            txt = lat(p.stdout)
            if not any(l.strip() and not l.startswith('>') for l in txt.split('\n')): return None
        if kind in ('delta', 'cluster') and not any(l[:1].isdigit() or l[:1] == ' ' for l in txt.split('\n')[3:]): return None
        return txt, files, [n for n, _ in refs], [n for n, _ in qrys]

    def mp_row(match, name, files, args):
        reset()
        open(WORK + '/' + name, 'w').write(match)
        for fn, t in files.items(): open(WORK + '/' + fn, 'w').write(t)
        keep = set(files) | {name}
        def go(env_extra):
            for f in os.listdir(WORK):
                if f not in keep:
                    fp = os.path.join(WORK, f)
                    shutil.rmtree(fp) if os.path.isdir(fp) else os.remove(fp)
            env = dict(os.environ, PATH=binp + ':' + os.environ['PATH'], GP_LOG=WORK + '/gp.log')
            env.update(env_extra)
            p = run(['mummerplot'] + args + [name], env=env, timeout=60)
            outs = {}
            for f in sorted(os.listdir(WORK)):
                if f in keep or f == 'gp.log': continue
                outs[f] = open(os.path.join(WORK, f), encoding='latin-1', newline='').read()
            cmd = open(WORK + '/gp.log').read().rstrip('\n').split('\n') if os.path.exists(WORK + '/gp.log') else None
            return p, outs, cmd
        p, outs, cmd = go({})
        row = {'kind': 'mp', 'match': match, 'name': name, 'files': files, 'args': args, 'out': lat(p.stdout), 'err': lat(p.stderr), 'rc': p.returncode, 'outfiles': outs, 'cmd': cmd}
        # perl randomises the order of hash keys per process: any difference between several runs (random seeds and two fixed ones) means the output depends on it
        base = (lat(p.stdout), lat(p.stderr), outs, cmd)
        envs = [{}, {}, {}, {}, {'PERL_HASH_SEED': '5000', 'PERL_PERTURB_KEYS': '0'}, {'PERL_HASH_SEED': '987654321', 'PERL_PERTURB_KEYS': '0'}]
        if '-l' in args or '--fat' in args: envs += [{}] * 10        # layouts are ordered by perl's hash order: more runs to see it
        for env in envs:
            q, o2, c2 = go(env)
            if (lat(q.stdout), lat(q.stderr), o2, c2) != base:
                row['ub'] = True
                break
        return row

    for k in range(n_rows):
        kind = rnd.choices(['delta', 'cluster', 'mummer', 'tiling'], [.6, .12, .16, .12])[0]
        got = make_inputs(kind)
        if got is None: continue
        match, files, refids, qids = got
        name = {'delta': 'x.delta', 'cluster': 'x.cluster', 'mummer': 'x.out', 'tiling': 'x.tiling'}[kind]
        args = []
        if rnd.random() < .3: args += ['-p', rnd.choice(['plot', 'x', 'sub_out'])]
        r = rnd.random()
        if r < .3: args.append('--postscript')
        elif r < .55: args.append('--png')
        elif r < .62: args += ['-t', rnd.choice(['postscript', 'png'])]
        elif r < .66: args += ['--png', '--postscript']
        r = rnd.random()
        if r < .15: args.append('--medium')
        elif r < .3: args.append('--large')
        elif r < .35: args += ['-s', 'large']
        elif r < .38: args += ['--large', '--small']
        if rnd.random() < .2: args.append('--color')
        elif rnd.random() < .08: args.append('--nocolor')
        if rnd.random() < .12: args.append('-c')
        elif rnd.random() < .05: args.append('--nocoverage')
        if rnd.random() < .1: args += ['-b', str(rnd.choice([0, 20, 100, 500]))]
        if rnd.random() < .1: args.append('-S' if kind == 'delta' else '-f')
        if rnd.random() < .12: args.append('-f')
        if rnd.random() < .08: args.append('-rv')
        if rnd.random() < .12: args += ['--title', rnd.choice(['MUMmer plot', 'a b', 'x'])]
        if rnd.random() < .12: args += ['-x', rnd.choice(['[1:1000]', '[10,2000]', '[0:5000]'])]
        if rnd.random() < .1: args += ['-y', rnd.choice(['[1:900]', '[5,3000]'])]
        if rnd.random() < .06: args += ['-rport', '7000']
        if rnd.random() < .06: args += ['-qport', '7001']
        if rnd.random() < .12 and refids: args += ['-r', rnd.choice(refids)]
        if rnd.random() < .12 and qids: args += ['-q', rnd.choice(qids)]
        # -R / -Q files: the multi-FASTA files themselves or lists of ids, lengths and directions
        if rnd.random() < .25:
            lines = []
            for nm in refids:
                ln = len([x for x in files['ref.fa'].split('>%s\n' % nm)[1].split('>')[0].split('\n') if x]) and sum(len(x) for x in files['ref.fa'].split('>%s\n' % nm)[1].split('>')[0].split('\n'))
                lines.append('%s %d %s\n' % (nm, ln, rnd.choice(['+', '-', ''])))
            files['rlist.txt'] = ''.join(lines)
            args += ['-R', 'rlist.txt' if rnd.random() < .6 else 'ref.fa']
        if rnd.random() < .25:
            args += ['-Q', rnd.choice(['qry.fa', 'qry.fa', 'qlist.txt'])]
            if args[-1] == 'qlist.txt':
                files['qlist.txt'] = ''.join('%s %d %s\n' % (nm, sum(len(x) for x in files['qry.fa'].split('>%s\n' % nm)[1].split('>')[0].split('\n')), rnd.choice(['+', '-', ''])) for nm in qids)
        if rnd.random() < .2: args += ['-l']
        elif rnd.random() < .05: args += ['--fat']
        rows.append(mp_row(match, name, files, args))
    for match, name in [('', 'e.delta'), ('foo bar\nbaz\n', 'f.delta'), ('> x\n  1 2 3 4\n', 'm.out')]:
        rows.append(mp_row(match, name, {}, []))
    shutil.rmtree(WORK, ignore_errors=True); shutil.rmtree(binp, ignore_errors=True)
    for r in rows: print(json.dumps(r))
    kinds = {}
    for r in rows:
        key = r['kind'] + (':rc%d' % r['rc'] if r['rc'] else '') + (':ub' if r.get('ub') else '')
        kinds[key] = kinds.get(key, 0) + 1
    sys.stderr.write('rows=%d %s\n' % (len(rows), kinds))


main()
