#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A2-14d, B09 F105) oracle: the alignment files of ClustalW 2.1 - CLUSTAL, GCG (MSF), GDE, PHYLIP, PIR, NEXUS, FASTA - byte for byte.

Seqeron API: ClustalWCommandLine.Parse(args) -> SequenceAligner.RunClustalW(...) -> ClustalWAlignmentWriter.Write(alignment, options, structure).

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7) `-INFILE=x.fa -TYPE=DNA|PROTEIN -ALIGN -OUTPUT=<format> -OUTFILE=o.out` with -OUTORDER=INPUT|ALIGNED (and `-CONVERT`: the unaligned input written in the format, 301 sets), -CASE=LOWER|UPPER, -SEQNOS=ON, -SEQNO_RANGE=ON, -RANGE=m,n,
sequence names of every length (1 - 40 characters, names that already carry /start-end), the alignment option variants of the earlier ClustalW oracles, and -GAPOPEN ... flags (they also test the command-line reader).
Output TSV (tests/.../TestData/Alignment/clustalw_output_oracle.tsv):  id  type  format  names(;)  seqs(;)  args(space separated)  file(\\n / \\t / \\\\ escaped)
Usage: python3 scripts/alignment/generate_clustalw_output_oracle.py OUTFILE [setsPerFormat] [seed]      (needs `clustalw` on PATH)
"""
import importlib.util, os, random, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


dna = load('generate_clustalw_msa_full_oracle')
prot = load('generate_clustalw_protein_msa_oracle')
FLAGS = dict(prot.FLAGS)
FLAGS['transweight'] = 'TRANSWEIGHT'
FORMATS = ['CLUSTAL', 'GCG', 'GDE', 'PHYLIP', 'PIR', 'NEXUS', 'FASTA']
EXT = {'CLUSTAL': 'aln', 'GCG': 'msf', 'GDE': 'gde', 'PHYLIP': 'phy', 'PIR': 'pir', 'NEXUS': 'nxs', 'FASTA': 'fasta'}


def esc(t):
    return t.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t')


def names_for(rnd, k):
    style = rnd.random()
    out, seen = [], set()
    while len(out) < k:
        if style < 0.35: n = f's{len(out):03d}'
        elif style < 0.55: n = ''.join(rnd.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-') for _ in range(rnd.randint(1, 12)))
        elif style < 0.8: n = ''.join(rnd.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_') for _ in range(rnd.randint(8, 40)))
        else:
            n = ''.join(rnd.choice('abcdefghijklmnopqrstuvwxyz0123456789_') for _ in range(rnd.randint(3, 10)))
            if rnd.random() < 0.7:
                a = rnd.randint(1, 50); n += f'/{a}-{a + rnd.randint(5, 60)}'
        if n[0] in '-.' or n in seen: continue
        seen.add(n); out.append(n)
    return out


def options(rnd, typ):
    m = dna if typ == 'dna' else prot
    m.rnd.seed(rnd.randint(0, 2 ** 31))
    o = m.options()
    for k in ('ktuple', 'window', 'topdiags', 'pairgap', 'novgap', 'endgaps', 'hgapresidues'): o.pop(k, None)
    return o


def main(outfile, per, seed, convert=301):
    rnd = random.Random(seed)
    rows = []
    i = 0
    plan = [(fmt, j, False) for fmt in FORMATS for j in range(per)] + [(FORMATS[j % 7], j, True) for j in range(convert)]
    for fmt, j, conv in plan:
        if True:
            typ = 'dna' if j % 2 == 0 else 'protein'
            dna.rnd.seed(rnd.randint(0, 2 ** 31)); prot.rnd.seed(rnd.randint(0, 2 ** 31))
            maxlen = rnd.choice([30, 50, 80, 130, 200])
            seqs = dna.make_set(maxlen, 'ACGT', 2, rnd.choice([3, 5, 8])) if typ == 'dna' else prot.make_set(maxlen, 2, rnd.choice([3, 5, 8]))
            if len(seqs) < 2: seqs = seqs + seqs[:1]
            names = names_for(rnd, len(seqs))
            o = options(rnd, typ)
            args = [f'-TYPE={typ.upper()}', '-CONVERT' if conv else '-ALIGN', f'-OUTPUT={fmt}']
            if not conv: args += [f'-{FLAGS[k]}' + ('' if v is None else f'={v}') for k, v in o.items()]
            if rnd.random() < 0.5: args.append('-OUTORDER=' + rnd.choice(['INPUT', 'ALIGNED']))
            if rnd.random() < 0.4: args.append('-CASE=' + rnd.choice(['LOWER', 'UPPER']))
            if rnd.random() < 0.4: args.append('-SEQNOS=' + rnd.choice(['ON', 'OFF']))
            if rnd.random() < 0.35: args.append('-SEQNO_RANGE=' + rnd.choice(['ON', 'OFF']))
            with tempfile.TemporaryDirectory() as d:
                with open(os.path.join(d, 'x.fa'), 'w') as f:
                    for n, s in zip(names, seqs): f.write(f'>{n}\n{s}\n')
                if rnd.random() < 0.3:
                    # -RANGE: from the alignment length (taken from a first run)
                    subprocess.run(['clustalw', '-INFILE=x.fa'] + args + ['-OUTFILE=probe.out'], capture_output=True, text=True, cwd=d, timeout=30)
                    probe = open(os.path.join(d, 'probe.out')).read() if os.path.exists(os.path.join(d, 'probe.out')) else ''
                    L = max((len(s) for s in seqs), default=10) + 10
                    a = rnd.randint(1, max(1, L // 2)); b = a + rnd.randint(0, L)
                    args.append(f'-RANGE={a},{b}')
                out = os.path.join(d, 'o.' + EXT[fmt])
                r = subprocess.run(['clustalw', '-INFILE=x.fa'] + args + [f'-OUTFILE=o.{EXT[fmt]}'], capture_output=True, text=True, cwd=d, timeout=30)
                text = open(out).read() if os.path.exists(out) and os.path.getsize(out) > 0 else ''
            if not text:
                continue
            rows.append((i, typ, fmt, names, seqs, args, text)); i += 1
    with open(outfile, 'w') as f:
        f.write('id\ttype\tformat\tnames\tseqs\targs\tfile\n')
        for (i, typ, fmt, names, seqs, args, text) in rows:
            f.write(f"{i}\t{typ}\t{fmt}\t{';'.join(names)}\t{';'.join(seqs)}\t{' '.join(args)}\t{esc(text)}\n")
    from collections import Counter
    print(len(rows), 'sets', Counter(r[2] for r in rows))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 330, int(sys.argv[3]) if len(sys.argv) > 3 else 20261122)
