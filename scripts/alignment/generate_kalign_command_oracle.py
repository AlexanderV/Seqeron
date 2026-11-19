#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/kalign_command_oracle.tsv (B09 audit A2-20): whole Kalign 3.3.5 command lines run with the installed `kalign`.

Columns: id, args (the argument list, '|' separated; the files in/in2 are given by the next columns), in, in2 (FASTA text, '\\n' escaped), rc, out (the written -o file, or stdout when no -o is named; escaped).
Only runs without a timestamp in their output are recorded (the serial recursion `-n 1` logs "Aligned ..." lines with the time of day to standard output, so -n 1 runs always use -o).
Run: python3 scripts/alignment/generate_kalign_command_oracle.py
"""
import os, random, subprocess, tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'kalign_command_oracle.tsv')
rng = random.Random(777)
DNA, PROT = 'ACGT', 'ACDEFGHIKLMNPQRSTVWY'


def fam(alpha, n, L):
    base = ''.join(rng.choice(alpha) for _ in range(L))
    out = []
    for _ in range(n):
        s = ''.join(c if rng.random() > 0.1 else rng.choice(alpha) for c in base)
        if rng.random() < 0.5:
            p = rng.randrange(len(s)); s = s[:p] + s[p + 2:]
        out.append(s)
    return out


def fasta(seqs, names=None, width=None):
    t = ''
    for i, s in enumerate(seqs):
        n = names[i] if names else 'seq%d' % (i + 1)
        t += '>' + n + '\n'
        if width:
            for p in range(0, len(s), width):
                t += s[p:p + width] + '\n'
        else:
            t += s + '\n'
    return t


cases = []


def add(args, a, b=''):
    cases.append((args, a, b))


dna = fasta(fam(DNA, 5, 70), ['seq one description', 'two', 'three|x', 'four', 'five'])
prot = fasta(fam(PROT, 6, 90), None, 25)
dna_long = fasta(fam(DNA, 3, 150))
for t in ['dna', 'rna', 'internal', '-', 'protein', 'divergent', 'RNA', 'xrnax']:
    for form in (['-i', 'in', '-o', 'out', '-q'], ['-q', '-i', 'in', '--outfile=out']):
        args = form + ([] if t == '-' else ['--type', t])
        add(args, dna)
        add(args, prot)
for g in (['--gpo', '3'], ['--gpe', '2'], ['--tgpe', '1.5'], ['--gpo', '4', '--gpe', '1', '--tgpe', '0'], ['-gpo', '12', '-gpe=3'], ['--gpo', 'abc'], ['--gpe', '1e1', '--tgpe', '.5x']):
    add(['-i', 'in', '-o', 'out', '-q'] + g, dna)
    add(['-i', 'in', '-o', 'out', '-q', '--type', 'protein'] + g, prot)
for n in ('1', '2', '3', '8', '0', '-2', 'x'):
    add(['-i', 'in', '-o', 'out', '-q', '-n', n], dna)
    add(['-i', 'in', '-o', 'out', '-q', '--nthreads=' + n], prot)
add(['-n4', '-i', 'in', '-o', 'out', '-q'], dna)
add(['-i', 'in', '-o', 'out', '-q', '-f', 'fasta'], dna)
add(['-i', 'in', '-o', 'out', '-q', '--format', 'fa'], dna)
add(['-i', 'in', '-o', 'out', '-q', '--format', 'xyz'], dna)
add(['-q', '-V'], '')
add(['-q', '--version'], '')
add(['-q', '-v', '-i', 'in'], dna)
add(['-i', 'in', '--ty', 'dna', '-o', 'out', '-q'], dna)       # unambiguous abbreviation of --type
add(['-i', 'in', '--gp', '3', '-o', 'out', '-q'], dna)         # ambiguous: --gpo / --gpe
add(['-i', 'in', '--bogus', '-o', 'out', '-q'], dna)
add(['-i', 'in', '-o', 'out', '-q', '--type'], dna)
add(['in', 'in2', '-o', 'out', '-q'], dna, dna_long)             # two inputs merged
add(['in', '-o', 'out', '-q'], prot)                             # input as a plain argument
add(['-i', 'missing', '-o', 'out', '-q'], dna)
add(['-i', 'in', '-o', 'out', '-q', '--set', '3'], dna)
add(['-i', 'in', '-o', 'out', '-q', '--showw'], dna)
add(['-i', 'in', '-o', 'out', '-q', '--type', 'dna'], prot)      # type / alphabet mismatch
add(['-i', 'in', '-o', 'out', '-q'], '>a\nACGT\n')              # one sequence only
add(['-i', 'in', '-o', 'out', '-q'], '>a\nACGTAC\n>b\n---\n>c\nACGAC\n')
add(['-i', 'in', '-o', 'out', '-q', '-n', '4'], fasta(fam(DNA, 4, 130), None, 60))
add(['-i', 'in', '-o', 'out', '-q'], fasta(fam(DNA, 3, 121), ['a b  c', 'd', 'e']))   # 60-column wrapping of the output, names with blanks
add(['-i', 'in', '-q'], dna)                                     # to standard output, 4 threads
add(['-i', 'in', '-q', '-n', '2'], prot)

rows = []
with tempfile.TemporaryDirectory() as d:
    for k, (args, a, b) in enumerate(cases):
        for f in ('in', 'in2', 'out'):
            if os.path.exists(os.path.join(d, f)):
                os.remove(os.path.join(d, f))
        open(os.path.join(d, 'in'), 'w', newline='\n').write(a)
        if b:
            open(os.path.join(d, 'in2'), 'w', newline='\n').write(b)
        p = subprocess.run(['kalign'] + args, cwd=d, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=60)
        out = ''
        if os.path.exists(os.path.join(d, 'out')):
            out = open(os.path.join(d, 'out')).read()
        else:
            out = p.stdout
        esc = lambda s: s.replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n')
        rows.append(['C%03d' % (k + 1), '|'.join(args), esc(a), esc(b), str(p.returncode), esc(out)])
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', newline='\n') as fh:
    fh.write('id\targs\tin\tin2\trc\tout\n')
    for r in rows:
        fh.write('\t'.join(r) + '\n')
print(len(rows), 'rows;', sum(1 for r in rows if r[4] != '0'), 'with rc != 0')
