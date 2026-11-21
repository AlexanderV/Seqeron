#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/kalign_formats_oracle.tsv (B09 audit A2-20b): whole Kalign 3.3.5 command lines with the formats msf / clu, MSF / Clustal inputs, several inputs and the non-quiet
standard output (banner, log lines, warnings), run with the installed `kalign` under the time() shim of faketime.c (the run date is a chosen epoch, TZ=UTC).

Columns: id, args ('|' separated; the files in / in2 / in3 are given by the next columns), in, in2, in3 (text, '\\n' escaped), epoch, rc, stdout, out (the -o file, '' when none).
The text after `CPU Time: ` of the stopwatch lines is replaced by `~` (it is the CPU time of the run).
Run: python3 scripts/alignment/generate_kalign_formats_oracle.py
"""
import os, random, shutil, subprocess, tempfile, re

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'kalign_formats_oracle.tsv')
rng = random.Random(20261301)
DNA, PROT = 'ACGT', 'ACDEFGHIKLMNPQRSTVWY'
EPOCHS = [1795279369, 1709629380, 1735689600, 1767225540, 1800000000, 1772582400]


def fam(alpha, n, L, indel=0.5):
    base = ''.join(rng.choice(alpha) for _ in range(L))
    out = []
    for _ in range(n):
        s = ''.join(c if rng.random() > 0.1 else rng.choice(alpha) for c in base)
        if rng.random() < indel:
            p = rng.randrange(len(s)); s = s[:p] + s[p + rng.randrange(1, 4):]
        if rng.random() < indel:
            p = rng.randrange(len(s)); s = s[:p] + ''.join(rng.choice(alpha) for _ in range(rng.randrange(1, 4))) + s[p:]
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


def add(args, a, b='', c='', epoch=None):
    cases.append((args, a, b, c, epoch if epoch is not None else rng.choice(EPOCHS)))


def names_for(n, spaces=True):
    pool = ['seq%d' % (i + 1) for i in range(n)]
    if rng.random() < 0.4:
        pool = [rng.choice(['alpha', 'b', 'a long sequence name' if spaces else 'a_long_sequence_name', 'x_1', 'gi|123|ref|NP_0.1|', 'tab']) + str(i) for i in range(n)]
    if rng.random() < 0.2:
        pool[0] = 'n' * rng.choice([100, 255, 256, 257, 300])
    return pool


tmp = tempfile.mkdtemp()
shim = os.path.join(tmp, 'faketime.so')
subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)


def run(args, files, epoch, cwd):
    shutil.rmtree(cwd)
    os.makedirs(os.path.join(cwd, 'sub'))
    for k, v in files.items():
        open(os.path.join(cwd, k), 'w', newline='\n').write(v)
    env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch), TZ='UTC')
    p = subprocess.run(['kalign'] + args, cwd=cwd, capture_output=True, text=True, stdin=subprocess.DEVNULL, env=env, timeout=60)
    out = ''
    for cand in ('out', os.path.join('sub', 'dir.out')):
        if os.path.exists(os.path.join(cwd, cand)):
            out = open(os.path.join(cwd, cand)).read()
    return p, out


work = tempfile.mkdtemp()
# 1. output formats on random families
for k in range(60):
    alpha = DNA if k % 2 == 0 else PROT
    n = rng.randrange(2, 9)
    L = rng.choice([20, 45, 60, 61, 90, 120, 130, 181])
    seqs = fam(alpha, n, L)
    text = fasta(seqs, names_for(n), rng.choice([None, 60, 25]))
    f = rng.choice(['msf', 'clu', 'fasta', 'fa', 'MSF', 'clustal', 'xmsfx', 'cluster'])
    extra = rng.choice([[], [], ['--type', 'dna' if alpha == DNA else 'protein'], ['-n', '1'], ['-n', '2'], ['--gpo', '5']])
    if k % 3 == 0:
        add(['-i', 'in', '-q', '-f', f] + extra, text)           # standard output
    elif k % 3 == 1:
        add(['-i', 'in', '-q', '--format', f, '-o', 'out'] + extra, text)
    else:
        add(['-q', '-f', f, '-i', 'in', '-o', 'sub/dir.out'] + extra, text)
# aligned width multiples of 60 and sequences of the same length
for L in (59, 60, 61, 119, 120, 121):
    s = fam(DNA, 3, L, 0.0)
    add(['-i', 'in', '-q', '-f', 'msf'], fasta(s))
    add(['-i', 'in', '-q', '-f', 'clu'], fasta(s))
# 2. non-quiet standard output
for k in range(40):
    alpha = DNA if k % 2 == 0 else PROT
    n = rng.randrange(2, 7)
    seqs = fam(alpha, n, rng.choice([30, 70, 100]))
    text = fasta(seqs, names_for(n))
    extra = rng.choice([[], [], ['-n', '1'], ['-n', '2'], ['-f', 'clu'], ['-f', 'msf'], ['--type', 'rna' if alpha == DNA else 'protein']])
    if k % 4 == 3:
        add(['-i', 'in', '-o', 'out'] + extra, text)
    else:
        add(['-i', 'in'] + extra, text)
add(['-i', 'in'], fasta(['ACGTACGT', 'ACGTACCT', 'ACGAACGT']))                            # same length, no gaps: the warning block
add(['-i', 'in'], fasta(['ACGT-ACGT', 'ACGTACCT', 'AC--ACGT']))                           # gaps, different lengths
add(['-i', 'in'], fasta(['ACGT-ACGT', 'ACGTACCT-', 'AC--ACGT-']))                         # aligned input (the same width)
add(['-i', 'in'], fasta(['AC.GT*ACGT', 'ACGT%ACCT', 'AC--ACGT1']))                        # punctuation, digits
add(['-i', 'in'], '>a\nACGT\n>b\n\n>c\nACGAC\n>d\nACGTT\n')                                # an empty sequence
add(['-i', 'in'], '>a\n\n>b\n\n>c\nACGAC\n>d\nACGTT\n>e\n')                                # three empty sequences
add(['-i', 'in'], '>a\nACGT\n>b\n---\n>c\nACGAC\n')                                       # gap-only sequence
add(['-i', 'in'], '>a\n>b\n')                                                              # no letters: alphabet tie
add(['-i', 'in'], '>a\n123\n>b\n456\n')
add(['-i', 'in'], 'hello\nthere\n')                                                         # no format detected
add(['-i', 'in'], '')                                                                       # empty input
add(['-i', 'in'], '>\n')
add(['-i', 'in'], '>a\nACGT\n')                                                             # one sequence
add(['-h'], '')
add(['--showw'], '')
add(['-n', '0', '-i', 'in'], fasta(['ACGT', 'ACGA']))
add(['-q', '-n', '0', '-i', 'in'], fasta(['ACGT', 'ACGA']))
add(['--set', '1', '-i', 'in'], fasta(['ACGTA', 'ACGAA']))
add(['-V'], '')
add([], '')
add(['-i', 'missing', '-q'], '')
add(['-q', '-i', 'in', 'missing'], fasta(['ACGTA', 'ACGAA']))
# 3. inputs in MSF / Clustal format: written by kalign itself, then perturbed
inputs = []
for k in range(40):
    alpha = DNA if k % 2 == 0 else PROT
    n = rng.randrange(2, 7)
    seqs = fam(alpha, n, rng.choice([30, 75, 130]))
    text = fasta(seqs, names_for(n, False))   # a name with a blank would be cut by the readers and its tail read as residues (letters outside the alphabet: uninitialised codes in the program)
    for fmt in ('msf', 'clu'):
        p, o = run(['-i', 'in', '-q', '-f', fmt, '-o', 'out'], {'in': text}, 1795279369, work)
        if p.returncode == 0:
            inputs.append((fmt, o))
for fmt, o in inputs:
    ex = rng.choice([[], [], ['-f', 'fasta'], ['-f', 'clu'], ['-f', 'msf'], ['-n', '1']])
    q = rng.choice([['-q'], []])
    if q == []:
        ex = [e for e in ex if e != '-n' and e != '1']
    add(['-i', 'in'] + q + ex, o)
# clustal files as other programs write them (conservation line, numbers, blank lines, tabs) and MSF without Name lines
clu = 'CLUSTAL W (1.83) multiple sequence alignment\n\n\nseqA          ACGT-ACGTAAGT\nseqB          ACGTAACGTAAGT\nseqC          ACG--ACGT-AGT\n               ***  ***  ***\n\n\nseqA          ACGTAC\nseqB          ACGTAC\nseqC          AC--AC\n\n'
add(['-i', 'in', '-q'], clu)
add(['-i', 'in'], clu)
add(['-i', 'in', '-q'], clu.replace('seqA          ', 'seqA   123  '))
add(['-i', 'in', '-q'], 'CLUSTAL multiple sequence alignment\n\nacgtacgt\ngattaca\n\nacgtacgt ACGT\ngattaca ACTT\n')
add(['-i', 'in', '-q'], 'CLUSTAL W\n\nseqA ACGTAC\nseqB ACGTTC\nseqC ACGAAC\n\nseqA ACGT\nseqB ACGT\nseqC ACGA\n')
add(['-i', 'in', '-q'], 'CLUSTAL W\n\nseqA\tACGTAC\nseqB ACGTTC\n')                       # a tab ends the line
add(['-i', 'in', '-q'], 'CLUSTAL W\n\n' + 'x' * 255 + ' ACGT\n' + 'y' * 260 + ' ACGA\n')   # names of 255 / 260 characters
add(['-i', 'in', '-q'], 'CLUSTAL W\n\n' + 'x' * 255 + 'ACGT ACGA\nsecond ACGTAT\n')
add(['-i', 'in', '-q'], 'CLUSTAL W\n\nacgtacgt\nseqB ACGTTC\nseqC ACCTTC\n')               # a name without residues: read as residues
msf = '!!NA_MULTIPLE_ALIGNMENT 1.0\n\n stdout  MSF: 13  Type: N  November 21, 2026 16:42  Check: 4658  ..\n\n Name: a         Len:     13  Check: 5668  Weight: 1.00\n Name: b         Len:     13  Check: 4528  Weight: 1.00\n Name: c_desc    Len:     13  Check: 4462  Weight: 1.00\n\n//\n\na            ACGTACGTAAGT-\nb            ACGTA-GTAAGTA\nc_desc       ACGAC-GTAAGT-\n\n\n'
add(['-i', 'in', '-q'], msf)
add(['-i', 'in'], msf)
add(['-i', 'in', '-q'], msf.replace('Name: ', 'Nam: '))                                      # no names: no sequences
add(['-i', 'in', '-q'], msf.replace('Len:', 'Length'))
add(['-i', 'in', '-q'], msf.replace('//\n', ''))                                            # no // line
add(['-i', 'in', '-q'], msf + 'a            ACGT\nb            ACGA\nc            ACGC\nextra        ACGT\n')
add(['-i', 'in', '-q'], msf.replace('a            ACGT', 'a ACGT'))
add(['-i', 'in', '-q'], msf.replace('b            ACGTA-GTAAGTA', 'ACGTA-GTAAGTAAAAAAAAAAA'))
# 4. several inputs
d1 = fasta(fam(DNA, 3, 60)); d2 = fasta(fam(DNA, 4, 70), ['m1', 'm2', 'm3', 'm4']); p1 = fasta(fam(PROT, 3, 60))
for q in (['-q'], []):
    add(['in', 'in2', '-o', 'out'] + q, d1, d2)
    add(['in', 'in2', 'in3', '-o', 'out'] + q, d1, d2, d1)
    add(['in', 'in2', '-o', 'out'] + q, d1, p1)                                             # different alphabets
    add(['in', 'in2', '-o', 'out'] + q, '>a\nACGTAC\n', d2)                                 # first file has one sequence
    add(['in', 'in2', '-o', 'out'] + q, d1, '>a\nACGTAC\n')
    add(['in', 'in2', '-o', 'out'] + q, d1, 'nothing here\n')
    add(['in', 'in2', '-o', 'out'] + q, d1, '')
    add(['in', 'in2', '-o', 'out'] + q, '', d1)
c1 = run(['-i', 'in', '-q', '-f', 'clu', '-o', 'out'], {'in': d1}, 1, work)[1]
m1 = run(['-i', 'in', '-q', '-f', 'msf', '-o', 'out'], {'in': d2}, 1, work)[1]
add(['in', 'in2', '-o', 'out', '-q'], c1, m1)
add(['in', 'in2', '-o', 'out', '-q', '-f', 'msf'], c1, d2)
add(['in', 'in2', '-o', 'out'], c1, m1)
add(['-i', 'in', 'in2', '-o', 'out', '-q'], d1, d2)

rows = []
bad_cpu = 0
for k, (args, a, b, c, epoch) in enumerate(cases):
    files = {'in': a}
    if b:
        files['in2'] = b
    if c:
        files['in3'] = c
    p, out = run(args, files, epoch, work)
    stdout = re.sub(r'CPU Time: .*', 'CPU Time: ~', p.stdout)   # the stopwatch lines are not reproducible (clock() of the run)
    esc = lambda s: s.replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n')
    rows.append(['F%03d' % (k + 1), '|'.join(args), esc(a), esc(b), esc(c), str(epoch), str(p.returncode), esc(stdout), esc(out)])
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', newline='\n') as fh:
    fh.write('id\targs\tin\tin2\tin3\tepoch\trc\tstdout\tout\n')
    for r in rows:
        fh.write('\t'.join(r) + '\n')
print(len(rows), 'rows;', sum(1 for r in rows if r[6] != '0'), 'with rc != 0;')
