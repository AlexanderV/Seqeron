#!/usr/bin/env python3
"""ALIGN-LOCAL-001 oracle (audit A4-9): the documented EDGE CASES of EMBOSS dotmatcher 6.6.0 vs DotMatcher.Compute - the cases the 600 rows of generate_dotmatcher_oracle.py never contain.

Case kinds (>= 100 rows each):
  lower      lower-case and mixed-case DNA / protein (explicit matrix)
  ambig      IUPAC ambiguity codes (DNA: N R Y K M S W B D H V U; protein: B Z X and the letters J O U that are not in the matrix, score 0)
  star       `*` (stop) inside protein sequences and inside DNA
  gap        gap characters - . ~ inside the sequences (removed by the reader before any coordinate; positions count the cleaned sequences)
  default    NO -matrixfile: the program picks EDNAFULL / EBLOSUM62 by the type of the first sequence (DNA/DNA, protein/protein, DNA/protein, protein/DNA)
  (a call the program refuses with "Sequence is not nucleic" - the first sequence is nucleic and the second is not - is recorded as ERR:notnucleic)
  overhang   a window as long as, or up to 8 residues longer than, one sequence (the program reads the zero terminator / allocation slack beyond the end; 9+ residues crash it)
TSV: case  matrix ('-' = default)  W  T  seq1  seq2  segments ('x1:y1:x2:y2' separated by spaces, '-' = none), like dotmatcher_oracle.tsv; a row whose call fails (a crash) is skipped.
Usage: python3 scripts/alignment/generate_dotmatcher_edge_oracle.py [per_case=120] > dotmatcher_edge_oracle.tsv
"""
import os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261201)
AA = 'ACDEFGHIKLMNPQRSTVWY'


def mutate(s, rate, alpha):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.6: out.append(rnd.choice(alpha))
        elif r < rate * 0.8: pass
        elif r < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def pair(alpha, lo=14, hi=90):
    base = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(lo, hi)))
    if rnd.random() < 0.75:
        return mutate(base, rnd.choice([0.05, 0.15, 0.3]), alpha), mutate(base, rnd.choice([0.05, 0.15, 0.3]), alpha)
    return base, ''.join(rnd.choice(alpha) for _ in range(rnd.randint(lo, hi)))


def sprinkle(s, chars, rate):
    return ''.join(rnd.choice(chars) if rnd.random() < rate else c for c in s)


def insert(s, chars, rate):
    out = []
    for c in s:
        if rnd.random() < rate: out.append(rnd.choice(chars))
        out.append(c)
    return ''.join(out)


def refold(s):
    mode = rnd.choice(['lower', 'mixed', 'mixed'])
    return s.lower() if mode == 'lower' else ''.join(c.lower() if rnd.random() < 0.5 else c for c in s)


def make(case):
    dna = rnd.random() < 0.5
    alpha = 'ACGT' if dna else AA
    a, b = pair(alpha)
    mat = rnd.choice(['EDNAFULL', 'EDNAFULL', 'EDNASIMPLE']) if dna else rnd.choice(['EBLOSUM62', 'EBLOSUM62', 'EPAM250', 'EBLOSUM30'])
    w = rnd.randint(3, min(15, len(a) - 1, len(b) - 1))
    if case == 'lower':
        a, b = refold(a), refold(b) if rnd.random() < 0.8 else b
    elif case == 'ambig':
        codes = 'NRYKMSWBDHVU' if dna else 'BZXJOU'
        a, b = sprinkle(a, codes, 0.12), sprinkle(b, codes, 0.12)
        if rnd.random() < 0.3: a, b = refold(a), refold(b)
    elif case == 'star':
        a, b = sprinkle(a, '*', 0.08), sprinkle(b, '*', 0.08)
    elif case == 'gap':
        a, b = insert(a, '-.~', 0.12), insert(b, '-.~', 0.12)
    elif case == 'default':
        kind = rnd.choice(['dd', 'pp', 'dp', 'pd'])
        a = pair('ACGT' if kind[0] == 'd' else AA)[0]
        b = pair('ACGT' if kind[1] == 'd' else AA)[0]
        mat = '-'
        w = rnd.randint(3, min(15, len(a) - 1, len(b) - 1))
    elif case == 'overhang':
        # the window is as long as the shorter-first sequence or up to 8 longer; the other sequence is longer than the window
        short = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(6, 20)))
        long_ = mutate(short, 0.1, alpha) + ''.join(rnd.choice(alpha) for _ in range(rnd.randint(10, 40)))
        if rnd.random() < 0.5: a, b = short, long_
        else: a, b = long_, short
        w = len(short) + rnd.randint(0, 8)
        if w >= max(len(a), len(b)): w = max(len(a), len(b)) - 1
    unit = 5
    t = rnd.randint(max(1, w * unit // 3), w * unit)
    return mat, w, t, a, b


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    rows = 0
    with tempfile.TemporaryDirectory() as d:
        for case in ['lower', 'ambig', 'star', 'gap', 'default', 'overhang']:
            made = 0
            attempts = 0
            while made < per and attempts < per * 4:
                attempts += 1
                mat, w, t, a, b = make(case)
                open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
                for f in os.listdir(d):
                    if f.startswith('dm'): os.remove(os.path.join(d, f))
                cmd = ['dotmatcher', d + '/a.fa', d + '/b.fa', '-windowsize', str(w), '-threshold', str(t), '-graph', 'data', '-goutfile', d + '/dm', '-auto']
                if mat != '-': cmd += ['-matrixfile', mat]
                r = subprocess.run(cmd, capture_output=True, text=True)
                files = [f for f in os.listdir(d) if f.startswith('dm') and f.endswith('.dat')]
                if r.returncode != 0 and 'Sequence is not nucleic' in (r.stdout + r.stderr):
                    # the second sequence must be nucleic when the first is (ACD: bsequence follows the type of asequence): an error row
                    print('\t'.join([case, mat, str(w), str(t), a, b, 'ERR:notnucleic']))
                    made += 1; rows += 1
                    continue
                if r.returncode != 0 or not files:
                    sys.stderr.write('skipped %s: rc=%d\n' % (case, r.returncode)); continue
                segs = []
                for line in open(os.path.join(d, files[0])):
                    if line.startswith('Rectangle'): break
                    m = re.match(r'Line x1 (\S+) y1 (\S+) x2 (\S+) y2 (\S+) colour', line)
                    if m: segs.append(':'.join(str(int(float(g))) for g in m.groups()))
                print('\t'.join([case, mat, str(w), str(t), a, b, ' '.join(segs) or '-']))
                made += 1; rows += 1
    sys.stderr.write('rows=%d\n' % rows)


main()
