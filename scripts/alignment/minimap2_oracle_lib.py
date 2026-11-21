#!/usr/bin/env python3
"""Shared helpers of the minimap2 2.26 oracle generators (A3-9c2 / A5-9b / A5-9c, B09 WP45): random references, mutated / chimeric / inverted / spliced queries, FASTA / FASTQ text and the run of the installed `minimap2` in a temporary
directory with the files named r.fa / q.fa / q1.fa / q2.fa (so that the @PG command line of a SAM header is reproducible).

Every oracle line is {"args": [...], "files": {name: text}, "stdout": "<the program's standard output>", "set": "<group>"}: the test feeds `args` and `files` to Minimap2Program.Run and compares the whole standard output as text. A reference text (r*.fa) longer than 800
characters that occurs in several cases is written once as a line {"pool": N, "text": "..."} and referred to as {"$": N} in "files".
"""
import json, os, random, subprocess, tempfile

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def revcomp(s):
    return s.translate(COMP)[::-1]


class Gen:
    def __init__(self, seed):
        self.rnd = random.Random(seed)

    def rseq(self, n, hp=False, gc=None):
        r = self.rnd
        s = []
        while len(s) < n:
            c = r.choice('ACGT') if gc is None else (r.choice('GC') if r.random() < gc else r.choice('AT'))
            s.extend([c] * (r.randint(3, 9) if hp and r.random() < 0.15 else 1))
        return ''.join(s[:n])

    def mutate(self, s, rate, big_indels=True):
        r = self.rnd
        out = []
        i = 0
        while i < len(s):
            c = s[i]
            x = r.random()
            if x < rate * 0.5:
                out.append(r.choice('ACGT'))
            elif x < rate * 0.75:
                if big_indels and r.random() < 0.05:
                    i += r.choice([8, 20, 60, 150])       # a long deletion
            elif x < rate:
                out.append(c)
                n = r.choice([1, 1, 1, 2, 3, 5]) if not (big_indels and r.random() < 0.05) else r.choice([10, 30, 80, 200])
                out.extend(r.choice('ACGT') for _ in range(n))
            else:
                out.append(c)
            i += 1
        return ''.join(out)

    def with_n(self, s, k):
        r = self.rnd
        c = list(s)
        for _ in range(k):
            c[r.randrange(len(c))] = 'N'
        return ''.join(c)

    def quals(self, n):
        r = self.rnd
        return ''.join(chr(33 + r.choice([2, 10, 20, 30, 38, 40])) for _ in range(n))


def fasta(names, seqs, width=None):
    out = []
    for n, s in zip(names, seqs):
        out.append('>%s\n' % n)
        if width:
            for i in range(0, len(s), width):
                out.append(s[i:i + width] + '\n')
        else:
            out.append(s + '\n')
    return ''.join(out)


def fastq(reads):
    return ''.join('@%s\n%s\n+\n%s\n' % (n, s, q) for n, s, q in reads)


def run(args, files, record_err=False, outfiles=()):
    """Run minimap2 with `args` in a fresh directory containing `files`; returns stdout (None on failure), or with record_err (stdout, stderr, exit code, {output file: text})."""
    with tempfile.TemporaryDirectory() as d:
        for n, t in files.items():
            open(os.path.join(d, n), 'w', newline='').write(t)
        p = subprocess.run(['minimap2'] + args, capture_output=True, cwd=d)
        p.stdout = p.stdout.decode('utf-8', 'surrogateescape')
        p.stderr = p.stderr.decode('utf-8', 'surrogateescape')
        if record_err:
            made = {}
            for f in outfiles:
                fp = os.path.join(d, f)
                if os.path.exists(fp):
                    made[f] = open(fp, newline='').read()
            return p.stdout, p.stderr, p.returncode, made
        if p.returncode != 0:
            return None
        return p.stdout


POOL = {}


def pooled(files, out):
    """Replace the long file texts by pool references, writing a pool line to `out` (stdout when None) the first time a text is used."""
    res = {}
    for n, t in files.items():
        if len(t) > 800 and n.startswith('r'):
            if t not in POOL:
                POOL[t] = len(POOL)
                line = json.dumps({'pool': POOL[t], 'text': t})
                if out is None:
                    print(line)
                else:
                    out.write(line + '\n')
            res[n] = {'$': POOL[t]}
        else:
            res[n] = t
    return res


def emit(group, args, files, out=None, record_err=False, outfiles=()):
    if record_err:
        stdout, stderr, rc, made = run(args, files, True, outfiles)
        if rc not in (0, 1):
            return False
        rec = {'set': group, 'args': args, 'files': pooled(files, out), 'stdout': stdout, 'stderr': stderr, 'exit': rc, 'outfiles': made}
        if any(a in args for a in ('--print-seeds', '--print-chains', '--print-qname', '--print-aln-seq', '--no-kalloc')):
            rec['stderr'] = None      # the debug log of the program is not reproduced
        if out is None:
            print(json.dumps(rec))
        else:
            out.write(json.dumps(rec) + '\n')
        return True
    stdout = run(args, files)
    if stdout is None:
        return False
    line = json.dumps({'set': group, 'args': args, 'files': pooled(files, out), 'stdout': stdout})
    if out is None:
        print(line)
    else:
        out.write(line + '\n')
    return True
