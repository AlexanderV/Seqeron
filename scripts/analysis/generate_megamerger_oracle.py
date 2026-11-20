#!/usr/bin/env python3
"""GENOMIC-COMMON-001 oracle (audit A3-13b): EMBOSS megamerger 6.6.0 vs GenomicAnalyzer.MegaMerge.

Pairs of related DNA sequences: a shared core (60-300 bp) with substitutions / indels / N / IUPAC letters (rate 0-8 %) between the two copies, flanks of random length (zero flanks give the overlaps at the
start / end of a sequence, two flanks give the "neither overlap is at the start / end" warnings), lower-case stretches; and some unrelated random pairs. Word size 2-16, -prefer N / Y.
JSON lines: {"seq1": ..., "seq2": ..., "wordsize": W, "prefer": bool, "merged": "<sequence>" | null, "report": "<text>" | null, "status": "ok"|"nomerge"|"abort"}
  nomerge = the program wrote nothing (no match left); abort = the program died (return code != 0).
Never aborts: a failing call is recorded as abort.
Usage: python3 scripts/analysis/generate_megamerger_oracle.py [n=500] > megamerger_oracle.jsonl
"""
import json, os, random, subprocess, sys, tempfile

rnd = random.Random(20261127)


def rseq(n, n_rate=0.0, amb_rate=0.0):
    out = []
    for _ in range(n):
        r = rnd.random()
        if r < n_rate: out.append('N')
        elif r < n_rate + amb_rate: out.append(rnd.choice('RYSWKM'))
        else: out.append(rnd.choice('ACGT'))
    return ''.join(out)


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice('ACGT'))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice('ACGT'))
        else: out.append(c)
    return ''.join(out)


def lowercase_stretch(s):
    if rnd.random() < 0.25 and len(s) > 10:
        a = rnd.randrange(len(s) - 5); b = a + rnd.randint(3, min(30, len(s) - a))
        return s[:a] + s[a:b].lower() + s[b:]
    return s


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    counts = {'ok': 0, 'nomerge': 0, 'abort': 0}
    with tempfile.TemporaryDirectory() as d:
        for _ in range(n):
            if rnd.random() < 0.85:
                core = rseq(rnd.randint(60, 300), rnd.choice([0, 0, 0.01]), rnd.choice([0, 0, 0.01]))
                f1 = rseq(rnd.choice([0, 0, rnd.randint(1, 80)])); f2 = rseq(rnd.choice([0, 0, rnd.randint(1, 80)]))
                g1 = rseq(rnd.choice([0, 0, rnd.randint(1, 80)])); g2 = rseq(rnd.choice([0, 0, rnd.randint(1, 80)]))
                rate = rnd.choice([0, 0.01, 0.03, 0.06, 0.08])
                s1 = f1 + mutate(core, rate / 2) + g1
                s2 = f2 + mutate(core, rate / 2) + g2
            else:
                s1 = rseq(rnd.randint(20, 200)); s2 = rseq(rnd.randint(20, 200))
            s1 = lowercase_stretch(s1); s2 = lowercase_stretch(s2)
            w = rnd.choice([2, 4, 6, 8, 10, 12, 16, 8, 10, 12])
            prefer = rnd.random() < 0.3
            fa, fb, so, sr = (os.path.join(d, x) for x in ('a.fa', 'b.fa', 'o.fa', 'r.txt'))
            open(fa, 'w').write('>seqA\n' + s1 + '\n'); open(fb, 'w').write('>seqB\n' + s2 + '\n')
            for f in (so, sr):
                if os.path.exists(f): os.remove(f)
            cmd = ['megamerger', fa, fb, '-wordsize', str(w), '-outseq', so, '-osformat', 'fasta', '-outfile', sr, '-auto'] + (['-prefer'] if prefer else [])
            r = subprocess.run(cmd, capture_output=True, text=True)
            merged = report = None
            if r.returncode != 0:
                status = 'abort'
            else:
                txt = open(sr).read() if os.path.exists(sr) else ''
                fa_txt = open(so).read() if os.path.exists(so) else ''
                if not txt.strip():
                    status = 'nomerge'
                else:
                    status = 'ok'; report = txt
                    merged = ''.join(l for l in fa_txt.split('\n') if l and not l.startswith('>'))
            counts[status] += 1
            print(json.dumps({'seq1': s1, 'seq2': s2, 'wordsize': w, 'prefer': prefer, 'merged': merged, 'report': report, 'status': status}))
    sys.stderr.write('rows=%d %s\n' % (n, counts))


main()
