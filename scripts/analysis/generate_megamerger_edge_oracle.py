#!/usr/bin/env python3
"""GENOMIC-COMMON-001 oracle (audit A4-9): the EDGE CASES of EMBOSS megamerger 6.6.0 vs GenomicAnalyzer.MegaMerge that generate_megamerger_oracle.py never produces (that script has lower-case STRETCHES and single N / IUPAC letters).

Case kinds (>= 100 rows each), pairs of related DNA sequences as in the main generator:
  lowerall   both sequences entirely lower case (the merged sequence keeps the case of the letters)
  lowerone   one sequence entirely lower case, the other upper case (which copy of a matched region is written, and in which case)
  nrun       runs of N longer than the word size (inside the shared core and in the flanks): the word matcher skips them (an N run never matches an N run)
  short      sequences shorter than / equal to the word size, identical sequences and one sequence inside the other (the program's "no match" / warning paths)
JSON lines: {"case": ..., "seq1": ..., "seq2": ..., "wordsize": W, "prefer": bool, "merged": "<sequence>" | null, "report": "<text>" | null, "status": "ok"|"nomerge"|"abort"}  (abort = the program died, return code != 0)
Usage: python3 scripts/analysis/generate_megamerger_edge_oracle.py [per_case=120] > megamerger_edge_oracle.jsonl
"""
import json, os, random, subprocess, sys, tempfile

rnd = random.Random(20261202)


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


def related(core_lo=60, core_hi=300):
    core = rseq(rnd.randint(core_lo, core_hi), rnd.choice([0, 0, 0.01]), rnd.choice([0, 0, 0.01]))
    f1 = rseq(rnd.choice([0, 0, rnd.randint(1, 80)])); f2 = rseq(rnd.choice([0, 0, rnd.randint(1, 80)]))
    g1 = rseq(rnd.choice([0, 0, rnd.randint(1, 80)])); g2 = rseq(rnd.choice([0, 0, rnd.randint(1, 80)]))
    rate = rnd.choice([0, 0.01, 0.03, 0.06, 0.08])
    return f1 + mutate(core, rate / 2) + g1, f2 + mutate(core, rate / 2) + g2


def make(case, w):
    if case == 'lowerall':
        s1, s2 = related()
        return s1.lower(), s2.lower()
    if case == 'lowerone':
        s1, s2 = related()
        return (s1.lower(), s2) if rnd.random() < 0.5 else (s1, s2.lower())
    if case == 'nrun':
        core = rseq(rnd.randint(80, 260))
        pos = sorted(rnd.sample(range(5, len(core) - 5), 2))
        for p in pos:
            ln = rnd.randint(w + 1, w + 12)
            core = core[:p] + 'N' * ln + core[p + ln:]
        f1 = rseq(rnd.choice([0, rnd.randint(1, 40)])) + ('N' * rnd.randint(w, w + 8) if rnd.random() < 0.3 else '')
        g2 = ('N' * rnd.randint(w, w + 8) if rnd.random() < 0.3 else '') + rseq(rnd.choice([0, rnd.randint(1, 40)]))
        rate = rnd.choice([0, 0.01, 0.03])
        s1 = f1 + mutate(core, rate / 2) + rseq(rnd.choice([0, rnd.randint(1, 40)]))
        s2 = rseq(rnd.choice([0, rnd.randint(1, 40)])) + mutate(core, rate / 2) + g2
        return s1, s2
    # short
    kind = rnd.choice(['tiny', 'equal', 'inside', 'identical'])
    if kind == 'tiny':
        return rseq(rnd.randint(1, w + 2)), rseq(rnd.randint(1, w + 20))
    if kind == 'equal':
        a = rseq(w)
        return a, a if rnd.random() < 0.5 else rseq(w)
    if kind == 'inside':
        big = rseq(rnd.randint(40, 120)); i = rnd.randrange(0, len(big) - 20)
        small = big[i:i + rnd.randint(w, 20)]
        return (big, small) if rnd.random() < 0.5 else (small, big)
    a = rseq(rnd.randint(30, 120))
    return a, a


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    counts = {}
    with tempfile.TemporaryDirectory() as d:
        for case in ['lowerall', 'lowerone', 'nrun', 'short']:
            for _ in range(per):
                w = rnd.choice([2, 4, 6, 8, 10, 12, 16, 8, 10, 12])
                s1, s2 = make(case, w)
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
                counts[(case, status)] = counts.get((case, status), 0) + 1
                print(json.dumps({'case': case, 'seq1': s1, 'seq2': s2, 'wordsize': w, 'prefer': prefer, 'merged': merged, 'report': report, 'status': status}))
    sys.stderr.write('rows %s\n' % sorted(counts.items()))


main()
