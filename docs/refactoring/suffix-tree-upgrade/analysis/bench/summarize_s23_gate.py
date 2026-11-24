# aggregates two raw files of run_1b9_gate.sh (base, new) into the table of results_s23_gate_rows*.txt: usage: summarize_s23_gate.py <base raw.tsv> <new raw.tsv>
import sys, collections, statistics
def load(path):
    d = collections.defaultdict(list)
    for l in open(path):
        p = l.rstrip('\n').split('\t')
        if p[0] == 'tool' or len(p) < 7: continue
        try: v = float(p[6])
        except ValueError: continue
        d[(p[2], p[4], p[5])].append(v)
    return d
def q(a, f):
    a = sorted(a); i = f * (len(a) - 1); lo = int(i); hi = min(lo + 1, len(a) - 1)
    return a[lo] + (a[hi] - a[lo]) * (i - lo)
a = load(sys.argv[1]); b = load(sys.argv[2])
print('input\toperation\tmetric\trounds\tbase_median\tbase_q25\tbase_q75\tnew_median\tnew_q25\tnew_q75\tnew/base')
for k in sorted(a):
    if k not in b or k[2] in ('input_chars',) : continue
    x, y = a[k], b[k]
    mx, my = statistics.median(x), statistics.median(y)
    print('%s\t%s\t%s\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%.3f' % (k[0], k[1], k[2], min(len(x), len(y)), mx, q(x, .25), q(x, .75), my, q(y, .25), q(y, .75), my / mx if mx else float('nan')))
