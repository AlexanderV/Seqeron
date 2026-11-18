#!/usr/bin/env python3
"""ALIGN-STATS-001 (audit A1-9) oracle: NCBI BLAST+ 2.12.0 blastn / blastp report statistics for random homologous pairs.

For every pair the best HSP is reported twice: tabular `-outfmt "6 qseq sseq length nident mismatch gapopen gaps positive pident"`
and the pairwise text format (`-outfmt 0`), whose line ` Identities = a/L (p%), Positives = b/L (q%), Gaps = g/L (r%)` carries the
INTEGER percentages (CAlignFormatUtil::GetPercentMatch, NCBI C++ toolkit align_format_util.cpp: round half up, 100 only if
numerator == denominator, else at most 99). Output TSV (embedded in tests/.../TestData/Alignment/blast_stats_oracle.tsv):
 program qseq sseq length nident mismatch gapopen gaps positive pident(3 decimals) idPct posPct gapPct   (posPct = NA for blastn)

Usage: python3 scripts/alignment/generate_blast_stats_oracle.py [n_per_program=150] > blast_stats_oracle.tsv
"""
import random, re, subprocess, sys, tempfile

rnd = random.Random(20261118)
DNA = 'ACGT'
PROT = 'ACDEFGHIKLMNPQRSTVWY'


def mutate(s, alpha, rate):
    out = []
    for c in s:
        x = rnd.random()
        if x < rate * 0.6: out.append(rnd.choice(alpha))
        elif x < rate * 0.8: pass
        elif x < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out)


def run(prog, d, q, s, fmt, extra):
    open(d + '/q.fa', 'w').write('>q\n' + q + '\n'); open(d + '/s.fa', 'w').write('>s\n' + s + '\n')
    cmd = [prog, '-query', d + '/q.fa', '-subject', d + '/s.fa', '-outfmt', fmt, '-max_hsps', '1', '-evalue', '1e-3'] + extra
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else ''


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    print('program\tqseq\tsseq\tlength\tnident\tmismatch\tgapopen\tgaps\tpositive\tpident\tidPct\tposPct\tgapPct')
    with tempfile.TemporaryDirectory() as d:
        for prog in ('blastn', 'blastp'):
            alpha = DNA if prog == 'blastn' else PROT
            made = 0
            while made < n:
                q = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(60, 400) if prog == 'blastn' else rnd.randint(40, 300)))
                s = mutate(q, alpha, rnd.choice([0.004, 0.006, 0.03, 0.08, 0.15, 0.25, 0.35] if prog == 'blastn' else [0.006, 0.01, 0.05, 0.15, 0.3, 0.45]))
                if len(s) < 30: continue
                extra = [] if prog == 'blastp' else ['-task', rnd.choice(['blastn', 'blastn', 'megablast'])]
                if prog == 'blastp': extra = ['-comp_based_stats', '0', '-seg', 'no']
                if prog == 'blastn': extra += ['-dust', 'no']
                tab = run(prog, d, q, s, '6 qseq sseq length nident mismatch gapopen gaps positive pident', extra).strip().split('\n')[0].split('\t')
                txt = run(prog, d, q, s, '0', extra)
                if len(tab) < 9: continue
                m = re.search(r'Identities = (\d+)/(\d+) \((\d+)%\)(?:, Positives = (\d+)/(\d+) \((\d+)%\))?, Gaps = (\d+)/(\d+) \((\d+)%\)', txt)
                if not m: continue
                qseq, sseq, length, nident, mism, gopen, gaps, pos, pident = tab[:9]
                assert int(m.group(1)) == int(nident) and int(m.group(2)) == int(length) and int(m.group(7)) == int(gaps)
                posPct = m.group(6) if m.group(6) else 'NA'
                print('\t'.join([prog, qseq, sseq, length, nident, mism, gopen, gaps, pos, pident, m.group(3), posPct, m.group(9)]))
                made += 1


main()
