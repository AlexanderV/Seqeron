#!/usr/bin/env python3
"""A5-8 (B09 F176) oracle: the feature-report texts of the EMBOSS 6.6.0 pattern / repeat programs - `-rformat gff` and `-rformat excel` of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem - as TEXT, byte for byte
(the `#!Date` line of the gff header is the date of the run: the test masks it).

TSV (tests/.../TestData/Alignment/emboss_report_oracle.tsv):  id  program  format  params  sequence  expected
  params = as generate_emboss_gap_oracle.py: fuzznuc 'pattern|mismatch|complement', fuzzpro 'pattern', fuzztran 'pattern|frame', etandem 'min,max,threshold,mismatch,uniform', equicktandem 'max,threshold'
  expected = the whole report (\\n, \\t, \\\\ escaped); an empty report (the program found nothing: it writes nothing) is the empty string.
Usage: python3 scripts/alignment/generate_emboss_report_oracle.py OUTFILE [per_program_format=120] [seed]      (needs the EMBOSS 6.6.0 programs on PATH)
"""
import os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_emboss_gap_oracle import Runner, mutate, NUC, AA, esc, tandem_seq  # noqa: E402


def fuzznuc_case(rnd):
    s = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(25, 90)))
    if rnd.random() < 0.2: s = s.lower()
    st = rnd.randint(0, len(s) - 6)
    pat = rnd.choice([s[st:st + rnd.randint(3, 6)].upper(), s[st:st + 4].upper().replace(s[st].upper(), 'N', 1), s[st:st + 3].upper() + '[AC]' + s[st + 4:st + 6].upper(), 'AC' + 'N' * rnd.randint(1, 2) + 'GT', 'RYR', s[st:st + 3].upper(), s[st].upper() + '(2,3)', 'AC(1,2)N', s[st:st + 2].upper() + 'N(1,2)' + s[st + 3].upper()])
    mm = rnd.choice([0, 0, 1, 1, 2]); comp = rnd.choice([0, 1])
    return '%s|%d|%d' % (pat, mm, comp), s, ['-pattern', pat, '-pmismatch', str(mm), '-complement', 'Y' if comp else 'N']


def fuzzpro_case(rnd):
    s = ''.join(rnd.choice(AA) for _ in range(rnd.randint(25, 90)))
    st = rnd.randint(0, len(s) - 5)
    pat = rnd.choice(['-'.join(s[st:st + 2]), '-'.join(s[st:st + 3]), s[st] + '-x-' + s[st + 2], '[' + s[st] + s[st + 1] + ']-' + s[st + 2], 'x(2)-' + s[st + 3], s[st] + '-x(1,2)-' + s[st + 3], 'x(1,3)-' + s[st + 3]])
    return pat, s, ['-pattern', pat]


def fuzztran_case(rnd):
    s = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(30, 100)))
    pat = rnd.choice(['L-x-S', 'S-S', 'G-R', 'x(2)-K', '[LV]-T', 'M', 'A-A', 'T-x(1,2)-G'])
    fr = rnd.choice(['F', 'F', '6', 'R', '1', '2', '3', '-1', '-2', '-3'])
    return pat + '|' + fr, s, ['-pattern', pat, '-frame', fr, '-table', '0']


def etandem_case(rnd):
    d = tandem_seq(rnd)
    mn = rnd.choice([2, 3]); mx = mn + rnd.randint(1, 5); th = rnd.choice([8, 12, 15, 20]); mm = rnd.choice([0, 1]); un = rnd.choice([0, 1])
    return '%d,%d,%d,%d,%d' % (mn, mx, th, mm, un), d, ['-minrepeat', str(mn), '-maxrepeat', str(mx), '-threshold', str(th)] + (['-mismatch'] if mm else []) + (['-uniform'] if un else [])


def equick_case(rnd):
    d = tandem_seq(rnd)
    mx = rnd.choice([4, 6, 8]); th = rnd.choice([5, 8, 12, 20])
    return '%d,%d' % (mx, th), d, ['-maxrepeat', str(mx), '-threshold', str(th)]


CASES = {'fuzznuc': fuzznuc_case, 'fuzzpro': fuzzpro_case, 'fuzztran': fuzztran_case, 'etandem': etandem_case, 'equicktandem': equick_case}


def main(outfile, per, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        rn = Runner(d)
        for prog, fn in CASES.items():
            for fmt in ('gff', 'excel'):
                made = 0; tries = 0; empties = 0
                while made < per and tries < per * 40:
                    tries += 1
                    params, seq, args = fn(rnd)
                    rn.w('s.fa', '>s\n' + seq + '\n'); rn.rm('o.rep', 'o.orig')
                    cmd = [prog, '-sequence', 's.fa', '-rformat', fmt, '-auto'] + args
                    cmd += ['-outfile', 'o.rep'] + (['-origfile', 'o.orig'] if prog in ('etandem', 'equicktandem') else [])
                    p = rn.run(cmd)
                    if p.returncode != 0: continue
                    text = rn.r('o.rep')
                    if text == '':                                              # keep a few empty reports (the program writes nothing), mostly rows with hits
                        if empties >= per // 12: continue
                        empties += 1
                    rows.append('\t'.join([str(len(rows)), prog, fmt, params, seq, esc(text)]))
                    made += 1
    with open(outfile, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('id\tprogram\tformat\tparams\tsequence\texpected\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 120, int(sys.argv[3]) if len(sys.argv) > 3 else 20261212)
