#!/usr/bin/env python3
"""GENOMIC-MOTIFS-001 oracle (A3-12a): EMBOSS fuzznuc 6.6.0 with an IUPAC ambiguity code IN THE SEQUENCE vs an independent Python reference of the measured rules.

Measured rule (fuzznuc 6.6.0, WP36 / WP38): sequence letters are compared LITERALLY.
  * a pattern IUPAC letter L stands for its bases AND for itself: R matches sequence A, G and R; sequence R does not match pattern A, M, Y or the class [AG];
  * pattern N matches every symbol;
  * a class [ACG] lists plain bases only (matches exactly those sequence letters);
  * an exclusion {AG} matches every symbol not listed (R and N included);
  * U is read as T; X is a symbol that only pattern N and an exclusion match (the reader removes '-' gaps, not part of the contract);
  * the reverse strand is the IUPAC reverse complement of the sequence (R<->Y, K<->M, B<->V, D<->H, S W N self) searched with the same pattern;
  * everything else (mismatch counting, repeats, ranges, anchors, the patOUBrute / regex search classes) is the one of generate_fuzznuc_pattern_oracle.py.
Random syntax-tree patterns (IUPAC letters, classes, exclusions, X(n), X(n,m), separators, anchors) x random sequences over A C G T plus all 11 ambiguity letters (14-40 symbols, with planted
instances), -pmismatch k, -complement Y. The two programs are asserted equal; the script reports every disagreement on stderr and writes only agreeing rows.
TSV (tests/.../TestData/Analysis/fuzznuc_iupac_sequence_oracle.tsv):  sequence  pattern  k  forward(start0:length:mismatches ...|-)  reverse(...|-)
Never aborts: a failing fuzznuc call is reported on stderr and the row skipped.
Usage: python3 scripts/analysis/generate_fuzznuc_iupac_sequence_oracle.py [n_patterns=1200] > fuzznuc_iupac_sequence_oracle.tsv
"""
import itertools, os, random, subprocess, sys, tempfile

rnd = random.Random(20261124)
BASES = 'ACGT'
AMB = 'RYSWKMBDHVN'
EXTRA = 'UX'   # U is read as T, X is a symbol that only N and an exclusion match (a '-' gap is removed by the sequence reader: outside the contract)
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
         'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'S': 'S', 'W': 'W', 'N': 'N', 'U': 'A', 'X': 'X'}


def revcomp(s):
    return ''.join(COMP[c] for c in reversed(s))


def matches(elem, sym):
    kind, val = elem[0], elem[1]
    if sym == 'U':
        sym = 'T'
    if kind == 'L':
        return True if val == 'N' else (set(IUPAC[sym]) <= set(IUPAC[val]) if sym in IUPAC else False)
    if kind == 'C':
        return sym in val
    return sym not in val   # exclusion


class Pat:
    def __init__(self, text, a_start, a_end, elems):
        self.text, self.a_start, self.a_end, self.elems = text, a_start, a_end, elems  # elems: (kind, val, lo, hi, is_range)

    @property
    def has_range(self):
        return any(e[4] for e in self.elems)


def gen_pattern():
    n = rnd.randint(1, 5)
    elems, parts = [], []
    for _ in range(n):
        kind = rnd.random()
        if kind < 0.5:
            c = rnd.choice(list('ACGT') * 2 + list(AMB))
            atom, e = c, ('L', c)
        elif kind < 0.75:
            k = ''.join(sorted(rnd.sample(BASES, rnd.randint(1, 4))))
            atom, e = '[' + k + ']', ('C', k)
        else:
            k = ''.join(sorted(rnd.sample(BASES, rnd.randint(1, 3))))
            atom, e = '{' + k + '}', ('X', k)
        r = rnd.random()
        is_range = False
        if r < 0.55:
            lo = hi = 1
        elif r < 0.75:
            lo = hi = rnd.randint(2, 3)
            atom += '(%d)' % lo
        else:
            lo = rnd.randint(0, 2)
            hi = lo + rnd.randint(1, 2)
            atom += '(%d,%d)' % (lo, hi)
            is_range = True
        elems.append(e + (lo, hi, is_range))
        parts.append(atom)
    if all(e[2] == 0 for e in elems):
        return gen_pattern()
    text = ('-' if rnd.random() < 0.3 else '').join(parts)
    a_start, a_end = rnd.random() < 0.15, rnd.random() < 0.15
    if a_start:
        text = '<' + text
    if a_end:
        text = text + '>'
    if rnd.random() < 0.1:
        text += '.'
    return Pat(text, a_start, a_end, elems)


def branches(elems):
    out = set()
    for counts in itertools.product(*[range(e[2], e[3] + 1) for e in elems]):
        row, extras = [], []
        for idx, (e, c) in enumerate(zip(elems, counts)):
            for j in range(c):
                if j >= e[2]:
                    extras.append(len(row))
                row.append(idx)
        if row:
            out.add((tuple(row), frozenset(extras), counts[-1] == 0))
    return out


def reference(text, pat, k):
    """hits on the strand `text` itself (the reverse strand is searched as the IUPAC reverse complement)."""
    n = len(text)
    brute = pat.has_range and k > 0
    best = {}
    for row, extras, ends_empty in branches(pat.elems):
        L = len(row)
        for i in range(n - L + 1):
            if pat.a_start and i != 0:
                continue
            if pat.a_end and i + L != n:
                continue
            bad = [j for j in range(L) if not matches(pat.elems[row[j]], text[i + j])]
            if len(bad) > k:
                continue
            if brute and (any(j in extras for j in bad) or (ends_empty and i + L == n)):
                continue
            if best.get((i, L), 99) > len(bad):
                best[(i, L)] = len(bad)
    if k == 0 and pat.has_range and pat.a_start and len(best) > 1:
        longest = max(L for (_, L) in best)
        best = {key: v for key, v in best.items() if key[1] == longest}
    return best


def fuzznuc(seqs, pattern, k):
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 's.fa')
        open(f, 'w').write(''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs)))
        r = subprocess.run(['fuzznuc', '-sequence', f, '-pattern', pattern, '-pmismatch', str(k), '-complement', 'Y',
                            '-rformat', 'gff', '-stdout', '-auto'], capture_output=True, text=True, timeout=60)
    hits = {}
    for ln in r.stdout.splitlines():
        if ln.startswith('#') or not ln.strip():
            continue
        c = ln.split('\t')
        sid, start, end, strand, attr = int(c[0][1:]), int(c[3]), int(c[4]), c[6], c[8]
        mm = 0
        for note in attr.split(';'):
            if note.startswith('note=*mismatch'):
                mm = int(note.split()[-1])
        d = hits.setdefault((sid, strand), {})
        key = (start - 1, end - start + 1)
        d[key] = min(mm, d.get(key, 99))
    return hits, r.stderr


def fmt(d):
    return ' '.join('%d:%d:%d' % (i, L, mm) for (i, L), mm in sorted(d.items())) or '-'


def rand_seq(n):
    alpha = rnd.choice(['ACGT', 'ACGT', 'AC'])
    amb_rate = rnd.choice([0.1, 0.2, 0.35])
    return ''.join(rnd.choice(AMB + (EXTRA if rnd.random() < 0.15 else '')) if rnd.random() < amb_rate else rnd.choice(alpha) for _ in range(n))


def main():
    n_patterns = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
    rows = bad = 0
    for _ in range(n_patterns):
        pat = gen_pattern()
        minlen = sum(e[2] for e in pat.elems)
        k = rnd.choice([kk for kk in (0, 0, 1, 1, 2) if kk < max(minlen, 1)] or [0])
        maxlen = sum(e[3] for e in pat.elems)
        seqs = []
        for _ in range(3):
            s = rand_seq(rnd.randint(max(14, maxlen + 2), 40))
            if rnd.random() < 0.7:
                def pick(e):
                    if e[0] == 'L':
                        return rnd.choice(IUPAC[e[1]] + e[1])
                    if e[0] == 'C':
                        return rnd.choice(e[1])
                    return rnd.choice([c for c in BASES + AMB if c not in e[1]])
                inst = ''.join(pick(e) for e in pat.elems for _ in range(rnd.randint(e[2], e[3])))
                if inst and len(inst) <= len(s):
                    if rnd.random() < 0.4:
                        inst = revcomp(inst)
                    p = rnd.choice([0, len(s) - len(inst), rnd.randint(0, len(s) - len(inst))])
                    s = s[:p] + inst + s[p + len(inst):]
            seqs.append(s)
        try:
            hits, err = fuzznuc(seqs, pat.text, k)
        except Exception as e:
            sys.stderr.write('ERR %s %s\n' % (pat.text, e)); continue
        if 'Illegal' in err or 'Error' in err or 'Died' in err:
            sys.stderr.write('FUZZNUC REJECTED %s: %s\n' % (pat.text, err.strip()[:120])); continue
        for sid, s in enumerate(seqs):
            f = reference(s, pat, k)
            rc = reference(revcomp(s), pat, k)
            r = {(len(s) - i - L, L): v for (i, L), v in rc.items()}
            ff, fr = hits.get((sid, '+'), {}), hits.get((sid, '-'), {})
            if ff != f or fr != r:
                bad += 1
                sys.stderr.write('DISAGREE %s k=%d seq=%s\n  fuzznuc +%s -%s\n  python  +%s -%s\n' % (pat.text, k, s, fmt(ff), fmt(fr), fmt(f), fmt(r)))
                continue
            print('%s\t%s\t%d\t%s\t%s' % (s, pat.text, k, fmt(f), fmt(r)))
            rows += 1
    sys.stderr.write('rows=%d disagreements=%d\n' % (rows, bad))


if __name__ == '__main__':
    main()
