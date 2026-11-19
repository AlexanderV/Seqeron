#!/usr/bin/env python3
"""GENOMIC-MOTIFS-001 oracle (A2-32): the EMBOSS fuzznuc pattern language, fuzznuc 6.6.0 vs an independent Python reference.

Random syntax-tree patterns over the documented language: IUPAC letters, classes [ACG], exclusions {AG}, repeats X(n), ranges X(n,m) (n may be 0), `-` separators,
'<' / '>' anchors, a trailing '.'; random DNA sequences (14-40 nt, alphabets ACGT / AC, with planted instances of the pattern on both strands);
-pmismatch k (k < the shortest hit length), -complement Y. For every (sequence, pattern, k):
  fuzznuc -sequence s.fa -pattern P -pmismatch k -complement Y -rformat gff   -> (strand, start, end, `*mismatch n` note)
  Python reference: expand every repeat choice to a fixed IUPAC string, scan windows by Hamming distance (a base agrees when it is in the letter's base set), keep the
                    smallest mismatch count per (start, end); reverse strand = reversed elements with complemented base sets and swapped anchors.
                    Two fuzznuc search classes (EMBOSS embpat.c, read on raw.githubusercontent.com/kimrutherford/EMBOSS and confirmed here) are part of the reference:
                    (1) with -pmismatch > 0 a pattern that has a range X(n,m) is searched by patOUBrute: mismatches only in the mandatory part, the extra repeats beyond n
                        must match exactly;
                    (2) the same recursive matcher fails when elements of the pattern are left over once the scanned sequence has ended, so a path whose last element is
                        an empty optional element cannot end on the last base of the scanned strand;
                    (3) an exact (k = 0) range pattern starting with '<' is searched by a regex that returns only the longest match at the anchor.
                    fuzznuc lists a window once per matching path (duplicates); the oracle keeps the smallest mismatch count per window.
The two are asserted equal; the script reports (stderr) every disagreement (expected: none) and writes only agreeing rows.
TSV (tests/.../TestData/Analysis/fuzznuc_pattern_oracle.tsv):  sequence  pattern  k  forward(start0:length:mismatches ...|-)  reverse(...|-)
Never aborts: a failing fuzznuc call is reported on stderr and the row skipped.
Usage: python3 scripts/analysis/generate_fuzznuc_pattern_oracle.py > fuzznuc_pattern_oracle.tsv
"""
import itertools, os, random, subprocess, sys, tempfile

rnd = random.Random(20261123)
BASES = 'ACGT'
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
         'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def rand_seq(n, alpha):
    return ''.join(rnd.choice(alpha) for _ in range(n))


class Pat:
    def __init__(self, text, a_start, a_end, elems):
        self.text, self.a_start, self.a_end, self.elems = text, a_start, a_end, elems  # elems: (set, lo, hi, is_range, is_n)

    @property
    def has_range(self):
        return any(e[3] for e in self.elems)

    @property
    def has_n(self):
        return any(e[4] for e in self.elems)


def gen_pattern():
    n = rnd.randint(1, 5)
    elems, parts = [], []
    for _ in range(n):
        kind = rnd.random()
        is_n = False
        if kind < 0.45:
            c = rnd.choice(list('ACGT') * 3 + list('RYSWKMBDHVN'))
            atom, s = c, set(IUPAC[c])
            is_n = c == 'N'
        elif kind < 0.72:
            k = rnd.sample(BASES, rnd.randint(1, 4))
            atom, s = '[' + ''.join(k) + ']', set(k)
        else:
            k = rnd.sample(BASES, rnd.randint(1, 3))
            atom, s = '{' + ''.join(k) + '}', set(BASES) - set(k)
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
        elems.append((frozenset(s), lo, hi, is_range, is_n))
        parts.append(atom)
    if all(e[1] == 0 for e in elems):
        return gen_pattern()
    sep = '-' if rnd.random() < 0.3 else ''
    text = sep.join(parts)
    a_start, a_end = rnd.random() < 0.15, rnd.random() < 0.15
    if a_start:
        text = '<' + text
    if a_end:
        text = text + '>'
    if rnd.random() < 0.1:
        text += '.'
    return Pat(text, a_start, a_end, elems)


def branches(elems):
    """distinct (tuple of base sets, frozenset of extra positions, ends_empty): extra = repeats of an X(n,m) beyond its first n; ends_empty = the last element chose 0 repeats"""
    out = set()
    for counts in itertools.product(*[range(e[1], e[2] + 1) for e in elems]):
        sets, extras = [], []
        for e, c in zip(elems, counts):
            for j in range(c):
                if j >= e[1]:
                    extras.append(len(sets))
                sets.append(e[0])
        if sets:
            out.add((tuple(sets), frozenset(extras), counts[-1] == 0))
    return out


def reference(seq, pat, k, reverse):
    """fuzznuc reverse-complements the sequence and runs the same pattern on it; here the window is read as its reverse complement instead."""
    n = len(seq)
    brute = pat.has_range and k > 0
    best = {}
    for sets, extras, ends_empty in branches(pat.elems):
        L = len(sets)
        for i in range(n - L + 1):
            win = seq[i:i + L]
            if reverse:
                win = ''.join(COMP[c] for c in reversed(win))
            starts_at_text_start = (i + L == n) if reverse else (i == 0)   # start of the scanned strand
            ends_at_text_end = (i == 0) if reverse else (i + L == n)       # end of the scanned strand
            if pat.a_start and not starts_at_text_start:
                continue
            if pat.a_end and not ends_at_text_end:
                continue
            bad = [j for j in range(L) if win[j] not in sets[j]]
            if len(bad) > k:
                continue
            if brute and (any(j in extras for j in bad) or (ends_empty and ends_at_text_end)):
                continue
            if best.get((i, L), 99) > len(bad):
                best[(i, L)] = len(bad)
    if k == 0 and pat.has_range and pat.a_start and len(best) > 1:  # regex search: only the longest match at the anchor
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


def main():
    rows = bad = 0
    for _ in range(1100):
        pat = gen_pattern()
        minlen = sum(e[1] for e in pat.elems)
        k = rnd.choice([kk for kk in (0, 0, 1, 1, 2) if kk < max(minlen, 1)] or [0])
        maxlen = sum(e[2] for e in pat.elems)
        seqs = []
        for _ in range(3):
            alpha = rnd.choice(['ACGT', 'ACGT', 'AC'])
            s = rand_seq(rnd.randint(max(14, maxlen + 2), 40), alpha)
            if rnd.random() < 0.7:
                inst = ''.join(rnd.choice(sorted(e[0])) for e in pat.elems if e[0] for _ in range(rnd.randint(e[1], e[2])))
                if inst and len(inst) <= len(s):
                    if rnd.random() < 0.4:
                        inst = ''.join(COMP[c] for c in reversed(inst))
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
            f = reference(s, pat, k, False)
            r = reference(s, pat, k, True)
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
