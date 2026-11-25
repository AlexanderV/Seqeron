"""ASSEMBLY-QC-001 oracle: CalculateLocalQuality / FindSuspiciousRegions.

Independent Python references:
  * windowed GC       - Bio.SeqUtils.gc_fraction on seq[i:i+w] (ACGT-only fixtures, so Biopython 'remove'
                        equals the QUAST (G+C)/(A+C+G+T) rule); 0 for a window without A/C/G/T
  * linguistic cmplx  - Orlov & Potapov 2004 / Rosalind LING (m=10, sum of distinct words / sum of maximum), written
                        from the formula, N counted as a fifth symbol
  * DUST intervals    - line-by-line port of lh3/sdust sdust.c (sdust_core, W=64, T=20 -> 'T' is 10x threshold),
                        read from raw.githubusercontent.com/lh3/sdust (not executed)
  * region assembly   - union of flagged windows / intervals, merged when overlapping or adjacent
Usage: python qc_oracle.py   (prints the fixture table used by the C# tests)
"""
import random
from collections import deque
from Bio.SeqUtils import gc_fraction

NT4 = {c: i for i, c in enumerate("ACGT")}
NT4.update({c.lower(): i for i, c in enumerate("ACGT")})


def sdust(seq, T=20, W=64):
    """Port of sdust_core; returns half-open [start, finish) intervals."""
    WL, WTOT = 3, 64
    WMSK = WTOT - 1
    P = []          # perfect intervals: dicts, sorted by descending start then ascending finish
    res = []        # result (start, finish)
    w = deque()
    state = dict(rv=0, rw=0, L=0)
    cv = [0] * WTOT
    cw = [0] * WTOT

    def shift_window(t):
        if len(w) >= W - WL + 1:
            s = w.popleft()
            cw[s] -= 1
            state['rw'] -= cw[s]
            if state['L'] > len(w):
                state['L'] -= 1
                cv[s] -= 1
                state['rv'] -= cv[s]
        w.append(t)
        state['L'] += 1
        state['rw'] += cw[t]; cw[t] += 1
        state['rv'] += cv[t]; cv[t] += 1
        if cv[t] * 10 > T << 1:
            while True:
                s = w[len(w) - state['L']]
                cv[s] -= 1
                state['rv'] -= cv[s]
                state['L'] -= 1
                if s == t:
                    break

    def save_masked(start):
        if not P or P[-1]['start'] >= start:
            return
        p = P[-1]
        saved = False
        if res:
            s, f = res[-1]
            if p['start'] <= f:
                saved = True
                res[-1] = (s, max(f, p['finish']))
        if not saved:
            res.append((p['start'], p['finish']))
        i = len(P) - 1
        while i >= 0 and P[i]['start'] < start:
            i -= 1
        del P[i + 1:]

    def find_perfect(start):
        c = cv[:]
        r = state['rv']
        max_r = max_l = 0
        i = len(w) - state['L'] - 1
        while i >= 0:
            t = w[i]
            r += c[t]; c[t] += 1
            new_r, new_l = r, len(w) - i - 1
            if new_r * 10 > T * new_l:
                j = 0
                while j < len(P) and P[j]['start'] >= i + start:
                    p = P[j]
                    if max_r == 0 or p['r'] * max_l > max_r * p['l']:
                        max_r, max_l = p['r'], p['l']
                    j += 1
                if max_r == 0 or new_r * max_l >= max_r * new_l:
                    max_r, max_l = new_r, new_l
                    P.insert(j, dict(start=i + start, finish=len(w) + (WL - 1) + start, r=new_r, l=new_l))
            i -= 1

    n = len(seq)
    l = t = 0
    for i in range(n + 1):
        b = NT4.get(seq[i], 4) if i < n else 4
        if b < 4:
            l += 1
            t = ((t << 2) | b) & WMSK
            if l >= WL:
                start = (l - W if l - W > 0 else 0) + (i + 1 - l)
                save_masked(start)
                shift_window(t)
                if state['rw'] * 10 > state['L'] * T:
                    find_perfect(start)
        else:
            start = (l - W + 1 if l - W + 1 > 0 else 0) + (i + 1 - l)
            while P:
                save_masked(start)
                start += 1
            l = t = 0
    return res


def dust_score(s):
    """Morgulis 2006: sum c_t(c_t-1)/2 / (l-1), l = number of triplets."""
    s = s.upper()
    words = [s[i:i + 3] for i in range(len(s) - 2)]
    if len(words) < 2:
        return 0.0
    from collections import Counter
    return sum(c * (c - 1) / 2 for c in Counter(words).values()) / (len(words) - 1)


def linguistic_complexity(s, m=10):
    """Orlov & Potapov 2004 / Rosalind LING: LC = sum_i U_i / sum_i Umax_i over word lengths i = 1..min(m, n),
    Umax_i = min(a^i, n - i + 1), a = 4 plus every extra symbol present (N -> a = 5), as the canonical unit."""
    n = len(s)
    if n == 0:
        return 0.0
    a = len(set(s) | set('ACGT'))      # N (or any other symbol) extends the alphabet to 5
    num = den = 0
    for i in range(1, min(m, n) + 1):
        num += len({s[j:j + i] for j in range(n - i + 1)})
        den += min(a ** i, n - i + 1)
    return num / den


def windows(n, w, step):
    return [(i, min(w, n - i)) for i in range(0, n, step)]


def local_quality(seq, w, step=None):
    step = step or max(1, w // 2)
    out = []
    for i, ln in windows(len(seq), w, step):
        win = seq[i:i + ln]
        valid = sum(c in "ACGTacgt" for c in win)
        gc = gc_fraction(win.replace("N", "").replace("n", "")) if valid else 0.0
        out.append((i, ln, gc, win.upper().count("N")))
    return out


def suspicious(seqs, gc_dev=0.15, dust_t=2.0, max_n=0.1, w=500):
    tot_gc = sum(sum(c in "GCgc" for c in s) for s in seqs.values())
    tot_valid = sum(sum(c in "ACGTacgt" for c in s) for s in seqs.values())
    g = tot_gc / tot_valid if tot_valid else None
    step = max(1, w // 2)
    result = []
    for sid, s in seqs.items():
        hits = []
        prev = -1
        for i, ln, gc, nn in local_quality(s, w):
            covered = ln < w and prev >= 0 and prev + w >= len(s)
            prev = i
            if covered:
                continue
            if nn / ln > max_n:
                hits.append((i, i + ln - 1, 1, nn / ln))
            valid = sum(c in "ACGTacgt" for c in s[i:i + ln])
            if g is not None and valid and abs(gc - g) > gc_dev:
                hits.append((i, i + ln - 1, 2, abs(gc - g)))
        for a, b in sdust(s, int(round(dust_t * 10)), 64):
            hits.append((a, b - 1, 4, dust_score(s[a:b])))
        hits.sort(key=lambda h: (h[0], h[1]))
        k = 0
        while k < len(hits):
            rs, re_, fl = hits[k][0], hits[k][1], 0
            mg = mn = md = 0.0
            while k < len(hits) and hits[k][0] <= re_ + 1:
                h = hits[k]
                re_ = max(re_, h[1]); fl |= h[2]
                if h[2] == 1: mn = max(mn, h[3])
                elif h[2] == 2: mg = max(mg, h[3])
                else: md = max(md, h[3])
                k += 1
            reasons = [n for b, n in ((1, "High N content"), (2, "GC deviation"), (4, "Low complexity")) if fl & b]
            result.append((sid, rs, re_, "; ".join(reasons), mg, mn, md))
    return result


def rand_dna(seed, n, gc=0.5):
    """Deterministic LCG generator mirrored by the C# tests (x = x*1103515245+12345 mod 2^31;
    u = (x >> 8 & 0xFFFF) / 65536 picks GC vs AT, the next draw's bit 16 picks the base)."""
    x = seed
    out = []
    for _ in range(n):
        x = (x * 1103515245 + 12345) & 0x7FFFFFFF
        u = ((x >> 8) & 0xFFFF) / 65536.0
        x = (x * 1103515245 + 12345) & 0x7FFFFFFF
        bit = (x >> 16) & 1
        out.append(("G" if bit else "C") if u < gc else ("A" if bit else "T"))
    return "".join(out)


FIXTURES = {
    "polyA_in_random": {"s": rand_dna(1, 2000) + "A" * 2000 + rand_dna(2, 2000)},
    "gc_island": {"s": rand_dna(3, 3000, 0.4) + rand_dna(4, 1000, 0.95) + rand_dna(5, 3000, 0.4)},
    "n_run": {"s": rand_dna(6, 1500) + "N" * 400 + rand_dna(7, 1500)},
    "polyA_small": {"s": rand_dna(15, 3000) + "A" * 300 + rand_dna(16, 3000)},
    "multi": {"a": rand_dna(8, 1200, 0.45), "b": rand_dna(9, 1100, 0.45) + "AC" * 150 + rand_dna(10, 900, 0.45)},
}

if __name__ == "__main__":
    for name, seqs in FIXTURES.items():
        print(name, suspicious(seqs))
    for s in ["ACGT" * 3, "A" * 40, rand_dna(11, 60)]:
        print(len(s), round(linguistic_complexity(s), 6), sdust(s), round(dust_score(s), 6))
    for row in local_quality(rand_dna(12, 1500), 1000):
        print(row)


def emit_csharp():
    """Prints the expected values pasted into GenomeAssemblyAnalyzer_QualityScreen_Tests."""
    for name, seqs in FIXTURES.items():
        print("//", name)
        for r in suspicious(seqs):
            print(f'  ("{r[0]}", {r[1]}, {r[2]}, "{r[3]}", {r[4]!r}, {r[5]!r}, {r[6]!r}),')
    for name, seq, w, step in [("lq_random1500", rand_dna(12, 1500), 1000, None), ("lq_n", rand_dna(13, 400) + "N" * 50 + rand_dna(14, 300), 300, 100)]:
        print("//", name)
        for i, ln, gc, nn in local_quality(seq, w, step):
            print(f'  ({i}, {ln}, {gc!r}, {nn}, {linguistic_complexity(seq[i:i+ln].upper())!r}),')
