"""Literal port of NCBI WindowMasker first pass: win_mask_gen_counts.cpp (counting),
seq_masker_window.cpp (CSeqMaskerWindow FillWindow/Advance), seq_masker_score_mean.cpp (sum/num),
seq_masker.cpp CSeqMasker::DoMask (no merge pass, trigger=mean, unit_step=1, window_step=1)."""
import json, random, sys
L = {'A': 1, 'C': 2, 'G': 3, 'T': 4}

def rc(u, size):
    r = 0
    for _ in range(size):
        r = (r << 2) | (3 - (u & 3)); u >>= 2
    return r

def count_units(seqs, k):
    mask = (1 << (2*k)) - 1
    counts = {}
    for s in seqs:
        unit = 0; count = 0
        for ch in s.upper():
            if ch not in 'ACGT':
                count = 0; unit = 0; continue
            unit = ((unit << 2) & mask) + L[ch] - 1
            if count >= k - 1:
                ru = rc(unit, k)
                if unit <= ru: counts[unit] = counts.get(unit, 0) + 1
                if ru <= unit: counts[ru] = counts.get(ru, 0) + 1
            count += 1
    return counts

def count_fwd(seqs, k):
    mask = (1 << (2*k)) - 1
    counts = {}
    for s in seqs:
        unit = 0; count = 0
        for ch in s.upper():
            if ch not in 'ACGT':
                count = 0; unit = 0; continue
            unit = ((unit << 2) & mask) + L[ch] - 1
            if count >= k - 1: counts[unit] = counts.get(unit, 0) + 1
            count += 1
    return counts

class Window:
    def __init__(s, data, k, W):
        s.data = data; s.k = k; s.W = W; s.mask = (1 << (2*k)) - 1
        s.units = [0]*(W-k+1); s.fill(0)
    def fill(s, ws):
        s.first = 0; unit = 0; it = 0; s.end = ws
        while it < s.W and s.end < len(s.data):
            ch = s.data[s.end]
            if ch not in L:
                it = 0; s.end += 1; continue   # C++: iter=-1 then for-loop ++iter
            unit = ((unit << 2) & s.mask) + L[ch] - 1
            if it >= s.k - 1: s.units[it + 1 - s.k] = unit
            it += 1; s.end += 1
        s.start = s.end - s.W; s.end -= 1
        s.state = (it == s.W)
    def advance(s):  # step 1, unit_step 1 : C++ fast path emulated by refill (identical window content)
        s.fill(s.start + 1) if False else None
        # literal fast path
        n = len(s.units); step = 1
        it = 0
        e = s.end + 1
        if e < len(s.data):
            ch = s.data[e]
            if ch not in L:
                s.fill(e); return
            # shift
            unit = ((s.units[-1] << 2) & s.mask) + L[ch] - 1
            s.units = s.units[1:] + [unit]
            s.end = e; s.start += 1
        else:
            s.state = False

def mask(seqs_all, seq, k, cutoff, W, textend, both):
    counts = count_units(seqs_all, k) if both else count_fwd(seqs_all, k)
    def cnt(u):
        if both: return counts.get(min(u, rc(u, k)), 0)
        return counts.get(u, 0)
    data = seq.upper()
    out = []
    if len(data) < W: return out
    w = Window(data, k, W)
    start = end = cend = 0
    while w.state:
        s = sum(cnt(u) for u in w.units) // len(w.units)
        ws, we = w.start, w.end
        if s < textend:
            if end > start and ws > cend:
                out.append((start, end)); start = end = cend = 0
        elif s < cutoff:
            if end > start:
                if ws > cend + 1:
                    out.append((start, end)); start = end = cend = 0
                else: cend = we
        else:
            if end > start:
                if ws > cend + 1:
                    out.append((start, end)); start = ws
            else: start = ws
            cend = end = we
        w.advance()
    if end > start: out.append((start, end))
    res = []
    for a, b in out:
        res.append((a, b, max(cnt(u) for u in unit_list(data, k, a, b))))
    return res

def unit_list(data, k, a, b):
    mask_ = (1 << (2*k)) - 1
    for j in range(a, b - k + 2):
        u = 0
        for ch in data[j:j+k]: u = (u << 2) + L[ch] - 1
        yield u

def gen(rng):
    k = rng.choice([3, 4, 5, 6, 7, 8, 15])
    W = rng.randint(k, k + rng.choice([2, 6, 20, 40]))
    t = rng.randint(1, 6); e = rng.randint(1, t)
    base = ''.join(rng.choice('ACGT') for _ in range(rng.randint(30, 150)))
    motif = ''.join(rng.choice('ACGT') for _ in range(rng.randint(k, k + 12)))
    nseq = rng.randint(1, 3)
    seqs = []
    for _ in range(nseq):
        s = ''.join(rng.choice('ACGT') for _ in range(rng.randint(20, 200)))
        for _ in range(rng.randint(0, 3)):
            p = rng.randint(0, len(s)); m = motif if rng.random() < .6 else motif[::-1].translate(str.maketrans('ACGT','TGCA'))
            s = s[:p] + m * rng.randint(1, 6) + s[p:]
        if rng.random() < .5:
            for _ in range(rng.randint(1, 4)):
                p = rng.randint(0, len(s)-1); s = s[:p] + rng.choice('NNRYn') + s[p+1:]
        if rng.random() < .3: s = s.lower()
        seqs.append(s)
    return dict(seqs=seqs, k=k, t=t, w=W, e=e, both=rng.random() < .7)

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1])); n = int(sys.argv[2])
    cases = [gen(rng) for _ in range(n)]
    json.dump(cases, open(sys.argv[3], 'w'))
    exp = []
    for c in cases:
        exp.append([[f's{i}', a, b, cp] for i, s in enumerate(c['seqs']) for a, b, cp in mask(c['seqs'], s, c['k'], c['t'], c['w'], c['e'], c['both'])])
    json.dump(exp, open(sys.argv[4], 'w'))
    print(sum(len(x) for x in exp), 'regions in', n, 'cases')
